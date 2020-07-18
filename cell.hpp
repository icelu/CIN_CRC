#ifndef CELL_HPP
#define CELL_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <string>

#include <assert.h>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include "util.hpp"


using namespace std;


struct Coord{
    double x;
    double y;
    double z;

    Coord() {}
    Coord(double x, double y, double z) : x(x), y(y), z(z) {}

    // sort on all axis
    // bool operator<(const Coord &o) const {
    //     if (x != o.x) {
    //      return x < o.x;
    //     }
    //     if (y != o.y) {
    //      return y < o.y;
    //     }
    //     return z < o.z;
    // }

    bool operator<(const Coord &o) const {
        return x < o.x;
    }
};


class Mutation
{
public:
    int mut_ID;

    int type;   // Mutation type. 0: SNV 1: CNV
    double time_occur;

    int cell_ID;
    int chr;  // chromosome on which the mutation occurs
    int arm;    // 0: no informaton, 1: p; 2: q
    int reciprocal;

    int start;  // start position
    int end;    // end position
    int size;

    int number;  // copy of this mutation

    ~Mutation() = default;
    Mutation(const Mutation& other) = default;
    Mutation(Mutation&& other) = default;
    Mutation& operator=(const Mutation& other) = default;
    Mutation& operator=(Mutation&& other) = default;

    Mutation(){
        mut_ID = 0;
        time_occur = 0;
        number = 1;
    }


    Mutation(int mut_ID, double time_occur){
        this->mut_ID = mut_ID;
        this->time_occur = time_occur;
        this->number = 1;
    }

    Mutation(int mut_ID, double time_occur, int chr, int arm, int type, int reciprocal){
        this->mut_ID = mut_ID;
        this->time_occur = time_occur;

        this->chr = chr;
        this->arm = arm;
        this->type = type;
        this->reciprocal = reciprocal;

        this->number = 1;
    }

    Mutation(int mut_ID, double time_occur, int chr, int start, int end, int arm, int type, int reciprocal){
        this->mut_ID = mut_ID;
        this->time_occur = time_occur;

        this->chr = chr;
        this->start = start;
        this->end = end;
        this->arm = arm;
        this->type = type;
        this->reciprocal = reciprocal;

        this->number = 1;
    }
};

class Cell;

// typedef shared_ptr<Cell> Cell_ptr;
typedef Cell* Cell_ptr;

class Cell
{
public:
    int cell_ID;
    int parent_ID;
    int clone_ID;

    double time_occur;

    int is_sampled;

    // parameters related to cell growth
    double birth_rate;
    double death_rate;

    Coord pos{0, 0, 0};

    // parameters related to mutation generations
    double mutation_rate;

    // parameters related to storing mutations
    int num_mut;    // number of mutations
    int loc_changes[NUM_LOC] = {0};


    ~Cell() = default;
    Cell(const Cell& other) = default;
    Cell(Cell&& other) = default;
    Cell& operator=(const Cell& other) = default;
    Cell& operator=(Cell&& other) = default;

    Cell() {
        cell_ID = 0;
        parent_ID = 0;
        clone_ID = 0;

        time_occur = 0;
        is_sampled = 0; // whether the cell is sampled or not

        birth_rate = log(2);
        death_rate = 0;

        mutation_rate = 0;
        num_mut = 0;
    }


    Cell(int cell_ID, int parent_ID, double time_occur) {
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;

        this->time_occur = time_occur;
        this->is_sampled = 0;

        this->birth_rate = birth_rate;
        this->death_rate = death_rate;

        this->mutation_rate = 0;
        this->num_mut = 0;
    }

    Cell(int cell_ID, int parent_ID, double birth_rate, double death_rate, double mutation_rate, double time_occur){
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;
        this->clone_ID = 0;

        this->time_occur = time_occur;
        this->is_sampled = 0;

        this->birth_rate = birth_rate;
        this->death_rate = death_rate;

        this->mutation_rate = mutation_rate;
        this->num_mut = 0;
    }


    bool operator<(const Cell &o) const {
        return pos < o.pos;
    }

    Cell_ptr get_parent(vector<Cell_ptr> cells){
        for(int i = 0; i < cells.size(); i++){
            Cell_ptr cell = cells[i];
            if(cell->cell_ID == parent_ID) return cell;
        }
        return NULL;
    }

    void copy_parent(const Cell& ncell){
        this->clone_ID = ncell.clone_ID;

        this->time_occur = ncell.time_occur;
        this->is_sampled = ncell.is_sampled;
        this->pos = ncell.pos;

        this->birth_rate = ncell.birth_rate;
        this->death_rate = ncell.death_rate;

        this->mutation_rate = ncell.mutation_rate;
        this->num_mut = ncell.num_mut;

        for(int i = 0; i < NUM_LOC; i++){
            loc_changes[i] = ncell.loc_changes[i];
        }
    }


    void update_loc(double offset_x, double offset_y, double offset_z) {
        pos.x += offset_x;
        pos.y += offset_y;
        pos.z += offset_z;
    }

    /*********************** functions related to mutation generations **************************/
    /*
    generate CNAs in a pseudo (abstract) way to save time
    */
    int generate_CNV_pseudo(int& mut_ID, double time_occur, int verbose = 0){
        int nu = gsl_ran_poisson(r, mutation_rate);
        if(nu <= 0) return nu;

        if(verbose > 1 && nu > 0) cout << "Generating " << nu << " CNAs under rate " << mutation_rate << " in cell " << cell_ID << endl;

        gsl_ran_discrete_t* dis_loc = gsl_ran_discrete_preproc(NUM_LOC, LOC_PROBS);
        double u = 0;
        int start = gsl_ran_discrete(r, dis_loc);
        double msize = 0;
        int len = 0;
        int end = 0;

        for (int j=0; j < nu; j++) {
            u = runiform(r, 0, 1);
            len = 1;    // add 1 to avoid 0 length
            if(u < 0.5){  // gain
                msize = MEAN_GAIN_SIZE;
                if(msize > 1) len = ceil(gsl_ran_exponential(r, msize));
                end = start + len;
                if(end > NUM_LOC) end = NUM_LOC;
                for(int i = start; i < end; i++){
                    loc_changes[i]++;
                }

            }else{
                msize = MEAN_LOSS_SIZE;
                if(msize > 1) len = ceil(gsl_ran_exponential(r, msize));
                end = start + len;
                if(end > NUM_LOC) end = NUM_LOC;
                for(int i = start; i < end; i++){
                    loc_changes[i]--;
                }
            }
            this->num_mut++;
        }

        return nu;
    }

    /*********************************************************************************/

    /*********************** function related to output generations **************************/

    void print_cell_info(){
        cout << "All infomation of cell " << this->cell_ID << endl;
        cout << "\t parent " << this->parent_ID << endl;
        cout << "\t in clone " << this->clone_ID << endl;

        cout << "\t is_sampled " << this->is_sampled << endl;

        cout << "\t position " << this->pos.x << endl;

        cout << "\t birth_rate " << this->birth_rate << endl;
        cout << "\t death_rate " << this->death_rate << endl;

        cout << "\t mutation_rate " << this->mutation_rate << endl;
    }
};


#endif
