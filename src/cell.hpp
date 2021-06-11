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

    // int flag;   // whether the cell is alive or not. 0:new, 1:divided, -1: death
    double time_occur;

    int is_sampled;

    // parameters related to cell growth
    double birth_rate;
    double death_rate;

    Coord pos{0, 0, 0};   // used in simulations of multi-region samples

    // parameters related to mutation generations
    double mutation_rate;

    // parameters related to storing mutations
    int num_mut;    // number of total mutations
    vector<string> muts;  // IDs of mutations in one cell, used to count #unique mutations in a random population of cells

    int loc_changes[NUM_LOC] = {0};
    int chr_changes[NUM_CHR] = {0};   // used when chr-level CNAs are simulated

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

        for(int i = 0; i < NUM_LOC; i++){
            loc_changes[i] = START_KARYOTYPE[i];
        }

        for(int i = 0; i < NUM_CHR; i++){
            chr_changes[i] = START_KARYOTYPE_CHR[i];
        }
    }


    // called when creating daughter cells
    // other properties will be copied from parent
    Cell(int cell_ID, int parent_ID, double time_occur) {
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;

        this->time_occur = time_occur;
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

        for(int i = 0; i < NUM_LOC; i++){
            loc_changes[i] = START_KARYOTYPE[i];
        }

        for(int i = 0; i < NUM_CHR; i++){
            chr_changes[i] = START_KARYOTYPE_CHR[i];
        }
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

        this->is_sampled = ncell.is_sampled;
        this->pos = ncell.pos;

        this->birth_rate = ncell.birth_rate;
        this->death_rate = ncell.death_rate;

        this->mutation_rate = ncell.mutation_rate;
        this->num_mut = ncell.num_mut;

        this->muts = ncell.muts;

        for(int i = 0; i < NUM_LOC; i++){
            loc_changes[i] = ncell.loc_changes[i];
        }

        for(int i = 0; i < NUM_CHR; i++){
            chr_changes[i] = ncell.chr_changes[i];
        }
    }


    void update_loc(double offset_x, double offset_y, double offset_z) {
        pos.x += offset_x;
        pos.y += offset_y;
        pos.z += offset_z;
    }


    // compute distance to optimum chr-level karyotype with potential weights
    double get_dist_chr_level(){
      double dist = 0;

      for(int i = 0; i < NUM_CHR; i++){
        double weight = 1;  // add weight for locations different from optimum karyotype
        double rcn = abs(this->chr_changes[i] - OPT_KARYOTYPE_CHR[i]);
        if((this->chr_changes[i] != NORM_PLOIDY && START_KARYOTYPE_CHR[i] == OPT_KARYOTYPE_CHR[i]) ||
            (this->chr_changes[i] != OPT_KARYOTYPE_CHR[i] && START_KARYOTYPE_CHR[i] != OPT_KARYOTYPE_CHR[i])){
          weight = WEIGHT_OPTIMUM;
        }
        dist += weight * rcn;
     }

      return dist;
    }


    // compute distance to optimum bin-level karyotype with potential weights
    double get_dist_bin_level(){
      double dist = 0;

      for(int i = 0; i < NUM_LOC; i++){
        double weight = 1;  // add weight for locations different from optimum karyotype
        double rcn = abs(this->loc_changes[i] - OPT_KARYOTYPE[i]);
        if((this->loc_changes[i] != NORM_PLOIDY && START_KARYOTYPE[i] == OPT_KARYOTYPE[i]) ||
            (this->loc_changes[i] != OPT_KARYOTYPE[i] && START_KARYOTYPE[i] != OPT_KARYOTYPE[i])){
          weight = WEIGHT_OPTIMUM;
        }
        dist += weight * rcn;
      }

      return dist;
    }


    /*********************** functions related to mutation generations **************************/
    /*
    generate CNAs in a pseudo (abstract) way to save time
    1. choose locations of mutations
    2. randomly select gain or loss
    not follow infinite site assumption
    */
    int generate_CNV_pseudo(int& mut_ID, double time_occur, int verbose = 0){
        int nu = gsl_ran_poisson(r, mutation_rate);

        if(nu <= 0) return nu;

        if(verbose > 1 && nu > 0) cout << "Generating " << nu << " CNAs under mutation rate " << mutation_rate << " in cell " << cell_ID << endl;

        gsl_ran_discrete_t* dis_loc;
        double u = 0;
        int start = 0;
        int len = 0;
        int end = 0;
        int chr = 0;

        // this chunk was outside loop previously, should have litte effect at bin level as the mutation rate is very low (0, 0.5)
        // most division will only introduce at most 1 mutation due to low rate
        // multiple mutations may overlap at some bins
        // save computation time when simulating a larger population
        // dis_loc = gsl_ran_discrete_preproc(NUM_LOC, LOC_PROBS);
        // start = gsl_ran_discrete(r, dis_loc);

        for (int j = 0; j < nu; j++) {
            dis_loc = gsl_ran_discrete_preproc(NUM_LOC, LOC_PROBS);
            start = gsl_ran_discrete(r, dis_loc);

            u = runiform(r, 0, 1);
            len = 1;    // add 1 to avoid 0 length
            if(u < 0.5){  // gain
                if(real_gain_sizes.size() > 0){
                  len = real_gain_sizes[myrng(real_gain_sizes.size())];
                }else{
                  len = ceil(gsl_ran_exponential(r, MEAN_GAIN_SIZE));
                }
                // cout << "gain size " << len << endl;
                end = start + len;

                if(end > NUM_LOC) end = NUM_LOC;

                // update loc_changes even when only chr-level CNAs are simulated to reuse some postprocessing functions
                for(int i = start; i < end; i++){
                    loc_changes[i]++;
                }
            }else{
                if(real_loss_sizes.size() > 0){
                  len = real_loss_sizes[myrng(real_loss_sizes.size())];
                }else{
                  len = ceil(gsl_ran_exponential(r, MEAN_LOSS_SIZE));
                }
                // cout << "loss size " << len << endl;
                end = start + len;

                if(end > NUM_LOC) end = NUM_LOC;

                for(int i = start; i < end; i++){
                    if(loc_changes[i] == 0){
                        // cout << "No more copies to lost!" << endl;
                        continue;
                    }
                    loc_changes[i]--;
                }
            }
            mut_ID = mut_ID + 1;
            this->muts.push_back(to_string(this->clone_ID) + "_" + to_string(mut_ID));
            this->num_mut++;
        }

        return nu;
    }


    /*
    generate chr-level CNAs in a pseudo (abstract) way to save time
    1. choose locations of mutations
    2. randomly select gain or loss
    not follow infinite site assumption
    */
    int generate_CNV_chr_level(int& mut_ID, double time_occur, int verbose = 0){
        int nu = gsl_ran_poisson(r, mutation_rate);

        if(nu <= 0) return nu;

        if(verbose > 1 && nu > 0) cout << "Generating " << nu << " chr-level CNAs under mutation rate " << mutation_rate << " in cell " << cell_ID << endl;

        gsl_ran_discrete_t* dis_loc;
        double u = 0;
        int start = 0;
        int end = 0;
        int chr = 0;

        for (int j = 0; j < nu; j++) {
            // only mutations at different chromosomes in a single division for large jumps
            // randomly select a chromosome
            dis_loc = gsl_ran_discrete_preproc(NUM_CHR, CHR_PROBS);
            chr = gsl_ran_discrete(r, dis_loc);
            // compute the start and end position (always change loc_changes for consistent with other computations)
            if(chr == 0){
              start = 0;
            }else{
              start = CHR_SSIZE[chr - 1];
            }
            end = CHR_SSIZE[chr];
            if(verbose > 1 && nu > 0) cout << " CNAs on chr " << chr << " start " << start << " end " << end << endl;

            u = runiform(r, 0, 1);
            if(u < 0.5){  // gain
              chr_changes[chr]++;
              if(end > NUM_LOC) end = NUM_LOC;
              // update loc_changes even when only chr-level CNAs are simulated to reuse some postprocessing functions
              for(int i = start; i < end; i++){
                  loc_changes[i]++;
              }
            }else{  // loss
              if(chr_changes[chr] > 0)  chr_changes[chr]--;
              if(end > NUM_LOC) end = NUM_LOC;
              for(int i = start; i < end; i++){
                  if(loc_changes[i] == 0){
                      // cout << "No more copies to lost!" << endl;
                      continue;
                  }
                  loc_changes[i]--;
              }
            }
            mut_ID = mut_ID + 1;
            this->muts.push_back(to_string(this->clone_ID) + "_" + to_string(mut_ID));
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
