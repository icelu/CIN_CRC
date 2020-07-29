#ifndef CLONE_HPP
#define CLONE_HPP

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
#include "cell.hpp"

using namespace std;


// The model of CNA evolution
class Model{
public:
    int model_ID;
    int genotype_diff;
    int growth_type;
    double fitness;
    double min_diff = 0;  // The mininal genotype differences that will cause fitness changes

    Model(){
        model_ID = 0;
        genotype_diff = 0;
        growth_type = 0;
        fitness=0;
    }

    Model(int model_ID, int genotype_diff, int growth_type, double fitness){
        this->model_ID = model_ID;
        this->genotype_diff = genotype_diff;
        this->growth_type = growth_type;
        this->fitness = fitness;
    }
};


/*
to represent a subpopulation of cells sampled from a clone
*/
class Sample{
public:
    int sample_ID;
    int clone_ID;

    Coord region_start;
    Coord region_end;
    vector<Cell_ptr> sample_cells;
    vector<string> sample_bps; // breakpoints
    map<string, int> sample_bps_count; // breakpoints and their frequency
    dpcn sample_pcn;        // using double
    double sample_loc_change[NUM_LOC]={0};

    ~Sample() = default;
    Sample(const Sample& other) = default;
    Sample& operator=(const Sample& other) = default;
    Sample& operator=(Sample&& other) = default;

    Sample(){

    }

    // create a sample with CNP
    Sample(int sample_ID, int clone_ID, Coord region_start, Coord region_end,
         vector<Cell_ptr>& sample_cells, dpcn& sample_pcn){
        this->sample_ID = sample_ID;
        this->clone_ID = clone_ID;
        this->region_start = region_start;
        this->region_end = region_end;
        this->sample_cells = sample_cells;
        this->sample_pcn = sample_pcn;
    }

    // create a sample with CNP
    Sample(int sample_ID, int clone_ID, Coord region_start, Coord region_end,
         vector<Cell_ptr>& sample_cells, double sample_loc_change[]){
        this->sample_ID = sample_ID;
        this->clone_ID = clone_ID;
        this->region_start = region_start;
        this->region_end = region_end;
        this->sample_cells = sample_cells;
        for(int i = 0; i < NUM_LOC; i++) this->sample_loc_change[i] = sample_loc_change[i];
    }

    // create a sample with breakpoints
    Sample(int sample_ID, int clone_ID, Coord region_start, Coord region_end,
         vector<Cell_ptr>& sample_cells, map<string, int>& sample_bps_count){
        this->sample_ID = sample_ID;
        this->clone_ID = clone_ID;
        this->region_start = region_start;
        this->region_end = region_end;
        this->sample_cells = sample_cells;
        this->sample_bps_count = sample_bps_count;
    }

    // Find all the breakpoints in the sample CNP
    void get_breakpoints(int loc_type = BIN, int verbose = 0){
        if(verbose > 0) cout << "\nFinding all the breakpoints in the sample CNP" << endl;
        if(loc_type == PSEUDO){
            loc2pcn(sample_loc_change, sample_pcn, verbose);
        }
        if(sample_pcn.size() == 0){
            if(verbose > 0) cout << "No bulk CNP detected" << endl;
            return;
        }
        // merge continuous regions
        dpcn::iterator it = sample_pcn.begin();
        pair<int, int> pos = (*it).first;
        int prev_chr = pos.first + 1;
        int prev_bin = pos.second;
        double prev_rcn = (*it).second;
        int start = prev_bin;
        ++it;

        for (;it != sample_pcn.end(); ++it){
            pos = (*it).first;
            int curr_chr = pos.first + 1;
            int curr_bin = pos.second;
            double curr_rcn = (*it).second;
            if(prev_chr != curr_chr || curr_bin - prev_bin != 1 || prev_rcn != curr_rcn){
                string bp = to_string(prev_chr) + "-" + to_string(prev_bin);
                sample_bps.push_back(bp);
                // start a new region
                start = curr_bin;
            }
            prev_chr = curr_chr;
            prev_bin = curr_bin;
            prev_rcn = curr_rcn;
        }
        // For regions ending at the last position
        string bp = to_string(prev_chr) + "-" + to_string(prev_bin);
        sample_bps.push_back(bp);

        sort(sample_bps.begin(), sample_bps.end());
    }

    // void get_breakpoints_from_pseudo(int verbose = 0){
    //     if(verbose > 1) cout << "Finding all the breakpoints in the sample CNP from just bins" << endl;
    //
    //     // merge continuous regions
    //     int prev_chr = 0;
    //     int prev_bin = 0;
    //     double prev_rcn = sample_loc_change[0];
    //     int start = prev_bin;
    //
    //     int curr_chr = 0;
    //     int curr_bin = 0;
    //     double curr_rcn = 0;
    //     string bp = "";
    //
    //     for(int loc = 1; loc < NUM_LOC; loc++){
    //         curr_rcn = sample_loc_change[loc];
    //         site2chr(loc, curr_chr, curr_bin, CHR_BIN_SIZE, verbose);
    //         if(prev_chr != curr_chr || curr_bin - prev_bin != 1 || prev_rcn != curr_rcn){
    //             bp = to_string(prev_chr) + "-" + to_string(prev_bin);
    //             sample_bps.push_back(bp);
    //             // start a new region
    //             start = curr_bin;
    //         }
    //         prev_chr = curr_chr;
    //         prev_bin = curr_bin;
    //         prev_rcn = curr_rcn;
    //     }
    //     // For regions ending at the last position
    //     bp = to_string(prev_chr) + "-" + to_string(prev_bin);
    //     sample_bps.push_back(bp);
    //
    //     sort(sample_bps.begin(), sample_bps.end());
    // }
    //

    // Get breakpoints without chromosome information
    void get_breakpoints_from_pseudo2(int verbose = 0){
        if(verbose > 1) cout << "Finding all the breakpoints in the sample CNP from just bins" << endl;

        // merge continuous regions
        int prev_bin = 0;
        double prev_rcn = sample_loc_change[0];
        int start = prev_bin;

        int curr_bin = 0;
        double curr_rcn = 0;
        string bp = "";

        for(int loc = 1; loc < NUM_LOC; loc++){
            curr_bin = loc;
            curr_rcn = sample_loc_change[loc];
            double diff = fabs(prev_rcn - curr_rcn);
            // cout << prev_rcn << "\t" << curr_rcn << "\t" << diff << endl;
            // cout << "cn difference is " << fabs(diff) << endl;
            if(curr_bin - prev_bin != 1 || diff >= BP_CUTOFF){
                bp = to_string(prev_bin);
                sample_bps.push_back(bp);
                // start a new region
                start = curr_bin;
            }
            prev_bin = curr_bin;
            prev_rcn = curr_rcn;
        }
        // For regions ending at the last position
        bp = to_string(prev_bin);
        sample_bps.push_back(bp);

        sort(sample_bps.begin(), sample_bps.end());
    }


    // Find all the breakpoints in the sample CNP from breakpoints with frequency higher than a cutoff
    void get_breakpoints_with_freq(double cutoff=0.1, int verbose = 0){
        if(verbose > 0) cout << "Finding all the breakpoints in the sample CNP" << endl;
        if(sample_bps_count.size() == 0){
            if(verbose > 0) cout << "No bulk CNP (breakpoints) detected" << endl;
            return;
        }
        // merge continuous regions
        for (auto bpos: sample_bps_count){
            double freq = bpos.second / sample_cells.size();
            if(freq > cutoff){
                sample_bps.push_back(bpos.first);
            }
        }
    }


    void print_breakpoints() {
        cout << "breakpoints:";
        for(auto bp: sample_bps)
            cout << "\t" << bp;
        cout << endl;
    }

};


/*
to represent a population of cells
*/
class Clone
{
public:
    int clone_ID;
    string name;

    // Used for tracking clone/deme expansion
    int parent_ID;

    vector<Cell_ptr> curr_cells;   // only available cells at present

    vector<int> id_curr_cells;      // Only store the IDs of cell to save space

    // double mutation_rate;       // assuming constant mutation rate for one clone
    // map<int, vector<int>> cnp_all[NUM_LOC];     // store all CNVs in a big array
    // vector<int> id_by_loc;      // store cell ID sorted by location to facilitate sampling
    // map<int, pair<double, double>> grates;  // birth/death rates for cells with fitness (dis)advantages

    vector<Sample*> samples;

    int ntot;   // total number of cells, used to increase cell ID
    int num_novel_mutation;
    double time_end;    // ending time of the simulation

    int model;  // the model of evolution

    Coord min_pos;
    Coord max_pos;

    // ~Clone() = default;
    Clone(const Clone& other) = default;
    Clone& operator=(const Clone& other) = default;
    Clone& operator=(Clone&& other) = default;

    Clone(){
        clone_ID = 0;
        parent_ID = -1;
        name = "";

        num_novel_mutation= 0;

        time_end = 0;
        ntot = 0;
    }


    Clone(int cID, string cname, double time_start){
        clone_ID = cID;
        parent_ID = -1;
        name = cname;
        time_end = time_start;

        num_novel_mutation= 0;

        time_end = 0;
        ntot = 0;
    }

    Clone(int cID, int pID){
        clone_ID = cID;
        parent_ID = pID;

        name = "";
        num_novel_mutation= 0;

        time_end = 0;
        ntot = 0;
    }

    ~Clone(){
        // cout << "Only delete available cells " << curr_cells.size() << endl;
        for(auto p : curr_cells){
            // cout << p->cell_ID << "\t" << p->clone_ID << endl;
            delete p;
        }

        for(auto s : samples){
            delete s;
        }
    }

    /*********************** functions related to simulating clonal growth **************************/
    void initialize(const Cell_ptr ncell, int verbose = 0){
        // this->cells.clear();
        this->id_curr_cells.clear();
        this->ntot = 1;

        this->id_curr_cells.push_back(ncell->cell_ID);
        this->time_end = ncell->time_occur;
    }

    /*
    Initialize the clone with clonal CNVs if exists
    nu: The number of new mutations
    */
    void initialize_with_cnv(const Cell_ptr ncell, Model model, int verbose = 0){
        this->curr_cells.clear();

        this->ntot = 1;
        this->num_novel_mutation = 0;

        this->model = model.model_ID;

        this->curr_cells.push_back(ncell);
        this->time_end = ncell->time_occur;
    }

    /*
    Initialize the clone with clonal CNVs if exists
    nu: The number of new mutations
    */
    void initialize_smpl(int Ns, Model model, int verbose = 0){
        this->id_curr_cells.clear();

        this->ntot = Ns;
        this->num_novel_mutation = 0;

        this->model = model.model_ID;
    }


    Cell_ptr get_cell_from_ID(vector<Cell_ptr>& cells, int cID){
        for(int i = 0; i < cells.size(); i++){
            Cell_ptr cell = cells[i];
            if(cell->cell_ID == cID) return cell;
        }
        return NULL;
    }


    int get_num_mut(){
        int sum=0;
        for (auto cell: curr_cells){
            sum += cell->num_mut;
        }
        return sum;
    }

    // Find  maximum birth rate (bmax) and maximum death rate (dmax) in the population
    double get_rmax(){
        double bmax = 0;
        double dmax = 0;
        for(unsigned int i=0; i<curr_cells.size(); i++) {
            Cell_ptr ci = curr_cells[i];
            double bi = ci->birth_rate;
            double di = ci->death_rate;
            if (bi > bmax) {
                bmax = bi;
            }
            if(di > dmax) {
                dmax = di;
            }
        }
        // cout << "Maximum birth rate: " << bmax << endl;
        // cout << "Maximum death rate: " << dmax << endl
        return bmax + dmax;
    }


    // Update the location of new-born daughter cells
    Coord get_neighbor_position(Coord p_pos, double sd){
        Coord d_pos;

        gsl_vector *mu = gsl_vector_alloc(3);
        gsl_vector *res = gsl_vector_alloc(3);
        gsl_matrix *lower_triangle = gsl_matrix_alloc(3, 3);

        gsl_vector_set(mu, 0, p_pos.x);
        gsl_vector_set(mu, 1, p_pos.y);
        gsl_vector_set(mu, 2, p_pos.z);

        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                if(i == j){
                    gsl_matrix_set(lower_triangle, i, j, sd * sd);
                }else{
                    gsl_matrix_set(lower_triangle, i, j, 0);
                }
            }
        }

        gsl_ran_multivariate_gaussian(r, mu, lower_triangle, res);

        d_pos.x = gsl_vector_get(res, 0);
        d_pos.y = gsl_vector_get(res, 1);
        d_pos.z = gsl_vector_get(res, 2);

        return d_pos;
    }


    // Update the birth/death rate of a cell (only occur when there are new mutations in daughter cells)
    // Different fitness values may apply at different locations of the genome
    void update_cell_growth(Cell_ptr dcell, const Cell_ptr ncell, const Model& model, int loc_type, int verbose = 0) {
        // Introduce selection to cells with CNAs using specified fitness
        double gdiff = 0;

        for(int i = 0; i < NUM_LOC; i++){
            gdiff += abs(dcell->loc_changes[i]);
        }
        // normalized by number of bins to avoid sharp change of net growth rate
        gdiff = gdiff / NUM_LOC;

        switch (model.growth_type) {
            case ONLY_BIRTH: {
                if(model.genotype_diff > 0){
                    dcell->birth_rate = ncell->birth_rate / (1 + gdiff * model.fitness);
                }
                else{
                    if(gdiff >= model.min_diff) dcell->birth_rate = ncell->birth_rate * (1 + model.fitness);
                }
                break;
            }
            case CHANGE_BIRTH: {   // decrease birth rate
                double new_grate = 0;
                if(model.genotype_diff > 0){
                    new_grate = (ncell->birth_rate - ncell->death_rate) / (1 + gdiff * model.fitness);
                    dcell->birth_rate = ncell->death_rate + new_grate;
                    dcell->death_rate = dcell->birth_rate - new_grate;
                }else{
                    if(gdiff >= model.min_diff){
                        new_grate = (ncell->birth_rate - ncell->death_rate) * (1 + model.fitness);
                        dcell->birth_rate = ncell->death_rate + new_grate;
                        dcell->death_rate = dcell->birth_rate - new_grate;
                    }
                }

                break;
            }
            case CHANGE_DEATH:{    // increase death rate
                double new_grate = 0;
                if(model.genotype_diff > 0){
                    new_grate = (ncell->birth_rate - ncell->death_rate) / (1 + gdiff * model.fitness);
                    double new_drate = (ncell->birth_rate - new_grate);
                    dcell->death_rate = new_drate > dcell->death_rate ? new_drate : dcell->death_rate;
                    dcell->birth_rate = dcell->death_rate + new_grate;
                }else{
                    if(gdiff >= model.min_diff){
                        new_grate = (ncell->birth_rate - ncell->death_rate) * (1 + model.fitness);
                        double new_drate = (ncell->birth_rate - new_grate);
                        dcell->death_rate = new_drate > dcell->death_rate ? new_drate : dcell->death_rate;
                        dcell->birth_rate = dcell->death_rate + new_grate;
                    }
                }

                break;
            }
            case CHANGE_BOTH:{
                double new_grate = 0;
                if(model.genotype_diff > 0){
                    new_grate = (ncell->birth_rate - ncell->death_rate) / (1 + gdiff * model.fitness);
                }else{
                    if(gdiff >= model.min_diff){
                        new_grate = (ncell->birth_rate - ncell->death_rate) * (1 + model.fitness);
                    }else{
                        new_grate = (ncell->birth_rate - ncell->death_rate);
                    }
                }
                double rn = runiform(r, 0, 1);
                if(rn < 0.5){
                    dcell->birth_rate = ncell->death_rate + new_grate;
                    dcell->death_rate = dcell->birth_rate - new_grate;
                }else{
                    double new_drate = (ncell->birth_rate - new_grate);
                    dcell->death_rate = new_drate > dcell->death_rate ? new_drate : dcell->death_rate;
                    dcell->birth_rate = dcell->death_rate + new_grate;
                }
                break;
            }
                 default: cout << "" << endl;
             }

        if(verbose > 1){
            cout << "new birth (death) rate for cell "<< dcell->cell_ID  << " is " << dcell->birth_rate << " (" << dcell->death_rate << ") with genotype difference " << gdiff << " and fitness " << model.fitness << endl;
        }
    }


     // Simulate the growth of demes. All cells have the same birth/death rate
     void grow(const Cell_ptr ncell, int size_deme, int verbose = 0, int restart = 1){
         if(restart == 1){
             initialize(ncell, verbose);
         }

         double t = 0;
         while(this->id_curr_cells.size() < size_deme) {
             // cout << this->id_curr_cells.size() << endl;

             if (this->id_curr_cells.size() == 0) {
                 t = 0;
                 initialize(ncell, verbose);
                 continue;
             }
             // Choose a random Cell from the current population
             int rindex = myrng(this->id_curr_cells.size());
             double rbrate = ncell->birth_rate;
             double rdrate = ncell->death_rate;

             double rmax = rbrate + rdrate;

             // increase time
             double tau = -log(runiform(r, 0, 1)); // an exponentially distributed random variable
             double deltaT = tau/(rmax * this->id_curr_cells.size());
             t += deltaT;

             // draw a random number
             double rb = runiform(r, 0, rmax);

             // cout << rbrate << "\t" << rdrate << "\t" << rb << "\n";
             if(rb < rbrate) {
                 // cout << "birth event" << endl;
                 this->id_curr_cells.push_back(++this->ntot);
                 this->id_curr_cells.push_back(++this->ntot);

                 this->id_curr_cells.erase(this->id_curr_cells.begin() + rindex);
             }
             // death event if b<r<b+d
             else if(rb >= rbrate && rb < rbrate + rdrate){
                 // cout << "death event" << endl;
                 this->id_curr_cells.erase(this->id_curr_cells.begin()+rindex);
             }
             else{

             }
         }
         if(verbose > 0) cout << "Generated " << this->id_curr_cells.size() << " cells before time " << t << " in clone  " << this->clone_ID  << endl;
     }


    /*
       This method simulates tumour growth with a rejection-kinetic Monte Carlo algorithm.
       intput:
        Nend -- the cell population size at the end
        ncell -- the starting cell-> Given a cell in another clone, it can mimick migragation
        model -- 0: neutral, 1: gradual, 2: punctuated
        restart -- 1: start with a new cell; 0: continue with current state
       output:
        a tree-like structure. For each Cell, its children, occurence time, birth rate, death rate
     */
     void grow_with_cnv(const Cell_ptr ncell, const Model& model, int Nend, int loc_type, double leap_size=0, int verbose = 0, int restart = 1, double tend = DBL_MAX){
         // Initialize the simulation with one cell
         if(restart == 1) initialize_with_cnv(ncell, model, verbose);

         double t = 0;  // starting time, relative to the time of last end
         int mut_ID = ncell->num_mut;
         int nu = 0; // The number of new mutations, relative to the time of last end
         int num_mut_event = 0; // count the number of times a CNA event is introduced

         if(verbose > 0) cout << "\nSimulating tumour growth with CNAs under model " << model.model_ID << " at time " << ncell->time_occur + t << endl;

         // && ncell->time_occur + t <= tend
         while(this->curr_cells.size() < Nend) {
             if (this->curr_cells.size() == 0) {
                 t = 0;
                 mut_ID = ncell->num_mut;
                 nu = 0;
                 initialize_with_cnv(ncell, model, verbose);
                 continue;
             }
             // print_all_cells(this->curr_cells, verbose);

             // Choose a random cell from the current population
             int rindex = myrng(this->curr_cells.size());
             Cell_ptr rcell = this->curr_cells[rindex];
             double rbrate = rcell->birth_rate;
             double rdrate = rcell->death_rate;
             int rID = rcell->cell_ID;

             double rmax = get_rmax();
             // increase time
             double tau = -log(runiform(r, 0, 1));  // an exponentially distributed random variable
             double deltaT = tau/(rmax * this->curr_cells.size());
             t += deltaT;

             if(ncell->time_occur + t > tend && this->curr_cells.size() >= MIN_NCELL){
                 t = tend - ncell->time_occur;
                 break;
             }

             // draw a random number
             double rb = runiform(r, 0, rmax);
             // cout << "random number " << rb << endl;
             if(rb < rbrate){
                 int rID = rcell->cell_ID;

                 Cell_ptr dcell1 = new Cell(++this->ntot, rID, ncell->time_occur + t);
                 dcell1->copy_parent((*rcell));

                 Cell_ptr dcell2 = new Cell(++this->ntot, rID, ncell->time_occur + t);
                 dcell2->copy_parent((*rcell));
                 dcell2->pos = get_neighbor_position(rcell->pos, POS_SIGMA);

                 // daughter cells aquire nu new mutations, where nu ~ Poisson(mutation_rate)
                 if (ncell->mutation_rate > 0) {
                     int nu1 = dcell1->generate_CNV_pseudo(mut_ID, t, verbose);
                     nu += nu1;
                     int nu2 = dcell2->generate_CNV_pseudo(mut_ID, t, verbose);
                     nu += nu2;

                     if(verbose > 1){
                         cout << "Number of mutations in cell " << dcell1->cell_ID << ": " << "\t" << nu1 << endl;
                         cout << "Number of mutations in cell " << dcell2->cell_ID << ": " << "\t" << nu2 << endl;
                     }

                     if(model.fitness != 0)
                     {
                         if(verbose > 1){
                             cout << "Update the grow parameters of daughter cells " << endl;
                         }
                         if(nu1 > 0) update_cell_growth(dcell1, ncell, model, loc_type, verbose);
                         if(nu2 > 0) update_cell_growth(dcell2, ncell, model, loc_type, verbose);
                     }
                 }
                 // Remove the parent cell from the list of current cells
                 if(rID != 1){
                     delete (this->curr_cells[rindex]);
                     this->curr_cells[rindex] = NULL;
                 }

                 this->curr_cells.erase(this->curr_cells.begin() + rindex);
                 this->curr_cells.push_back(dcell1);
                 this->curr_cells.push_back(dcell2);
             }
             // death event if b<r<b+d
             else if(rb >= rbrate && rb < rbrate + rdrate) {
                 // cout << " death event" << endl;
                 if(rID != 1){
                     delete (this->curr_cells[rindex]);
                     this->curr_cells[rindex] = NULL;
                 }
                 this->curr_cells.erase(this->curr_cells.begin() + rindex);
             }else{

             }
         }

         this->time_end += t;
         this->num_novel_mutation += nu;

         if(verbose > 0){
             cout << "Generated " << nu << " mutations during time " << t << endl;
             if(verbose > 1) print_all_cells(this->curr_cells, verbose);
         }
     }


     // /*
     // generate CNAs in a pseudo (abstract) way and store in a global data structure
     // */
     // int generate_CNV_global(int cell_ID, int verbose = 0){
     //     int nu = gsl_ran_poisson(r, mutation_rate);
     //     if(nu <= 0) return nu;
     //
     //     if(verbose > 1 && nu > 0) cout << "Generating " << nu << " CNAs under rate " << mutation_rate << endl;
     //
     //     gsl_ran_discrete_t* dis_loc = gsl_ran_discrete_preproc(NUM_LOC, LOC_PROBS);
     //     double u = 0;
     //     int start = gsl_ran_discrete(r, dis_loc);
     //     double msize = 0;
     //     int len = 0;
     //     int end = 0;
     //
     //     for (int j=0; j < nu; j++) {
     //         u = runiform(r, 0, 1);
     //         len = 1;    // add 1 to avoid 0 length
     //         if(u < 0.5){  // gain
     //             msize = MEAN_GAIN_SIZE;
     //             if(msize > 1) len = ceil(gsl_ran_exponential(r, msize));
     //             end = start + len;
     //             if(end > NUM_LOC) end = NUM_LOC;
     //             for(int i = start; i < end; i++){
     //                 cna[i][cell_ID]++;
     //             }
     //
     //         }else{
     //             msize = MEAN_LOSS_SIZE;
     //             if(msize > 1) len = ceil(gsl_ran_exponential(r, msize));
     //             end = start + len;
     //             if(end > NUM_LOC) end = NUM_LOC;
     //             for(int i = start; i < end; i++){
     //                 cna[i][cell_ID]--;
     //             }
     //         }
     //         this->num_mut++;
     //     }
     //
     //     return nu;
     // }


     /*
        This method simulates tumour growth with a rejection-kinetic Monte Carlo algorithm, simplified to simulate a very large population (assuming starting from a large population).
        TODO: seem not necessary to simulate to a large number
        intput:
         Nend -- the cell population size at the end
         Ns -- the number of starting cells
         model -- 0: neutral, 1: gradual, 2: punctuated
         restart -- 1: start with a new cell; 0: continue with current state
        output:
         a tree-like structure. For each Cell, its children, occurence time, birth rate, death rate
      */
//       void grow_with_cnv_smpl(int Ns, const Model& model, int Nend, int loc_type, double leap_size=0, int verbose = 0, int restart = 1, double tend = DBL_MAX){
//           // Initialize the simulation with one cell
//           if(restart == 1) initialize_smpl(Ns, model, verbose);
//
//           double t = 0;  // starting time, relative to the time of last end
//           int mut_ID = ncell->num_mut;
//           int nu = 0; // The number of new mutations, relative to the time of last end
//           int num_mut_event = 0; // count the number of times a CNA event is introduced
//
//           if(verbose > 0) cout << "\nSimulating tumour growth with CNAs under model " << model.model_ID << " at time " << this->time_occur + t << endl;
//
//           // && ncell->time_occur + t <= tend
//           while(this->ntot < Nend) {
//               if (this->ntot == 0) {
//                   t = 0;
//                   mut_ID = ncell->num_mut;
//                   nu = 0;
//                   initialize_smpl(Ns, model, verbose);
//                   continue;
//               }
//               // print_all_cells(this->curr_cells, verbose);
//
//               // Choose a random cell from the current population
//               int rID = myrng(this->ntot);
//               double rbrate = grates[rID][0];
//               double rdrate = grates[rID][1];
//
//               double rmax = get_rmax();
//               // increase time
//               double tau = -log(runiform(r, 0, 1));  // an exponentially distributed random variable
//               double deltaT = tau/(rmax * this->curr_cells.size());
//               t += deltaT;
//
//               if(this->time_occur + t > tend && this->ntot >= MIN_NCELL){
//                   t = tend - ncell->time_occur;
//                   break;
//               }
//
//               // draw a random number
//               double rb = runiform(r, 0, rmax);
//               // cout << "random number " << rb << endl;
//               if(rb < rbrate){
//                   // this->id_curr_cells.push_back(++this->ntot);
//                   //
//                   // this->id_curr_cells.push_back(++this->ntot);
//                   // ++this->ntot;
//
//                   // daughter cells aquire nu new (different) mutations, where nu ~ Poisson(mutation_rate)
//                   if (this->mutation_rate > 0) {
//                       int nu1 = generate_CNV_global(rID, verbose);  // One of the daughter cells get the ID of parent
//                       nu += nu1;
//                       int cell_ID2 = ++this->ntot;
//                       // copy the original mutations to daughter cells
//                       for(int i = 0; i < NUM_LOC; i++){
//                           cnp_all[i][cell_ID2] = cna[i][rID];
//                       }
//                       int nu2 = generate_CNV_global(cell_ID2, verbose);
//                       nu += nu2;
//
//                       if(verbose > 1){
//                           cout << "Number of mutations in cell " << dcell1->cell_ID << ": " << "\t" << nu1 << endl;
//                           cout << "Number of mutations in cell " << dcell2->cell_ID << ": " << "\t" << nu2 << endl;
//                       }
//
//                       if(model.fitness != 0)
//                       {
//                           if(verbose > 1){
//                               cout << "Update the grow parameters of daughter cells " << endl;
//                           }
//                       //     double gdiff;
//                       //     if(nu1 > 0){
//                       //         gdiff = 0;
//                       //         // update_cell_growth(dcell1, ncell, model, loc_type, verbose);
//                       //         for(int i = 0; i < NUM_LOC; i++){
//                       //             gdiff += abs(dcell->loc_changes[i]);
//                       //         }
//                       //         // normalized by number of bins to avoid sharp change of net growth rate
//                       //         gdiff = gdiff / NUM_LOC;
//                       //     }
//                       //     if(nu2 > 0) update_cell_growth(dcell2, ncell, model, loc_type, verbose);
//                       // }
//                   }
//               }
//               // death event if b<r<b+d
//               else if(rb >= rbrate && rb < rbrate + rdrate) {
//                   // cout << " death event" << endl;
//                   this->curr_cells.erase(this->curr_cells.begin() + rindex);
//
//               }else{
//
//               }
//           }
//
//           this->time_end += t;
//           this->num_novel_mutation += nu;
//
//           if(verbose > 0){
//               cout << "Generated " << nu << " mutations during time " << t << endl;
//               if(verbose > 1) print_all_cells(this->curr_cells, verbose);
//           }
//       }
//
//
     /***************************************************************************************************************************/


     /*********************** functions related to extracting information from simulated data **************************/
     // Get the CNPs of bulk sample (no noise)
     void set_bulk_cnp_from_pseudo(double bulk_cnp[NUM_LOC], const vector<Cell_ptr>& cells) {
         for(int i = 0; i < NUM_LOC; i++){
             for(auto cell : cells) {
                 bulk_cnp[i] += cell->loc_changes[i];
             }
             double avg_cn = bulk_cnp[i] / cells.size();
             double rcn =  rint(avg_cn * DEC_PLACE) / DEC_PLACE;  // 1 decimal place
             bulk_cnp[i] = rcn;
         }
     }


    /********************************************************************************/


    /*********************** functions related to output **************************/
    void get_cell_loc(const vector<Cell_ptr>& cells, int verbose = 0){
        // Check the location of all cells
        if(verbose > 1) cout << "Printing the location of all cells" << endl;
        // Find the min and max of each dimension
        double min_x=INT_MAX, max_x=INT_MIN;
        double min_y=INT_MAX, max_y=INT_MIN;
        double min_z=INT_MAX, max_z=INT_MIN;
        for(auto cell : cells){
            Coord pos = cell->pos;
            if(pos.x < min_x){
                min_x = pos.x;
            }
            if(pos.y < min_y){
                min_y = pos.y;
            }
            if(pos.z < min_z){
                min_z = pos.z;
            }
            if(pos.x > max_x){
                max_x = pos.x;
            }
            if(pos.y > max_y){
                max_y = pos.y;
            }
            if(pos.z > max_z){
                max_z = pos.z;
            }
            if(verbose > 1) cout << "position for cell " << cell->cell_ID << " is: " << pos.x << "\t" << pos.y << "\t" << pos.z << "\n";
        }
        min_pos = Coord{min_x, min_y, min_z};
        max_pos = Coord{max_x, max_y, max_z};
        if(verbose > 0){
            cout << "The minimal x,y,z positions of all cells are: " << min_x << "\t" << min_y << "\t" << min_z << "\n";
            cout << "The maximal x,y,z positions of all cells are: " << max_x << "\t" << max_y << "\t" << max_z << "\n";
        }
    }

    /*
       This method prints out the copy numbers of each final cell in a clone
     */
    void print_all_cells(vector<Cell_ptr> cells, int verbose = 0){
        int num_cell = cells.size();
        if(verbose > 0) cout << "Printing " << num_cell << " cells" << endl;

        for(unsigned int i = 0; i < num_cell; i++) {
                Cell_ptr cell = cells[i];
                cout << cell->cell_ID << "\t" << cell->parent_ID << endl;
        }
    }


    /*
    This function prints summary informaton about the simulated clone
     */
    void print_summary(string outfile) {
        ofstream out;
        out.setf(ios::fixed);
        out.setf(ios::showpoint);
        out.precision(9);
        out.open(outfile);

        double lambda = log(2);
        out << "Information for host population:" << endl;
        for(auto cell : curr_cells){
            if(cell->clone_ID == 0){
                lambda = cell->birth_rate - cell->death_rate;
                out << "\tCell ID: " << cell->cell_ID << endl;
                out << "\tMutation rate: " << cell->mutation_rate << endl;

                out << "\tBirth rate: " << cell->birth_rate << endl;
                out << "\tDeath rate: " << cell->death_rate << endl;

                out << "\tEffective mutation rate (μ/β): " << cell->mutation_rate / ((cell->birth_rate-cell->death_rate)/cell->birth_rate) << endl;
                out << endl;
                if(model == 0) break;     // same rates under neutral evolution
            }
        }
        out << endl;
        out << "\tNumber of de novo mutations: "<< num_novel_mutation << endl;
        out << "\tEnd time of simulation: "<< time_end << endl;

        out.close();
    }


     // Print out all informaton for one clone
     void print_single_clone(Cell_ptr start_cell, Model model, string outdir, string suffix, int loc_type=BIN, int verbose = 1){
         double lambda = start_cell->birth_rate - start_cell->death_rate;
         double tend = log(curr_cells.size())/(lambda); // The latest time that a subclone occurs
         cout << "\nPrinting informaton for clone " << clone_ID << endl;

         cout << "Model of evolution: " << model.model_ID << endl;
         cout << "Initial Net growth rate: " << lambda << endl;
         cout << "Mutation rate: " << start_cell->mutation_rate << endl;
         cout << "Estimated simulation finish time (tumor doublings): " << tend << endl;

         string outfile = "";
         outfile = outdir + "summary" + suffix;
         cout << "Printing summary" << endl;
         print_summary(outfile);
     }


     /*********************** functions related to multi-region sampling **************************/
      // Select multiple regions from a clone according to positions
      void sample_multi_regions(int nsample, const vector<Coord>& region_start, const vector<Coord>& region_end, int loc_type=BIN, int verbose = 0){
          for(int i = 0; i < nsample; i++){
              Coord pos_s = region_start[i];
              Coord pos_e = region_end[i];

              // Extract all cells in this region
              vector<Cell_ptr> sample_cells;
              for(auto cell : curr_cells){
                  if(cell->is_sampled == 0){
                      Coord pos_curr = cell->pos;
                      if((pos_curr.x >= pos_s.x && pos_curr.x <= pos_e.x) &&
                         (pos_curr.y >= pos_s.y && pos_curr.y <= pos_e.y) &&
                         (pos_curr.z >= pos_s.z && pos_curr.z <= pos_e.z)){
                          cell->is_sampled = 1;
                          sample_cells.push_back(cell);
                      }
                  }
              }

              // Find the average CNP of these cells
              double bulk_loc_change[NUM_LOC] = {0};
              set_bulk_cnp_from_pseudo(bulk_loc_change, sample_cells);
              Sample* s = new Sample(i, this->clone_ID, pos_s, pos_e, sample_cells, bulk_loc_change);
              // cout << "Generating samples from just bins" << endl << endl;
              samples.push_back(s);

              if(verbose > 0){
                  cout << "Generating sample " << i + 1 << " in clone " << this->clone_ID << '\n';
                  cout << "start position: " << pos_s.x << "\t" << pos_s.y << "\t" << pos_s.z << "\n";
                  cout << "end position: " << pos_e.x << "\t" << pos_e.y << "\t" << pos_e.z << "\n";
                  cout << "There are " << sample_cells.size() << " cells" << endl;

                  if(verbose > 1){
                      cout << "Cells in this clone: ";
                      for(auto cell : sample_cells){
                          cout << cell->cell_ID << '\t';
                      }
                      cout << '\n';
                  }
              }
          }
      }

};


#endif
