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



// Node of a binary lineage tree, used for sampling
struct node {
    int data;
    node* left;
    node* right;
    int flag; // whether the node is available (dead) or not

    node(int data)
    {
        this->data = data;
        this->flag = 0;
        left = NULL;
        right = NULL;
    }

    node(int data, int flag)
    {
        this->data = data;
        this->flag = flag;
        left = NULL;
        right = NULL;
    }
};

node* search_node(int key, node* p)
{
    // cout << "search node " << key << endl;
    if(p == NULL) return NULL;

    if(key == p->data){
        return p;
    }
    else{
        node* t = search_node(key, p->left);
        if(t == NULL){
            return search_node(key, p->right);
        }
        else{
            return t;
        }
    }
}

// find all leaves below an internal node
void get_leaves_below(node* p, vector<int>& nodes){
    if(p == NULL) return;
    if(p->left == NULL && p->right == NULL){
        nodes.push_back(p->data);
    }else{
        get_leaves_below(p->left, nodes);
        get_leaves_below(p->right, nodes);
    }
}

// find all leaves below an internal node which are available
void get_leaves_below_avail(node* p, vector<int>& nodes){
    if(p == NULL) return;
    if(p->left == NULL && p->right == NULL && p->flag >= 0){
        nodes.push_back(p->data);
    }else{
        get_leaves_below_avail(p->left, nodes);
        get_leaves_below_avail(p->right, nodes);
    }
}

void destroy_tree(node* root){
    if(root != NULL)
    {
      destroy_tree(root->left);
      destroy_tree(root->right);
      delete root;
    }
}

// The model of CNA evolution
class Model{
public:
    int model_ID;
    int genotype_diff;    // whether or not to consider genotype differences
    int growth_type;   // how cells grow under selection
    double fitness;    // the strength of selection
    int use_alpha;    // use the model based on alpha, in which fitness is alpha; Otherwise fitness is classic selection coefficient
    double min_diff = 0;  // The mininal genotype differences that will cause fitness changes
    int norm_by_bin = 1;  // normalized genotype differences by the number of bins

    Model(){
        model_ID = 0;
        genotype_diff = 0;
        growth_type = 0;
        use_alpha = 0;
        fitness = 0;
    }

    Model(int model_ID, int genotype_diff, int growth_type, double fitness){
        this->model_ID = model_ID;
        this->genotype_diff = genotype_diff;
        this->growth_type = growth_type;
        this->fitness = fitness;
    }

    Model(int model_ID, int genotype_diff, int growth_type, double fitness, int use_alpha, int norm_by_bin = 1){
        this->model_ID = model_ID;
        this->genotype_diff = genotype_diff;
        this->growth_type = growth_type;
        this->fitness = fitness;
        this->use_alpha = use_alpha;
        this->norm_by_bin = norm_by_bin;
    }
};


/*
A Sample represents a subpopulation of cells sampled from a clone
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
A Clone represents a population of cells
*/
class Clone
{
public:
    int clone_ID;   // the default initial clone has ID 0
    string name;

    // Used for tracking clone/deme expansion
    int parent_ID;

    vector<Cell_ptr> curr_cells;   // only available cells at present
    map<int, vector<map<int, pcn>>> cnp_by_time;   // CNPs at different population sizes
    map<int, map<int, double*>> lchange_by_time;   // CNAs at different population sizes (format used for computing summary statistics)
    map<int, double> meanfit_by_time ;
    map<int, double> maxfit_by_time;
    map<int, double> minfit_by_time;
    // map<int, vector<double>> fits_by_time;
    map<int, double> dist_by_time;
    // map<int, double> dm_by_time;

    // added for simulating glands
    vector<Cell_ptr> cells;   // all cells at present
    // vector<string> bps; // breakpoints

    vector<int> id_curr_cells;      // Only store the IDs of cell to save space

    // double mutation_rate;       // assuming constant mutation rate for one clone
    // map<int, vector<int>> cnp_all[NUM_LOC];     // store all CNVs in a big array
    // vector<int> id_by_loc;      // store cell ID sorted by location to facilitate sampling
    // map<int, pair<double, double>> grates;  // birth/death rates for cells with fitness (dis)advantages

    vector<Sample*> samples;

    int ntot;   // total number of cells, used to increase cell ID
    int num_novel_mutation;
    // set<int> muts;  // unique set of mutations present in all cells of the clone
    map<string, double> muts_freq;

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
        if(cells.size() > 0){
          // when all cells are recorded
          for(auto p : cells){
              // cout << p->cell_ID << "\t" << p->clone_ID << endl;
              delete p;
          }
        }else{
          for(auto p : curr_cells){
              // cout << p->cell_ID << "\t" << p->clone_ID << endl;
              delete p;
          }
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
    void initialize_with_cnv_cmpl(const Cell_ptr ncell, Model model, int verbose = 0){
        this->curr_cells.clear();
        this->cells.clear();

        this->ntot = 1;
        this->num_novel_mutation = 0;

        this->model = model.model_ID;

        this->curr_cells.push_back(ncell);
        this->cells.push_back(ncell);
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

    // Find mean copy number for each chromosome
    void get_mean_karyotype(const vector<Cell_ptr> cells, vector<double>& mean_karyotype, int verbose = 0){
      int nloc = NUM_LOC;
      if(CHR_CNA == 1){
        nloc = NUM_CHR;
        for(int i = 0; i < NUM_CHR; i++){
          double mean_cn = 0;
          for(auto c : cells){
            mean_cn += c->chr_changes[i];
          }
          mean_cn = mean_cn / cells.size();
          mean_karyotype.push_back(mean_cn);
        }
      }else{
        for(int i = 0; i < NUM_LOC; i++){
          double mean_cn = 0;
          for(auto c : cells){
            mean_cn += c->loc_changes[i];
          }
          mean_cn = mean_cn / cells.size();
          mean_karyotype.push_back(mean_cn);
        }
      }

      if(verbose > 1){
        cout << "Mean karyotype:";
        for(int i = 0; i < nloc; i++){
          cout << "\t" << mean_karyotype[i];
        }
        cout << endl;
      }
    }


    // get fitness statistics of a clone
    // here fitness is actually growth rate, computed over all available cells in the clone
    void get_fitness_stats() {
      double avg_fitness = 0.0;
      double max_fitness = 0.0;
      double min_fitness = 1000000.0;
      double avg_dist = 0.0;
      // double dm = 0.0;
      // vector<double> fits;
      for(auto c : this->curr_cells){
        double s = c->birth_rate - c->death_rate;
        // fits.push_back(s);
        avg_fitness += s;
        if(s < min_fitness){
          min_fitness = s;
        }
        if(s > max_fitness){
          max_fitness = s;
        }
        if(CHR_CNA == 1){
          avg_dist += c->get_dist_chr_level();
        }else{
          avg_dist += c->get_dist_bin_level();
        }
        // fitnesses.push_back(s);
      }


      double currsize = this->curr_cells.size();
      // this->fits_by_time[currsize] = fits;
      avg_fitness = avg_fitness / currsize;
      avg_dist = avg_dist / currsize;
      // // if(verbose > 0){

      //   vector<double> mean_karyotype;
      //   get_mean_karyotype(this->curr_cells, mean_karyotype, verbose);
      //   if(CHR_CNA == 1){
      //     for(int i = 0; i < NUM_CHR; i++){
      //       dm += abs(mean_karyotype[i] - OPT_KARYOTYPE_CHR[i]);
      //     }
      //   }else{
      //     for(int i = 0; i < NUM_LOC; i++){
      //       dm += abs(mean_karyotype[i] - OPT_KARYOTYPE[i]);
      //     }
      //   // }

      // cout << "Record CNPs at size " << currsize << " with " << samples.size() << " cells, average fitness " << avg_fitness << ", average distance (bin level) to optimum " << avg_dist << endl;
      // << ", distance (chr level) of mean karyotype to optimum " << dm << endl;;

      this->meanfit_by_time[currsize] = avg_fitness;
      // double max_fitness = *max_element(fits.begin(), fits.end());
      this->maxfit_by_time[currsize] = max_fitness;
      // double min_fitness = *min_element(fits.begin(), fits.end());
      this->minfit_by_time[currsize] = min_fitness;
      this->dist_by_time[currsize] = avg_dist;
      // this->dm_by_time[currsize] = dm;
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


    /*
    * Key function when introducing selection into the stochastic branching process
    * Update the birth/death rate of a cell (only occur when there are new mutations in daughter cells) according to that of a reference cell
    * TODO: Different fitness values may apply at different locations of the genome
    * start_cell: the baseline for fitness changes (should be the cell with optimum karyotype)
    * Two models considered here:
    * When genotype_diff > 0, selection model is based on alpha (negative selection -- positive alpha)
    * Otherwise, selection model is based on classic selection coefficient (negative selection -- negative s)
    * Genotype differences are considered to impose different selective advantages to different cells
    */
    void update_cell_growth(Cell_ptr dcell, const Cell_ptr start_cell, const Model& model, int loc_type, int verbose = 0) {
        // Introduce selection to cells with CNAs using specified fitness
        double gdiff = 0;
        // assume start cell has optimum karyotype by default
        double b0 = start_cell->birth_rate;
        double d0 = start_cell->death_rate;
        if(START_WITH_OPT == 0){   // when starting from diploid, rate changes relative to optimum
          b0 = OPT_BIRTHRATE;
          d0 = OPT_DEATHRATE;
        }
        // cout << "optimum birth rate " << b0 << ", optimum death rate " << d0 << endl;

        if(model.use_alpha > 0 || model.genotype_diff > 0){
          if(model.genotype_diff == 3){
            if(verbose > 1) cout << "Using number of mutations to measure genotype difference" << endl;
            gdiff += dcell->num_mut - start_cell->num_mut;
          }else if(model.genotype_diff == 2){
            if(verbose > 1) cout << "Using relative CN to measure genotype difference" << endl;
            for(int i = 0; i < NUM_LOC; i++){
              double rcn = abs(dcell->loc_changes[i] - OPT_KARYOTYPE[i]);
              gdiff += rcn;
            }
          }else if(model.genotype_diff == 1){  // use PGA to measure genotype difference
            if(verbose > 1) cout << "Using percentage of genome altered to measure genotype difference" << endl;
            for(int i = 0; i < NUM_LOC; i++){
              double rcn = abs(dcell->loc_changes[i] - OPT_KARYOTYPE[i]);
              if(rcn > 0) gdiff += 1;   // just count #bins to have a fixed range
            }
          }else{  // used for initialization from diploid, chr-level changes allowable
            assert(model.genotype_diff == 4);
            if(CHR_CNA == 1){
              gdiff += dcell->get_dist_chr_level();

              if(verbose > 1){
                // output cell karyotype
                cout << "current cell karyotype: ";
                for(int i = 0; i < NUM_CHR; i++){
                  cout << dcell->chr_changes[i] << ", ";
                }
                cout << endl;

                cout << "distance between them: " << gdiff << endl;
              }
            }else{
              gdiff += dcell->get_dist_bin_level();
            }
          }
          // normalized by number of bins to avoid sharp change of net growth rate, will be in (0,1) when PGA is used to compute genotype difference
          if(model.norm_by_bin == 1){
            if(verbose > 1) cout << "Nomalize genotype differences by number of bins" << endl;
            gdiff = gdiff / NUM_LOC;
          }
          assert((1 + gdiff * model.fitness) > 0);
        }else{
          assert(1 + model.fitness > 0);
        }


        switch (model.growth_type) {
            case ONLY_BIRTH: {
                if(model.use_alpha > 0){
                    dcell->birth_rate = b0 / (1 + gdiff * model.fitness);
                }else{
                    // To ignore the effect of genotype difference, set gdiff = 1
                    if(model.genotype_diff > 0 && gdiff >= model.min_diff){
                       dcell->birth_rate = b0 * (1 + gdiff * model.fitness);
                    }else{
                       dcell->birth_rate = b0 * (1 + model.fitness);
                    }
                }
                break;
            }
            case CHANGE_BIRTH: {   // decrease birth rate only
                double new_grate = 0;
                if(model.use_alpha > 0){
                    new_grate = (b0 - d0) / (1 + gdiff * model.fitness);
                    dcell->birth_rate = d0 + new_grate;
                    dcell->death_rate = dcell->birth_rate - new_grate;
                }else{
                    if(model.genotype_diff > 0 && gdiff >= model.min_diff){
                        new_grate = (b0 - d0) * (1 + gdiff * model.fitness);
                    }else{
                        new_grate = (b0 - d0) * (1 + model.fitness);
                    }
                    dcell->birth_rate = d0 + new_grate;
                    dcell->death_rate = dcell->birth_rate - new_grate;
                }
                break;
            }
            case CHANGE_DEATH:{    // increase death rate only
                double new_grate = 0;
                if(model.use_alpha > 0){
                    new_grate = (b0 - d0) / (1 + gdiff * model.fitness);
                }else{
                    if(model.genotype_diff > 0 && gdiff >= model.min_diff){
                        new_grate = (b0 - d0) * (1 + gdiff * model.fitness);
                    }else{
                        new_grate = (b0 - d0) * (1 + model.fitness);
                    }
                }
                // cout << "new growth rate " << new_grate << endl;
                double new_drate = (b0 - new_grate);
                dcell->death_rate = new_drate > dcell->death_rate ? new_drate : dcell->death_rate;
                dcell->birth_rate = dcell->death_rate + new_grate;
                break;
            }
            case CHANGE_BOTH:{
                double new_grate = 0;
                if(model.use_alpha > 0){
                    new_grate = (b0 - d0) / (1 + gdiff * model.fitness);
                }else{
                    if(model.genotype_diff > 0 && gdiff >= model.min_diff){
                        new_grate = (b0 - d0) * (1 + gdiff * model.fitness);
                    }else{
                        new_grate = (b0 - d0) * (1 + model.fitness);
                    }
                }
                double rn = runiform(r, 0, 1);
                if(rn < 0.5){
                    dcell->birth_rate = d0 + new_grate;
                    dcell->death_rate = dcell->birth_rate - new_grate;
                }else{
                    double new_drate = (b0 - new_grate);
                    dcell->death_rate = new_drate > dcell->death_rate ? new_drate : dcell->death_rate;
                    dcell->birth_rate = dcell->death_rate + new_grate;
                }
                break;
            }
          default: cout << "" << endl;
        }

        if(verbose > 1){
            cout << "new birth-death rate for cell "<< dcell->cell_ID << " is " << dcell->birth_rate << "-" << dcell->death_rate << " with fitness " << model.fitness << " and genotype difference " << gdiff << endl;
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
         int mut_ID = ncell->num_mut;  // used to distinguish different mutations
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
                 // if(verbose > 1){
                 //   cout << "Select cell " << rID << " with " << rcell->num_mut << " mutations and birth rate " << rcell->birth_rate << ", death rate " << rcell->death_rate << " to divide" << endl;
                 // }

                 Cell_ptr dcell1 = new Cell(++this->ntot, rID, ncell->time_occur + t);
                 dcell1->copy_parent((*rcell));

                 Cell_ptr dcell2 = new Cell(++this->ntot, rID, ncell->time_occur + t);
                 dcell2->copy_parent((*rcell));
                 dcell2->pos = get_neighbor_position(rcell->pos, POS_SIGMA);

                 // daughter cells aquire nu new mutations, where nu ~ Poisson(mutation_rate)
                 if (ncell->mutation_rate > 0) {
                     int nu1 = 0;
                     int nu2 = 0;

                     if(CHR_CNA == 1){
                       nu1 = dcell1->generate_CNV_chr_level(mut_ID, t, verbose);
                       nu2 = dcell2->generate_CNV_chr_level(mut_ID, t, verbose);
                     }else{
                       nu1 = dcell1->generate_CNV_pseudo(mut_ID, t, verbose);
                       nu2 = dcell2->generate_CNV_pseudo(mut_ID, t, verbose);
                     }

                     nu += nu1;
                     nu += nu2;

                     if(verbose > 1){
                         cout << "Number of new mutations in cell " << dcell1->cell_ID << ": " << "\t" << nu1 << endl;
                         cout << "Number of new mutations in cell " << dcell2->cell_ID << ": " << "\t" << nu2 << endl;
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
                 // if(verbose > 1){
                 //   cout << "Select cell " << rID << " with " << rcell->num_mut << " mutations and birth rate " << rcell->birth_rate << ", death rate " << rcell->death_rate << " to disappear" << endl;
                 // }
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



     /*
        This method simulates tumour growth with a rejection-kinetic Monte Carlo algorithm.
        Record all cells in the lineages, only suitable for a small number of cells (e.g. simulating gland as cell)
        Allow deleting cells in case a larger population is simulated
        intput:
         Nend -- the cell population size at the end
         start_cell -- the starting cell. Given a cell in another clone, it can mimick migragation
         model -- 0: neutral, 1: gradual, 2: punctuated
         restart -- 1: start with a new cell; 0: continue with current state
        output:
         a tree-like structure. For each Cell, its children, occurence time, birth rate, death rate
      */
      void grow_with_cnv_cmpl(const Cell_ptr start_cell, const Model& model, int Nend, node* root, int loc_type, double leap_size=0, int track_lineage = 0, int multiple_output = 0, int verbose = 0, int restart = 1, double tend = DBL_MAX){
          // Initialize the simulation with one cell
          if(restart == 1){
              initialize_with_cnv_cmpl(start_cell, model, verbose);
          }

          double t = 0;  // starting time, relative to the time of last end
          int mut_ID = start_cell->num_mut;
          int nu = 0; // The number of new mutations, relative to the time of last end
          int num_mut_event = 0; // count the number of times a CNA event is introduced

          if(verbose > 0) cout << "\nSimulating tumour growth with CNAs under model " << model.model_ID << " at time " << start_cell->time_occur + t << endl;

          // && start_cell->time_occur + t <= tend
          while(this->curr_cells.size() < Nend) {
              if (this->curr_cells.size() == 0) {
                  cout << "restart from first cell" << endl;
                  t = 0;
                  mut_ID = start_cell->num_mut;
                  nu = 0;
                  initialize_with_cnv_cmpl(start_cell, model, verbose);
                  restart = 1;
                  continue;
              }
              // print_all_cells(this->curr_cells, verbose);

              // The population may be shrinked at some time, only sample at the first apperance
              if(multiple_output == 1){
                int currsize = this->curr_cells.size();
                if (find(NTIMES.begin(), NTIMES.end(), currsize) != NTIMES.end() &&    // current population size is to be sampled
                    meanfit_by_time .find(currsize) == meanfit_by_time .end())    // not sampled yet
                {
                  vector<map<int, pcn>> pcn1;
                  map<int, double*> avg_loc_changes;  // used for summary statistics by bin
                  // map<int, double*> avg_chr_changes;  // used for summary statistics by chr
                  unordered_set<int> samples;
                  if(TOSAMPLE < currsize){
                    samples = BobFloydAlgo(TOSAMPLE, currsize);
                  }

                  if(verbose > 0){
                    cout << "sampled cell IDs at " << currsize << ":";
                    for(auto cid : samples){
                      cout << " " << cid;
                    }
                    cout << endl;
                  }
                  // randomly sample 100 cells and store their CNPs
                  int scount = 0;
                  for(auto s : this->curr_cells){
                    if(samples.find(scount) != samples.end() || TOSAMPLE == currsize){
                      // string sname = "G" + to_string(s->cell_ID);
                        map<int, pcn> cell_pcn;
                        pcn cnp;
                        loc2pcn(s->loc_changes, cnp, verbose);
                        cell_pcn[s->cell_ID] = cnp;
                        pcn1.push_back(cell_pcn);

                        avg_loc_changes[s->cell_ID] = new double[NUM_LOC];
                        for(int i = 0; i < NUM_LOC; i++){
                            avg_loc_changes[s->cell_ID][i] = s->loc_changes[i];
                        }
                    }
                    scount++;
                  }

                  cnp_by_time[currsize] = pcn1;
                  lchange_by_time[currsize] = avg_loc_changes;

                  // compute average fitness of the current population
                  // vector<double> fitnesses;
                  this->get_fitness_stats();
                }
              }

              // Choose a random cell from the current population
              int rindex = myrng(this->curr_cells.size());
              Cell_ptr rcell = this->curr_cells[rindex];
              double rbrate = rcell->birth_rate;
              double rdrate = rcell->death_rate;
              int rID = rcell->cell_ID;

              // cout << "root has data " << root->data << endl;
              node* parent = NULL;
              // only record tree structure when starting from a single cell
              if(restart == 1){
                parent = search_node(rID, root);
                // cout << "finish seraching" << endl;
                // cout << parent->data << endl;
                assert(parent->data == rID);
              }

              double rmax = get_rmax();
              // increase time
              double tau = -log(runiform(r, 0, 1));  // an exponentially distributed random variable
              double deltaT = tau/(rmax * this->curr_cells.size());
              t += deltaT;

              // draw a random number
              double rb = runiform(r, 0, rmax);
              // cout << "random number " << rb << endl;
              if(rb < rbrate){
                  if(verbose > 1){
                    cout << "Select cell " << rID << " at index " << rindex << " with " << rcell->num_mut << " mutations and birth rate " << rcell->birth_rate << ", death rate " << rcell->death_rate << " to divide" << endl;
                    // output the number of mutations and birth rates for all the current cells
                    cout << "#mutations for all the current cells:";
                    vector<int> nmuts;
                    for(auto c : this->curr_cells){
                      nmuts.push_back(c->num_mut);
                      // cout << " " << c->num_mut << ";";
                      // << "," << c->birth_rate << "," << c->death_rate << ";";
                    }
                    sort(nmuts.begin(), nmuts.end());
                    for(auto n : nmuts){
                      cout << " " << n << ";";
                    }
                    cout << endl;
                  }

                  int cID1 = this->ntot + 1;
                  Cell_ptr dcell1 = new Cell(cID1, rID, start_cell->time_occur + t);
                  dcell1->copy_parent((*rcell));

                  int cID2 = this->ntot + 2;
                  Cell_ptr dcell2 = new Cell(cID2, rID, start_cell->time_occur + t);
                  dcell2->copy_parent((*rcell));
                  // dcell2->pos = get_neighbor_position(rcell->pos, POS_SIGMA);

                  if(restart == 1){
                    node* node1 = new node(dcell1->cell_ID);
                    parent->left = node1;
                    node* node2 = new node(dcell2->cell_ID);
                    parent->right = node2;
                  }

                  this->ntot = this->ntot + 2;

                  // daughter cells aquire nu new mutations, where nu ~ Poisson(mutation_rate)
                  if (start_cell->mutation_rate > 0) {
                      int nu1 = 0;
                      int nu2 = 0;
                      if(CHR_CNA == 1){
                        nu1 = dcell1->generate_CNV_chr_level(mut_ID, t, verbose);
                        nu2 = dcell2->generate_CNV_chr_level(mut_ID, t, verbose);
                      }else{
                        nu1 = dcell1->generate_CNV_pseudo(mut_ID, t, verbose);
                        nu2 = dcell2->generate_CNV_pseudo(mut_ID, t, verbose);
                      }

                      if(verbose > 1){
                          cout << "Number of mutations in cell " << dcell1->cell_ID << ": " << "\t" << nu1 << endl;
                          cout << "Number of mutations in cell " << dcell2->cell_ID << ": " << "\t" << nu2 << endl;
                      }

                      if(model.fitness != 0)
                      {
                          if(verbose > 1){
                              cout << "Update the growth parameters of daughter cells " << endl;
                          }
                          if(nu1 > 0) update_cell_growth(dcell1, start_cell, model, loc_type, verbose);
                          if(nu2 > 0) update_cell_growth(dcell2, start_cell, model, loc_type, verbose);
                      }
                  }

                  if(track_lineage == 1 && restart == 1){ // all created cells are stored in "cells"
                    if(verbose > 1) cout << "pseudo-remove parent cell " << rID << " at index " << rindex << endl;
                    parent->flag = -1;
                    this->cells.push_back(dcell1);
                    this->cells.push_back(dcell2);
                  }else{
                    if(rID != 1){  // avoid destoying start cell for restarting
                        if(verbose > 1) cout << "delete parent cell " << rID << " at index " << rindex << endl;
                        delete (this->curr_cells[rindex]);
                        this->curr_cells[rindex] = NULL;
                    }
                  }

                  // Remove the parent cell from the list of current cells, not really delete it when tracing lineage
                  this->curr_cells.erase(this->curr_cells.begin() + rindex);

                  this->curr_cells.push_back(dcell1);
                  this->curr_cells.push_back(dcell2);
              }
              // death event if b<r<b+d
              else if(rb >= rbrate && rb < rbrate + rdrate) {
                  if(verbose > 1){
                    cout << "Select cell " << rID << " at index " << rindex << " with " << rcell->num_mut << " mutations and birth rate " << rcell->birth_rate << ", death rate " << rcell->death_rate << " to disappear" << endl;
                  }
                  if(track_lineage == 1 && restart == 1){
                    if(verbose > 1) cout << "pseudo-remove parent cell " << rID << " at index " << rindex << endl;
                    parent->flag = -1;
                  }else{
                    if(rID != 1){
                      if(verbose > 1) cout << "delete parent cell " << rID << " at index " << rindex << endl;
                        delete (this->curr_cells[rindex]);
                        this->curr_cells[rindex] = NULL;
                    }
                  }
                  // Remove the parent cell from the list of current cells, not really delete it for tracing lineage
                  this->curr_cells.erase(this->curr_cells.begin() + rindex);
              }else{

              }
          }

          this->time_end += t;
          this->num_novel_mutation += nu;

          if(verbose > 0){
              cout << "Generated " << nu << " mutations during time " << t << endl;
              if(verbose > 1 && track_lineage == 1) print_all_cells(this->cells, verbose);
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
     void set_bulk_cnp_from_pseudo(double avg_loc_change[NUM_LOC], const vector<Cell_ptr>& cells) {
         for(int i = 0; i < NUM_LOC; i++){
             for(auto cell : cells) {
                 avg_loc_change[i] += cell->loc_changes[i];
             }
             double avg_cn = avg_loc_change[i] / cells.size();
             double rcn =  rint(avg_cn * DEC_PLACE) / DEC_PLACE;  // 1 decimal place
             avg_loc_change[i] = rcn;
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
