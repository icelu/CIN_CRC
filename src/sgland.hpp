#ifndef GLAND_HPP
#define GLAND_HPP


#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <queue>
#include <random>
#include <string>

#include <assert.h>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include "util.hpp"


using namespace std;

// record GMutation ID to track clone expansion
class GMutation
{
public:
    int mut_ID;
    int cell_ID;
    int gen_ID;    
    int chr;  // chromosome on which the GMutation occurs
    int arm;    // 0: whole chr, 1: p; 2: q

    ~GMutation() = default;
    GMutation(const GMutation& other) = default;
    GMutation(GMutation&& other) = default;
    GMutation& operator=(const GMutation& other) = default;
    GMutation& operator=(GMutation&& other) = default;

    GMutation(){

    }

    GMutation(int mut_ID, int cell_ID, int gen_ID, int chr, int arm){
        this->mut_ID = mut_ID;
        this->cell_ID = cell_ID;
        this->gen_ID = gen_ID;

        this->chr = chr;
        this->arm = arm;
    }

    bool operator<(const GMutation& mut) const{
        if(mut.mut_ID < this->mut_ID){
            return true;
        }else{
            return false;
        }
            
    }

    void print(ofstream& fout) const{
        fout << mut_ID << "\t" << cell_ID << "\t" << gen_ID << "\t" << chr << "\t" << arm;
    }

};


// cells in a cancer gland
class GCell{
public:
    int cell_ID;
    int parent_ID;
    int gen_ID;  // generation ID to track lineages
    int num_divisions; // the number of generations that this cell exists since 1st division
    int status;  // 0: stem, 1: mitotic, 2: differentiated

    // copy number for chr, arm; relative so that it is safe to ignore clonal events
    map<pair<int, int>, int> cn_profile;        // only applied to each arm
    map<pair<int, int>, int> obs_cn_profile;    // grouped by chr-/arm- level
    vector<GMutation> muts; 

    ~GCell() = default;
    GCell(const GCell& other) = default;
    GCell(GCell&& other) = default;
    GCell& operator=(const GCell& other) = default;
    GCell& operator=(GCell&& other) = default;   


    GCell(int cell_ID, int parent_ID, int gen_ID, int num_divisions, int status){
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;
        this->gen_ID = gen_ID;
        this->num_divisions = num_divisions;
        this->status = status;
    }


    void copy_parent(GCell* parent){
        this->cn_profile = parent->cn_profile;
        this->muts = parent->muts;    
    }


    // generate a random number of CNAs
    void generate_CNA(int& mut_ID, double arm_rate, double chr_rate, int verbose = 0){
        assert(arm_rate >= 0 && chr_rate >= 0);
        // A chr is not applicable to GMutation if it has 0 copy. Since only relative changes are considered, assume there are sufficient copies
        gsl_ran_discrete_t* dis_chr = gsl_ran_discrete_preproc(NUM_CHR, CHR_PROBS);
        gsl_ran_discrete_t* dis_arm = gsl_ran_discrete_preproc(2, ARM_PROBS);
        int chr, arm;
        double u = 0;

        int nu1 = gsl_ran_poisson(r, chr_rate);
        if(verbose > 1 && nu1 > 0) cout << "Generating " << nu1 << " chr-level mutations under rate " << chr_rate << " in cell " << cell_ID << endl;
        if(nu1 > 0){
            for(int j = 0; j < nu1; j++){
                int chr = gsl_ran_discrete(r, dis_chr);
                pair<int, int> pos1(chr, 1);
                pair<int, int> pos2(chr, 2);   
                if(verbose > 1) cout << "\tmutation on cell " << cell_ID << " chr " << chr+1 << endl;  

                u = runiform(r, 0, 1);
                if(u < 0.5){  // gain
                    cn_profile[pos1]++;       
                    cn_profile[pos2]++;              
                }else{
                    cn_profile[pos1]--;
                    cn_profile[pos2]--;
                }
                GMutation mut(++mut_ID, cell_ID, gen_ID, chr, 0);
                this->muts.push_back(mut);
            }
        }

        int nu2 = gsl_ran_poisson(r, arm_rate);
        if(verbose > 1 && nu2 > 0) cout << "Generating " << nu2 << " arm-level mutations under rate " << arm_rate << " in cell " << cell_ID << endl;
        if(nu2 > 0){
            for(int j = 0; j < nu2; j++){
                int chr = gsl_ran_discrete(r, dis_chr);
                int arm = gsl_ran_discrete(r, dis_arm) + 1;
                pair<int, int> pos(chr, arm); 
                if(verbose > 1) cout << "\tmutation on cell " << cell_ID << " chr " << chr+1 << " arm " << arm+1 << endl;   

                u = runiform(r, 0, 1);
                if(u < 0.5){  // gain
                    cn_profile[pos]++;
                }else{
                    cn_profile[pos]--;
                }  
                GMutation mut(++mut_ID, cell_ID, gen_ID, chr, arm);
                this->muts.push_back(mut);                           
            }
        }

    }


    // only output CNAs
    void set_obs_cn(){
        obs_cn_profile.clear();
        for(int c = 0; c < NUM_CHR; c++){
            bool has_cna = false;
            if(cn_profile[pair<int, int>(c,1)] == cn_profile[pair<int, int>(c,2)] && cn_profile[pair<int, int>(c,1)] != 0){
                obs_cn_profile[pair<int, int>(c,0)] = cn_profile[pair<int, int>(c,1)];
                has_cna = true;
                continue;
            }
            if(cn_profile[pair<int, int>(c,1)] != 0){
                obs_cn_profile[pair<int, int>(c,1)] = cn_profile[pair<int, int>(c,1)];
                has_cna = true;
            }
            if(cn_profile[pair<int, int>(c,2)] != 0){
                obs_cn_profile[pair<int, int>(c,2)] = cn_profile[pair<int, int>(c,2)];
                has_cna = true;
            }
            if(!has_cna){  // no cn changes
                obs_cn_profile[pair<int, int>(c,0)] = 0;
            }
        }
    }


    // output the copy number profiles of a cell
    void print_cn_profile(ofstream& fout){
        // fout << "\tCNAs in cell " << cell_ID << " with flag " << flag << endl;
        for(auto cn : cn_profile){
            if(cn.second != 0){
                fout << cell_ID << "\t" << cn.first.first + 1 << "\t" << cn.first.second << "\t" << cn.second << endl;
            }
        }
    }


    // output the mutations in a cell
    void print_muts(ofstream& fout){
        // fout << "\tCNAs in cell " << cell_ID << " with flag " << flag << endl;
        for(auto mut : muts){           
            fout << mut.mut_ID << "\t" << mut.cell_ID << "\t" << mut.gen_ID << "\t" << mut.chr << "\t" << mut.arm << endl;
        }
    }


    void print_cell_info(){
        cout << "cell " << this->cell_ID << ", parent " << this->parent_ID << ", generation " << this->gen_ID  << ", number of divisions " << this->num_divisions << ", status " << this->status << endl;
    }


};


// a gland with N stem cells and n-N non-stem cells, constant size model
class SGland{
public:
    int num_stem;  // N stem cells starting with core karyotype
    int gsize;   // total number of cells in the gland 
    int num_gen; 
    int num_gen_nonstem;
    int num_gen_steady;
    double prob_d1;      // p (prob_d1) = 1: immortal stem cells; 0 <= p < 1: stem cell niche

    // assume all cells in the gland have the same GMutation rate
    double arm_rate;
    double chr_rate;

    vector<GCell*>  stem_cells;  // only available cells at present, use object for easy implementation
    queue<GCell*>  non_stem_cells;  // only available cells at present, use object for easy implementation

    set<tuple<int, int, int>> uniq_cns;
    map<pair<int, int>, int> locus_diversity;

    map<pair<int, int>, double> prob_children;

    ~SGland() = default;
    SGland(const SGland& other) = default;
    SGland(SGland&& other) = default;
    SGland& operator=(const SGland& other) = default;
    SGland& operator=(SGland&& other) = default;   

    SGland(){

    }

    SGland(int num_stem, int gsize, int num_gen, int num_gen_nonstem, int num_gen_steady, double prob_d1, double arm_rate, double chr_rate){
        this->num_stem = num_stem;
        this->gsize = gsize;
        this->num_gen = num_gen;
        this->num_gen_nonstem = num_gen_nonstem;
        this->num_gen_steady = num_gen_steady;
        this->prob_d1 = prob_d1;
        this->arm_rate = arm_rate;
        this->chr_rate = chr_rate;
    }


    void initialize(){
        // randomly assign num_stem cells to be stem
        int i = 0;
        for(; i < num_stem; i++){
            GCell* c = new GCell(i, -1, 1, 0, 0);
            stem_cells.push_back(c);
        }
        for(int j = 0; j < gsize - num_stem; j++){
            GCell* c = new GCell(j + i, -1, 1, 0, 1);
            non_stem_cells.push(c);
        }
    }


    void initialize_stem_cells(){
        // randomly assign num_stem cells to be stem
        int i = 0;
        for(; i < num_stem; i++){
            GCell* c = new GCell(i, -1, 1, 0, 0);
            stem_cells.push_back(c);
        }
    }


    // compute convolution probabilities P(l, i), i = 0, 1, ..., N
    void get_rate_children(int num_stem, double prob_d1, int verbose = 0){
        double q = (1 - prob_d1) / 2;
        prob_children[pair<int,int>(1, 0)] = q;
        prob_children[pair<int,int>(1, 2)] = q;
        prob_children[pair<int,int>(1, 1)] = prob_d1;
        for(int i = 3; i <= num_stem; i++){
            prob_children[pair<int,int>(1, i)] = 0;
        }

        for(int j = 1; j < num_stem; j++){
            for(int i = 0; i <= num_stem; i++){
                prob_children[pair<int,int>(j + 1, i)] = q * prob_children[pair<int,int>(j, i)];
                if(i > 1) prob_children[pair<int,int>(j + 1, i)] += prob_d1 * prob_children[pair<int,int>(j, i - 1)];
                if(i > 2) prob_children[pair<int,int>(j + 1, i)] += q * prob_children[pair<int,int>(j, i - 2)];
            }           
        }

        if(verbose > 1){
            for(int j = 1; j <= num_stem; j++){
                cout << "P(" << j << ",i):";
                double sum = 0;
                for(int i = 0; i <= num_stem; i++){
                    cout << "\t" << prob_children[pair<int,int>(j, i)]; 
                    sum += prob_children[pair<int,int>(j, i)];
                }
                cout << "\t" << sum;
                cout << endl;
            }
        }
    }


    void update_stem_cell(GCell* rcell, int& mut_ID, int rindex, int& n_new_stem, int& n_new_nonstem, int& index_old_stem, int& cell_count, int gen_count, int gen_nonstem, int verbose = 0){
        // Each stem cell produces two cells, of which either 0, 1, or 2 are stem cells with probabilites q, p, and q, where p + 2q = 1
        // sample an observation of (v1, ..., vN), the number of daughter stem cells after each division, from the joint distribution
        double q = prob_children[pair<int,int>(1, 0)];
        int n_child_stem = 0;    
        
        if(index_old_stem == 1){                  
            double probs[3] = {q, prob_d1, q};
            gsl_ran_discrete_t* dis = gsl_ran_discrete_preproc(3, probs);
            n_child_stem = gsl_ran_discrete(r, dis);
            n_new_stem += n_child_stem; 
            if(verbose > 1) cout << "The distribution of P(" << index_old_stem << ") " << probs[0] << "\t" << probs[1] << "\t" << probs[2] << endl;                      
        }else{
            // compute conditional probability dependent on index of stem cell 
            int j = index_old_stem - 1;    // assume current one is v_(j+1)       
            assert(j >= 1);
            int N = num_stem;
            if(j < N - 1){
                double p0 = q * prob_children[pair<int,int>(N - j - 1, N - n_new_stem)] / prob_children[pair<int,int>(N - j,  N - n_new_stem)];
                double p1 = prob_d1 * prob_children[pair<int,int>(N - j - 1, N - n_new_stem - 1)] / prob_children[pair<int,int>(N - j,  N - n_new_stem)];
                double p2 = q * prob_children[pair<int,int>(N - j - 1, N - n_new_stem - 2)] / prob_children[pair<int,int>(N - j,  N - n_new_stem)];
                double probs[3] = {p0, p1, p2};
                if(verbose > 1) cout << "The conditional distribution of P(" << j + 1 << ") " << probs[0] << "\t" << probs[1] << "\t" << probs[2] << endl;
                gsl_ran_discrete_t* dis = gsl_ran_discrete_preproc(3, probs);
                n_child_stem = gsl_ran_discrete(r, dis);
                n_new_stem += n_child_stem;
            }else{
                assert(j == N - 1);
                n_child_stem = N - n_new_stem;
                n_new_stem += n_child_stem;
            }           
        }
        index_old_stem++;

        if(verbose > 1) cout << "v" << index_old_stem << " = " << n_child_stem << ", new stem cells so far " << n_new_stem << endl;

        // daughter cells of a stem cell always have division count 0
        if(n_child_stem == 0){
            // check the life time of the non-stem cells, no need to store if they are very early 
            if(gen_count > gen_nonstem){               
                GCell* dcell1 = new GCell(++cell_count, rcell->cell_ID, gen_count, 0, 1);  
                dcell1->copy_parent(rcell);             
                dcell1->generate_CNA(mut_ID, arm_rate, chr_rate, verbose); 

                GCell* dcell2 = new GCell(++cell_count, rcell->cell_ID, gen_count, 0, 1);   
                dcell2->copy_parent(rcell);              
                dcell2->generate_CNA(mut_ID, arm_rate, chr_rate, verbose);  

                this->non_stem_cells.push(dcell1);
                this->non_stem_cells.push(dcell2); 
                n_new_nonstem += 2;
            }
        }else if(n_child_stem == 1){
            GCell* dcell1 = new GCell(++cell_count, rcell->cell_ID, gen_count, 0, 0);  
            dcell1->copy_parent(rcell);          
            dcell1->generate_CNA(mut_ID, arm_rate, chr_rate, verbose);               
            this->stem_cells.push_back(dcell1);
            if(gen_count > gen_nonstem){
                GCell* dcell2 = new GCell(++cell_count, rcell->cell_ID, gen_count, 0, 1);  
                dcell2->copy_parent(rcell);             
                dcell2->generate_CNA(mut_ID, arm_rate, chr_rate, verbose);   
                dcell2->num_divisions = 0;
                this->non_stem_cells.push(dcell2);
                n_new_nonstem += 1;
            }
        }else{
            assert(n_child_stem == 2);   

            GCell* dcell1 = new GCell(++cell_count, rcell->cell_ID, gen_count, 0, 0);
            dcell1->copy_parent(rcell);                          
            dcell1->generate_CNA(mut_ID, arm_rate, chr_rate, verbose);   

            GCell* dcell2 = new GCell(++cell_count, rcell->cell_ID, gen_count, 0, 0);
            dcell2->copy_parent(rcell);               
            dcell2->generate_CNA(mut_ID, arm_rate, chr_rate, verbose);  

            this->stem_cells.push_back(dcell1);
            this->stem_cells.push_back(dcell2); 
        }
    }


    // assume N stem cells in a gland base after each division
    // stop when reaching n cells (not all cells are stem-like), use 5000
    // N initial cancer stem cells in the gland with core karyotype, new CNAs generated during dynamic process of the cancer gland  
    void grow_with_cna(int verbose = 0){       
        int gen_count = 0;
        int mut_ID = 0;
        int num_nonstem = gsize - num_stem;   
        int gen_nonstem = num_gen - num_gen_nonstem - num_gen_steady;   // the generation ID when the non-stem cells can last to present
        if(verbose > 0){
            cout << "Start simulating non-stem cells at generation " << gen_nonstem << endl;
            cout << "Maximum number of non-stem cells " << num_nonstem << endl;
        }

        while(gen_count < num_gen){            
            assert(this->stem_cells.size() == num_stem); 

            int cell_count = 0;  // count for each generation, to track the total number of cells in history
            int n_new_nonstem = 0;                       
            int n_new_stem = 0;            
            int index_old_stem = 1; 
                               
            gen_count++;
            if(verbose > 0) cout << "\nGeneration " << gen_count << endl;

            // update non-stem cells first to avoid dividing new non-stem cells
            int num_nonstem = gsize - num_stem;  // the maximum number of non-stem cells
            int num_nonstem_curr = this->non_stem_cells.size();
            if(verbose > 0){
                cout << "Total number of non-stem cells so far " << num_nonstem_curr << endl;
            }

            if(gen_count > gen_nonstem){  // start from previous stem-cell daughters
                if(verbose > 0){
                    cout << "Simulating non-stem cells at generation " << gen_count << endl;
                }
                int n_dc = 0;
                // no need to grow new non-stem cells
                for(int i = 0; i < num_nonstem_curr; i++){
                    GCell* rcell = this->non_stem_cells.front();
                    // cout << "Dividing non-stem cell "  << i << endl;
                    // rcell->print_cell_info();
                    assert(rcell->status != 0);
                    if(rcell->num_divisions >= num_gen_nonstem){
                        rcell->status = 2;
                        n_dc += 1;
                        // cout << "differentiated cell" << endl;
                        // rcell->print_cell_info();
                        this->non_stem_cells.pop(); 
                        // move this cell from head to tail
                        this->non_stem_cells.push(rcell);                       
                        continue;
                    }
                    // a non-stem cell always generates two cells
                    GCell* dcell1 = new GCell(++cell_count, rcell->cell_ID, gen_count, rcell->num_divisions+1, 1);
                    dcell1->copy_parent(rcell);               
                    dcell1->generate_CNA(mut_ID, arm_rate, chr_rate, verbose);     
                    this->non_stem_cells.push(dcell1);          

                    GCell* dcell2 = new GCell(++cell_count, rcell->cell_ID, gen_count, rcell->num_divisions+1, 1); 
                    dcell2->copy_parent(rcell);              
                    dcell2->generate_CNA(mut_ID, arm_rate, chr_rate, verbose);  
                    this->non_stem_cells.push(dcell2);
                    
                    n_new_nonstem += 2; 

                    delete rcell;
                    this->non_stem_cells.pop();
                }
                if(verbose > 0){
                    cout << "Number of new non-stem cells generated during non-stem cell divisions " << n_new_nonstem << endl;                   
                    cout << "Number of differentiated cells " << n_dc << endl;
                }
            }                              
       
            // Choose a random Cell from the current population
            if(verbose > 0) cout << "Updating stem cells" << endl;            
            for(int i = 0; i < num_stem; i++){
                GCell* rcell = this->stem_cells[i]; 
                assert(rcell->status == 0);
                update_stem_cell(rcell, mut_ID, i, n_new_stem, n_new_nonstem, index_old_stem, cell_count, gen_count, gen_nonstem, verbose);
            }           
            // remove old stem cells at the end to avoid breaking index
            for(int i = 0; i < num_stem; i++){ 
                delete (this->stem_cells[i]);
                this->stem_cells[i] = NULL;
                this->stem_cells.erase(this->stem_cells.begin() + i);
            }
            assert(stem_cells.size() == num_stem); 
            if(verbose > 0) cout << "Number of new non-stem cells generated during stem cell divisions " << n_new_nonstem << endl;

        }   
        // only check for the number of non-stem cells in the end
        assert(this->non_stem_cells.size() == num_nonstem);     
    } 



    // Take a random sample and compute summary statistics
    void sample_cell_from_gland(vector<GCell*>& sampled_cells, int nsample, int verbose = 0){
        // merge all cells into one vector
        vector<GCell*> cells_all;
        for(auto c : stem_cells){
            cells_all.push_back(c);
        }
        queue<GCell*> tmp_q = non_stem_cells;
        while(!tmp_q.empty()){
           GCell* c = tmp_q.front();
           cells_all.push_back(c);
           tmp_q.pop();
        }
        assert(cells_all.size() == gsize);

        vector<int> index(gsize);
        iota(index.begin(), index.end(), 0);
        vector<int>::iterator sindex = random_unique(index.begin(), index.end(), num_stem);

        for(int i = 0; i < nsample; i++){
            sampled_cells.push_back(cells_all[i]);
            if(verbose > 1){
                cout << "sample cell at " << i << endl;
                cells_all[i]->print_cell_info();
            }
        }
        assert(sampled_cells.size() == nsample);
    }


    // get unique CNVs
    void set_uniq_cn(vector<GCell*> cells){
        for(int c = 0; c < cells.size(); c++){
            cells[c]->set_obs_cn();
            map<pair<int, int>, int> obs_cns = cells[c]->obs_cn_profile;
            for (auto cp : obs_cns){
                tuple<int, int, int> cnvec(cp.first.first, cp.first.second, cp.second);
                uniq_cns.insert(cnvec);
            }
        }

        // cout << "There are " << uniq_cns.size() << " unique CNA events in the clone " << endl;
        // for(auto cv : uniq_cns){
        //    cout << get<0>(cv) << "\t"  << get<1>(cv) << "\t"  << get<2>(cv) << endl;
        // }
    }


    void get_avg_reciprocal_cn(const vector<GCell*> cells, vector<int>& avg_cn){
         map<pair<int, int>, int> cp;

         set_uniq_cn(cells);

         for(int i = 0; i < NUM_CHR; i++){
             for(int j = 0; j < 3; j++){
                 pair<int, int> pos(i,j);
                 cp[pos] = 0;
             }
         }

         for(int j = 0; j < 2; j++){
             avg_cn[j] = 0;
         }

         for(auto cv: uniq_cns){
             pair<int, int> pos(get<0>(cv), get<1>(cv));
             cp[pos] += abs(get<2>(cv));
         }

         for(auto cv: cp){
             // cout << cv.first.first + 1 << "\t"  << cv.first.second << "\t" << cv.second << endl;
             if(cv.first.second==0){
                 avg_cn[0] += cv.second;
             }else{
                 avg_cn[1] += cv.second;
             }
         }   
         // return avg_cn;
     }


     // percentage of cells with CNA
     void get_locus_diversity(){
        for(auto c : stem_cells){
             for(int i = 0; i < NUM_CHR; i++){
                for(int j = 1; j < 3; j++){
                    pair<int, int> pos(i,j);
                    if(c->cn_profile[pos] != 0){
                        locus_diversity[pos]++;
                    }
                }
             }
        }

        queue<GCell*> non_stem_cells2(non_stem_cells);
        while(!non_stem_cells2.empty()){
            GCell* c = non_stem_cells2.front();
            for(int i = 0; i < NUM_CHR; i++){
                for(int j = 1; j < 3; j++){
                    pair<int, int> pos(i,j);
                    if(c->cn_profile[pos] != 0){
                        locus_diversity[pos]++;
                    }
                }
            }
            non_stem_cells2.pop();             
        }
     }


     void write_cna(ofstream& fout){
        for(auto c : stem_cells){
            c->print_cn_profile(fout); 
        }

        queue<GCell*> non_stem_cells2(non_stem_cells);
        while(!non_stem_cells2.empty()){
            GCell* c = non_stem_cells2.front();
            c->print_cn_profile(fout);
            non_stem_cells2.pop();             
        }
     }


     // write mutation frequency   
     void write_mut(ofstream& fout){
        map<GMutation, int> mut_count;
        for(auto c : stem_cells){
            // string mu = to_string(c->mut.mut_ID) << "\t" << to_string(mut.cell_ID) << "\t" << to_string(mut.gen_ID) << "\t" << to_string(mut.chr << "\t" << to_string(mut.arm
            // c->print_muts(fout);
            for(auto mut : c->muts){  
                mut_count[mut]++;
            }
        }

        queue<GCell*> non_stem_cells2(non_stem_cells);
        while(!non_stem_cells2.empty()){
            GCell* c = non_stem_cells2.front();
            for(auto mut : c->muts){  
                mut_count[mut]++;
            }
            non_stem_cells2.pop();             
        }

        for(auto mc : mut_count){
            mc.first.print(fout);
            fout << "\t" << to_string(mc.second) << endl;
        }
     }     
};


#endif