#ifndef Gland_HPP
#define Gland_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <string>

#include <assert.h>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>


#include "util.hpp"
#include "cell.hpp"
#include "clone.hpp"


using namespace std;


typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance>> VAR;
// typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::count>> NBIN;

bool compare_ptr2cell(Cell_ptr a, Cell_ptr b) { return (a->cell_ID < b->cell_ID); }


// Glands may have several clones sampled
class Glands{
public:
    int glands_ID;

    vector<Clone*> clones;      // For all glands
    vector<Clone*> samples;     // For selected glands
    vector<Cell_ptr> sampled_cells;     // For selected glands
    vector<int> sample_IDs;
    vector<int> sample_IDs_side1;
    vector<int> sample_IDs_side2;

    node* root;     // root node of gland lineage tree (storing integer IDs of each node)

    set<string> bps_all;
    set<string> muts;

    Glands(const Glands& other) = default;
    Glands& operator=(const Glands& other) = default;
    Glands& operator=(Glands&& other) = default;

    Glands(){
        glands_ID = 0;
        root = NULL;
    }

    ~Glands(){
        // cout << "Release " << clones.size() << " clones " << endl;
        for(auto cl: clones){
            // cout << "clone at " << cl << endl;
            delete cl;
        }
        if(root != NULL)
            destroy_tree(root);
    }


/*********************** Functions related to gland growth *************************************/
// Sample demes (glands) from all the current demes. Output the side, piece, gland ID, for each deme.
void sample_demes(int nregion){
    // for(int i = 0; i < nglands.size(); i++){
    //     int num_s1 = nglands[i];
    //     // vector<Clone*> glands1;     // sampled glands at one side
    //     for(int j = 0; j < num_s1; j++){
    //         int gID = myrng(100);
    //         samples.push_back(this->clones[gID]);
    //         // Print out CNPs of sampled gland
    //     }
    // }
    // if(verbose > 0)   print_sample_cnp(outdir, suffix, verbose);
}

// Simulate CNP for each gland.
// nglands: number of glands for each sides
// ndeme: total number of glands
void simulate_gland_growth(const Cell_ptr start_cell, int ndeme, int max_deme_size, const Model& start_model, vector<string>& lineages, int store_lineage = 0, int loc_type=BIN, double leap_size=0, int verbose = 0){
    int num_clone = 1;
    Clone* s = new Clone(num_clone, 0);
    // s->grow(start_cell, max_deme_size, verbose, 1);
    s->grow_with_cnv(start_cell, start_model, max_deme_size, loc_type, leap_size, verbose);
    this->clones.push_back(s);

    while (this->clones.size() < ndeme) {
        // randomly select a deme to grow and split
        // cout << "Deme number " << this->clones.size() << endl;
        int rindex = myrng(this->clones.size());
        Clone* s0 = this->clones[rindex];

        // start growing from current population
        while (s0->curr_cells.size() < max_deme_size) {
            // s0->grow(start_cell, max_deme_size, verbose, 0);
            s0->grow_with_cnv(start_cell, start_model, max_deme_size, loc_type, leap_size, verbose, 0);
        }

        // split deme
        Clone* s1 = new Clone(num_clone+1, s0->clone_ID);
        Clone* s2 = new Clone(num_clone+2, s0->clone_ID);
        num_clone += 2;
        for(auto c : s0->curr_cells){
            double rc = runiform(r, 0, 1);
            if(rc < 0.5){
                c->clone_ID = s1->clone_ID;
                s1->curr_cells.push_back(c);
                // s1->muts.insert(c->muts.begin(), c->muts.end());
                for(auto mut : c->muts){
                  s1->muts_freq[mut]++;
                }
            }else{
                c->clone_ID = s2->clone_ID;
                s2->curr_cells.push_back(c);
                // s2->muts.insert(c->muts.begin(), c->muts.end());
                for(auto mut : c->muts){
                  s2->muts_freq[mut]++;
                }
            }
            // this->muts.insert(c->muts.begin(), c->muts.end());
        }
        s0->curr_cells.clear();
        // keep growing the current demes
        this->clones.push_back(s1);
        this->clones.push_back(s2);

        if(store_lineage != 0 && s0->clone_ID != 0){
            string lstr = to_string(s0->clone_ID) + "\t" + to_string(s1->clone_ID) + "\t" + to_string(s1->curr_cells.size())  + "\t" + to_string(s1->muts_freq.size())  + "\t" + to_string(s0->num_novel_mutation);
            lineages.push_back(lstr);
            lstr = to_string(s0->clone_ID) + "\t" + to_string(s2->clone_ID) + "\t" + to_string(s2->curr_cells.size()) + "\t" + to_string(s2->muts_freq.size()) + "\t" + to_string(s0->num_novel_mutation);
            lineages.push_back(lstr);
        }

        delete (this->clones[rindex]);
        this->clones[rindex] = NULL;
        this->clones.erase(this->clones.begin() + rindex);
    }
}


void print_gland_lineage(string fdeme, string header, const vector<string>& lineages) {
    ofstream fout(fdeme);
    fout << header << endl;

    for(auto s : lineages){
        fout << s << endl;
    }

    fout.close();
}



// Simulate CNP for each gland (cell).
// nglands: number of glands for each sides
// ndeme: total number of glands
void simulate_gland_as_cell(const Cell_ptr start_cell, int ndeme, const Model& start_model, vector<string>& lineages, int store_lineage = 0, int loc_type=BIN, double leap_size=0, int verbose = 0){
    int num_clone = 1;
    Clone* s = new Clone(num_clone, 0);
    // s->grow(start_cell, max_deme_size, verbose, 1);
    assert(root != NULL);
    s->grow_with_cnv_cmpl(start_cell, start_model, ndeme, root, loc_type, leap_size, verbose);
    this->clones.push_back(s);

    if(store_lineage != 0){
        for(auto d : s->cells){
            string lstr = to_string(d->parent_ID) +  "\t" + to_string(d->cell_ID) + "\t" + to_string(d->time_occur) + "\t" + to_string(d->num_mut);
            if(d->parent_ID != 0) lineages.push_back(lstr);
        }
    }
}


// Sample from a larger population of cells to intimate spatial sampling and increase diversity among cells
// Sometimes, the tree starting from the root is very skewed such that there is not enough glands at one side.
// In that case, ignore the lineage with less than n tips, where n is the minimal number of glands to sample
void sample_gland_from_cell(const vector<int>& nglands, int verbose = 0){
    sample_IDs.clear();
    sample_IDs_side1.clear();
    sample_IDs_side2.clear();

    int nside = nglands.size();
    // if(nside == 0) return;
    assert(nside = 2);

    int ng1 = nglands[0];
    int ng2 = nglands[1];

    // cout << ng1 << " glands taken at side 1" << endl;
    // cout << ng2 << " glands taken at side 2" << endl;

    int s1 = 0; // track number of sampled glands at each side

    // find all leaves at one side
    vector<int> nodes1;
    get_leaves_below_avail(root->left, nodes1);
    vector<int> nodes2;
    get_leaves_below_avail(root->right, nodes2);

    int nsize1 = nodes1.size();
    int nsize2 = nodes2.size();

    int min_nsize = nsize1 < nsize2? nsize1 : nsize2;
    int min_ng = ng1 < ng2? ng1: ng2;

    node* r = root;
    while(min_nsize < min_ng){
        if(nsize1 < min_ng){
          r = r->right;
        }else{
          r = r->left;
        }
        if(verbose > 0) cout << "sampling from node " << r->data << endl;
        nodes1.clear();
        get_leaves_below_avail(r->left, nodes1);
        nodes2.clear();
        get_leaves_below_avail(r->right, nodes2);

        nsize1 = nodes1.size();
        nsize2 = nodes2.size();

        min_nsize = nsize1 < nsize2? nsize1 : nsize2;
    }
    assert(min_nsize >= min_ng);

    random_shuffle(nodes1.begin(), nodes1.end());
    random_shuffle(nodes2.begin(), nodes2.end());

    if(nsize1 < nsize2){
        for(int i = 0; i < ng1; i++){
            sample_IDs_side1.push_back(nodes1[i]);
            sample_IDs.push_back(nodes1[i]);
        }
        for(int i = 0; i < ng2; i++){
            sample_IDs_side2.push_back(nodes2[i]);
            sample_IDs.push_back(nodes2[i]);
        }
    }else{
        for(int i = 0; i < ng1; i++){
            sample_IDs_side1.push_back(nodes2[i]);
            sample_IDs.push_back(nodes2[i]);
        }
        for(int i = 0; i < ng2; i++){
            sample_IDs_side2.push_back(nodes1[i]);
            sample_IDs.push_back(nodes1[i]);
        }
    }

    sort(sample_IDs.begin(), sample_IDs.end());

    if(verbose > 0){
        cout << sample_IDs_side1.size() << " glands at side 1" << endl;
        for(auto id : sample_IDs_side1){
          cout << id << endl;
        }
        cout << sample_IDs_side2.size() << " glands at side 2" << endl;
        for(auto id : sample_IDs_side2){
          cout << id << endl;
        }
        cout << sample_IDs.size() << " glands total taken" << endl;
    }

}


void  print_sstat(const vector<int>& ids, map<int, double*>& avg_loc_changes, int stat_type, int loc_type, double min_freq=0, double max_freq=0.25, double delta = 0.001, int use_std = 0, int verbose = 0) {
    if(verbose > 0){
        cout << "Printing summary statistics" << endl;
    }
    // LVAR, ADIFF, AVG, CMPL, DIFF, BP, ALTBIN, ALTBIN_SEP, BP_BIN, ALL
    switch (stat_type) {
        case LVAR:{ // 0
            // print_variance(avg_loc_changes, loc_type);
            print_mut_freq(min_freq, max_freq, delta);
            break;
        }
        case ADIFF: { // 1
            // cout << "Using pairwise difference across locations" << endl;
            print_bin_subclonal_stat_by_type(avg_loc_changes, 1, verbose);
            print_pairwise_divergence(ids, avg_loc_changes, 1, verbose);
            break;
        }
        case AVG: { // 2
            // cout << "Using pairwise difference across locations" << endl;
            print_bin_subclonal_stat_by_type(avg_loc_changes, 0, verbose);
            print_pairwise_divergence(ids, avg_loc_changes, 1,  verbose);
            break;
        }
        case DIFF: {  // 4
            // cout << "Using pairwise difference" << endl;
            print_bin_subclonal_stat_by_type(avg_loc_changes, 0, verbose);
            print_pairwise_divergence(ids, avg_loc_changes, 0, verbose);
            collect_private_subclonal_bps(avg_loc_changes, verbose);
            print_num_uniq_mut(verbose);
            print_pairwise_mismatch(verbose);
            break;
        }
        case BP: {  // 5
            // collect_private_subclonal_bps(avg_loc_changes, verbose);
            print_nmut_per_gland();
            break;
        }
        case ALTBIN: {  // 6
            print_bin_subclonal_stat_by_type(avg_loc_changes, 0, verbose);
            print_variance(avg_loc_changes, loc_type, 0, use_std);
            break;
        }
        case ALTBIN_SEP: { // 7
            print_bin_subclonal_stat_by_type(avg_loc_changes, 1, verbose);
            break;
        }
        case BP_BIN: { //8
            collect_private_subclonal_bps(avg_loc_changes, verbose);
            // print_bin_subclonal_stat_by_type(avg_loc_changes, verbose);
            break;
        }
        case ALL: {  // 9
            if(verbose > 0) cout << "All summary statistics" << endl;
            // collect_private_subclonal_bps(avg_loc_changes, verbose);
            print_bin_subclonal_stat_by_type(avg_loc_changes, 1, verbose);
            print_variance(avg_loc_changes, loc_type, 0);
            print_pairwise_divergence(ids, avg_loc_changes, verbose);
            break;
        }
        default: cout << "" << endl;
    }
}


// Get breakpoints without chromosome information
void get_breakpoints_from_pseudo(const double* avg_loc_change, vector<string>& bps, int verbose = 0){
    if(verbose > 1) cout << "Finding all the breakpoints in the CNP from just bins" << endl;

    // merge continuous regions
    int prev_bin = 0;
    double prev_rcn = avg_loc_change[0];

    int curr_bin = 0;  // pointer to current location
    double curr_rcn = 0;
    string bp = "";

    for(int loc = 1; loc < NUM_LOC; loc++){
        curr_bin = loc;
        curr_rcn = avg_loc_change[loc];
        // when starting with a specified genotype, values in avg_loc_change are absolute CNs
        double diff = fabs(prev_rcn - curr_rcn);
        // cout << prev_rcn << "\t" << curr_rcn << "\t" << diff << "\t" << fabs(diff) << endl;
        if(diff >= BP_CUTOFF){
            // cout << prev_rcn << "\t" << curr_rcn << "\t" << diff << "\t" << fabs(diff) << endl;
            bp = to_string(prev_bin);
            bps.push_back(bp);
        }
        prev_bin = curr_bin;
        prev_rcn = curr_rcn;
    }
    // For regions ending at the last position
    bp = to_string(prev_bin);
    bps.push_back(bp);

    sort(bps.begin(), bps.end());
}


// only consider subclonal events to exclude the effect of starting cell
void collect_private_subclonal_bps(const map<int, double*>& avg_loc_changes, int verbose=0){
    vector<string> bps_common;
    vector<string>::iterator it_common;
    set<string> bps_sep;

    // int nsample = samples.size();
    // cout << "There are " << nsample << " samples" << endl;

    // Unique breakpoints for each clone
    int i = 0;
    for(auto s1 : avg_loc_changes){
        // Clone* s1 = samples[i];
        vector<string> s1_bps;
        get_breakpoints_from_pseudo(s1.second, s1_bps, verbose);

        bps_sep.insert(s1_bps.begin(), s1_bps.end());

        if(i==0){
            bps_common.assign(s1_bps.begin(), s1_bps.end());
        }
        it_common = set_intersection(s1_bps.begin(), s1_bps.end(), bps_common.begin(), bps_common.end(), bps_common.begin());
        i++;
    }
    // #unique breakpoints
    cout << bps_sep.size() << endl;

    if(verbose > 0){
      bps_common.resize(it_common - bps_common.begin());
      cout << bps_common.size() << endl;

      // remove common elements from private elements
      vector<string> bps_private(bps_sep.size());
      vector<string>::iterator it_private;
      it_private = set_difference(bps_sep.begin(), bps_sep.end(), bps_common.begin(), bps_common.end(), bps_private.begin());
      bps_private.resize(it_private - bps_private.begin());
      cout << bps_private.size() << endl;

      if(verbose >= 1){
          cout << "Common breakpoints:";
          for(auto bp : bps_common){
              cout << "\t" << bp;
          }
          cout << "\ncommon breakpoint number " << bps_common.size() << endl;

          cout << "union of breakpoints:";
          for(auto bp: bps_sep){
              cout << "\t" << bp;
          }
          cout << "\nunion breakpoint number " << bps_sep.size() << endl;

          cout << "private breakpoints:";
          for(auto bp: bps_private){
              cout << "\t" << bp;
          }
          cout << endl;
      }
    }
}


// distinguish gain/loss to account for different size distribution
void print_bin_subclonal_stat_by_type(const map<int, double*>& avg_loc_changes, int by_type = 1, int verbose=0){
    int nsample = avg_loc_changes.size();
    // cout << "There are " << nsample << " samples" << endl;

    int alter_indicator_sep_gain[NUM_LOC] = {0};
    int alter_indicator_sep_loss[NUM_LOC] = {0};
    int alter_indicator_sep[NUM_LOC] = {0};
    double avg_nalter_sep_gain = 0;
    double avg_nalter_sep_loss = 0;
    double avg_nalter_sep = 0;

    for(auto s : avg_loc_changes){
        for(int i = 0; i < NUM_LOC; i++){
            if(round(s.second[i]) - START_GENOTYPE[i] > 0){
                alter_indicator_sep_gain[i] += 1;
                alter_indicator_sep[i] += 1;
            }
            else if(round(s.second[i]) - START_GENOTYPE[i] < 0){
                alter_indicator_sep_loss[i] += 1;
                alter_indicator_sep[i] += 1;
            }else{

            }
        }
    }

    // exclude clonal regions
    for(int i = 0; i < NUM_LOC; i++){
        if(alter_indicator_sep[i] > 0 && alter_indicator_sep[i] < nsample)
            avg_nalter_sep++;
    }

    for(int i = 0; i < NUM_LOC; i++){
        if(alter_indicator_sep_gain[i] > 0 && alter_indicator_sep_gain[i] < nsample)
            avg_nalter_sep_gain++;
    }

    for(int i = 0; i < NUM_LOC; i++){
        if(alter_indicator_sep_loss[i] > 0 && alter_indicator_sep_loss[i] < nsample)
            avg_nalter_sep_loss++;
    }

    // normalized by NUM_LOC to account for different choices of bin size (in real data)
    avg_nalter_sep = (double) avg_nalter_sep / NUM_LOC;
    cout << avg_nalter_sep << endl;

    if(by_type == 1){
        avg_nalter_sep_gain = (double) avg_nalter_sep_gain / NUM_LOC;
        cout << avg_nalter_sep_gain << endl;

        avg_nalter_sep_loss = (double) avg_nalter_sep_loss / NUM_LOC;
        cout << avg_nalter_sep_loss << endl;
    }
}


// count number of mutations in each gland and output their frequency. E.g how many glands with n mutations
void print_nmut_per_gland(int max_mut = 10) {
  map<int, double> nmut_by_gland;
  map<int, double> gland_by_nmut;
  for(auto c : this->clones[0]->curr_cells){
    // only count sampled cells
    if(this->sample_IDs.size() > 0 && find(this->sample_IDs.begin(), this->sample_IDs.end(), c->cell_ID) == this->sample_IDs.end())
      continue;
    nmut_by_gland[c->cell_ID] = c->muts.size();
  }

  for(auto ng: nmut_by_gland){
    // cout << ng.first << "\t" << ng.second << endl;
    gland_by_nmut[ng.second] += 1;
  }

  // cout << "number of glands for certain number of mutations" << endl;
  // for(auto gm: gland_by_nmut){
  //   cout << gm.first << "\t" << gm.second << endl;
  // }

  for(auto gm: gland_by_nmut){
    // cout << gm.first << "\t" << gm.second << endl;
    double freq = (double) gm.second / this->sample_IDs.size();
    gland_by_nmut[gm.first] = freq;
  }

  // cout << "for frequency" << endl;
  // for(auto gm: gland_by_nmut){
  //   cout << gm.first << "\t" << gm.second << endl;
  // }

  // Output cumulative distribution
  double sum = 0;
  for(int i = 0; i < max_mut; i++){
    sum += gland_by_nmut[i];
    cout << sum << endl;
  }
}


// count number of each unique mutation across all glands and output their frequency
void print_mut_freq(double min_freq=0, double max_freq=0.25, double delta = 0.001) {
  this->clones[0]->muts_freq.clear();

  for(auto c : this->clones[0]->curr_cells){
    // only count sampled cells
    if(this->sample_IDs.size() > 0 && find(this->sample_IDs.begin(), this->sample_IDs.end(), c->cell_ID) == this->sample_IDs.end())
      continue;
    for(auto mut : c->muts){
      this->clones[0]->muts_freq[mut]++;
    }
  }

  // Output mutation frequency
  vector<double> freqs;
  for(auto f : this->clones[0]->muts_freq){
    double freq = f.second / this->sample_IDs.size();
    this->clones[0]->muts_freq[f.first] = freq;
    // cout << freq << '\n';
    freqs.push_back(freq);
  }

  // Output cumulative distribution
  sort(freqs.begin(), freqs.end());
  // for(auto f: freqs){
  //   cout << f << endl;
  // }

  // make Histogram of frequencies (number of mutation in a bin (f1, f2))
  vector<int> bins;
  int start = 0;
  double f1 = min_freq;
  for(double f = min_freq + delta; f <= max_freq; f += delta){
    int count = 0;  // number of mutations in this bin
    int i = start;
    for(; i < freqs.size(); i++){
      if(freqs[i] > f){
        start = i;
        // cout << "exit at index " << i << endl;
        break;
      }
      else count += 1;
    }
    if(i == freqs.size()){
      start = i;
      // cout << "exit at index " << i << endl;
    }
    // cout << f1 << "\t" << f << "\t" << start << "\t" << count << endl;

    f1 = f;
    bins.push_back(count);
  }

  // cout << "Histogram of frequencies: " << endl;
  // for(auto b: bins){
  //   cout << b << endl;
  // }

  // compute cumulative sum
  vector<double> psum(bins.size());
  partial_sum(bins.begin(), bins.end(), psum.begin());
  for(auto s: psum){
    cout << s / freqs.size() << endl;
  }
}


// variance of CNs at each bin across all glands
void print_variance(const map<int, double*>& avg_loc_changes, int loc_type=BIN, int normalize = 1, int use_std = 0){
    // CNs at each bin across all samples
    map<int, VAR> bin_pcn;
    // map<int, NBIN> bin_cout;
    set<double> vars;   // only record unique value to ignore the effect of CNA size
    // vector<double> vars;

    int nsample = avg_loc_changes.size();
    // cout << "There are " << nsample << " samples" << endl;

    for(auto s : avg_loc_changes){
        for(int i = 0; i < NUM_LOC; i++){
            // if (s.second[i] > START_GENOTYPE[i] || s.second[i] < START_GENOTYPE[i])
            // bin_pcn[i](s.second[i]);
            bin_pcn[i](s.second[i] - START_GENOTYPE[i]);
        }
    }

    // double sum_var_orig = 0;     // use sum to avoid stochascity of events on different locations
    for(int i = 0; i < bin_pcn.size(); i++){
        // assert(boost::accumulators::count(bin_pcn[i])==nsample);
        // Compute sample variance
        double lvar = boost::accumulators::variance(bin_pcn[i]) * nsample / (nsample - 1);
        // cout << "var of loc " << i << " is " << lvar << endl;
        // only consider locations with different copy numbers
        if(lvar > 0){
          // vars.push_back(lvar);
          if(use_std == 1){
            vars.insert(sqrt(lvar));
          }else{
            vars.insert(lvar);
          }
          // sum_var_orig += lvar;
        }
    }

    // normalize variance to (0, 1)
    double sum_var = 0;
    if(vars.size() > 0){
      if (normalize == 1){
        double var_min = *min_element(vars.begin(), vars.end());
        double var_max = *max_element(vars.begin(), vars.end());
        if(var_max  - var_min > 0){
          for(auto v : vars){
              // cout << "var " << i << " is " << vars[i] << endl;
              double norm_var = (v - var_min) / (var_max - var_min);
              sum_var += norm_var;
          }
        }
      }else{
        for(auto v : vars){
            sum_var += v;
        }
      }
    }

    // double avg_var = sum_var / NUM_LOC;
    double avg_var = sum_var / vars.size();
    // cout << "original average var is " << sum_var_orig / vars.size() << '\n';
    cout << avg_var << '\n';
}


/*
Divergence of CNAs across samples within each patient was quantified by computing the
pairwise divergence between samples.
Specifically, this was the proportion of altered bins
(copy number not equal to 2 in either or both samples) that had different copy number in
each sample.
ids: The ID of each clone (cell), sorted
*/
void print_pairwise_divergence(const vector<int>& ids, map<int, double*>& avg_loc_changes, int use_cdf = 1, int verbose = 0){
    vector<double> alters;  // proportion of altered bins that are different for each pair of glands
    int ntotal = 0;
    double avg_alter;
    VAR ppalters;

    for(int i = 0; i < ids.size(); i++){
        double* lchanges1 = avg_loc_changes[ids[i]];

        for(int j = i+1; j < ids.size(); j++){
            ntotal++;

            double* lchanges2 = avg_loc_changes[ids[j]];

            int num_alter = 0;
            int num_diff = 0;

            for(int k = 0; k < NUM_LOC; k++){
                if(round(lchanges1[k]) != NORM_PLOIDY || round(lchanges2[k]) != NORM_PLOIDY){
                    num_alter++;
                    if(round(lchanges2[k]) != round(lchanges1[k])){
                        num_diff++;
                    }
                }
            }

            double prop_alter = 0;
            if(num_alter > 0) prop_alter = (double) num_diff / num_alter;
            alters.push_back(prop_alter);
            ppalters(prop_alter);
            // cout << prop_alter << endl;
        }
    }

    for(int i = 0; i < alters.size(); i++){
      avg_alter += alters[i];
    }
    avg_alter = avg_alter / alters.size();
    cout << avg_alter << endl;
    double var_alter = boost::accumulators::variance(ppalters);
    cout << var_alter << endl;
    // count number of pairs smaller than a threshold
    // map<double, int> num_pairs;
    if(verbose > 0){
        cout << "total number of pairs " << ntotal << endl;
        cout << "count number of pairs smaller than a threshold" << endl;
    }

    if(use_cdf == 1){
      for(double i = 0.1; i <= 0.5; i += 0.1){
          // cout << i << endl;
          int npair = 0;
          for(auto p : alters){
              if(p < i) npair++;
          }
          cout << (double) npair / ntotal << endl;
      }
    }
}


// Extract the sampled cells from the given IDs
void get_sampled_cells(){
  if(this->sample_IDs.size() == 0)  return;
  sampled_cells.clear();
  // cells are stored in a vector of pointer, hard to access by int IDs directly
  for(auto c1 : this->clones[0]->curr_cells){
    // only count sampled cells
    if(find(this->sample_IDs.begin(), this->sample_IDs.end(), c1->cell_ID) == this->sample_IDs.end()){
      continue;
    }
    sampled_cells.push_back(c1);
  }
  assert(sampled_cells.size() == this->sample_IDs.size());
}

/*
Print number of segragting mutations (mutations present in at least one sample)
Under infinite site model, each site corresponds to a mutation
*/
void print_num_uniq_mut(int verbose = 0){
  get_sampled_cells();
  set<string> all_muts;

  for(auto c : sampled_cells){
    copy(c->muts.begin(), c->muts.end(), inserter(all_muts,all_muts.end()));
  }

  if(verbose > 0) cout << "Number of unique mutation is " << endl;
  cout << all_muts.size() << endl;
}



/*
Mismatch of CNAs across samples within each patient was quantified by computing the
pairwise mismatch between samples. (similar to Tajima's estimator)
In Tajima's estimator, each mutation occurs on a different site
Here, only unique mutations are counted by their IDs
*/
void print_pairwise_mismatch(int verbose = 0){
    vector<int> ndiffs;  // proportion of altered bins that are different for each pair of glands
    double avg_ndiff = 0;

    get_sampled_cells();
    // cout << this->sample_IDs.size() << endl;

    sort(sampled_cells.begin(), sampled_cells.end(), compare_ptr2cell);
    // for(auto c1 : sampled_cells){
    //   cout << c1->cell_ID << endl;
    // }

    // cells are stored in a vector of pointer, hard to access by int IDs directly
    for(auto c1 : sampled_cells){
        for(auto c2 : sampled_cells){
          // only count sampled cells
          if (c2->cell_ID <= c1->cell_ID) continue;
          // Find the number of different mutations between two samples (cells)
          vector<string> diff_muts;
          set_difference(c1->muts.begin(), c1->muts.end(), c2->muts.begin(), c2->muts.end(), inserter(diff_muts, diff_muts.begin()));
          ndiffs.push_back(diff_muts.size());
        }
    }

    for(int i = 0; i < ndiffs.size(); i++){
      avg_ndiff += ndiffs[i];
    }
    avg_ndiff = avg_ndiff / ndiffs.size();
    cout << avg_ndiff << '\n';

    // count number of pairs smaller than a threshold
    // map<double, int> num_pairs;
    if(verbose > 0){
        cout << "total number of pairs " << ndiffs.size() << endl;
        for(int i = 0; i < ndiffs.size(); i++){
          cout << ndiffs[i] << endl;
        }
    }
}



// Print CNP of selected glands in a single file
void print_sample_cnp(const map<int, double*>& avg_loc_changes, string outdir, string suffix, int verbose = 0){
    string fname = outdir + "gland_cn"  + suffix + ".txt";
    ofstream fcn;
    fcn.open(fname, ofstream::trunc | ofstream::out);

    for(auto s : avg_loc_changes){
        // if(verbose > 0) cout << "Print CNP of gland " << s->clone_ID << endl;
        if(verbose > 0) cout << "Print CNP of gland " << s.first << endl;
        // Find the average CNP of all cells
        // double avg_loc_change[NUM_LOC] = {0};
        // s->set_bulk_cnp_from_pseudo(avg_loc_change, s->curr_cells);

        dpcn cnp;
        loc2pcn(s.second, cnp, verbose);

        // string fname = outdir + "bulk_cn"  + suffix  + "_clone" + to_string(clone_ID) + "_sample" + to_string(sample_ID) + ".txt";
        string sname = "g" + to_string(s.first);
        print_bulk_cn(sname, cnp, fcn, verbose);
    }

    fcn.close();
}


// Print CNP of selected glands in a single file
void print_cell_cnp(Clone* gland, string outdir, string suffix, int verbose = 0){
    string fname = outdir + "gland_cn"  + suffix + ".txt";
    ofstream fcn;
    fcn.open(fname, ofstream::trunc | ofstream::out);
    // vector<pcn> cnp_s1, cnp_s2;
    // vector<string> sname_s1, sname_s2;

    for(auto s : gland->curr_cells){
        if(sample_IDs.size() > 0 && find(sample_IDs.begin(), sample_IDs.end(), s->cell_ID) == sample_IDs.end()) continue;

        if(verbose > 0) cout << "Print CNP of gland " << s->cell_ID << endl;

        pcn cnp;
        loc2pcn(s->loc_changes, cnp, verbose);

        // string fname = outdir + "bulk_cn"  + suffix  + "_clone" + to_string(clone_ID) + "_sample" + to_string(sample_ID) + ".txt";
        string sname = "G" + to_string(s->cell_ID);
        if(sample_IDs_side1.size() > 0 && find(sample_IDs_side1.begin(), sample_IDs_side1.end(), s->cell_ID) != sample_IDs_side1.end()){
            sname = "E1_" + sname;
            // cnp_s1.push_back(cnp);
            // sname_s1.push_back(sname);
        }
        else if(sample_IDs_side2.size() > 0 && find(sample_IDs_side2.begin(), sample_IDs_side2.end(), s->cell_ID) != sample_IDs_side2.end()){
            sname = "E2_" + sname;
            // cnp_s2.push_back(cnp);
            // sname_s2.push_back(sname);
        }else{

        }

        print_bulk_cn(sname, cnp, fcn, verbose);
    }
    // for(int i = 0; i < cnp_s1.size(); i++){
    //     print_bulk_cn(sname_s1[i], cnp_s1[i], fcn, verbose);
    // }
    // for(int i = 0; i < cnp_s2.size(); i++){
    //     print_bulk_cn(sname_s2[i], cnp_s2[i], fcn, verbose);
    // }

    fcn.close();
}


};




#endif
