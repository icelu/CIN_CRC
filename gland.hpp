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
#include "clone.hpp"

using namespace std;


typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance>> VAR;
// typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::count>> NBIN;


// Glands may have several clones sampled
class Glands{
public:
    int glands_ID;

    vector<Clone*> clones;      // For all glands
    vector<Clone*> samples;     // For selected glands
    vector<int> sample_IDs;
    vector<int> sample_IDs_side1;
    vector<int> sample_IDs_side2;

    node* root;     // root node of gland lineage tree

    set<string> bps_all;

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
            }else{
                c->clone_ID = s2->clone_ID;
                s2->curr_cells.push_back(c);
            }
        }
        s0->curr_cells.clear();
        // keep grow the current demes
        this->clones.push_back(s1);
        this->clones.push_back(s2);

        if(store_lineage != 0){
            string lstr = to_string(s0->clone_ID) +  "\t" + to_string(s1->clone_ID) + "\t" + to_string(s1->curr_cells.size());
            lineages.push_back(lstr);
            lstr = to_string(s0->clone_ID) +  "\t" + to_string(s2->clone_ID) + "\t" + to_string(s2->curr_cells.size());
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
            string lstr = to_string(d->parent_ID) +  "\t" + to_string(d->cell_ID) + "\t" + to_string(d->time_occur);
            lineages.push_back(lstr);
        }
    }
}


// sample from a larger population of cells to intimate spatial sampling and increase diversity among cells
void sample_gland_from_cell(const vector<int>& nglands, int verbose = 0){
    sample_IDs.clear();
    sample_IDs_side1.clear();
    sample_IDs_side2.clear();

    int nside = nglands.size();
    assert(nside = 2);

    int ng1 = nglands[0];
    int ng2 = nglands[1];

    // cout << ng1 << " glands taken at side 1" << endl;
    // cout << ng2 << " glands taken at side 2" << endl;

    int s1 = 0; // track number of sampled glands at each side

    // find all leaves at one side
    vector<int> nodes1;
    get_leaves_below(root->left, nodes1);
    vector<int> nodes2;
    get_leaves_below(root->right, nodes2);

    int nsize1 = nodes1.size();
    int nsize2 = nodes2.size();

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

    if(verbose > 0){
        cout << nsize1 << " glands at side 1" << endl;
        cout << nsize2 << " glands at side 2" << endl;
        cout << sample_IDs.size() << " glands total taken" << endl;
    }

}


void  print_sstat(const vector<int>& ids, map<int, double*>& avg_loc_changes, int stat_type, int loc_type, int verbose = 0) {
    if(verbose > 0){
        cout << "Printing summary statistics" << endl;
    }
    switch (stat_type) {
        case LVAR:{ // 0
            print_variance(avg_loc_changes, loc_type);
            break;
        }
        case ADIFF: {
            // cout << "Using pairwise difference across locations" << endl;
            // print_diversity_by_site();
            break;
        }
        case DIFF: {
            // cout << "Using pairwise difference" << endl;
            // print_diversity_by_site();
            break;
        }
        case BP: {
            collect_private_subclonal_bps(avg_loc_changes, verbose);
            break;
        }
        case ALTBIN_SEP: {
            print_bin_subclonal_stat_by_type(avg_loc_changes, verbose);
            break;
        }
        case BP_BIN: {
            collect_private_subclonal_bps(avg_loc_changes, verbose);
            print_bin_subclonal_stat_by_type(avg_loc_changes, verbose);
            break;
        }
        case ALL: {
            if(verbose > 0) cout << "All summary statistics" << endl;
            // collect_private_subclonal_bps(avg_loc_changes, verbose);
            print_bin_subclonal_stat_by_type(avg_loc_changes, verbose);
            print_variance(avg_loc_changes, loc_type);
            print_pairwise_divergence(ids, avg_loc_changes, verbose);
            // print_diversity_by_site();
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
    int start = prev_bin;

    int curr_bin = 0;
    double curr_rcn = 0;
    string bp = "";

    for(int loc = 1; loc < NUM_LOC; loc++){
        curr_bin = loc;
        curr_rcn = avg_loc_change[loc];
        double diff = fabs(prev_rcn - curr_rcn);
        // cout << prev_rcn << "\t" << curr_rcn << "\t" << diff << endl;
        // cout << "cn difference is " << fabs(diff) << endl;
        if(curr_bin - prev_bin != 1 || diff >= BP_CUTOFF){
            bp = to_string(prev_bin);
            bps.push_back(bp);
            // start a new region
            start = curr_bin;
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
    bps_common.resize(it_common - bps_common.begin());

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


// distinguish gain/loss to account for different size distribution
void print_bin_subclonal_stat_by_type(const map<int, double*>& avg_loc_changes, int verbose=0){
    int nsample = avg_loc_changes.size();
    // cout << "There are " << nsample << " samples" << endl;
    // map<int, double> avg_nalter_sep_gain;
    // map<int, double> avg_nalter_sep_loss;
    // count number of times a bin is altered in a clone. May separate different sides?
    // map<int, int*> alter_indicator_sep_gain;
    // map<int, int*> alter_indicator_sep_loss;
    int alter_indicator_sep_gain[NUM_LOC] = {0};
    int alter_indicator_sep_loss[NUM_LOC] = {0};
    int avg_nalter_sep_gain = 0;
    int avg_nalter_sep_loss = 0;

    for(auto s : avg_loc_changes){
        for(int i = 0; i < NUM_LOC; i++){
            if(s.second[i] - START_GENOTYPE[i] >= BIN_CUOFF){
                alter_indicator_sep_gain[i] += 1;
            }
            else if(s.second[i] - START_GENOTYPE[i] <= -BIN_CUOFF){
                alter_indicator_sep_loss[i] += 1;
            }else{

            }
        }
    }

    // exclude clonal regions
    // for(auto cn: alter_indicator_sep_gain){
    for(int i = 0; i < NUM_LOC; i++){
        // cout << cn.first << "\t" << cn.second[i] << endl;
        if(alter_indicator_sep_gain[i] > 0 && alter_indicator_sep_gain[i] < nsample)
            avg_nalter_sep_gain++;
    }
    // }

    // for(auto cn: alter_indicator_sep_loss){
    for(int i = 0; i < NUM_LOC; i++){
        // cout << cn.first << "\t" << cn.second[i] << endl;
        if(alter_indicator_sep_loss[i] > 0 && alter_indicator_sep_loss[i] < nsample)
            avg_nalter_sep_loss++;
    }
    // }

    // normalized by NUM_LOC to account for different choices of bin size (in real data)
    avg_nalter_sep_gain = avg_nalter_sep_gain * 100 / NUM_LOC;
    cout << avg_nalter_sep_gain << endl;

    avg_nalter_sep_loss = avg_nalter_sep_loss * 100 / NUM_LOC;
    cout << avg_nalter_sep_loss << endl;
}


// variance of CNs at each bin across all glands
void print_variance(const map<int, double*>& avg_loc_changes, int loc_type=BIN){
    // CNs at each bin across all samples
    map<int, VAR> bin_pcn;
    // map<int, NBIN> bin_cout;

    int nsample = avg_loc_changes.size();
    // cout << "There are " << nsample << " samples" << endl;

    for(auto s : avg_loc_changes){
        for(int i = 0; i < NUM_LOC; i++){
            bin_pcn[i](s.second[i]);
        }
    }

    // int nvar = 0;   // number of locations with variance != 0
    double sum_var = 0;     // use sum to avoid stochascity of events on different locations
    for(int i = 0; i < NUM_LOC; i++){
        // assert(boost::accumulators::count(bin_pcn[i])==nsample);
        // Compute sample variance
        double lvar = boost::accumulators::variance(bin_pcn[i]) * nsample / (nsample - 1);
        // if(lvar > 0) nvar++;
        // cout << boost::accumulators::variance(bin_pcn[i]) * nsample / (nsample - 1) << endl;
        sum_var += lvar;
    }

    double avg_var = sum_var * 100 / NUM_LOC;
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
void print_pairwise_divergence(const vector<int>& ids, map<int, double*>& avg_loc_changes, int verbose = 0){
    vector<double> alters;
    int ntotal = 0;
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
            // cout << prop_alter << endl;
        }
    }

    // count number of pairs smaller than a threshold
    // map<double, int> num_pairs;
    if(verbose > 0){
        cout << "total number of pairs " << ntotal << endl;
        cout << "count number of pairs smaller than a threshold" << endl;
    }
    for(double i = 0.1; i <= 0.5; i += 0.1){
        // cout << i << endl;
        int npair = 0;
        for(auto p : alters){
            if(p < i) npair++;
        }
        cout << (double) npair / ntotal << endl;
    }
}


// Print CNP of selected glands in a single file
void print_sample_cnp(const map<int, double*>& avg_loc_changes, string outdir, string suffix, int verbose = 0){
    string fname = outdir + "gland_cn"  + suffix + ".txt";
    ofstream fcn;
    fcn.open(fname, std::ofstream::trunc | std::ofstream::out);

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
    fcn.open(fname, std::ofstream::trunc | std::ofstream::out);
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
        if(sample_IDs_side2.size() > 0 && find(sample_IDs_side2.begin(), sample_IDs_side2.end(), s->cell_ID) != sample_IDs_side2.end()){
            sname = "E2_" + sname;
            // cnp_s2.push_back(cnp);
            // sname_s2.push_back(sname);
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
