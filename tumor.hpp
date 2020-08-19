#ifndef TUMOR_HPP
#define TUMOR_HPP

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


// a tumor may have several clones sampled
class Tumor{
public:
    int tumor_ID;

    vector<Clone*> clones;   // For all clones
    vector<Clone*> metastasis;
    vector<Clone*> primary;

    vector<Sample*> samples;

    set<string> bps_all;
    set<string> bps_prim;
    set<int> id_mets;
    map<int, set<string>> bps_met_sep;
    map<string, int> bps_prim_count;
    map<int, map<string, int>> bps_met_sep_count;

    Tumor(const Tumor& other) = default;
    Tumor& operator=(const Tumor& other) = default;
    Tumor& operator=(Tumor&& other) = default;

    Tumor(){
        tumor_ID = 0;
    }

    ~Tumor(){
        // cout << "Release " << clones.size() << " clones " << endl;
        for(auto cl: clones){
            // cout << "clone at " << cl << endl;
            delete cl;
        }
    }


    // Simulate the growth of demes (crypts). Used for tracking deme relationships
    // A deme represents a specific population of cells, say cells in a crypt
    // The first deme is generated via birth-and-death process, beginning with a single transformed founding tumor cell
    // Once a deme exceeds the maximum size (10,000 cells), it splits into two offspring demes via random sampling of cells
    void simulate_deme_partition(const Cell_ptr start_cell, int ndeme = 100, int max_deme_size = 10000, string fname = "", int verbose = 0){
        ofstream fout(fname);
        string header = "ID\tParent\tNcell";
        fout << header << endl;

        int num_clone = 1;
        Clone* s = new Clone(num_clone, 0);
        s->grow(start_cell, max_deme_size, verbose, 1);
        this->clones.push_back(s);

        while (this->clones.size() < ndeme) {
            // randomly select a deme to grow and split
            // cout << "Deme number " << this->clones.size() << endl;
            int rindex = myrng(this->clones.size());
            Clone* s0 = this->clones[rindex];
            while (s0->id_curr_cells.size() < max_deme_size) {
                s0->grow(start_cell, max_deme_size, verbose, 0);
            }
            // split deme
            Clone* s1 = new Clone(num_clone+1, s0->clone_ID);
            Clone* s2 = new Clone(num_clone+2, s0->clone_ID);
            num_clone += 2;
            for(auto c : s0->id_curr_cells){
                double rc = runiform(r, 0, 1);
                if(rc < 0.5){
                    s1->id_curr_cells.push_back(c);
                }else{
                    s2->id_curr_cells.push_back(c);
                }
            }
            // keep grow the current demes
            this->clones.push_back(s1);
            this->clones.push_back(s2);

            fout << s0->clone_ID << "\t"  << s1->clone_ID << "\t" << s1->id_curr_cells.size() << endl;
            fout << s0->clone_ID << "\t"  << s2->clone_ID << "\t" << s2->id_curr_cells.size() << endl;

            delete (this->clones[rindex]);
            this->clones[rindex] = NULL;
            this->clones.erase(this->clones.begin() + rindex);
        }
        if(verbose > 0){
            for(auto d : this->clones){
                cout << d->clone_ID << "\t"  << d->parent_ID << "\t"  << d->id_curr_cells.size() << endl;
            }
        }
    }


    /*********************** Functions related to primary-metastasis samples *************************************/
    /*
    All clones end at the same time
    */
    void simulate_metastasis(int nmeta, const vector<int>& time_migration, const vector<double>& tend_migration, const vector<Model>& model_migration, const Cell_ptr start_cell, const Model& start_model, int size_primary, int loc_type=BIN, double leap_size=0, int verbose = 0) {
        // grow the primary clone (there may be subclones)
        int num_clone = 0;
        Clone* s0 = new Clone(num_clone++, "Primary", 0);
        s0->grow_with_cnv(start_cell, start_model, time_migration[0], loc_type, leap_size, verbose);
        if(verbose > 0) cout << "Primary clone stops at time " << s0->time_end << " with " << time_migration[0] << " cells at address " << s0 << endl;

        // Store the starting cells of metastasis clones
        vector<Cell_ptr> start_rcells_meta;
        // randomly pick a cell to start a metastasis
        cout << "#CNAs in starting cell of migration:";
        for(int c = 1; c <= nmeta; c++){
            if(verbose > 0){
                cout << "\nSimulating metastasis clone " << c << " when primary clone has cell number " << time_migration[c-1] << endl;
            }

            // randomly pick a cell in the clone which appear at the stop time
            int sel = runiform(r, 0, s0->curr_cells.size());
            // int sel = s0->curr_cells.size() - 1;
            // int sel = 0;
            // // Find cells with most CN changes
            // double max_cn_change = 0;
            // for(int i = 0; i < s0->curr_cells.size(); i++){
            //     Cell_ptr cell = s0->curr_cells[i];
            //     double gdiff = 0;
            //     for(int i = 0; i < NUM_LOC; i++){
            //         gdiff += abs(cell->loc_changes[i]);
            //     }
            //     if(gdiff > max_cn_change){
            //         max_cn_change = gdiff;
            //         sel = i;
            //     }
            // }

            Cell_ptr rcell = s0->curr_cells[sel];
            // Remove this cell from primary clone
            s0->curr_cells.erase(s0->curr_cells.begin()+sel);

            assert(s0->curr_cells.size() == time_migration[c-1] - 1);
            if(verbose > 0){
                cout << "Selected cell " << rcell->cell_ID << " with " << rcell->num_mut << " CNAs " << endl;
                // rcell->print_cell_info();
                // cout << "Primary clone has " << s0->curr_cells.size() << " cells after migration" << endl;
            }
            cout << "\t" << rcell->num_mut;

            // start a new clone with selected cell
            double offset = MIGRATE_OFFSET * c;
            rcell->update_loc(offset, offset, offset);
            start_rcells_meta.push_back(rcell);

            // continue the growth of primary clone
            if(c < nmeta){
                s0->grow_with_cnv(start_cell, start_model, time_migration[c], loc_type, leap_size, verbose, 0);
                if(verbose > 0) cout << "Keep growing primary clone until time " << s0->time_end << " with " << s0->curr_cells.size() << " cells" << endl;

            }
        }
        cout << endl;

        if(s0->curr_cells.size() < size_primary){
            s0->grow_with_cnv(start_cell, start_model, size_primary, loc_type, leap_size, verbose, 0);
            if(verbose > 0){
                cout << "Simulated primary clone has " << s0->num_novel_mutation << " de novo CNAs in " << s0->curr_cells.size() << " cells, growth stopped at time " << s0->time_end << endl;
            }
        }
        primary.push_back(s0);
        clones.push_back(s0);

        int c = 1;
        for(auto rcell : start_rcells_meta){
            string name = "Metastasis " + to_string(c);
            Clone* s1 = new Clone(num_clone++, name, s0->time_end);
            rcell->clone_ID = s1->clone_ID;
            rcell->cell_ID = 1;
            // allowing to grow to at most size_primary cells to save time
            s1->grow_with_cnv(rcell, model_migration[c-1], size_primary, loc_type, leap_size, verbose, 1, s0->time_end + tend_migration[c-1]);
            metastasis.push_back(s1);
            clones.push_back(s1);
            if(verbose > 0){
                // , stored at address " << s1
                cout << "Simulated metastasis clone " << c << " has " << s1->num_novel_mutation << " de novo CNAs in " << s1->curr_cells.size() << " cells, growth stopped at time " << s1->time_end << endl;
            }
            c++;
        }

        cout << "End time of each clone:";
        for(auto c : clones){
            cout << "\t" << c->time_end;
        }
        cout << endl;

        cout << "Size of each clone:";
        for(auto c : clones){
            cout << "\t" << c->curr_cells.size();
        }
        cout << endl;
    }



        /*
        Assume the growing process starts at fitness peak
        */
        void simulate_metastasis_from_peak(int Nprim, int nmeta, const vector<int>& time_migration, const vector<double>& tend_migration, const vector<Model>& model_migration, const Cell_ptr start_cell, const Model& start_model, int size_primary, int loc_type=BIN, double leap_size=0, int verbose = 0) {
            // grow the primary clone (there may be subclones)
            int num_clone = 0;
            Clone* s0 = new Clone(num_clone++, "Primary", 0);
            // insert available cells
            for(int i = 0; i < Nprim; i++){
                s0->id_curr_cells.push_back(i+1);
            }
            s0->grow_with_cnv(start_cell, start_model, time_migration[0], loc_type, leap_size, verbose, 0);
            if(verbose > 0) cout << "Primary clone stops at time " << s0->time_end << " with " << time_migration[0] << " cells at address " << s0 << endl;

            // Store the starting cells of metastasis clones
            vector<Cell_ptr> start_rcells_meta;
            // randomly pick a cell to start a metastasis
            cout << "#CNAs in starting cell of migration:";
            for(int c = 1; c <= nmeta; c++){
                if(verbose > 0){
                    cout << "\nSimulating metastasis clone " << c << " when primary clone has cell number " << time_migration[c-1] << endl;
                }

                // randomly pick a cell in the clone which appear at the stop time
                int sel = runiform(r, 0, s0->curr_cells.size());
                // int sel = s0->curr_cells.size() - 1;
                // int sel = 0;
                // // Find cells with most CN changes
                // double max_cn_change = 0;
                // for(int i = 0; i < s0->curr_cells.size(); i++){
                //     Cell_ptr cell = s0->curr_cells[i];
                //     double gdiff = 0;
                //     for(int i = 0; i < NUM_LOC; i++){
                //         gdiff += abs(cell->loc_changes[i]);
                //     }
                //     if(gdiff > max_cn_change){
                //         max_cn_change = gdiff;
                //         sel = i;
                //     }
                // }

                Cell_ptr rcell = s0->curr_cells[sel];
                // Remove this cell from primary clone
                s0->curr_cells.erase(s0->curr_cells.begin()+sel);

                assert(s0->curr_cells.size() == time_migration[c-1] - 1);
                if(verbose > 0){
                    cout << "Selected cell " << rcell->cell_ID << " with " << rcell->num_mut << " CNAs " << endl;
                    // rcell->print_cell_info();
                    // cout << "Primary clone has " << s0->curr_cells.size() << " cells after migration" << endl;
                }
                cout << "\t" << rcell->num_mut;

                // start a new clone with selected cell
                double offset = MIGRATE_OFFSET * c;
                rcell->update_loc(offset, offset, offset);
                start_rcells_meta.push_back(rcell);

                // continue the growth of primary clone
                if(c < nmeta){
                    s0->grow_with_cnv(start_cell, start_model, time_migration[c], loc_type, leap_size, verbose, 0);
                    if(verbose > 0) cout << "Keep growing primary clone until time " << s0->time_end << " with " << s0->curr_cells.size() << " cells" << endl;

                }
            }
            cout << endl;

            if(s0->curr_cells.size() < size_primary){
                s0->grow_with_cnv(start_cell, start_model, size_primary, loc_type, leap_size, verbose, 0);
                if(verbose > 0){
                    cout << "Simulated primary clone has " << s0->num_novel_mutation << " de novo CNAs in " << s0->curr_cells.size() << " cells, growth stopped at time " << s0->time_end << endl;
                }
            }
            primary.push_back(s0);
            clones.push_back(s0);

            int c = 1;
            for(auto rcell : start_rcells_meta){
                string name = "Metastasis " + to_string(c);
                Clone* s1 = new Clone(num_clone++, name, s0->time_end);
                rcell->clone_ID = s1->clone_ID;
                rcell->cell_ID = 1;
                // allowing to grow to at most size_primary cells to save time
                s1->grow_with_cnv(rcell, model_migration[c-1], size_primary, loc_type, leap_size, verbose, 1, s0->time_end + tend_migration[c-1]);
                metastasis.push_back(s1);
                clones.push_back(s1);
                if(verbose > 0){
                    // , stored at address " << s1
                    cout << "Simulated metastasis clone " << c << " has " << s1->num_novel_mutation << " de novo CNAs in " << s1->curr_cells.size() << " cells, growth stopped at time " << s1->time_end << endl;
                }
                c++;
            }

            cout << "End time of each clone:";
            for(auto c : clones){
                cout << "\t" << c->time_end;
            }
            cout << endl;

            cout << "Size of each clone:";
            for(auto c : clones){
                cout << "\t" << c->curr_cells.size();
            }
            cout << endl;
        }



    // sort the cells by distance to the origin and take same number of cells for each region
    void get_sample_regions_even(int num_sample, const Clone* clone, vector<Coord>& region_start, vector<Coord>& region_end, int verbose = 0){
        vector<Coord> all_pos;
        for(auto c : clone->curr_cells){
            all_pos.push_back(c->pos);
        }
        sort(all_pos.begin(), all_pos.end());

        double min_x = clone->min_pos.x;
        double min_y = clone->min_pos.y;
        double min_z = clone->min_pos.z;
        double max_x = clone->max_pos.x;
        double max_y = clone->max_pos.y;
        double max_z = clone->max_pos.z;

        vector<double> bp_x;
        bp_x.push_back(min_x);
        int num_cell = clone->curr_cells.size();
        for(int i = 1; i < num_sample; i++){
            int b = floor(i * num_cell / num_sample) - 1;
            double bp = all_pos[b].x;
            assert(bp >= bp_x[i-1]);
            bp_x.push_back(bp);
        }
        assert(max_x >= bp_x[num_sample-1]);
        bp_x.push_back(max_x);

        double x1, y1, z1, x2, y2, z2;
        for(int i = 0; i < num_sample; i++){
            x1 = bp_x[i];
            y1 = min_y;
            z1 = min_z;
            x2 = bp_x[i+1];
            y2 = max_y;
            z2 = max_z;
            assert(x1<=x2 && y1<=y2 && z1<=z2);
            Coord start{x1, y1, z1};
            Coord end{x2, y2, z2};
            region_start.push_back(start);
            region_end.push_back(end);
        }
    }


    void simulate_multiple_samples(const vector<Clone*>& clones, const vector<int>& nsamples, string outdir, string suffix, double frac_cutoff=0.5, int stat_type = 0, int loc_type=BIN, int verbose=0){
        int nclone = clones.size();
        if(verbose > 0) cout << "\nSimulating multi-sample bulk CNP from " << nclone << " clones" << endl;
        assert(clones.size() == nsamples.size());
        // cout << "Number of mutations:";
        for(int k = 0; k < nclone; k++){
            int num_sample = nsamples[k];
            if(verbose > 0){
                cout << "\nSimulating " << num_sample << " samples for clone " << clones[k]->clone_ID << " " << clones[k]->name << " with mutations " << clones[k]->get_num_mut() << endl;
            }

            clones[k]->get_cell_loc(clones[k]->curr_cells, verbose);

            vector<Coord> region_start;
            vector<Coord> region_end;
            if(verbose > 0){
                std::cout << "randomly pick non-overlapping regions" << '\n';
            }
            get_sample_regions_even(num_sample, clones[k], region_start, region_end, verbose);

            clones[k]->sample_multi_regions(num_sample, region_start, region_end, loc_type, verbose);
            for(auto s : clones[k]->samples)  samples.push_back(s);
        }
        // cout << "Printing pairwise differences of CNs between samples" << endl;
        // cout << "\n";
        cout << "Sample size:";
        for(auto s : samples){
            cout << "\t" << s->sample_cells.size();
        }
        cout << "\n";
        switch (stat_type) {
            case LVAR:{ // 0
                // print_variance(samples, loc_type);
                print_variance_by_site(loc_type);
                break;
            }
            case ADIFF: {
                // cout << "Using pairwise difference across locations" << endl;
                print_diversity_by_site();
                break;
            }
            case AVG: {
                // cout << "Using average CNP" << endl;
                print_avg_cnp();
                break;
            }
            case CMPL: {
                print_cmpl_cnp(loc_type);
                break;
            }
            case DIFF: {
                // cout << "Using pairwise difference" << endl;
                print_diversity_by_site();
                break;
            }
            case BP: {
                collect_private_subclonal_bps(verbose);
                break;
            }
            // print_breakpoints_stat(nsamples, loc_type, verbose); break;
            case ALTBIN: {
                print_bin_subclonal_stat(verbose);
                break;
            }
            case ALTBIN_SEP: {
                print_bin_subclonal_stat_by_type(verbose);
                break;
            }
            case BP_BIN: {
                collect_private_subclonal_bps(verbose);
                print_bin_subclonal_stat_by_type(verbose);
                break;
            }
            case ALL: {
                collect_private_subclonal_bps(verbose);
                print_bin_subclonal_stat_by_type(verbose);
                // print_variance(samples, loc_type);
                print_variance_by_site(loc_type);
                print_diversity_by_site();
                break;
            }
            default: cout << "" << endl;
        }
        if(verbose > 0)     print_sample_cnp(outdir, suffix, verbose);
        if(verbose > 1)     print_sample_cnp_orig(outdir, suffix, verbose);
    }

    // Print the complete CNP of all samples, separated by ; so that it can be read by Julia into a matrix
    void print_cmpl_cnp(int loc_type=BIN){
        int nsample = samples.size();
        // cout << "There are " << nsample << " samples" << endl;
        for(int i = 0; i < nsample; i++){
                Sample* s1 = samples[i];
                // cout << s1->sample_ID << "," << s1->clone_ID << "\t" << s2->sample_ID << "," << s2->clone_ID << "\n";
                int k = 0;
                for(int i = 0; i < NUM_LOC; i++){
                    cout << s1->sample_loc_change[i] << " ";
                }
                cout << ";";
        }
    }


    // Print the average CNP of all samples
    void print_avg_cnp(){
        // compute pairwise differences of CNP between samples
        vector<double> avg_pcn(NUM_LOC, 0.0);

        int nsample = samples.size();
        // cout << "There are " << nsample << " samples" << endl;
        for(int i = 0; i < nsample; i++){
                Sample* s1 = samples[i];
                // cout << s1->sample_ID << "," << s1->clone_ID << "\t" << s2->sample_ID << "," << s2->clone_ID << "\n";
                int k = 0;
                for(int i = 0; i < NUM_CHR; i++){
                    for(int j = 0; j < CHR_BIN_SIZE[i]; j++){
                        pair<int, int> pos(i,j);
                        int cn1 = 0;
                        if(s1->sample_pcn.find(pos) != s1->sample_pcn.end()) cn1 = s1->sample_pcn[pos];
                        avg_pcn[k++] += abs(cn1);
                    }
                }
                assert(k == NUM_LOC);
        }

        for(int i = 0; i < avg_pcn.size(); i++){
            // if(diff_pcn[i]!=0) cout << i << "\t" << diff_pcn[i] << "\n";
            cout << avg_pcn[i]/nsample << endl;
        }
    }


    // Print the diversity vector between samples
    void print_diversity(const vector<Sample*>& samples){
        // compute pairwise differences of CNP between samples
        // vector<double> diff_pcn(NUM_LOC, 0.0);

        double avg_diff_all = 0;
        int npair = 0;
        int nsample = samples.size();
        // cout << "There are " << nsample << " samples" << endl;
        for(int i = 0; i < nsample - 1; i++){
            for(int j = i + 1; j < nsample; j++){
                Sample* s1 = samples[i];
                Sample* s2 = samples[j];
                // cout << s1->sample_ID << "," << s1->clone_ID << "\t" << s2->sample_ID << "," << s2->clone_ID << "\n";
                // int k = 0;
                // for(int i = 0; i < NUM_CHR; i++){
                //     for(int j = 0; j < CHR_BIN_SIZE[i]; j++){
                //         pair<int, int> pos(i,j);
                //         double cn1 = 0;
                //         double cn2 = 0;
                //         if(s1->sample_pcn.find(pos) != s1->sample_pcn.end()) cn1 = s1->sample_pcn[pos];
                //         if(s2->sample_pcn.find(pos) != s2->sample_pcn.end()) cn2 = s2->sample_pcn[pos];
                //         double d = abs(cn1 - cn2);
                //         // cout << i << "\t" << j << "\t" << d << "\n";
                //         diff_pcn[k++] += d;
                //     }
                // }
                // assert(k == NUM_LOC);
                double avg_diff = 0;
                npair++;
                for(int i = 0; i < NUM_LOC; i++){
                    double cn1 = s1->sample_loc_change[i];
                    double cn2 = s2->sample_loc_change[i];
                    avg_diff += abs(cn1 - cn2);
                }
                // avg_diff = avg_diff / NUM_LOC;
                // cout << avg_diff << '\n';
                avg_diff_all += avg_diff;
            }
        }
        avg_diff_all = (double) avg_diff_all / npair;
        cout << avg_diff_all << endl;
        // for(int i = 0; i < diff_pcn.size(); i++){
        //     // if(diff_pcn[i]!=0) cout << i << "\t" << diff_pcn[i] << "\n";
        //     cout << diff_pcn[i] << endl;
        // }
    }


    // Print the average of diversity between samples taken at a site
    void print_diversity_by_site(){
        for(auto c : clones){
            // cout << "pairwise diversity for clone " << c->clone_ID << endl;
            print_diversity(c->samples);
        }
    }

    void collect_bps(vector<string>& bps_prim_private, map<int, vector<string>>& bps_met_privates, int verbose=0){
        set<string> bps_mets;
        int nsample = samples.size();
        // cout << "There are " << nsample << " samples" << endl;

        // Unique breakpoints for each clone
        for(int i = 0; i < nsample; i++){
            Sample* s1 = samples[i];
            if(s1->clone_ID==0){
                bps_prim.insert(s1->sample_bps.begin(), s1->sample_bps.end());
                bps_all.insert(s1->sample_bps.begin(), s1->sample_bps.end());
                for(auto bp : s1->sample_bps){
                    bps_prim_count[bp]++;
                }
            }else{
                bps_mets.insert(s1->sample_bps.begin(), s1->sample_bps.end());
                bps_all.insert(s1->sample_bps.begin(), s1->sample_bps.end());
                id_mets.insert(s1->clone_ID);
                bps_met_sep[s1->clone_ID].insert(s1->sample_bps.begin(), s1->sample_bps.end());
                for(auto bp : s1->sample_bps){
                    bps_met_sep_count[s1->clone_ID][bp]++;
                }
            }
        }

        if(verbose >= 1){
            cout << "union of primary breakpoints of size " << bps_prim.size() << ":";
            for(auto bp: bps_prim){
                cout << "\t" << bp;
            }
            cout << endl;

            cout << "union of metastasis breakpoints of size " << bps_mets.size() << ":";
            for(auto bp: bps_mets){
                cout << "\t" << bp;
            }
            cout << endl;
        }

        // for metastasis samples
        for(auto k : id_mets){
            set<string> met = bps_met_sep[k];
            vector<string>::iterator it_met_private;
            vector<string> bps_met_private(met.size());
            it_met_private = set_difference(met.begin(), met.end(), bps_prim.begin(), bps_prim.end(), bps_met_private.begin());
            bps_met_private.resize(it_met_private - bps_met_private.begin());
            bps_met_privates[k] = bps_met_private;

            if(verbose >= 1){
                cout << "\nunion of breakpoints for clone " << k << " of size " << met.size() << ":";
                for(auto bp: met){
                    cout << "\t" << bp;
                }
                cout << endl;

                cout << "metastasis-private breakpoints for clone " << k << " of size " << bps_met_private.size() << ":";
                for(auto bp: bps_met_private){
                    cout << "\t" << bp;
                }
                cout << endl;
            }
        }
    }

    // only consider subclonal events to exclude the effect of starting cell
    void collect_private_subclonal_bps(int verbose=0){
        map<int, set<string>> bps_sep;
        map<int, map<string, int>> bps_sep_count;
        map<int, vector<Sample*>> samples_sep;

        int nsample = samples.size();
        // cout << "There are " << nsample << " samples" << endl;

        // Unique breakpoints for each clone
        for(int i = 0; i < nsample; i++){
            Sample* s1 = samples[i];

            if(s1->sample_bps.size()==0){
                s1->get_breakpoints_from_pseudo2(verbose);
            }

            samples_sep[s1->clone_ID].push_back(s1);

            bps_sep[s1->clone_ID].insert(s1->sample_bps.begin(), s1->sample_bps.end());

            for(auto bp : s1->sample_bps){
                bps_sep_count[s1->clone_ID][bp]++;
            }
        }

        if(verbose >= 1){
            for(auto bps : bps_sep){
                cout << "union of breakpoints for clone " << bps.first << " of size " << bps.second.size() << ":";
                for(auto bp: bps.second){
                    cout << "\t" << bp;
                }
                cout << endl;
            }
        }

        // remove common elements from private elements
        map<int, vector<string>> bps_private_subclonal;
        for(auto c: samples_sep){
            vector<string> bps_common;  // number of common breakpoints present in all samples
            find_common_bps(bps_common, c.second, verbose);
            if(verbose >= 1) cout << "\nFind common breakpoints for clone " << c.first << " with " << c.second.size() << " samples and common breakpoint number " << bps_common.size() << endl;

            vector<string> bps_private(bps_sep[c.first].size());
            // if(bps_common.size() > 0){
                vector<string>::iterator it_private;
                it_private = set_difference(bps_sep[c.first].begin(), bps_sep[c.first].end(), bps_common.begin(), bps_common.end(), bps_private.begin());
                bps_private.resize(it_private - bps_private.begin());
            // }else{
            //     copy(bps_sep[c.first].begin(), bps_sep[c.first].end(), bps_private.begin());
            // }

            cout << bps_private.size() << endl;

            bps_private_subclonal[c.first] = bps_private;

            // if(verbose >= 1) cout << "Number of private breakpoints at different frequency cutoff" << endl;
            // double frac_delta = 1.0/c.second.size();
            // for(double frac_cutoff = frac_delta; frac_cutoff < 1; frac_cutoff += frac_delta){
            //     int num_private_shared = 0;
            //     for(auto bp: bps_private){
            //         double frac = (double) bps_sep_count[c.first][bp] / c.second.size();
            //         if(bps_sep_count[c.first][bp] > 0 && verbose >= 1) cout << bp << ": " << bps_met_sep_count[c.first][bp] << "\t" << frac << endl;
            //         if(frac <= frac_cutoff) num_private_shared++;
            //     }
            //     cout << num_private_shared << endl;
            // }
        }

        if(verbose >= 1){
            for(auto bps : bps_private_subclonal){
                cout << "private subclonal breakpoints for clone" << bps.first << " of size " << bps.second.size() << ":";
                for(auto bp: bps.second){
                    cout << "\t" << bp;
                }
                cout << endl;
            }
        }
    }


    // Find common breakpoints from all regions taken at a site
    void find_common_bps(vector<string>& bps_common, const vector<Sample*>& samples, int verbose=0) {
        vector<string>::iterator it_common;

        int nsample = samples.size();

        for(int i = 0; i < nsample; i++){
            Sample* s1 = samples[i];

            if(verbose > 1){
                string sname = "c" + to_string(s1->clone_ID) + "_s" + to_string(s1->sample_ID);
                cout << sname << endl;
                s1->print_breakpoints();
            }

            if(i==0){
                bps_common.assign(s1->sample_bps.begin(), s1->sample_bps.end());
            }
            it_common = set_intersection(s1->sample_bps.begin(), s1->sample_bps.end(), bps_common.begin(), bps_common.end(), bps_common.begin());
        }
        bps_common.resize(it_common - bps_common.begin());
        if(verbose > 1){
            cout << "Common breakpoints:";
            for(auto bp : bps_common){
                cout << "\t" << bp;
            }
            cout << "\n";
        }
    }


    // Print number of all/primary/metastasis common breakpoints (early CNVs)
    void print_breakpoints_stat_common(int verbose=0) {
        // number of common breakpoints present in all samples (not useful when comparing to real data)
        vector<string> bps_common;
        find_common_bps(bps_common, samples, verbose);
        cout << "\n" << bps_common.size() << "\n";

        map<int, vector<Sample*>> samples_sep;
        // cout << "There are " << nsample << " samples" << endl;
        for(int i = 0; i < samples.size(); i++){
            Sample* s1 = samples[i];
            samples_sep[s1->clone_ID].push_back(s1);
        }
        for(auto c: samples_sep){
            if(verbose>1) cout << "Find common breakpoints for clone " << c.first << endl;
            vector<string> bps_common;  // number of common breakpoints present in all samples
            find_common_bps(bps_common, c.second, verbose);
            cout << bps_common.size() << "\n";
        }
    }


    // Only count X-private breakpoints present in certain fraction of samples
    void print_breakpoints_stat_private2(const vector<int>& nsamples, const vector<string>& bps_prim_private, map<int, vector<string>>& bps_met_privates, double frac_cutoff = 0.5, int loc_type=BIN, int verbose=0) {
        int num_prim_private_shared = 0;
        for(auto bp: bps_prim_private){
            double frac = (double) bps_prim_count[bp] / nsamples[0];
            if(bps_prim_count[bp] > 0 && verbose >= 1) cout << bp << ": " << bps_prim_count[bp] << "\t" << frac << endl;
            if(frac <= frac_cutoff) num_prim_private_shared++;
        }
        // << "\n"  << bps_prim_private.size()
        // cout << num_prim_private_shared;
         // / bps_prim.size()
        double frac_prim = (double)num_prim_private_shared;
        cout << frac_prim;

        // for metastasis samples, excluding clonal events which are strongly affected by starting cell
        if(frac_cutoff==1) return;
        assert(id_mets.size() == nsamples.size() - 1);
        for(auto k : id_mets){
            int num_met_private_shared = 0;
            vector<string> bps_met_private = bps_met_privates[k];
            for(auto bp: bps_met_private){
                double frac = (double) bps_met_sep_count[k][bp] / nsamples[k];
                if(bps_met_sep_count[k][bp] > 0 && verbose >= 1) cout << bp << ": " << bps_met_sep_count[k][bp] << "\t" << frac << endl;
                if(frac <= frac_cutoff) num_met_private_shared++;
            }
            // << "\n" << bps_met_private.size()
            // cout << "\n" << num_met_private_shared;
            //  / bps_met_sep[k].size()
            double frac_met = (double)num_met_private_shared;
            cout << "\n" << frac_met;
        }
        cout << endl;
    }

    // Print proportion of altered bins
    void print_bin_stat(int verbose=0){
        int nsample = samples.size();

        // cout << "There are " << nsample << " samples" << endl;
        double avg_nalter = 0;
        map<int, double> avg_nalter_sep;

        int alter_indicator[NUM_LOC]={0};
        map<int, int*> alter_indicator_sep;

        for(auto c : this->clones){
            alter_indicator_sep[c->clone_ID] = new int[NUM_LOC];
            avg_nalter_sep[c->clone_ID] = 0;
        }

        for(int i = 0; i < nsample; i++){
            Sample* s1 = samples[i];
            for(int i = 0; i < NUM_LOC; i++){
                if(abs(s1->sample_loc_change[i]) > BIN_CUOFF){
                    alter_indicator_sep[s1->clone_ID][i] = 1;
                    alter_indicator[i] = 1;
                }else{
                    alter_indicator_sep[s1->clone_ID][i] = 0;
                }
            }
        }

        for(int i = 0; i < NUM_LOC; i++){
            if(alter_indicator[i] == 1) avg_nalter++;
        }

        for(auto cn: alter_indicator_sep){
            for(int i = 0; i < NUM_LOC; i++){
                if(cn.second[i] == 1) avg_nalter_sep[cn.first]++;
            }
        }

        avg_nalter = avg_nalter / NUM_LOC;
        if(verbose > 0) cout << "Average number of altered bins" << endl;
        cout << avg_nalter << endl;;

        for(auto na : avg_nalter_sep){
            if(verbose > 0) cout << "Average number of altered bins for clone " << na.first << endl;
            na.second = na.second / NUM_LOC;
            cout << na.second << endl;;
        }
    }


    // Print proportion of altered bins (subclonal)
    void print_bin_subclonal_stat(int verbose=0){
        int nsample = samples.size();
        // cout << "There are " << nsample << " samples" << endl;
        map<int, double> avg_nalter_sep;
        map<int, int*> alter_indicator_sep;     // count number of times a bin is altered in a clone
        map<int, vector<Sample*>> samples_sep;

        for(auto c : this->clones){
            alter_indicator_sep[c->clone_ID] = new int[NUM_LOC];
            for(int i = 0; i < NUM_LOC; i++){
                alter_indicator_sep[c->clone_ID][i] = 0;
            }
        }

        for(int i = 0; i < nsample; i++){
            Sample* s1 = samples[i];
            samples_sep[s1->clone_ID].push_back(s1);

            for(int i = 0; i < NUM_LOC; i++){
                if(fabs(s1->sample_loc_change[i]) >= BIN_CUOFF){
                    alter_indicator_sep[s1->clone_ID][i] += 1;
                }
            }
        }

        // exclude clonal regions
        for(auto cn: alter_indicator_sep){
            avg_nalter_sep[cn.first] = 0;
            for(int i = 0; i < NUM_LOC; i++){
                // cout << cn.first << "\t" << cn.second[i] << endl;
                if(cn.second[i] > 0 && cn.second[i] < samples_sep[cn.first].size())
                    avg_nalter_sep[cn.first]++;
            }
        }

        // normalized by NUM_LOC to account for different choices of bin size (in real data)
        for(auto na : avg_nalter_sep){
            if(verbose > 0) cout << "Average number of altered bins for clone " << na.first << endl;
            na.second = na.second * 100 / NUM_LOC;
            cout << na.second << endl;
        }
    }

    // distinguish gain/loss to account for different size distribution
    void print_bin_subclonal_stat_by_type(int verbose=0){
        int nsample = samples.size();
        // cout << "There are " << nsample << " samples" << endl;
        map<int, double> avg_nalter_sep_gain;
        map<int, double> avg_nalter_sep_loss;
        map<int, int*> alter_indicator_sep_gain;     // count number of times a bin is altered in a clone
        map<int, int*> alter_indicator_sep_loss;
        map<int, vector<Sample*>> samples_sep;

        for(auto c : this->clones){
            alter_indicator_sep_gain[c->clone_ID] = new int[NUM_LOC];
            alter_indicator_sep_loss[c->clone_ID] = new int[NUM_LOC];
            for(int i = 0; i < NUM_LOC; i++){
                alter_indicator_sep_gain[c->clone_ID][i] = 0;
                alter_indicator_sep_loss[c->clone_ID][i] = 0;
            }
        }

        for(int i = 0; i < nsample; i++){
            Sample* s1 = samples[i];
            samples_sep[s1->clone_ID].push_back(s1);

            for(int i = 0; i < NUM_LOC; i++){
                // cout << s1->sample_loc_change[i] << endl;
                if(s1->sample_loc_change[i] >= BIN_CUOFF){
                    alter_indicator_sep_gain[s1->clone_ID][i] += 1;
                }
                else if(s1->sample_loc_change[i] <= -BIN_CUOFF){
                    alter_indicator_sep_loss[s1->clone_ID][i] += 1;
                }else{

                }
            }
        }

        // exclude clonal regions
        for(auto cn: alter_indicator_sep_gain){
            avg_nalter_sep_gain[cn.first] = 0;
            for(int i = 0; i < NUM_LOC; i++){
                // cout << cn.first << "\t" << cn.second[i] << endl;
                if(cn.second[i] > 0 && cn.second[i] < samples_sep[cn.first].size())
                    avg_nalter_sep_gain[cn.first]++;
            }
        }

        for(auto cn: alter_indicator_sep_loss){
            avg_nalter_sep_loss[cn.first] = 0;
            for(int i = 0; i < NUM_LOC; i++){
                // cout << cn.first << "\t" << cn.second[i] << endl;
                if(cn.second[i] > 0 && cn.second[i] < samples_sep[cn.first].size())
                    avg_nalter_sep_loss[cn.first]++;
            }
        }

        // normalized by NUM_LOC to account for different choices of bin size (in real data)
        for(auto na : avg_nalter_sep_gain){
            if(verbose > 0) cout << "Average number of gained bins for clone " << na.first << endl;
            na.second = na.second * 100 / NUM_LOC;
            cout << na.second << endl;
        }
        for(auto na : avg_nalter_sep_loss){
            if(verbose > 0) cout << "Average number of lost bins for clone " << na.first << endl;
            na.second = na.second * 100 / NUM_LOC;
            cout << na.second << endl;
        }
    }


    void print_breakpoints_stat(const vector<int>& nsamples, int loc_type=BIN, int verbose=0){
        int nsample = samples.size();
        // cout << "There are " << nsample << " samples" << endl;
        for(int i = 0; i < nsample; i++){
            Sample* s1 = samples[i];
            if(s1->sample_bps.size()==0){
                s1->get_breakpoints_from_pseudo2(verbose);
            }
        }

        vector<string> bps_prim_private(NUM_LOC);
        map<int, vector<string>> bps_met_privates;
        collect_bps(bps_prim_private, bps_met_privates, verbose);

        // if(verbose > 1) cout << "\nNumber of unique breakpoints" << endl;
        // cout << bps_all.size() << endl;
        // cout << bps_prim.size();
        // for(auto k : id_mets){
        //     cout << "\n" << bps_met_sep[k].size();
        // }

        // if(verbose > 1) cout << "\nNumber of common breakpoints" << endl;
        // print_breakpoints_stat_common(verbose);
        vector<string> bps_common;
        find_common_bps(bps_common, samples, verbose);
        cout << "\n" << bps_common.size() << "\n";

        map<int, vector<Sample*>> samples_sep;
        // cout << "There are " << nsample << " samples" << endl;
        for(int i = 0; i < samples.size(); i++){
            Sample* s1 = samples[i];
            samples_sep[s1->clone_ID].push_back(s1);
        }
        for(auto c: samples_sep){
            if(verbose>1) cout << "Find common breakpoints for clone " << c.first << endl;
            vector<string> bps_common;  // number of common breakpoints present in all samples
            find_common_bps(bps_common, c.second, verbose);
            cout << bps_common.size() << "\n";
        }

        // get clone-private subclonal events


        if(verbose > 1) cout << "\nNumber of private breakpoints at different frequency cutoff" << endl;
        int max_nsample = *max_element(nsamples.begin(), nsamples.end());
        double frac_delta = 1.0/max_nsample;
        for(double frac_cutoff = frac_delta; frac_cutoff <= 1; frac_cutoff += frac_delta){
            if(verbose > 1) cout << "freq cutoff " << frac_cutoff << endl;
            print_breakpoints_stat_private2(nsamples, bps_prim_private, bps_met_privates, frac_cutoff, loc_type, verbose);
        }
    }


    // variance of CNs at each bin across all samples
    void print_variance(const vector<Sample*>& samples, int loc_type=BIN){
        // CNs at each bin across all samples
        map<int, VAR> bin_pcn;
        // map<int, NBIN> bin_cout;

        int nsample = samples.size();
        // cout << "There are " << nsample << " samples" << endl;

        for(int i = 0; i < nsample; i++){
            Sample* s1 = samples[i];
            // cout << s1->sample_ID << "," << s1->clone_ID << "\n";
            for(int i = 0; i < NUM_LOC; i++){
                bin_pcn[i](s1->sample_loc_change[i]);
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
        cout << sum_var << '\n';
    }

    // variance of CNs at each bin across all samples taken at one site
    void print_variance_by_site(int loc_type=BIN){
        for(auto c : clones){
            print_variance(c->samples, loc_type);
        }
    }


    // Print CNP in a single file
    void print_sample_cnp(string outdir, string suffix, int verbose = 0){
        int i = 0;
        string fname = outdir + "bulk_cn"  + suffix + ".txt";
        ofstream fcn;
        fcn.open(fname, std::ofstream::trunc | std::ofstream::out);

        for(auto s : samples){
            i++;
            if(verbose > 1) cout << "Print CNP of sample " << i << endl;

            int clone_ID = s->clone_ID;
            int sample_ID = s->sample_ID;
            if(s->sample_pcn.size()==0){
                loc2pcn(s->sample_loc_change, s->sample_pcn, verbose);
            }
            dpcn cnp = s->sample_pcn;
            // string fname = outdir + "bulk_cn"  + suffix  + "_clone" + to_string(clone_ID) + "_sample" + to_string(sample_ID) + ".txt";
            string sname = "c" + to_string(clone_ID) + "_s" + to_string(sample_ID);
            print_bulk_cn(sname, cnp, fcn, verbose);
        }

        fcn.close();
    }


    // Print CNP in a single file
    void print_sample_cnp_orig(string outdir, string suffix, int verbose = 0){
        int i = 0;
        string fname = outdir + "bulk_cn_orig"  + suffix + ".txt";
        ofstream fcn_orig;
        fcn_orig.open(fname, std::ofstream::trunc | std::ofstream::out);

        for(auto s : samples){
            i++;
            int clone_ID = s->clone_ID;
            int sample_ID = s->sample_ID;
            dpcn cnp = s->sample_pcn;

            if(verbose > 1) cout << "Print CNP of sample " << i << endl;
            // fname = outdir + "bulk_cn_orig"  + suffix  + "clone" + to_string(clone_ID) + "_sample" + to_string(sample_ID) + ".txt";
            string sname = "c" + to_string(clone_ID) + "_s" + to_string(sample_ID);
            print_bulk_cn_orig(sname, cnp, fcn_orig, verbose);
        }

        fcn_orig.close();
    }


    /*
     print out the lineages of multiple clones (for Trevor's data), sample tree
    */
    void print_lineage(){

    }
};


#endif
