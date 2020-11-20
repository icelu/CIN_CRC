#include <cstdlib>
#include <sstream>
#include <string>
#include <climits>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "gland.hpp"


using namespace std;


// Simulate bulk CNPs of individual glands taken from one patient


int main(int argc, char const *argv[]) {
    int mode;

    int ndeme;
    int max_deme_size;
    string fdeme;

    string fmut;

    string fgenotype;     // starting genotype
    string flprob;     // location mutation probability
    string fgain;
    string floss;

    double min_freq;
    double max_freq;
    double delta;
    int use_std;

    double leap_size;

    string num_glands; // number of samples to take at each side

    int Nend;
    double birth_rate, death_rate;
    double mutation_rate;
    // parameters for mean of dup/del size distributions
    // int mean_gain_size, mean_loss_size;
    int loc_type;

    int model_ID;   // not making differences, fitness is used to introduce selection
    int use_alpha;  // when alpha is used, genotype differences are considered
    double fitness;
    int genotype_diff;
    int growth_type;
    int norm_by_bin;

    int stat_type;
    double frac_cutoff;
    // int num_chr;

    string outdir, suffix; // output

    unsigned long seed;

    int verbose;

    namespace po = boost::program_options;

    po::options_description generic("Generic options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ;

    po::options_description required("Required parameters");
    required.add_options()
      ("odir,o", po::value<string>(&outdir)->required()->default_value("./"), "output directory")
       ;

    po::options_description optional("Optional parameters");
    optional.add_options()
      ("mode", po::value<int>(&mode)->default_value(0), "mode of simulation. 0: gland as cell; 1: gland fission")

      ("ndeme", po::value<int>(&ndeme)->default_value(0), "number of demes in the final tumor")
      ("max_deme_size", po::value<int>(&max_deme_size)->default_value(10000), "maximum number of cells in a deme")
      ("num_glands", po::value<string>(&num_glands)->default_value("2 1 1"), "the number of sides followed by the number of glands taken at different sides of a crypt, separated by space")
      ("fdeme", po::value<string>(&fdeme)->default_value(""), "file with deme (gland) relationships")
      ("fmut", po::value<string>(&fmut)->default_value(""), "file with mutation counting informaton")
      ("fgenotype", po::value<string>(&fgenotype)->default_value(""), "TSV file with starting genotype of the first cell")
      ("flprob", po::value<string>(&flprob)->default_value(""), "TSV file with probability of mutation at each location")
      ("fgain", po::value<string>(&fgain)->default_value(""), "TSV file with size of all copy number gains in the real data")
      ("floss", po::value<string>(&floss)->default_value(""), "TSV file with size of all copy number losses in the real data")

      ("stat_type", po::value<int>(&stat_type)->default_value(0), "type of summary statistics. 0: variance; 1: clone-pairwise differences; 2: average CNP; 3: complete CNP; 4: sample-pairwise differences")
      ("loc_type", po::value<int>(&loc_type)->default_value(0), "type of CNA unit. 0: BIN-level at one chromsome; 1: arm-level; 2: BIN-level across multiple chromsomes")

      ("min_freq", po::value<double>(&min_freq)->default_value(0.0), "minimal mutation frequency to consider")
      ("max_freq", po::value<double>(&max_freq)->default_value(0.25), "maximal mutation frequency to consider")
      ("delta", po::value<double>(&delta)->default_value(0.025), "step size of mutation frequency histogram")

      ("use_std", po::value<int>(&use_std)->default_value(0), "whether or not to use standard deviation")

      ("leap_size", po::value<double>(&leap_size)->default_value(0.0), "step size of tau-leaping to accelerate the simulations")
      ("frac_cutoff", po::value<double>(&frac_cutoff)->default_value(0.5), "cutoff of counting breakpoints")

      ("birth_rate", po::value<double>(&birth_rate)->default_value(1), "birth rate")
      ("death_rate", po::value<double>(&death_rate)->default_value(0), "death rate")
      ("Nend,e", po::value<int>(&Nend)->default_value(100), "size of final cell populations")
      ("mutation_rate", po::value<double>(&mutation_rate)->default_value(0), "mutation rate of CNAs")
      ("mean_gain_size", po::value<int>(&MEAN_GAIN_SIZE)->default_value(0), "mean size of segment gain (in terms of bins)")
      ("mean_loss_size", po::value<int>(&MEAN_LOSS_SIZE)->default_value(0), "mean size of segment loss (in terms of bins)")
      // ("num_chr", po::value<int>(&num_chr)->default_value(22), "number of chromosomes to consider")

      ("bp_cutoff", po::value<double>(&BP_CUTOFF)->default_value(0.1), "threshold to determine whether a breakpoint can be detected or not")
      ("bin_cutoff", po::value<double>(&BIN_CUOFF)->default_value(0.1), "threshold to determine whether a bin can be detected as altered or not")

      ("model", po::value<int>(&model_ID)->default_value(0), "model of evolution. 0: neutral; 1: selection")
      ("use_alpha", po::value<int>(&use_alpha)->default_value(1), "whether or not to use alpha in selection model. 0: use selection coefficient; 1: use alpha")
      ("fitness", po::value<double>(&fitness)->default_value(0), "fitness values of mutatants")
      ("genotype_diff", po::value<int>(&genotype_diff)->default_value(3), "type of genotype difference (L1 distance) in simulating selection. 0: no; 1: L1 distance of CN; 2: Hamming distance of CN; 3: number of new mutations")
      ("norm_by_bin", po::value<int>(&norm_by_bin)->default_value(0), "whether or not to normalize genotype difference by number of bins (segments) in the genome. 0: no; 1: yes")
      ("growth_type", po::value<int>(&growth_type)->default_value(0), "Type of growth. 0: only birth; 1: change birth rate; 2: change death rate; 3: change both birth and death rate")

      ("suffix,s", po::value<string>(&suffix)->default_value(""), "suffix of output file")

      ("seed", po::value<unsigned long>(&seed)->default_value(0), "seed used for generating random numbers")
      ("verbose", po::value<int>(&verbose)->default_value(0), "verbose level (0: default, 1: print information of final cells; 2: print information of all cells)")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(required).add(optional);
    po::variables_map vm;

    try {
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
        if(vm.count("help")){
            cout << cmdline_options << endl;
            return 1;
        }
        if(vm.count("version")){
            cout << "sim_gland [version 0.1], a program to simulate copy number variations along a stochastic branching tree" << endl;
            return 1;
        }
        po::notify(vm);
    }
    catch (const exception& e) {
          cerr << e.what() << endl;
          return 1;
    }

    unsigned long rseed = setup_rng(seed);

    if(verbose > 0){
        cout << "Random seed: " << rseed << endl;
        cout << "Simulating bulk WGS and WES samples from glands sampled from a patient crypt" << endl;
        cout << "Fitness:\t" << fitness << endl;
        cout << "Mutation rate:\t" << mutation_rate << endl;
    }

    // cout << "Using Boost "
    //   << BOOST_VERSION / 100000     << "."  // major version
    //   << BOOST_VERSION / 100 % 1000 << "."  // minor version
    //   << BOOST_VERSION % 100                // patch level
    //   << std::endl;

    string size;

    ifstream fcn_gain(fgain);
    while(getline(fcn_gain, size)){
        real_gain_sizes.push_back(atof(size.c_str()));
    }
    // cout << real_gain_sizes.size() << endl;

    ifstream fcn_loss(floss);
    while(getline(fcn_loss, size)){
        real_loss_sizes.push_back(atof(size.c_str()));
    }
    // cout << real_loss_sizes.size() << endl;


    if(fgenotype != ""){
        if(verbose > 0)
            cout << "Assigning a genotype to the starting cell" << endl;
        ifstream infile(fgenotype);

        // cout << "Skip header line" << endl;
        string line;
        getline(infile, line); // skip the first line
        // cout << line << endl;

        while(getline(infile, line))
        {
          istringstream iss(line);
          string token;
          vector<int> vals;
          while (std::getline(iss, token, '\t'))
          {
              // process each token
              // std::cout << token << " ";
              vals.push_back(atoi(token.c_str()));
          }
          START_GENOTYPE[vals[0] - 1] = vals[1];
          // std::cout << std::endl;
        }
        infile.close();
    }
    // for(auto cn : start_genotype){
    //     cout << cn << endl;
    // }

    // cout << "LOC_PROBS before" << endl;
    // for(int i = 0; i < NUM_LOC; i++){
    //     cout << LOC_PROBS[i] << endl;
    // }
    if(flprob != ""){
        if(verbose > 0)
            cout << "Assigning a mutation probability to each location" << endl;
        ifstream infile(flprob);

        // cout << "Skip header line" << endl;
        string line;
        getline(infile, line); // skip the first line
        // cout << line << endl;

        while(getline(infile, line))
        {
          istringstream iss(line);
          string token;
          vector<int> vals;
          while (std::getline(iss, token, '\t'))
          {
              // process each token
              // std::cout << token << " ";
              vals.push_back(atof(token.c_str()));
          }
          LOC_PROBS[vals[0] - 1] = (double) vals[1] / NUM_LOC;
          // std::cout << std::endl;
        }
        infile.close();
    }
    // cout << "LOC_PROBS after" << endl;
    // for(int i = 0; i < NUM_LOC; i++){
    //     cout << LOC_PROBS[i] << endl;
    // }

    Glands glands;
    Cell_ptr start_cell = new Cell(1, 0, birth_rate, death_rate, mutation_rate, 0);
    Model start_model(model_ID, genotype_diff, growth_type, fitness, use_alpha, norm_by_bin);
    // if(verbose > 0) cout << "Simulating deme partition history only" << endl;
    // tumor.simulate_deme_partition(start_cell, ndeme, MAX_DEME_SIZE, fdeme, verbose);

    int num_sampled_glands = 0;
    vector<int> nglands;    // only do two-side sampling if nglands is not empty
    if(num_glands != "")
    {
        get_vals_from_str(nglands, num_glands, 0);
        sort(nglands.begin(), nglands.end());
    }

    for(int i = 0; i < nglands.size(); i++){
        num_sampled_glands += nglands[i];
    }

    if(verbose > 0){
        for(int i = 0; i < nglands.size(); i++){
            cout << "number of glands at side " << i + 1 << " is " << nglands[i] << endl;
        }
        cout << "total number of glands to sample is " << num_sampled_glands << endl;
    }

    // If the total number of glands is not specified, it is the sum of given glands
    if(ndeme == 0){
        for(int i = 0; i < nglands.size(); i++){
            ndeme += nglands[i];
        }
    }
    assert(ndeme > 0);

    if(verbose > 0){
        cout << "\nSimulating gland growth with (allele-specific) CNAs, with " <<  ndeme << " glands in total" << endl;
    }

    set_outdir(outdir, verbose);
    if(suffix == ""){
        string sep = "-";
        suffix = sep + to_string(ndeme) + sep + to_string(int(birth_rate*10)) + sep + to_string(int(death_rate*10)) + sep + to_string(int(mutation_rate)) + sep + to_string(int(rseed));
    }

    if(mode == 1){
        if(verbose > 0) cout << "Simulating gland fission" << endl;

        int store_lineage = 0;
        vector<string> lineages;
        if(fdeme != ""){
            store_lineage = 1;
        }

        glands.simulate_gland_growth(start_cell, ndeme, max_deme_size, start_model, lineages, store_lineage, loc_type, leap_size, verbose);

        if(verbose > 0 && fdeme != ""){
            cout << "Print gland lineages" << endl;
            string header = "Parent\tID\tNcell\tNmut\tNmut_parent";
            glands.print_gland_lineage(fdeme, header, lineages);
        }

        vector<int> gids;
        // Get average CNP of all cells
        map<int, double*> avg_loc_changes;
        // Initialize avg_loc_changes to be global
        for(auto s : glands.clones){
            gids.push_back(s->clone_ID);
            avg_loc_changes[s->clone_ID] = new double[NUM_LOC];
            for(int i = 0; i < NUM_LOC; i++){
                avg_loc_changes[s->clone_ID][i] = 0;
            }
            s->set_bulk_cnp_from_pseudo(avg_loc_changes[s->clone_ID], s->curr_cells);
        }
        sort(gids.begin(), gids.end());
        glands.print_sstat(gids, avg_loc_changes, stat_type, loc_type, min_freq, max_freq, delta, use_std, verbose);

        if(verbose > 0){
          glands.print_sample_cnp(avg_loc_changes, outdir, suffix, verbose);
        }
        if(fmut != ""){
          ofstream fcf(fmut);
          for(auto g : glands.clones){
            for(auto cf : g->muts_freq){
              // Compute frequency of mutation in each gland
              double freq = g->muts_freq[cf.first] / (g->curr_cells.size() - 1);
              // g->muts_freq[cf.first] = freq;
              if(freq >= 0.5) glands.muts.insert(cf.first);
              // fcf << g->clone_ID << "\t" << cf.first << "\t" << cf.second << "\t" << g->curr_cells.size() << "\t" << freq << endl;
            }
          }
          // fcf << "Total number of unique (clonal) mutations in all glands is: " << glands.muts.size() << endl;
          double fission_rate = (double) glands.muts.size() / (glands.clones.size() - 1);
          fcf << "0\t0\t" << glands.muts.size() << "\t" << glands.clones.size() << "\t" << fission_rate << endl;
        }
    }else{
        if(verbose > 0) cout << "Simulating gland as cell" << endl;
        int store_lineage = 0;
        vector<string> lineages;
        if(fdeme != ""){
            store_lineage = 1;
        }

        glands.root = new node(start_cell->cell_ID);
        glands.simulate_gland_as_cell(start_cell, ndeme, start_model, lineages, store_lineage, loc_type, leap_size, verbose);

        if(nglands.size() > 0){
            glands.sample_gland_from_cell(nglands, verbose);
            // assert(glands.sample_IDs.size() == num_sampled_glands);
            if(verbose > 0) cout << "Sampling " << glands.sample_IDs.size() << " glands from different sides" << endl;
        }

        if(verbose > 0 && fdeme != ""){
            string header = "Parent\tID\tTime\tNmut";
            glands.print_gland_lineage(fdeme, header, lineages);
        }

        vector<int> gids;
        // Get average absolute CNP of all cells
        map<int, double*> avg_loc_changes;
        // assert(glands.clones[0]->curr_cells.size() == ndeme);
        // Initialize avg_loc_changes to be global
        vector<int> curr_cIDs;
        for(auto c : glands.clones[0]->curr_cells){
            // cout << c->cell_ID << endl;
            curr_cIDs.push_back(c->cell_ID);
            if(glands.sample_IDs.size() > 0 && find(glands.sample_IDs.begin(), glands.sample_IDs.end(), c->cell_ID) == glands.sample_IDs.end()) continue;
            // cout << " adding into list" << endl;
            gids.push_back(c->cell_ID);
            avg_loc_changes[c->cell_ID] = new double[NUM_LOC];
            for(int i = 0; i < NUM_LOC; i++){
                avg_loc_changes[c->cell_ID][i] = c->loc_changes[i];
            }
        }
        // cout << curr_cIDs.size() << " cells in the end" << endl;
        // for(auto sid : glands.sample_IDs){
        //   if(find(curr_cIDs.begin(), curr_cIDs.end(), sid) == curr_cIDs.end()){
        //     cout << "Weird! cell not in curr_cells " << sid << endl;
        //   }
        // }
        // cout << "glands used for summary statistics " << gids.size() << endl;
        // assert(gids.size() == num_sampled_glands);

        sort(gids.begin(), gids.end());
        glands.print_sstat(gids, avg_loc_changes, stat_type, loc_type, min_freq, max_freq, delta, use_std, verbose);

        if(verbose > 0) glands.print_cell_cnp(glands.clones[0], outdir, suffix, verbose);

        // output mutations in each gland (cell)
        if(fmut != ""){
          ofstream fcm(fmut);
          for(auto c : glands.clones[0]->curr_cells){
            int is_sampled = 1;
            if(glands.sample_IDs.size() > 0 && find(glands.sample_IDs.begin(), glands.sample_IDs.end(), c->cell_ID) == glands.sample_IDs.end()) is_sampled = 0;
            for(auto mid : c->muts){
              fcm << c->cell_ID << "\t" << mid << "\t" << is_sampled << endl;
            }
          }
        }
    }

    return 0;
}
