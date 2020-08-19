#include <cstdlib>
#include <sstream>
#include <string>
#include <climits>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "tumor.hpp"


using namespace std;

const int MIN_NCELL_MET = 10;

// TODO: pass parameters
int main(int argc, char const *argv[]) {
    int mode;

    int ndeme;
    string fdeme;

    double leap_size;

    int num_migration;
    string time_migration_vals;
    string size_migration_vals;
    string tend_migration_vals;

    string nsample_vals;    // number of samples to take at each clone

    int Nend;
    double birth_rate, death_rate;
    double mutation_rate;
    // parameters for mean of dup/del size distributions
    // int mean_gain_size, mean_loss_size;
    int loc_type;

    int model_ID;
    double fitness;
    int genotype_diff;
    int growth_type;

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
      ("mode", po::value<int>(&mode)->default_value(0), "mode of simulation. 0: single clone; 1: multi samples")

      ("ndeme", po::value<int>(&ndeme)->default_value(0), "number of demes in the final tumor")
      ("fdeme", po::value<string>(&fdeme)->default_value(""), "file with deme (gland) relationships")

      ("nsample", po::value<string>(&nsample_vals)->default_value(""), "number of samples to take at each clone (NSAMPLE1 ... NSAMPLEn)")

      ("stat_type", po::value<int>(&stat_type)->default_value(0), "type of summary statistics. 0: variance; 1: clone-pairwise differences; 2: average CNP; 3: complete CNP; 4: sample-pairwise differences")
      ("loc_type", po::value<int>(&loc_type)->default_value(0), "type of CNA unit. 0: BIN-level at one chromsome; 1: arm-level; 2: BIN-level across multiple chromsomes")
      // ("use_diff", po::value<int>(&use_diff)->default_value(0), "using  as summary statistics")
      // ("use_avg", po::value<int>(&use_avg)->default_value(0), "using  as summary statistics")
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

      // parameters related to metastasis clones
      ("num_migration", po::value<int>(&num_migration)->default_value(0), "number of metastasis/migration clones")
      ("time_migration", po::value<string>(&time_migration_vals)->default_value(""), "time of metastasis/migration clones appearing in term of cell number in primary clones")
      ("size_migration", po::value<string>(&size_migration_vals)->default_value(""), "number of cells in metastasis clones")
      ("tend_migration", po::value<string>(&tend_migration_vals)->default_value(""), "time when metastasis clones stop growing (when samples were taken)")

      ("model", po::value<int>(&model_ID)->default_value(0), "model of evolution. 0: neutral; 1: gradual; 2: punctuated; 3: positive")
      ("fitness", po::value<double>(&fitness)->default_value(0), "fitness values of mutatants")
      ("genotype_diff", po::value<int>(&genotype_diff)->default_value(0), "whether or not to use genotype difference (L1 distance) in simulating selection. 0: no; 1: yes")
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
            cout << "sim_crc [version 0.1], a program to simulate copy number variations in primary-metastasis multi-region samples" << endl;
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
        // cout << "Loc probability:";
        // for(int i = 0; i < NUM_LOC; i++){
        //     cout<< "\t" << LOC_PROBS[i];
        // }
        // cout << endl;
        // cout << "Chr probability:";
        // for(int i = 0; i < NUM_CHR; i++){
        //     cout<< "\t" << CHR_PROBS[i];
        // }
        // cout << endl;
    }

    set_outdir(outdir, verbose);
    if(suffix==""){
        string sep = "-";
        suffix = sep + to_string(Nend) + sep + to_string(int(birth_rate*10)) + sep + to_string(int(death_rate*10)) + sep + to_string(int(mutation_rate)) + sep + to_string(int(rseed));
    }

    Tumor tumor;
    Cell_ptr start_cell = new Cell(1, 0, birth_rate, death_rate, mutation_rate, 0);

    if(ndeme > 0){
        if(verbose > 0) cout << "Simulating deme partition history" << endl;
        tumor.simulate_deme_partition(start_cell, ndeme, MAX_DEME_SIZE, fdeme, verbose);
        return 0;
    }

    Model start_model(model_ID, genotype_diff, growth_type, fitness);

    if(verbose > 0) cout << "Simulating multi-region samples from a growing tumor with metastasis clones" << endl;

    // time from early to late
    vector<int> time_migration;
    if(time_migration_vals != ""){
        get_vals_from_str(time_migration, time_migration_vals, num_migration);
    }else{
        // Randomly pick a metastasis time
        for(int c = 1; c <= num_migration; c++){
            int n = runiform(r, MIN_NCELL_MET, Nend/2);
            time_migration.push_back(n);
        }
        sort(time_migration.begin(), time_migration.end());
    }

    cout << "Fitness:\t" << fitness << endl;
    cout << "Mutation rate:\t" << mutation_rate << endl;

    cout << "Time of migration:";
    for(int c = 0; c < num_migration; c++){
        cout << "\t" << time_migration[c];
    }
    cout << endl;

    vector<double> tend_migration;
    if(tend_migration_vals != ""){
        get_vals_from_str(tend_migration, tend_migration_vals, num_migration);
    }else{
        // same time as when primary clones stop growing
        for(int c = 1; c <= num_migration; c++){
            tend_migration.push_back(0);
        }
    }

    vector<Model> model_migration;
    for(int i = 0; i < num_migration; i++){
        Model model(model_ID, genotype_diff, growth_type, fitness);
        model_migration.push_back(model);
    }
    if(verbose > 0) cout << "\nSimulating metastasis clones" << endl;
    tumor.simulate_metastasis(num_migration, time_migration, tend_migration, model_migration, start_cell, start_model, Nend, loc_type, leap_size, verbose);


    if(nsample_vals=="") return 0;

    if(verbose > 0) cout << "\nSimulating multi-region samples" << endl;
    vector<int> num_samples;
    get_vals_from_str(num_samples, nsample_vals, tumor.clones.size());
    int num_tot_samples = 0;
    for(int i = 0; i < num_samples.size(); i++){
        num_tot_samples += num_samples[i];
    }

    tumor.simulate_multiple_samples(tumor.clones, num_samples, outdir, suffix, frac_cutoff, stat_type, loc_type, verbose);

    return 0;
}
