#include <cstdlib>
#include <sstream>
#include <string>
#include <climits>
#include <time.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "gland.hpp"


using namespace std;


// Simulate bulk CNPs of individual glands taken from one patient


// Start simulations from diploid, to show negative selection suppress karyotype diversity when close to fitness peak (optimal karyotype)
void simulate_from_diploid(unsigned long seed, Cell_ptr start_cell, Glands glands, Model start_model, int ndeme, double WEIGHT_OPTIMUM,
  int bottleneck, int time_bottleneck, int track_lineage, int loc_type, double leap_size, int multiple_output,
  string fdeme, string fgenotype_chr, string outdir, string suffix, string fstat, string fdiv, string ffit, clock_t tStart, int verbose = 0){
  // simulate cell (gland) growth, allowing output at different population sizes and metastasis from some cells
  if(verbose > 0){
    cout << "\nSimulating gland time-series growth" << endl;
    if(CHR_CNA == 1){
      cout << "\tSimulating chr-level CNAs for quick convergence" << endl;
    }
  }
  int store_lineage = 0;
  vector<string> lineages;
  if(fdeme != ""){
      store_lineage = 1;
  }

  if(fgenotype_chr != ""){
      ifstream infile(fgenotype_chr);

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
        OPT_KARYOTYPE_CHR[vals[0] - 1] = vals[1];
        // std::cout << std::endl;
      }
      infile.close();

      // Increase the mutation likelihood of CNAs at optimum karyotype
      for(int i = 0; i < NUM_CHR; i++){
        if(OPT_KARYOTYPE_CHR[i] != NORM_PLOIDY){
          CHR_PROBS[i] = CHR_PROBS[i] * WEIGHT_OPTIMUM;
        }
      }

      if(verbose > 1){
        cout << "probability of chr-level CNAs:";
        for(int i = 0; i < NUM_CHR; i++){
          cout << " " << CHR_PROBS[i];
        }
        cout << endl;

        cout << "start karyotype: ";
        for(int i = 0; i < NUM_CHR; i++){
          cout << START_KARYOTYPE_CHR[i] << ", ";
        }
        cout << endl;

        // output optimal karyotype
        cout << "optimal karyotype: ";
        for(int i = 0; i < NUM_CHR; i++){
          cout << OPT_KARYOTYPE_CHR[i] << ", ";
        }
        cout << endl;
      }
  }

  glands.root = new node(start_cell->cell_ID);
  if(bottleneck > 0){
    cout << "Simulate cell growth with " << bottleneck <<  "-cell bottleneck" << endl;
    glands.simulate_gland_as_cell(start_cell, time_bottleneck, start_model, lineages, store_lineage, loc_type, leap_size, 0, 0, verbose);
    cout << "Grow population until " << glands.clones[0]->curr_cells.size() << endl;
    // sample bottleneck cells from current population
    // cout << "Simulate bottleneck by deleting cells" << endl;
    cout << "Simulate bottleneck by selecting cells" << endl;
    int i = 0;
    // create a new clone with selected cells since deletion is likely to cause pointer errors
    Clone* s = new Clone(1, 0);
    while(i < bottleneck){
    // while(i < time_bottleneck - bottleneck){
      // cout << "delete " << i+1 << "th cell" << endl;
      int rindex = myrng(glands.clones[0]->curr_cells.size());
      Cell_ptr cell = glands.clones[0]->curr_cells[rindex];
      // cout << "select surviving cell " << cell->cell_ID << " as cell " << i+1 << endl;
      // restart counting cell ID to be consistent with grow function
      cell->cell_ID = i + 1;
      s->curr_cells.push_back(cell);
      // remove cells from current population
      // if(glands.clones[0]->curr_cells[rindex]->cell_ID != 1){ // keep start cell
      //   delete (glands.clones[0]->curr_cells[rindex]);
      //   glands.clones[0]->curr_cells[rindex] = NULL;
      // }
      glands.clones[0]->curr_cells.erase(glands.clones[0]->curr_cells.begin() + rindex);
      i++;
    }
    // insert new clone at the top
    glands.clones.insert(glands.clones.begin(), s);
    // update ntot which is baseline for new cell ID for cases when the population starts from several cells
    glands.clones[0]->ntot = glands.clones[0]->curr_cells.size();
    // start from remaining cells
    cout << "Regrow from the bottleneck population " << glands.clones[0]->curr_cells.size() << endl;
    cout << "old population now has size " << glands.clones[1]->curr_cells.size() << endl;
    glands.simulate_gland_as_cell(start_cell, ndeme, start_model, lineages, store_lineage, loc_type, leap_size, track_lineage, multiple_output, verbose, 0);
  }else{
    glands.simulate_gland_as_cell(start_cell, ndeme, start_model, lineages, store_lineage, loc_type, leap_size, track_lineage, multiple_output, verbose);
  }
  cout << "Simulation stops and output samples" << endl;
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

  if(verbose > 0 && fdeme != ""){
      string header = "Parent\tID\tTime\tNmut";
      glands.print_gland_lineage(fdeme, header, lineages);
  }

  map<int, vector<double>> sstat_by_time;
  map<int, vector<double>> div_by_time;  // divergence by time
  // print CNPs taken at different times
  for(auto cnp : glands.clones[0]->cnp_by_time){
      int nsize = cnp.first;
      vector<int> gids;
      string fname = outdir + "gland_cn_size" + to_string(nsize) + "_" + suffix + ".txt";
      ofstream fcn;
      fcn.open(fname, ofstream::trunc | ofstream::out);

      for(auto pcn : cnp.second){  // a vector of all CNPs
        for(auto cell_pcn : pcn){
          gids.push_back(cell_pcn.first);
          string sname = to_string(cell_pcn.first);
          print_bulk_cn(sname, cell_pcn.second, fcn, verbose);
        }
      }

      // compute summary statistics
      if(verbose > 0) cout << "sampled cells " << gids.size() << endl;
      sort(gids.begin(), gids.end());
      map<int, double*> avg_loc_changes = glands.clones[0]->lchange_by_time[nsize];
      double pga = glands.print_bin_subclonal_stat_by_type(avg_loc_changes, 0, verbose);
      vector<double> alters;
      double avg_var = glands.print_pairwise_divergence(gids, avg_loc_changes, alters, 0, 1, verbose);
      vector<double> sstats{pga, avg_var};
      sstat_by_time[nsize] = sstats;
      div_by_time[nsize] = alters;
      // vector<double> svars = print_pairwise_divergence(ids, avg_loc_changes, 0, verbose);
      // double nbp = collect_private_subclonal_bps(avg_loc_changes, verbose);
      // double nmut = print_num_uniq_mut(verbose);
      // double avg_ndiff = print_pairwise_mismatch(verbose);
  }

  vector<int> gids;
  // Get average absolute CNP of all glands (cells)
  map<int, double*> avg_loc_changes;
  // assert(glands.clones[0]->curr_cells.size() == ndeme);
  // Initialize avg_loc_changes to be global
  int nsampled = 0;
  // randomly sample n glands
  unordered_set<int> samples = BobFloydAlgo(TOSAMPLE, glands.clones[0]->curr_cells.size());
  if(verbose > 0){
    cout << "sampled cell IDs at " << glands.clones[0]->curr_cells.size() << ":";
    for(auto cid : samples){
      cout << " " << cid;
    }
    cout << endl;
  }
  int scount = 0;
  for(auto c : glands.clones[0]->curr_cells){
      // cout << c->cell_ID << endl;
      // only store sampled glands if there is sampling
      // if(nsampled >= TOSAMPLE) break;
      // double rnum = runiform(r, 0, 1);
      // if(rnum > 0.5){
      //   nsampled += 1;
      if(samples.find(scount) != samples.end()){
        glands.sample_IDs.push_back(c->cell_ID);  // used for plotting

        gids.push_back(c->cell_ID);
        avg_loc_changes[c->cell_ID] = new double[NUM_LOC];
        for(int i = 0; i < NUM_LOC; i++){
            avg_loc_changes[c->cell_ID][i] = c->loc_changes[i];
            // if(c->loc_changes[i] != NORM_PLOIDY)
            //   cout << i << "\t" << avg_loc_changes[c->cell_ID][i] << endl;
        }
      }
      scount++;
  }
  if(verbose > 0) cout << "sampled cells " << gids.size() << endl;
  sort(gids.begin(), gids.end());
  double pga = glands.print_bin_subclonal_stat_by_type(avg_loc_changes, 0, verbose);
  vector<double> alters;
  double avg_var = glands.print_pairwise_divergence(gids, avg_loc_changes, alters, 0, 1, verbose);
  vector<double> sstats{pga, avg_var};
  int currsize = glands.clones[0]->curr_cells.size();
  sstat_by_time[currsize] = sstats;
  div_by_time[currsize] = alters;

  glands.clones[0]->get_fitness_stats();

  glands.print_cell_cnp(glands.clones[0], outdir, suffix, verbose);

  // seed size PGA divergence mean_fitness max_fitness min_fitness distance_to_optimum
  if(fstat != ""){
    ofstream fout;
    fout.open(fstat, ofstream::app);
    for(auto ss : sstat_by_time){
      string msg = to_string(seed) + "\t" + to_string(ss.first) + "\t" + to_string(ss.second[0]) + "\t" + to_string(ss.second[1]) + "\t" + to_string(glands.clones[0]->meanfit_by_time[ss.first]) + "\t" + to_string(glands.clones[0]->maxfit_by_time[ss.first]) + "\t" + to_string(glands.clones[0]->minfit_by_time[ss.first]) + "\t" + to_string(glands.clones[0]->dist_by_time[ss.first]);
      // + "\t" + to_string(glands.clones[0]->dm_by_time[ss.first]);
      cout << msg << '\n';
      fout << msg << endl;
    }
    fout.close();
  }

  // Print all pairwise divergences
  if(fdiv != ""){
    ofstream fout;
    fout.open(fdiv, ofstream::trunc | ofstream::out);
    for(auto ss : div_by_time){
      string msg = to_string(seed) + "\t" + to_string(ss.first);
      for(auto div : ss.second){
        msg += "\t" + to_string(div);
      }
      // cout << msg << '\n';
      fout << msg << endl;
    }
    fout.close();
  }

  // Print all individual fitnesses (time and memory consuming, not much informaton revealed)
  // if(ffit != ""){
  //   ofstream fout;
  //   fout.open(ffit, ofstream::trunc | ofstream::out);
  //   for(auto ss : glands.clones[0]->fits_by_time){
  //     string msg = to_string(seed) + "\t" + to_string(ss.first);
  //     for(auto fit : ss.second){
  //       msg += "\t" + to_string(fit);
  //     }
  //     // cout << msg << '\n';
  //     fout << msg << endl;
  //   }
  //   fout.close();
  // }
}


int main(int argc, char const *argv[]) {
    clock_t tStart = clock();

    int mode;

    int ndeme;
    int max_deme_size;
    string fdeme;
    int track_lineage;

    string fmut;

    string fgenotype;     // starting karyotype
    string flprob;     // location mutation probability
    string fgain;
    string floss;
    string fgenotype_chr;

    double leap_size;

    string num_glands; // number of samples to take at each side

    double birth_rate, death_rate;
    double mutation_rate;
    // parameters for mean of dup/del size distributions
    // int mean_gain_size, mean_loss_size;
    int loc_type;
    int bottleneck;  // bottleneck size
    int time_bottleneck;

    int model_ID;   // not making differences, fitness is used to introduce selection
    int use_alpha;  // when alpha is used, karyotype differences are considered
    double fitness;
    int genotype_diff;
    int growth_type;
    int norm_by_bin;

    int stat_type;
    double frac_cutoff;
    double min_freq;
    double max_freq;
    double delta;
    int use_std;

    string outdir, suffix; // output
    string fstat, fdiv, ffit;
    int multiple_output;

    unsigned long seed;

    int verbose;

    namespace po = boost::program_options;

    po::options_description generic("Generic options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ;

    // parameters that are essential to simulations
    po::options_description required("Required parameters");
    required.add_options()
      ("mode,m", po::value<int>(&mode)->default_value(0), "mode of simulation. 0: gland as cell; 1: gland fission when reaching certain size; 2: with potential metastasis")

      ("birth_rate,b", po::value<double>(&birth_rate)->default_value(1), "birth rate")
      ("death_rate,d", po::value<double>(&death_rate)->default_value(0), "death rate")

      ("mutation_rate,r", po::value<double>(&mutation_rate)->default_value(0.2), "copy number alteration (CNA) rate per division")
      ("mean_gain_size", po::value<int>(&MEAN_GAIN_SIZE)->default_value(MEAN_GAIN_SIZE), "mean size of segment gain (in terms of bins)")
      ("mean_loss_size", po::value<int>(&MEAN_LOSS_SIZE)->default_value(MEAN_LOSS_SIZE), "mean size of segment loss (in terms of bins)")
      // by default, not consider chromosome boundaries to save computation
      ("loc_type", po::value<int>(&loc_type)->default_value(2), "type of CNA unit. 0: BIN-level at one chromosome; 1: arm-level; 2: BIN-level across multiple chromsomes")

      ("ndeme,n", po::value<int>(&ndeme)->default_value(1000), "number of glands when simulation stops")
      ("max_deme_size", po::value<int>(&max_deme_size)->default_value(10000), "maximum number of cells in a gland before gland fission (used in mode 1)")
      ("num_glands,g", po::value<string>(&num_glands)->default_value("2 30 20"), "the number of sides followed by the number of glands taken at different sides of a crypt, separated by space")

      ("odir,o", po::value<string>(&outdir)->required()->default_value("./"), "output directory")
       ;

    po::options_description optional("Optional parameters");
    optional.add_options()
      // input files which specifies #sampled glands and CNA informaton
      ("fgenotype", po::value<string>(&fgenotype)->default_value(""), "TSV file with starting karyotype of the first cell")
      ("flprob", po::value<string>(&flprob)->default_value(""), "TSV file with probability of mutation at each location")
      ("fgain", po::value<string>(&fgain)->default_value(""), "TSV file with size of all copy number gains in the real data")
      ("floss", po::value<string>(&floss)->default_value(""), "TSV file with size of all copy number losses in the real data")
      ("fgenotype_chr", po::value<string>(&fgenotype_chr)->default_value(""), "TSV file with starting karyotype of the first cell at chromosome level")

      // options related to model of evolution
      ("model", po::value<int>(&model_ID)->default_value(0), "model of evolution. 0: neutral; 1: selection")
      ("use_alpha", po::value<int>(&use_alpha)->default_value(1), "whether or not to use alpha in selection model. 0: use selection coefficient; 1: use alpha")
      ("fitness,f", po::value<double>(&fitness)->default_value(0), "fitness values of mutatants")
      ("genotype_diff", po::value<int>(&genotype_diff)->default_value(3), "type of karyotype difference (L1 distance) in simulating selection. 0: not considered; 1: number of altered bins; 2: L1 distance of CN vectors at bin level; 3: number of new mutations (only known in simulated data); 4: distance of CN vectors computed when starting from diploid")
      ("norm_by_bin", po::value<int>(&norm_by_bin)->default_value(0), "whether or not to normalize karyotype difference by number of bins (segments) in the genome. 0: no; 1: yes")
      ("growth_type,t", po::value<int>(&growth_type)->default_value(0), "Type of growth when adding selection. 0: only birth; 1: change birth rate; 2: change death rate; 3: change both birth or death rate")

      // options for simulation starting from diploid
      ("start_with_opt", po::value<int>(&START_WITH_OPT)->default_value(1), "the start karyotype of first node. 0: diploid. 1: optimum karyotype specified by input file if specified or using diploid as optimum")
      ("bottleneck", po::value<int>(&bottleneck)->default_value(0), "bottleneck size during cell growth")
      ("time_bottleneck", po::value<int>(&time_bottleneck)->default_value(1000), "population size when bottleneck is introduced")
      ("multiple_output", po::value<int>(&multiple_output)->default_value(0), "output CNP at different sampling time points")
      ("cna_level", po::value<int>(&CHR_CNA)->default_value(0), "the level of CNA. 0: bin-level. 1: chr-level")
      ("nsample", po::value<int>(&TOSAMPLE)->default_value(100), "number of samples taken at each time point")
      ("weight_optimum", po::value<double>(&WEIGHT_OPTIMUM)->default_value(1.0), "controlling probability of mutation at locations in optimum karyotype")

      // options related to summary statistics
      ("stat_type,s", po::value<int>(&stat_type)->default_value(4), "type of summary statistics. 0: variance; 1: clone-pairwise differences; 2: average CNP; 3: complete CNP; 4: sample-pairwise differences")
      ("min_freq", po::value<double>(&min_freq)->default_value(0.0), "minimal mutation frequency to consider")
      ("max_freq", po::value<double>(&max_freq)->default_value(0.25), "maximal mutation frequency to consider")
      ("delta", po::value<double>(&delta)->default_value(0.025), "step size of mutation frequency histogram")
      ("use_std", po::value<int>(&use_std)->default_value(0), "whether or not to use standard deviation of pairwise divergence. If not, the variance is output")
      ("bp_cutoff", po::value<double>(&BP_CUTOFF)->default_value(BP_CUTOFF), "threshold to determine whether a breakpoint can be detected or not")
      ("bin_cutoff", po::value<double>(&BIN_CUOFF)->default_value(BIN_CUOFF), "threshold to determine whether a bin can be detected as altered or not")
      ("frac_cutoff", po::value<double>(&frac_cutoff)->default_value(0.5), "cutoff of counting breakpoints")

      // options related to output
      ("suffix", po::value<string>(&suffix)->default_value(""), "suffix of output file")
      ("track_lineage", po::value<int>(&track_lineage)->default_value(0), "whether or not to trace lineages. 0: No (only current population are kept). 1: Yes")
      ("fdeme", po::value<string>(&fdeme)->default_value(""), "file with deme (gland) relationships")
      ("fmut", po::value<string>(&fmut)->default_value(""), "file with mutation counting informaton")
      ("fstat", po::value<string>(&fstat)->default_value(""), "file with summary statistics at different population size")
      ("fdiv", po::value<string>(&fdiv)->default_value(""), "file with pairwise divergences at different population size")
      ("ffit", po::value<string>(&ffit)->default_value(""), "file with individual fitnesses at different population size")

      ("leap_size", po::value<double>(&leap_size)->default_value(0.0), "step size of tau-leaping to accelerate the simulations (not used for now)")

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
            cout << "simgland [version 0.1], a program to simulate copy number variations along a stochastic branching tree" << endl;
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
        cout << "Simulating copy numbers called from bulk WGS and WES samples of glands sampled from a patient crypt" << endl;
        cout << "Fitness:\t" << fitness << endl;
        cout << "Mutation rate:\t" << mutation_rate << endl;
        cout << "Mutation level: \t" << CHR_CNA << endl;
    }

    // cout << "Using Boost "
    //   << BOOST_VERSION / 100000     << "."  // major version
    //   << BOOST_VERSION / 100 % 1000 << "."  // minor version
    //   << BOOST_VERSION % 100                // patch level
    //   << std::endl;

    // read real CNA sizes from file
    string size;
    if(fgain !=  ""){
      ifstream fcn_gain(fgain);
      while(getline(fcn_gain, size)){
          real_gain_sizes.push_back(atof(size.c_str()));
      }
      // cout << real_gain_sizes.size() << endl;
    }

    if(floss !=  ""){
      ifstream fcn_loss(floss);
      while(getline(fcn_loss, size)){
          real_loss_sizes.push_back(atof(size.c_str()));
      }
      // cout << real_loss_sizes.size() << endl;
    }

    // read real mode karyotype from file, used to get optimum karyotype (bin level) and probability of CNAs at different bins
    if(fgenotype != ""){
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
          OPT_KARYOTYPE[vals[0] - 1] = vals[1];
          // std::cout << std::endl;
        }
        infile.close();

        for(int i = 1; i < NUM_LOC; i++){
          // only at boundaries of CNAs
          if(OPT_KARYOTYPE[i] != NORM_PLOIDY && OPT_KARYOTYPE[i-1] == NORM_PLOIDY){
            LOC_PROBS[i] = LOC_PROBS[i] * WEIGHT_OPTIMUM;
          }
        }
    }

    // when no karyotype file is specified, it is normal diploid
    if(START_WITH_OPT == 1){
      // assert(fgenotype != "");
      if(verbose > 1){
        cout << "Assigning optimum karyotype to the starting node" << endl;
      }
      for(int i = 0; i < NUM_LOC; i++){
        START_KARYOTYPE[i] = OPT_KARYOTYPE[i];
      }
    }else{
      if(verbose > 1){
        cout << "Assigning diploid karyotype to the starting node" << endl;
      }
      if(CHR_CNA == 1){
        for(int i = 0; i < NUM_CHR; i++){
          START_KARYOTYPE_CHR[i] = NORM_PLOIDY;
        }
      }

      // always initialize bin level changes (used for plotting, computation of summary statistics)
      for(int i = 0; i < NUM_LOC; i++){
        START_KARYOTYPE[i] = NORM_PLOIDY;
      }

      OPT_BIRTHRATE = birth_rate;
      OPT_DEATHRATE = death_rate;
    }
    // for(auto cn : start_genotype){
    //     cout << cn << endl;
    // }

    // cout << "LOC_PROBS before" << endl;
    // for(int i = 0; i < NUM_LOC; i++){
    //     cout << LOC_PROBS[i] << endl;
    // }

    // read real #CNAs per bin from file
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
    Model start_model(model_ID, genotype_diff, growth_type, fitness, use_alpha, norm_by_bin);
    Cell_ptr start_cell = new Cell(1, 0, birth_rate, death_rate, mutation_rate, 0);

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
        cout << "\nSimulating gland growth with total CNAs, with " <<  ndeme << " glands in total" << endl;
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
        vector<double> alters;
        glands.print_sstat(gids, avg_loc_changes, alters, stat_type, loc_type, min_freq, max_freq, delta, use_std, verbose);

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
          // first 2 columns have no meanings
          fcf << "0\t0\t" << glands.muts.size() << "\t" << glands.clones.size() << "\t" << fission_rate << endl;
        }
    }else if(mode == 0){
        if(verbose > 0) cout << "Simulating gland as cell" << endl;
        int store_lineage = 0;
        vector<string> lineages;
        if(fdeme != ""){
            store_lineage = 1;
        }

        glands.root = new node(start_cell->cell_ID);
        glands.simulate_gland_as_cell(start_cell, ndeme, start_model, lineages, store_lineage, loc_type, leap_size, track_lineage, 0, verbose);

        // sample from two sides of the lineage tree if the numbers are specified
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
        // Get average absolute CNP of all glands (cells)
        map<int, double*> avg_loc_changes;
        // assert(glands.clones[0]->curr_cells.size() == ndeme);
        // Initialize avg_loc_changes to be global
        vector<int> curr_cIDs;
        for(auto c : glands.clones[0]->curr_cells){
            // cout << c->cell_ID << endl;
            curr_cIDs.push_back(c->cell_ID);
            // only store sampled glands if there is sampling
            if(glands.sample_IDs.size() > 0 && find(glands.sample_IDs.begin(), glands.sample_IDs.end(), c->cell_ID) == glands.sample_IDs.end()) continue;
            // cout << " adding into list" << endl;
            gids.push_back(c->cell_ID);
            avg_loc_changes[c->cell_ID] = new double[NUM_LOC];
            for(int i = 0; i < NUM_LOC; i++){
                avg_loc_changes[c->cell_ID][i] = c->loc_changes[i];
                // if(c->loc_changes[i] != NORM_PLOIDY)
                //   cout << i << "\t" << avg_loc_changes[c->cell_ID][i] << endl;
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
        vector<double> alters;
        glands.print_sstat(gids, avg_loc_changes, alters, stat_type, loc_type, min_freq, max_freq, delta, use_std, verbose);

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
    }else{
      assert(mode == 2);
      simulate_from_diploid(seed, start_cell, glands, start_model, ndeme, WEIGHT_OPTIMUM,
         bottleneck, time_bottleneck, track_lineage, loc_type, leap_size, multiple_output,
         fdeme, fgenotype_chr, outdir, suffix, fstat, fdiv, ffit, tStart, verbose);
    }

    return 0;
}
