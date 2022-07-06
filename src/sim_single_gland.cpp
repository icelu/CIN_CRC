#include <cstdlib>
#include <sstream>
#include <string>
#include <climits>
#include <time.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "sgland.hpp"


using namespace std;


const int DAY_IN_YEAR = 365;


// number of divisions, patient age: 75, the age of the MRCA cell of the cancer for patient C274 (11 or 7), use 9
// originally treat N, p as variables, mutation (methylation) rate fixed, #divisions determined by age with fixed division rate one division per day, 2048 cells in total
// For C274: the age of the most recent common ancestor (MRCA) cell of the cancer
// observed data: CNAs of 112 single cells with one bulk


// Modeling single cell expansion in a gland/crypt to fit sampled single cell CNA data
// aim: estimate mutation rate from observed CNAs in sampled cells
// assumption: discrete and non-overlapping generations
int main(int argc, char const *argv[]){
    string outdir, prefix, suffix, flocus, fcna, fmut;

    /******************* paramters obtained from literature *******************/
    int num_stem;  // N stem cells starting with core karyotype, range
    int gsize;   // total number of cells in the gland, range 2000 to 10000
    double birth_rate;    // cell division rate, use 0.5 per day? , range
    double prob_d1;   // asymmetric division rate, range
    int num_gen_nonstem;  // non-stem cell lineages have a maximal life time of h generations, depending on crypt size
    int num_gen_steady;   // steady-state process where oldest non-stem cells no longer divide but stay around before exiting the crypt

    /******************* paramters derived from real data *******************/
    int age_cancer;  // in year, estimated based on SNV, range 7 to 11
    int num_sample;   // taking random sample to fit real data

    /******************* paramters to estimate *******************/
    // CNA rate as variable
    // separate arm- and chr- level to better fit real data
    double chr_rate;
    double arm_rate;

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
      ("age_cancer", po::value<int>(&age_cancer)->default_value(9), "the age of the MRCA cell of the cancer for patient")
      ("birth_rate,r", po::value<double>(&birth_rate)->default_value(0.4), "cell division rate")
      // ("death_rate", po::value<double>(&death_rate)->default_value(0), "death rate")
      ("prob_d1,d", po::value<double>(&prob_d1)->default_value(0.8), "asymmetric division probability")
      ("num_stem,n", po::value<int>(&num_stem)->default_value(100), "number of stem cells in a gland")
      ("num_gen_nonstem", po::value<int>(&num_gen_nonstem)->default_value(3), "number of generations where non-stem cells keep dividing")
      ("num_gen_steady", po::value<int>(&num_gen_steady)->default_value(3), "number of generations where non-stem cells keep steady")
      ("gsize,g", po::value<int>(&gsize)->default_value(5000), "size of cell populations in a gland")

      // ("file_cmut", po::value<string>(&file_cmut)->default_value(""), "the TSV file which contains clonal CNVs, with three columns (chr, arm, cn)")

      ("arm_rate,a", po::value<double>(&arm_rate)->default_value(0.1), "arm-level CNA rate per cell division")
      ("chr_rate,c", po::value<double>(&chr_rate)->default_value(0.1), "chr-level CNA rate per cell division")
      // ("multi_rate, m", po::value<double>(&multi_rate)->default_value(0), "multiple simultaneous chr-level CNA rate per cell division")

      ("num_sample", po::value<int>(&num_sample)->default_value(112), "number of cells sampled from the gland")

      ("prefix,p", po::value<string>(&prefix)->default_value(""), "prefix of output file")
      ("suffix,s", po::value<string>(&suffix)->default_value(""), "suffix of output file")
      ("flocus,s", po::value<string>(&flocus)->default_value(""), "output file recording locus-specific diversity")
      ("fcna", po::value<string>(&fcna)->default_value(""), "output file recording copy number profile")
      ("fmut", po::value<string>(&fmut)->default_value(""), "output file recording mutations")

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
            cout << "sim_gland [version 0.1], a program to simulate copy number variations during cancer gland growth" << endl;
            return 1;
        }
        po::notify(vm);
    }
    catch (const exception& e) {
          cerr << e.what() << endl;
          return 1;
    }


    unsigned long rseed = setup_rng(seed);

    int num_gen = age_cancer * DAY_IN_YEAR * birth_rate;
    if(verbose > 0){
      cout << "patient age " << age_cancer << endl;
      cout << "number of stem cells " << num_stem << endl;
      cout << "number of total cells " << gsize << endl;
      cout << "probability of asymmetric division " << prob_d1 << endl;
      cout << "number of generations " << num_gen << endl;
    }

    SGland g(num_stem, gsize, num_gen, num_gen_nonstem, num_gen_steady, prob_d1, arm_rate, chr_rate);
    g.initialize_stem_cells();
    g.get_rate_children(num_stem, prob_d1, verbose);
    g.grow_with_cna(verbose);

    vector<GCell*> sampled_cells;
    g.sample_cell_from_gland(sampled_cells, num_sample, verbose);

    cout.precision(9);
    vector<int> avg_cn(2, 0);
    g.get_avg_reciprocal_cn(sampled_cells, avg_cn);
    vector<double> ac_scaled(2, 0.0);
    for(int i = 0; i < avg_cn.size(); i++){
        double ac = (double) avg_cn[i] / sampled_cells.size();
        // sum_stats.push_back(ac);
        ac_scaled[i] = ac;
    }
    cout << ac_scaled[0] << "\t" << ac_scaled[1] << "\n";

    if(flocus != ""){
      g.get_locus_diversity();
      cout << "writing locus-specific diversity" << endl;
      ofstream fout(flocus);
      for(auto ld : g.locus_diversity){
          // chr, arm, #cells
          fout << ld.first.first + 1 << "\t" << ld.first.second << "\t"  << ld.second << "\n";
      }
      fout.close();
    }

    // output CNAs in each cell for visualization
    if(fcna != ""){
      cout << "writing copy numbers in each cell" << endl;
      ofstream fout(fcna);
      g.write_cna(fout);
      fout.close();
    }

    // check the distribution of CNA frequency, similar to VAF distribution
    if(fmut != ""){
      cout << "writing mutation information" << endl;
      ofstream fout(fmut);
      g.write_mut(fout);
      fout.close();
    }

}
