#ifndef UTIL_HPP
#define UTIL_HPP


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>

#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <set>
#include <climits>
#include <cstdlib>

#include <unistd.h>

using namespace std;


typedef map<pair<int, int>, int> pcn;   // copy number at a position
typedef map<pair<int, int>, double> dpcn;   // copy number at a position


const int N_MUT_TYPE = 5;
// 22 for CRC data, 23 for PDO data
// const int NUM_CHR = 23;
const int NUM_CHR = 22;
const int NORM_PLOIDY = 2;

// Number of chromsomes affected in a multipolar event
const int MULTI_NCHR = 16;

// The number of bins for each chromosome. Each bin corresponds to a window of size 500,000 bp. 313,115 for chr X, Y
// const vector<int> CHR_BIN_SIZE{498,485,397,381,364,342,319,291,277,268,271,267,229,215,204,181,167,161,118,129,94,102};
// const int NUM_LOC = 5760;
// const int MEAN_GAIN_SIZE = 62;
// const int MEAN_LOSS_SIZE = 110;
// The number of bins for each chromosome. Each bin corresponds to a window of size 1,500,000 bp.
const vector<int> CHR_BIN_SIZE{166,162,133,127,122,114,107,97,93,90,91,89,77,72,68,61,56,54,40,43,32,34};
const int NUM_LOC = 1928;
const int MEAN_GAIN_SIZE = 21;
const int MEAN_LOSS_SIZE = 37;
const int MAX_DEME_SIZE = 10000;

const vector<double> LOC_PROBS_vec(NUM_LOC, 1.0/NUM_LOC);
const double* LOC_PROBS = &LOC_PROBS_vec[0];


enum CNA_type{BIN, ARM, PSEUDO, BPOINT};
const double POS_SIGMA = 1;
const int MIN_NCELL = 100;
const int DEC_PLACE = 10;

// Position offset for migration clones
const double MIGRATE_OFFSET = 100.0;

enum SStat_type{LVAR, ADIFF, AVG, CMPL, DIFF, BP, ALTBIN};
enum Growth_type{ONLY_BIRTH, CHANGE_BIRTH, CHANGE_DEATH, CHANGE_BOTH};

double BP_CUTOFF = 0.1;
double BIN_CUOFF = 0.1;

gsl_rng * r;
std::mt19937 eng;

unsigned long setup_rng(unsigned long set_seed){
  gsl_rng_env_setup();

  const gsl_rng_type* T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  if( set_seed != 0 ){
    gsl_rng_set(r, set_seed);
    return(set_seed);
  }else{
    int t = time(NULL);
    int pid = getpid();
    long s = t*pid;
    //cout << "pid:" << "\t" << getpid() << endl;
    // cout << "seed:" << "\t" << t << "\t" << pid << "\t" << abs(s) << endl;
    std::cout << "seed:" << "\t" << abs(s) << std::endl;
    gsl_rng_set (r, abs(s));
    return(abs(s));
  }
}

// unary function and pointer to unary function
// allows use of gsl rng for standard template algorithms
long unsigned myrng(long unsigned n)
{
  return gsl_rng_uniform_int(r, n);
}

long unsigned (*fp)(long unsigned) = myrng;

// factorial for choose(n,k)
int fact(int n){
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}

// wrapper for uniform
double runiform(gsl_rng* r, double a, double b){
  double myrandom = a + (b-a)*gsl_rng_uniform (r);

  while(myrandom==0){
    myrandom = a + (b-a)*gsl_rng_uniform (r);
  }
  return myrandom;
}


// sample an element according to a probability vector
// gsl_ran_multinomial?
int rchoose(gsl_rng* r, const std::vector<double>& rates){
  //cout << "rchoose, rates:";
  //for(int i=0; i<rates.size(); ++i) cout << "\t" << rates[i];
  //cout << endl;

  std::vector<double> p;
  double s = accumulate(rates.begin(), rates.end(), 0.0);

  for(int i=0; i<rates.size(); ++i){
    p.push_back(rates[i]/s);
  }

  std::vector<double> psum(p.size(),0.0);
  partial_sum(p.begin(),p.end(),psum.begin());

  double u = gsl_rng_uniform(r);
  int ret = -1;
  for(int i=0; i<rates.size(); ++i){
    if( u < psum[i] ){
      ret = i;
      break;
    }
  }

  //cout << "u=\t" << u << "\t" << ret << endl;
  return ret;
}



void set_outdir(string outdir, int verbose = 0){
    const char* path = outdir.c_str();
    boost::filesystem::path dir(path);
    if(boost::filesystem::create_directory(dir))
    {
        if(verbose > 0) cerr << "Directory Created: " << outdir <<endl;
    }
}


// Read strings separated by space into a vector
template <typename T>
void get_vals_from_str(vector<T>& vals, string str_vals, int num=0){
    assert(str_vals != "");
    stringstream ss(str_vals);
    if(num==0)   ss >> num;
    for(int i = 0; i < num; i++){
        T d1;
        ss >> d1;
        vals.push_back(d1);
    }
    // cout << "vector from " << str_vals << ": ";
    // for(auto v : vals){
    //     cout << "\t" << v;
    // }
    // cout << endl;
}



// Note: the sites are counted for all the haplotypes. The positions on different haplotypes are counted as one site.
// The starting point of sites/chr/seg is 0
void site2chr(int site, int& chr, int& seg, const vector<int>& chr_lengths, int verbose = 0){
    // cout << "called here" << endl;
    // if(verbose > 1) cout << "Site to convert is " << site << endl;
    if(site >= 0 && site < chr_lengths[0]){
        chr = 0;
        seg = site;
    }
    else{
        // if(verbose > 1) cout << "There are " << chr_lengths.size() << " chromosomes" << endl;
        for(int i=0; i<chr_lengths.size(); i++){
            int sum1 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 1, 0);
            int sum2 = accumulate(chr_lengths.begin(), chr_lengths.begin() + i + 2, 0);
            // if(debug){
            //     cout << "Size for chr " << i+1 << " is " << chr_lengths[i] << endl;
            //     cout << "sum until chromosome " << i + 1 << " is " << sum1 << endl;
            //     cout << "sum until chromosome " << i + 2 << " is " << sum2 << endl;
            // }
            if(site >= sum1 && site < sum2){
                chr = i + 1;
                seg = site - sum1;
                break;
            }
        }
    }
    assert(seg < chr_lengths[chr]);
}


template <typename T>
void loc2pcn(T loc_changes[], map<pair<int, int>, T>& bulk_cnp, int verbose = 0){
    bulk_cnp.clear();
    for(int loc = 0; loc < NUM_LOC; loc++){
        int chr = 0;
        int start = 0;
        site2chr(loc, chr, start, CHR_BIN_SIZE, verbose);
        pair<int, int> pos(chr, start);
        bulk_cnp[pos] = loc_changes[loc];
    }
    // if(verbose > 0){
    //     cout << "CNP converted from simple bins" << endl;
    //     for(auto cp : bulk_cnp){
    //         cout << cp.first.first << "\t" << cp.first.second  << "\t"<< cp.second << endl;
    //     }
    // }
}

/*
   This method prints out the average copy numbers of final cells in a clone
 */
template <typename T>
void print_bulk_cn_orig(string sample, const map<pair<int, int>, T>& bulk_cnp, ofstream& fcn, int verbose = 0){
    if(verbose > 1) cout << "Printing average copy numbers of all current cells (mimicking bulk sample)" << endl;

    for(auto cp : bulk_cnp){
        // assert(cp.second != 0);
        pair<int, int> pos = cp.first;
        fcn << sample << "\t" << pos.first + 1 << "\t" << pos.second << "\t" << cp.second << endl;
    }
}


// Cannot use const when using iterator
template <typename T>
void print_bulk_cn(string sample, map<pair<int, int>, T>& bulk_cnp, ofstream& fcn, int verbose = 0){
    if(verbose > 1) cout << "Printing merged average copy numbers of all current cells (mimicking bulk sample)" << endl;
    if(bulk_cnp.size() == 0){
        if(verbose > 0) cout << "No bulk CNP detected" << endl;
        return;
    }
    // merge continuous regions
    bool has_change = true;
    typename map<pair<int, int>, T>::iterator it = bulk_cnp.begin();
    pair<int, int> pos = (*it).first;
    int prev_chr = pos.first + 1;
    int prev_bin = pos.second;
    T prev_rcn = (*it).second;
    int start = prev_bin;
    ++it;

    for (;it != bulk_cnp.end(); ++it){
        pos = (*it).first;
        int curr_chr = pos.first + 1;
        int curr_bin = pos.second;
        T curr_rcn = (*it).second;
        if(prev_chr != curr_chr || curr_bin - prev_bin != 1 || prev_rcn != curr_rcn){
            fcn << sample << "\t" << prev_chr << "\t" << start << "\t" << prev_bin << "\t" << prev_rcn << endl;
            // start a new region
            start = curr_bin;
        }
        prev_chr = curr_chr;
        prev_bin = curr_bin;
        prev_rcn = curr_rcn;
    }
    // For regions ending at the last position
    fcn << sample << "\t" << prev_chr << "\t" << start << "\t" << prev_bin << "\t" << prev_rcn << endl;
}


#endif
