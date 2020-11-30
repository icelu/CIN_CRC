# Introduction
This program mainly simulates the growth of individual glands with a stochastic birth-death branching process to learn mutation rates and selection forces.

The major output are summary statistics extracted from copy number profiles (CNPs) of multiple glands sampled from a patient.
Only copy number gain and loss are simulated.
These summary statistics can be used for inference with ABC (approximate Bayesian computation) approach.

There is also a program to simulate copy numbers obtained from multi-region sequencing data of a patient (not user friendly for now).


# Install
## Required library
* [GSL](https://www.gnu.org/software/gsl/)
* [boost](https://www.boost.org)

## Compilation
Download the source code and then run `make` in folder "code" to generate the simulation program 'simgland'.

Comment out relevant line to generate the simulation program 'simcrc'.


# Usage of simgland
Please run `./simgland -h` for all options.

There are two modes of simulation.
* Mode 0 (default): simulating each gland as a node.
* Mode 1: simulating each cell as a node and a gland is composed of a group of cells. The gland will divide into two glands when reaching a certain number of cells.

In each mode, there are two models of evolution.
* Model 0 (default): neutral.
* Model 1: selection. Under this model, the selection strengths can be specified by option -f.


Option "verbose" controls how much to output.

* When verbose = 0, a vector of summary statistics is written to the standard output.
The summary statistics represents
percentage of genome altered (PGA) relative to start karyotype, average pairwise divergence (proportion of different bins that are altered), variance of pairwise divergences, number of unique breakpoints across all sampled glands.
We also computed the number of unique mutations across all sampled glands and average pairwise different number of mutations.
The first three values were used for ABC inference of mutation rates and selection strengths.

* When verbose > 0, more information will be written to the standard output and there is a file containing the copy numbers for each sampled gland. Using '--fdeme FILENAME' will generate the lineages of glands.

<!-- using --fmut FILENAME will generate mutation information of the glands. -->


## Example
### Simulate multiple glands sampled from a single patient

The example commands show the basic usage of the program.
They will simulate 1000 glands and then sample 30 glands from each side of the lineage tree under mutation rate 0.1.

* Assume neutral evolution:
```
./simgland -o ./ -n 1000 -g "2 30 30" --seed $RANDOM --verbose 0 -r 0.1
```

* Assume (negative) selection:
```
./simgland -o ./  -n 1000 -g "2 30 30" --seed $RANDOM --verbose 0 -r 0.1 --model 1 -f 2
```
Here, Option '--model 1' indicates selection is imposed.
Option '--fitness 2' specifies the strengths of selection.
