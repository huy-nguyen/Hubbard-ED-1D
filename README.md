# Exact Diagonalization of the Hubbard Model in 1-D

This repository contains the `MATLAB` code to perform exact calculations of the *imaginary-time correlation functions* of the *Hubbard* model in one dimension.

The Hubbard model is widely believed to be the model that describes [high-temperature superconductivity](https://en.wikipedia.org/wiki/High-temperature_superconductivity). I presented the theory behind this model in a manner accessible to senior-year physics majors in Chapter 2 of my [undergraduate thesis](http://www.huy-nguyen.com/Reed-thesis).

The git history of this repo represents the work I did over a 14-day period in 2014. Between the very first commit (208e393) and the last commit (f29f165), the efficiency of the code improved by a factor of over a thousand due to my extensive profiling and optimization.

The code is designed to run in parallel within a single compute node on multi-core Linux high-performance clusters using `MATLAB`'s [Parallel Computing Toolbox](http://www.mathworks.com/products/parallel-computing/). Furthermore, `MATLAB` will automatically distributes the work load over multiple compute nodes if you have a license for their [Distributed Computing Server](http://www.mathworks.com/products/distriben/).

To run the code:

1. Modify the parameters of the model in `ED_01_sparse.m`
2. Open the Linux command line and `cd` to the directory containing that `m` file.
3. Execute the `matlab` command as shown in the example in the `EDScript01_sparse` shell script. If you are on a cluster that uses the [PBS job scheduler](http://www.arc.ox.ac.uk/content/pbs), you may find the scheduler configuration options included as comments in the beginning of that file helpful.

The output is a `mat` file containing 2 square matrices whose entries are the values of the correlation functions for the spin-up and spin-down sectors.

This repository includes a full suite of unit tests. All files whose names start with *test_* are unit test files that can be run run with [xUnit Test Framework](http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework/content/matlab_xunit_3_1_1/xunit/assertEqual.m). For a quick look at what the output look like for a given input, see the `assert` statements in the `test_unequal_time_gf_6_sites.m` file.


