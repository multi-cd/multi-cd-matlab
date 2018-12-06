# Multi-CD-Matlab #

**Multi-CD: Multi-scale discovery of Chromatin Domains**

Matlab code for a flexible algorithm to determine chromosome domains solutions that capture the correlation patterns in the Hi-C data.

The code is intended to accompany a manuscript under review. The link to the reference paper will be added later, as well as a supplementary material with details of the algorithm. For inquiries in the meantime, please send an email to Ji Hyun Bak (jhbak@kias.re.kr).

## Installation

You can do one of the following to obtain the latest code package.

* **Download**: click to download a zipped archive  [multi-cd-matlab-master.zip](https://github.com/multi-cd/multi-cd-matlab/archive/master.zip)
* **Clone**: clone the repository by typing this to the command line: 
```git clone https://github.com/multi-cd/multi-cd-matlab.git```

To run the code, you need Matlab with the Statistics toolbox.
Code was developed in Matlab R2016b, but we believe it should run in earlier releases of Matlab, as well. Please feel free to let us know if you encounter any problem.


## Documentation

We provide two example scripts to demonstrate how Multi-CD can be applied, with step-by-step tutorials.


### Pre-processing (Hi-C to correlation matrix)

`demo1_prepCorr_HiC.m` - demonstrates the pre-processing, going from the Hi-C matrix to the correlation matrix.

- Data files are not included in this repository. You need to provide the original Hi-C data yourself to run this script -- we are sorry for this inconvenience.
Hi-C datasets can be easily downloaded from the public databases; 
alternatively, you can input your own Hi-C data matrix.


### Identification of domain solutions (Multi-CD)

`demo2_MultiCD.m` - demonstrates how Multi-CD works at a fixed lambda (which can be changed by user), by showing each step of the simulated annealing process.

- By default, this demo script runs with a model correlation matrix, generated by a multi-scale group model. Alternatively, if you have already run `demo1` to get a "real" correlation matrix from Hi-C, you can provide its output (the correlation matrix) as the input here.

