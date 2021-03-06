# Multi-CD-Matlab #

**Multi-CD: Multi-scale discovery of Chromatin Domains**

Matlab code for a flexible algorithm to determine chromosome domains solutions that capture the correlation patterns in the Hi-C data.

The code is intended to accompany a manuscript under review. The link to the reference paper will be added later, as well as a supplementary material with details of the algorithm. For inquiries in the meantime, please send an email to Ji Hyun Bak (jhbak@kias.re.kr), or [create an issue](https://github.com/multi-cd/multi-cd-matlab/issues/new).

## System Requirements

**Software Requirements:**
To run the package, you need Matlab with the Statistics toolbox.
Code was tested in Matlab versions R2016b-R2017b, on Mac OS (OSX 10.12) and Linux (Centos 7). 
We tried to use basic Matlab functionalities whenever possible, so that the code could be used in other versions of Matlab.
Please feel free to let us know if you encounter any problem in other versions of Matlab.

**Hardware Requirements:**
The package requires only a standard computer, with enough RAM to support the handling of the data matrix.


## Installation Guide

You can do one of the following to obtain the latest code package.
This is a lightweight package, so download should be quick in normal internet conditions. 

* **Download**: click to download a zipped archive  [multi-cd-matlab-master.zip](https://github.com/multi-cd/multi-cd-matlab/archive/master.zip)
* **Clone**: clone the repository by typing this to the command line: 
```git clone https://github.com/multi-cd/multi-cd-matlab.git```


To run the provided scripts, set your Matlab's working directory to the location where you downloaded this code package (the directory where the demo scripts are in). 
You can do this by typing 

```
cd [your-path]/multi-cd-matlab
```

on the Matlab command line. The relative path to the rest of the code is automatically set in each demo script, by calling `setpaths.m` in the beginning.


## Instructions for Use

We provide two example scripts to demonstrate how Multi-CD can be applied, with step-by-step tutorials.
Simply open a demo script in Matlab and run; preferably, run by sections to see the intermediate outcome of each step.
Example output plots from the demo scripts are included under [Data](Data).

All custom functions (that do the real job) are in the [Code](Code) folder.
If you find a function call in the demo and want to see a quick documentation of what it does, you can type 

```
help [function-name]
```

on the Matlab command line. For example, `help HS_calculation_all`. Note that `setpaths.m` should have been run in advance, if not already called in the beginning of the demo script.


### demo1: Pre-processing (Hi-C to correlation matrix)

`demo1_prepCorr_HiC.m` - demonstrates the pre-processing, going from the Hi-C matrix to the correlation matrix.

- If you are only interested a demo of how our Multi-CD algorithm works, you can skip this step and proceed directly to `demo2`. 
- Pre-processing is necessary if you want to start from a real Hi-C dataset. 
No Hi-C dataset is included in this package, but you can download one from a public database. 
See [Data/readme](Data/readme.md) for more information.

<!--- ![Example output from demo1](Data/example-output-demo1.png) --->
<img src="Data/example-output-demo1.png" alt="Example output from demo1" width="600">


### demo2: Identification of domain solutions (Multi-CD)

`demo2_MultiCD.m` - demonstrates how Multi-CD works at a fixed lambda (which can be changed by user), by showing each step of the simulated annealing process.

- By default, this demo script runs with a simulated correlation matrix, generated by a multi-scale group model (see a separate [demo script](Code/gen_model/demo_gen_corrMat.m) for how it is generated). On a typical desktop computer, it takes only a few seconds to run this demo as provided. Feel free to try different values of `lambda`, to find domains at different scales.

<!--- ![Example output from demo2](Data/example-output-demo2.png) --->
<img src="Data/example-output-demo2.png" alt="Example output from demo2" width="600">


### demo3: Testing Multi-CD on synthetic data

`demo3_test_synth.m` - demonstrates how Multi-CD can be applied to a synthetic dataset (a model correlation matrix) with varying lambda, to recover the underlying domains that generates the correlation pattern in the data.

- Again, see a separate [demo script](Code/gen_model/demo_gen_corrMat.m) for more details of model correlation matrix generation.

<!--- ![Example output from demo3](Data/example-output-demo3.png) --->
<img src="Data/example-output-demo3.png" alt="Example output from demo3" width="600">


### Application to real data

The demo scripts can be directly used to perform a Multi-CD analysis on real data.

- Prepare Hi-C data in a compatible `.mat` file format (see 
[Data/readme](Data/readme.md) for more information). Provide it as the input to `demo1`, and run to get a correlation matrix. Save the resulting correlation matrix to a file. This can be done by uncommenting the last lines in `demo1`. 
- The correlation matrix is now used as the input to `demo2`. To do this, change the data import option to `useRealData=true` in `demo2`, and make sure the file path/name is correctly specified.
- Options and parameters for the Multi-CD algorithm can be adjusted easily by changing the numbers in `demo2`. We tried to keep all parameter-setting explicit at the script level, so that the user does not need to look into the contents of the package.
- Run `demo2` repeatedly at different values of the Multi-CD scale parameter `lambda` (or modify the script to add a loop over `lambda`), to find the multi-scale domain solution.
    Feel free to adapt from `demo3` to set up the loop.

The run time now depends on the size of your data, on your choice of temperature schedule, and/or how many different values of `lambda` you wish to try.
See the [Configs](Configs) directory for a full list of parameters that was used for the main analysis.
