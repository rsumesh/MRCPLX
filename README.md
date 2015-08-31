This file contains information related to installation and execution of routines that are part of the MRCPLX github repository.
MRCPLX repository contains MATLAB scripts and functions and some data related to various techniques that estimate fundamental MRI parameters, namely, T1, T2 and T2*. Some of the routines implemented here work with magnitude data and others work with complex-valued data. The details regarding these methods have been presented in the article
Improved estimation of MR relaxation parameters using complex-valued data, S Umesh Rudrapatna, C J G Bakker, M A Viergever, A van der Toorn and R M Dijkhuizen

All the matlab files distributed under MRCPLX were developed and tested under MATLAB 8.3.0.532 (R2014a), running under Linux (Ubuntu 14.04), on a HP Z420 workstation with 4-core x86_64 architecture and 12 GB RAM.
These files were used to generate the data (in case of simulations) and obtain results reported in the above referenced article.

#Installation and use:
Download the matlab_code and data directories in the MRCPLX repository to your location of choice. All MATLAB files (in matlab_code directory) that begin with name sc_* are script files and can be run from the MATLAB editor or command window. The rest *.m files are functions called by these scripts. A brief overview of the functionality of the various MATLAB scripts is given below.
1. sc_generate_simulated_data.m: Generates simulated data on which various estimators can be tested (sample data generated from this script have been provided in the data directory). This script generates simulated data files with names t2s_data_sim.mat (simulated T2* data), t2_data_sim.mat (T2 data) and t1_data_sim.mat (T1 data).
2. sc_compare_r2s_fits_simulations.m: Script that is used to generate results from various T2* estimators.
3. sc_compare_r2_fits_simulations.m: Script that is used to generate results from various T2 estimators.
4. sc_compare_r1_fits_simulations.m: Script that is used to generate results from various T1 estimators.
See each of these individual files (2-4) for details regarding the various techniques available for fitting respective data.
5. sc_benchmark_real_t1.m: This script estimates the execution times of various T1 estimation techniques reported in the article. This uses t1_benchmarking_data.mat from the data folder for the analysis.
6. sc_benchmark_real_t2.m: This script estimates the execution times of various T2 estimation techniques reported in the article. This uses t2_benchmarking_data.mat from the data folder for the analysis.
7. sc_benchmark_real_t2s.m: This script estimates the execution times of various T2* estimation techniques reported in the article. This uses t2s_benchmarking_data.mat from the data folder for the analysis.
8. sc_benchmark_t2s_multi_rx.m: This script estimates the execution times of two T2* estimation techniques that mimic usage of multi-receive datasets. This uses t2s_mrx_benchmarking_data.mat from the data folder for the analysis.
9. sc_benchmark_t1_multi_rx.m: This script estimates the execution times of two T1 estimation techniques that mimic usage of multi-receive datasets. This uses t1_mrx_benchmarking_data.mat from the data folder for the analysis.

Please modify the paths in these files as you feel necessary (especially the '/' and '\' operators according to your operating system).

#License:
Please read the associated GPL V3 license.

#Contact information:
S. Umesh Rudrapatna
e-mail: umeshrs at gmail.com

