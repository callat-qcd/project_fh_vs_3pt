# This directory exists to calculate the ground-state energy and the $g\_{A00}$ and $g\_{V00}$ values of a proton based on nonlinear least-squares fits to simulation data obtained from Lattice QCD at a lattice spacing 0.12 fm
## Initial pion mass used: 310 MeV, initial lattice length 24 (both of those are intended to be adjustable)

Data files
- a12m310_a_avg.h5: the primary HDF5 file with the 2-point and 3-point correlation function data
- a12m310_a_fh.h5: the HDF5 file with the Feynman-Hellman data

Libraries
- h5data_methods.py: containes functions for getting data from files and plotting data and fitfunctions

Code-files
- Model_comparator: compares the results for ground-state energy, $g\_{A00}$, and $g\_{V00}$ for the different fit-models made by the files below
- Parameter_charts: For each fit-model used below, creates an energy-spectrum of the prior and posterior and compares it to the pion-excitation energies
- 2point_final: for ascertaining that 2-point data from the a12m310_a_avg.h5 and the a12m310_a_fh.h5 files are the same, doing fits on 2-point data alone, and seeing what happens to ground state energy and fit quality when number of excited states and minimum and maximum time of data is changed
- 2pt_FH_fits: does a simulaneous 2-point and sum-subtracted fit
- chained_2pt_FH: does a chained fit of 2-point and sum-subtracted data
All data files below compare results for fits where the 3-point data has all points between the creation and annihilation point and where the points neighboring the creation and annihilation point are also removed, results where the Feynman-Hellman data is not included, and where it is included. Also, all except All_fits_combined (due to time constraints) test what happens to fit quality, ground-state energy,  g\_{A00}, and g\_{V00} as the number of states used, the amount of points cut from the $\tau$ edges, and the width of the first excited state energy prior are vaired
- All_fits_combined: for doing a simultaneous fit of 2-point, sum-subtracted, and 3-point data
- chainedfits_2pt_FH_3pt: for doing a chained fit of 2-point, sum-subtracted, and 3-point data (in that order)
- chainedfits_first: for doing a chained fit of 2-point, 3-point, and sum-subtracted data (in that order)
- 2pt_3pt_fits_combined: for doing simultaneous fits of 2-point and 3-point data
- chained fit_2pt_3pt: for doing chained fits of 2-point and 3-point data
- singlefits_3pt_FH.ipynb: for doing 3-point and Feynman-Hellman fits individually

Folders with other data
- energy_level_calculator: a program that uses a brute-force algorithm to find all possible energy levels based on pion excitation up to a certain number of pions and a certain maximum particle momentum (warning: it does not eliminate unphysical zero-momentum states with odd numbers of pions)
- old: files that are not currently used. Includes (but is not necessarily limited to) Energy_Levels (for visually representing the energy levels of energy_level_calculator as an energy-level diagram), covar_matrices (for examining the covariance matrices of priors and posteriors in fits to the data), Higher_Order_Effective_States (an abandoned attempt at prioring excited state energies by subtracting out lower order energies and finding the effective mass), First_axial_dataset (a first attempt at fitting the 2-point data, using almost every single energy level up to about 1.68 obtained by energy_level_calculator), 3_pt_sum_subtracted_trial1 (a trial run at plotting both 2-point and 3-point data), 2point_model_comparison (a version of 2point_final but with no data used from a12m310_a_fh.h5 and several abandoned projects), and 3pt_data_A3 and 3pt_data_V4 (test files for the fits to the 3-point data)



