This repository contains MATLAB code to reproduce the results in:
   Vattiato G, Maclaren O, Lustig A, Binny RN, Hendy SC, Plank MJ (2022) An assessment of the potential impact of the Omicron variant of SARS-CoV-2 in Aotearoa New Zealand, Infectious Disease Modelling, Volume 7, Issue 2, Pages 94-105, https://doi.org/10.1016/j.idm.2022.04.002.

README for SimLeaky_Omicron_waning

This folder contains:
- main.m: the main Matlab run file, contains the different “scenario” parameters
- getParOmiWane.m: a Matlab function generating the remaining parameters, called by main.m
- runSimWaning.m: Matlab file that runs the branching process model, called by main.m
- README.txt: this README

- data: folder containing spreadsheets called by the model to generate the contact matrix and the vaccination coverage time series
- dependencies: folder containing Matlab function dependencies of the main files
- plots: folder where plots are automatically saved by main
- summaries: folder where results summaries are automatically saved by main
- timeseries: folder where the timeseries for infections, cases, hospitalisations, cumulative deaths, and hospital beds occupied are automatically saved by main

To run:
1. In the main.m file, change the variable “savefolder”, if required, with the name of the folder where the timeseries will be saved. Make sure a folder with that name exists within the timeseries folder
2. In the main.m file, change the variable “nReps” with the number of simulation repetitions required
3. In the main.m file, change any other parameters as required. If left as it is, the program will run 100 simulation repetitions of the “baseline scenario” described in the paper
4. Run the main.m file
5. Results will be saved in the summaries, timeseries and plots folders
