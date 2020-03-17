COVID-19
This repository contains the necessary information to re-create Figure 1 and Tables 3, 4, and 5 from the pre-print publication Estimating the scale of COVID-19 Epidemic in the United States: Simulations Based on Air Traffic directly from Wuhan, China doi: https://doi.org/10.1101/2020.03.06.20031880.

Getting Started

Step 1
Clone this source repository and enter the project directory.

git pull dalinlee/COVID19_US

Step 2

Modify parameters of interest in the 1_model_and_sim_params.R file. 

Step 3

Modify parameters in run_file.R and run run_file.R.  Please modify the mc.cores parameter to match the number of cores and key parameters you wish to assign for this analysis. You might need to change the permission of the file prior to running it so it can be executable.



#################

License
Creative Commons License
The analytical work is licensed under a Creative Commons Attribution 4.0 International License. The source code is licensed under a BSD 3-Clause License. Copyright (c) 2020, Dalin Li & Gregory Botwin
