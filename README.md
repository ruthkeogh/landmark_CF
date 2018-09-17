This repository provides R code for implementation of methods referred to in the following paper:

"Dynamic prediction of survival in cystic fibrosis: A landmarking analysis using UK patient registry data" 
Ruth Keogh, Shaun Seaman, Jessica Barrett, David Taylor-Robinson, Rhonda Szczesniak. Epidemiology 2018. In Press.

mixoutsamp_v2.R: R code for obtaining out-of-sample predictions from a mixed model fitted using lme.

mixed_models_outofsample_prediction: Illustrating the use of mixoutsamp on a freely-available exmple data set. Rmd file and corresponding pdf file. 

estimated_survival_probabilities: R code for obtaining estimated survival probabilities using the dynamic prediction model developed in the above manuscript. Rmd file and corresponding pdf file. 

baseline_cumulative_hazards, times, log_hazard_ratios: csv and RData files containing estimated baseline cumulative hzards, event/censoring time and log hazard ratios from the dynamic prediction model developed in the above manuscript. These are used to obtain the estimated survival probabilities. 
