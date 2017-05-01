Overall_Script: Organizes the estimation step and maximization step functions.
Initial_PH_Est_Step: Performs estimation of the initial proportional hazard coefficients under the assumption that K=1.
PH_Est_Step_K_5: Takes as input the parameter values estimated from Initial_PH_Est_Step and finds the updated values for the proportional hazards coefficients given current cluster parameters, Z, assuming K=5.
Multinom_Est_Step: Takes the current cluster parameters, Z, and updates the multinomial logistic regression coefficients.
Estep: Performs the expectation step of the EM algorithm given current values of the multinomial and proportional hazards parameters.

More info can be found in Supplemental Material H.