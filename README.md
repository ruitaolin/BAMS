# BAMS
R codes for implementing the Bayesian Adaptive Model Selection (BAMS) Design.
# Implementation details:
The function is developed using Rcpp:
* The file "AMS-cal2.cpp" is used to generate prior samples under each model. This file is should be included in both R files.
* The file "BAMS-fixed.r" includes the implementation code for the fixed-threshold BAMS design.
* The file "BAMS-adaptive.r" includes the implementation code for the adaptive-threshold BAMS design.
