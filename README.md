# Design for validation studies

This page contains all functions used in the paper "Two-Phase Sampling Designs for Data Validation in Settings with Covariate Measurement Error and Continuous Outcome". The paper compares different ways to select observations to have their records validated, aiming to reduce the variance of some parameter of interest while keepting the size of the validation study fixed.

We discussed sampling schemes designed for both likelihood and design-based estimators. For the former we applied multiple imputation (MI) and semiparametric maximum likelihood estimators (SPMLE) and for the latter we applied inverse probability weighting (IPW) and generalized raking. 

It includes:

  - R code for simulations
    - Generalized raking and IPW
    - MI and SPMLE
 
  - R code for case study:
    - Vanderbilt Comprehensive Care Clinic (VCCC)
