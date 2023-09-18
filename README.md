# Flagging-unusual-clusters
Repository of code and examples for flagging unusual clusters manuscript
This repository contains examples of code, datasets and analyses in support of the manuscript "Flagging unusual clusters in linear mixed
models using weighted and self-calibrated predictors" by Charles E. McCulloch, John M. Neuhaus and Ross D. Boylan. 
The following files are contained herein:
1. CABG-Surgeon-Summary-Table-2015-2016-Final.xlsx:  An Excel file containing the overall risk-adjusted operative mortality rates for 
coronary artery bypass surgery for each of 276 surgeons in the state of California in 2017-2018.
2. A Stata .do file (analyze surgeon data.do) that reproduces the descriptive statistics noted in the Introduction using the Excel file.
3. A Stata .do file (Simulate incorrect flag rate calibrated predictors tau 1645 alpha 5 v3 gh.do) the gives example code for calculating the incorrect flagging rates.
4. A Stata .do file (Simulate correct flag rate calibrated predictors tau 1645 alpha 5 v3 gh.do) the gives example code for calculating the incorrect flagging rates.
5. A Stata .do file (Calibrate lambda for z_AB general v3 gh.do) that illustrates how the values of \lambda^* in Table 3 were derived.  This code is for \lambda^*_{AB} for the predictor \tilde{z}_{AB}.
6. A Stata data file (PIPA data.dta) that gives the data used in the Example in Section 10 of the manuscript.
7. A Stata .do file (PIPA flagging analysis v4 gh.do) that produces the results in Table 4 based on the PIPA datafile.
8. A Stata .do file (Gauss Hermite calibrate lambda for each R z_SQ v2.do) that solves for the value of lambda to self-calibrate the absolute weighted predictor using semi-infinite Gaussian quadrature. 
