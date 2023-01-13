*Kaiser data analysis
*v3 corrected z_FX_dc
*Set to 10% flagging

clear all
set more off
*Set cutoff for validation data at 6 months into the year
local mth_ct=6


cd "C:\Users\CMcCulloch\Documents\Grants\BLUP R01 application\Suni Kaiser data"
local rnname="PIPA flagging analysis v4"


capture log close
log using "`rnname'", replace
use "PIPA data", replace


*Transform LOS
gen log_los=log(los+1)

*Check distribution
*histogram log_los
*graph export "log_los histogram.pdf", as(pdf) replace

*Calculate sample sizes
bysort hospid_in: gen n_i=_N
*In early period
gen cnt_e=1
replace cnt_e=0 if mth_in>`mth_ct'
replace cnt_e=0 if los==.
bysort hospid_in: egen n_i_erly=total(cnt_e)
*and late period
gen cnt_l=1
replace cnt_l=0 if mth_in<=`mth_ct'
replace cnt_l=0 if los==.
bysort hospid_in: egen n_i_late=total(cnt_l)
*Drop if there is no data either early or late
drop if n_i_erly==0|n_i_late==0

*Mixed model fit using first 6 months of data
mixed log_los if mth_in<=`mth_ct' || hospid_in:
predict prd_BP, reffects
*Extract variance components for simple fit
local sig2u=exp(2*_b[lns1_1_1:_cons])
local sig2e=exp(2*_b[lnsig_e:_cons])
*Calculate intraclass correlation coefficient
local ICC=`sig2u'/(`sig2u'+`sig2e')
di "ICC estimated to be `ICC'"
*Calculate values of R
gen R_i=`sig2e'/n_i_erly/`sig2u'

*Fixed effects fits based on first 6 months
gen y_f=log_los
replace y_f=. if mth_in>`mth_ct'
bysort hospid_in: egen mean_y=mean(y_f)
*Calculate fixed effects predictions
gen temp=mean_y
*Make sum to 0 with equal weighting per cluster
bysort hospid_in: replace temp=. if _n>1
qui summ temp
local tmean=r(mean)
gen prd_FX=mean_y-`tmean'

*Check on reproducing the z_BP from z_FX
gen prd_BP_frm_FX=prd_FX/(1+R_i)

*Mean value in first and second half of year
*Second half
gen log_los_late=log_los
replace log_los_late=. if mth_in<=`mth_ct'
*First half
gen log_los_erly=log_los
replace log_los_erly=. if mth_in>`mth_ct'

*Make early data missing for late length of stay data
gen los_late=los
replace los_late=. if mth_in<=`mth_ct'
*and the converse
gen los_erly=los
replace los_erly=. if mth_in>`mth_ct'

summ los_erly, detail

*Collapse down to one observation per subject
collapse (first) n_i_erly n_i_late n_i R_i (mean) prd* log_los* mean_los_late=los_late mean_los_erly=los_erly (median) median_los_late=los_late median_los_erly=los_erly  , by(hospid_in)

foreach prd in FX BP {
    gen z_`prd'=prd_`prd'/sqrt(`sig2u')
}

*From calibration at \tau=1.28, Table 3 of manuscript 
local lam_SQ=0.25
local lam_AB=0.72
local lam_CT=0.61

*Calculate predictors
gen z_SQ=z_FX/(1+R_i*(1-2*`lam_SQ'))

*z_CT calculation
gen A_minus_CT=(-`lam_CT'-z_BP)/sqrt(R_i/(1+R_i))
gen A_plus_CT =(`lam_CT'-z_BP)/sqrt(R_i/(1+R_i))
gen z_CT=z_BP+sqrt(R_i/(1+R_i))/sqrt(2*_pi)*(exp(-A_plus_CT*A_plus_CT/2)-exp(-A_minus_CT*A_minus_CT/2))/(normal(A_minus_CT)+1-normal(A_plus_CT))

*z_AB calculation
gen del=sqrt(R_i/(1+R_i))
gen A_plus_AB=-(z_BP-`lam_AB'*del*del)/del
gen A_minus_AB=-(z_BP+`lam_AB'*del*del)/del
#delimit ;
gen z_AB=
	((-del*normalden(A_plus_AB)+(z_BP-`lam_AB'*del*del)*normal(A_plus_AB))*exp(-z_BP*`lam_AB')+
	(del*normalden(A_minus_AB)+(z_BP+`lam_AB'*del*del)*(1-normal(A_minus_AB)))*exp(z_BP*`lam_AB'))
	/
	(normal(A_plus_AB)*exp(-z_BP*`lam_AB')+
	(1-normal(A_minus_AB))*exp(z_BP*`lam_AB'))
	;
#delimit cr;

*Lower bounds of prediction intervals
gen LB_BP=z_BP
gen LB_SQ=z_SQ
gen LB_CT=z_CT
gen LB_AB=z_AB

gen LB_FX_uc=z_FX-invnormal(0.05)*sqrt(R_i)

gen scMSEP_FX=R_i/(1+R_i)+z_FX^2*R_i/(1+R_i)*R_i/(1+R_i)
gen LB_FX_sc=z_FX-invnormal(.186)*sqrt(scMSEP_FX)

*Calculate dcMSEP for z_FX
*Note - need to flip all the z_FX to get smaller MSEP for negative values
gen mz_FX=-z_FX
****NOTE:  tau hard-coded
local tauc=1.28
gen argM=(`tau'-mz_FX/(1+R_i))/sqrt(R_i/(1+R_i))
gen M=normalden(argM)/(1-normal(argM))
*Computations break down for very large argM
*e.g., argM=8.  This is E[z|z>argM] and hence must be larger than argM
replace M=argM if argM>7
gen p3=R_i/(1+R_i)*(1+`tauc'*M*sqrt((1+R_i)/R_i))
gen p1_FX=z_FX^2*4*0.5^2*R_i^2/(1+R_i*(1-2*0.5))^2/(1+R_i)^2
gen p2_FX=-mz_FX*M*(1+R_i*(1+2*0.5))/(1+R_i*(1-2*0.5))/(1+R_i)*sqrt(R_i/(1+R_i))

gen dcMSEP_FX=p1_FX+p2_FX+p3

gen LB_FX_dc=z_FX-invnormal(0.068)*sqrt(dcMSEP_FX)

*Check out using z_FX by itself
gen LB_FX_nc=z_FX

*Calibrated for bottom 10%, so flag bottom 7 by mean LOS late
gen mean_los_late_low=mean_los_late<=25

summ R_i, detail

*Calculate flagging for all predictors and various cutoffs
foreach tau100 in 200 150 128 100 {
    local tau=-`tau100'/100
    foreach prd in BP FX_uc FX_sc FX_dc FX_nc CT AB SQ {
	    gen flg_`tau100'_`prd'=(LB_`prd'<`tau')&(LB_`prd'<.)
	}
}

drop if mean_los_late==.

*Generate values reported in Table 4
foreach tau100 in 128 {
    foreach prd in BP FX_uc FX_sc FX_dc CT AB SQ {
	    table flg_`tau100'_`prd' mean_los_late_low
}
}

save "BPs for `rnname'", replace

log close
