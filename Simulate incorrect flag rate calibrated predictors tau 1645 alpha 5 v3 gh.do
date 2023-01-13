*Simulation of calibrated predictors

clear all
set seed 31415

*Set cutoff for flagging, alpha value for incorrect flagging rate (times 100) and number of simulations
local tau=1.645
local alpha=5
local nreps=5000000

*Set values of delta from Table 2 of manuscript
local del_FX_uc=`alpha'/100 
local del_FX_dc=0.075
*Set values of \lambda^* from Table 3 of manuscript
local lam_SQ=0.35
local lam_AB=1.25
local lam_CT=1.11

local tau1000=1000*`tau'

local run_nme="Incorrect flagging rates for calibrated predictors tau `tau1000' alpha `alpha' v3"

set more off

cd "C:\Users\cmcculloch\Documents\Grants\BLUP R01 application\flagging methods\simulations"
capture log close
log using "`run_nme'", replace

capture program drop flag_prob

*Define program for calculating flagging probability
#delimit;
program define flag_prob, rclass;
        version 16;
        syntax [, 
		nreps(integer 500) rval(real 1) tau(real 2) 
		lam_AB(real 1) lam_CT(real 1) lam_SQ(real 1)
		del_FX_uc(real 0.05) del_FX_sc(real 0.05) del_FX_dc(real 0.05)
		 ];
#delimit cr;

clear
set obs `nreps'

*True values of random effects
gen z=rnormal(0,1)
*Only interested in incorrect flagging rates, so only keep those < \tau
keep if z<`tau'

*Error part
*Recall that R=\sigma^2_e/n/\sigma^2_u
*This is the model to estimate the standardized
*portion of the random effect
local sqrt_R=sqrt(`rval')
gen err=rnormal(0,1)*`sqrt_R'

*Calculate predicted values
gen z_FX=z+err
gen z_BP=z_FX/(1+`rval')
gen z_SQ=z_FX/(1+`rval'*(1-2*`lam_SQ'))

*z_CT calculation
gen A_minus_CT=(-`lam_CT'-z_BP)/sqrt(`rval'/(1+`rval'))
gen A_plus_CT =(`lam_CT'-z_BP)/sqrt(`rval'/(1+`rval'))
gen z_CT=z_BP+sqrt(`rval'/(1+`rval'))/sqrt(2*_pi)*(exp(-A_plus_CT*A_plus_CT/2)-exp(-A_minus_CT*A_minus_CT/2))/(normal(A_minus_CT)+1-normal(A_plus_CT))

*z_AB calculation
gen del=sqrt(`rval'/(1+`rval'))
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

*Flagging rules
*Self flagging
foreach typ in BP SQ AB CT {
    gen flg_`typ'=z_`typ'>`tau'
}
*z_FX rules
gen ucMSEP_FX=`rval'
gen flg_FXuc=z_FX+invnormal(`del_FX_uc')*sqrt(ucMSEP_FX)>`tau'

gen scMSEP_FX=`rval'/(1+`rval')+z_FX^2*`rval'/(1+`rval')*`rval'/(1+`rval')
*delta is fixed for the singly conditioned
gen flg_FXsc=z_FX-1.001*sqrt(scMSEP_FX)>`tau'

*Calculate dcMSEP for z_FX
gen argM=(`tau'-z_FX/(1+`rval'))/sqrt(`rval'/(1+`rval'))
gen M=normalden(argM)/(1-normal(argM))
*Computations break down for very large argM
*e.g., argM=8.  This is E[z|z>argM] and hence must be larger than argM
*Asymptotically, the value is argM.
replace M=argM if argM>7
gen p3=`rval'/(1+`rval')*(1+`tau'*M*sqrt((1+`rval')/`rval'))
gen p1_FX=z_FX^2*4*0.5^2*`rval'^2/(1+`rval'*(1-2*0.5))^2/(1+`rval')^2
gen p2_FX=-z_FX*M*(1+`rval'*(1+2*0.5))/(1+`rval'*(1-2*0.5))/(1+`rval')*sqrt(`rval'/(1+`rval'))

gen dcMSEP_FX=p1_FX+p2_FX+p3
gen flg_FXdc=z_FX+invnormal(`del_FX_dc')*sqrt(dcMSEP_FX)>`tau'


end

*Start dataset for storing values
set obs 0
gen tau=.
foreach typ in BP SQ AB CT FXuc FXsc FXdc {
		gen flg_prob_`typ'=.	
}

save "Results `run_nme'", replace

*Evaluate for a grid of values of R
forvalues R=0.1(0.1)4 {
qui: flag_prob, tau(`tau') rval(`R') nreps(`nreps') lam_AB(`lam_AB') lam_CT(`lam_CT') lam_SQ(`lam_SQ') del_FX_dc(`del_FX_dc') del_FX_uc(`del_FX_uc')
di "Working on R=`R'"
foreach typ in BP SQ AB CT FXuc FXsc FXdc {
		summ flg_`typ'
		gen flg_prob_`typ'=r(mean)
}
gen tau=`tau'
gen R=`R'
keep tau R flg_prob* 
keep in 1
append using "Results `run_nme'"
save "Results `run_nme'", replace
}


log close
