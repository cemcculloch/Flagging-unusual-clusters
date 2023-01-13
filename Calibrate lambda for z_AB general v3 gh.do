*Simulation to find lambda such that the incorrect flagging rate is 0.05
clear all
set more off

local run_nme="Calibrate lambda for z_AB general v3"

*Set number of replications
local nreps=1000000
*Set minimum and maximum values of R to search
local Rmin=0.1
local Rmax=10
local Rstep=0.1
*After a coarse grid search, these are the ones that bracket the needed values of \lamba^2 (times 1000)
local lam1000_values="700 710 720 730 740 750 1190 1200 1210 1220 1230 1240 1250 1630 1640 1650 1660 1670 1680 1700 2130 2140 2150 2170 2200 2450 2460 2470 2500"

cd "C:\Users\cmcculloch\Documents\Grants\BLUP R01 application\flagging methods\simulations"
capture log close
log using "`run_nme'", replace

set obs `nreps'
*Set up blank columns to store results
gen R=.
gen tau=.
gen lam=.
gen flag_rate=.
gen flag=.
gen z_AB=.
gen z_BP=.
gen A_plus_AB=.
gen A_minus_AB=.

*Ingredients for calculating predictor
gen z_FX=.

*Generate true z values
gen z=rnormal(0,1)
*And error portion
gen z0=rnormal(0,1)

local i=1

foreach tau in 1.28 1.645 2.33  {
*Calculate incorrect flagging rates
qui {
forvalues R=`Rmin'(`Rstep')`Rmax' {
	noisily di "Working on tau=`tau' and R=`R'"
	replace z_FX=(z+sqrt(`R')*z0)
	local del=sqrt(`R'/(1+`R'))
	foreach lam1000 in `lam1000_values' {
	local lam_AB=`lam1000'/1000
	replace z_BP=z_FX/(1+`R')
    replace A_plus_AB=-(z_BP-`lam_AB'*`del'*`del')/`del'
    replace A_minus_AB=-(z_BP+`lam_AB'*`del'*`del')/`del'
	#delimit ;
    replace z_AB=
	((-`del'*normalden(A_plus_AB)+(z_BP-`lam_AB'*`del'*`del')*normal(A_plus_AB))*exp(-z_BP*`lam_AB')+
	(`del'*normalden(A_minus_AB)+(z_BP+`lam_AB'*`del'*`del')*(1-normal(A_minus_AB)))*exp(z_BP*`lam_AB'))
	/
	(normal(A_plus_AB)*exp(-z_BP*`lam_AB')+
	(1-normal(A_minus_AB))*exp(z_BP*`lam_AB'))
	;
#delimit cr;
		*One tailed flagging
		replace flag=z_AB>`tau'
		summ flag if z<`tau'
		replace flag_rate=r(mean) in `i'
		replace R=`R' in `i'
		replace lam=`lam1000'/1000 in `i'
		replace tau=`tau' in `i'
		local i=`i'+1
	}
}
}
}

keep if R~=.
keep R lam tau flag_rate

*Find maximum value of incorrect flagging rate by lambda
bysort lam tau: egen max_flag_rate=max(flag_rate)

*Only keep the value with the maximum
keep if abs(flag_rate-max_flag_rate)<0.0001

di "Maximum (across values of R) of incorrect flagging rate for tau=`tau' for various lambda"

sort tau lam

*Generate warning if max is at boundary
gen warning=""
replace warning ="WARNING" if abs(R-`Rmin')<0.2|abs(R-`Rmax')<0.2

list R tau lam max_flag_rate warning

log close
