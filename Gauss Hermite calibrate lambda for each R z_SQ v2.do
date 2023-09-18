*Given inputs of \alpha and \tau, attempt to self-calibrate for a given R

clear all

local rnname="Compare given vs all R z_SQ v2"

*Bounded Gauss-Hermite with 15 quadrature points
*Read in quadrature weights and abscissae
import excel "C:\Users\CMcCulloch\Documents\Grants\BLUP R01 application\Flagging methods\GH quadrature\bounded GH quadrature points.xlsx", sheet("Pages from gauss hermite bounde") firstrow

*Rename weights and abscissae
forvalues i=1/15 {
	global gh_w_`i'=weight[`i']
	global gh_x_`i'=abscissae[`i']
}

di "Global first weight is $gh_w_1"

*Make space to save results
set obs 100

*Write function to return incorrect flagging probability
*as a function of lambda, tau, and R

capture program drop inc_flag_SQ

*Define a program to calculate the incorrect flagging rate for z_SQ
*using bounded GH quadrature
#delimit;
program define inc_flag_SQ, rclass;
        version 16;
        syntax [, 
		lam_SQ(real .25) rval(real 1) tau(real 1.645)
			 ];
#delimit cr;

*Do calculation for each quadrature point
forvalues i=1/15 {
local inc_flag_`i'=${gh_w_`i'}*normal((`rval'*(1-2*`lam_SQ')*`tau'+sqrt(2)*${gh_x_`i'})/sqrt(`rval'))*exp(-`tau'^2/2+sqrt(2)*`tau'*${gh_x_`i'})/sqrt(_pi)
}

local integral=`inc_flag_1'+`inc_flag_2'+`inc_flag_3'+`inc_flag_4'+`inc_flag_5'+`inc_flag_6'+`inc_flag_7'+`inc_flag_8'+`inc_flag_9'+`inc_flag_10'+`inc_flag_11'+`inc_flag_12'+`inc_flag_13'+`inc_flag_14'+`inc_flag_15'

return scalar prob_inc_flag=1-`integral'/normal(`tau')

end

*Try function once
inc_flag_SQ, lam_SQ(0.35) rval(2) tau(1.645)
return list


capture program drop cor_flag_SQ

*Define a program to calculate the correct flagging rate for z_SQ
*using bounded GH quadrature
#delimit;
program define cor_flag_SQ, rclass;
        version 16;
        syntax [, 
		lam_SQ(real .25) rval(real 1) tau(real 1.645)
			 ];
#delimit cr;

*Do calculation for each quadrature point
forvalues i=1/15 {
local cor_flag_`i'=${gh_w_`i'}*normal((`rval'*(1-2*`lam_SQ')*`tau'-sqrt(2)*${gh_x_`i'})/sqrt(`rval'))*exp(-`tau'^2/2-sqrt(2)*`tau'*${gh_x_`i'})/sqrt(_pi)
}

local integral=`cor_flag_1'+`cor_flag_2'+`cor_flag_3'+`cor_flag_4'+`cor_flag_5'+`cor_flag_6'+`cor_flag_7'+`cor_flag_8'+`cor_flag_9'+`cor_flag_10'+`cor_flag_11'+`cor_flag_12'+`cor_flag_13'+`cor_flag_14'+`cor_flag_15'

return scalar prob_cor_flag=1-`integral'/(1-normal(`tau'))

end

*Try function once
cor_flag_SQ, lam_SQ(0.35) rval(1) tau(1.645)
return list

*Set tau and alpha
local tau=1.645
local alpha=0.05

*Do just once for now
local R=1

capture program drop self_calib_R

*Define a program to calculate the value of lambda that gives an incorrect flagging rate of alpha, 
*as a function of \alpha and \tau
#delimit;
program define self_calib_R, rclass;
        version 16;
        syntax [, 
		alpha(real 0.05) tau(real 1.645) rval(real 1)
			 ];
#delimit cr;

*Start with wide bounds for \lambda_SQ 
local lamsqlow=0
local lamsqhigh=0.5
*Set tolerance to high value
local tol=1

*Check to make sure range is feasible
inc_flag_SQ, lam_SQ(`lamsqlow') rval(`rval') tau(`tau')
local incflaglow=r(prob_inc_flag)
if 	`incflaglow'>`alpha' {
		di "Not feasible incflaglow"
		exit
	}

inc_flag_SQ, lam_SQ(`lamsqhigh') rval(`rval') tau(`tau')
local incflaghigh=r(prob_inc_flag)
di "High inc flagging is `incflaghigh'"
if 	`incflaghigh'<`alpha' {
		di "Not feasible incflaghigh"
		return scalar lam_SQ_calib=`lamsqhigh' 
	}
	
	
*Bisection loop to find root
while `tol'>1e-5 {
	local lamsqmid=(`lamsqhigh'+`lamsqlow')/2
	inc_flag_SQ, lam_SQ(`lamsqmid') rval(`rval') tau(`tau')
	local incflagmid=r(prob_inc_flag)
	if `incflagmid' > `alpha' {
		local lamsqhigh=`lamsqmid'
	}
	else {
		local lamsqlow=`lamsqmid'
	}
	local tol=abs(`lamsqhigh'-`lamsqlow')
}	

return scalar lam_SQ_calib=`lamsqlow'
end

local tau=1.645
local alpha=0.05

set obs 300
gen R=.
gen self_lam=.
gen self_inc_flag=.
gen self_cor_flag=.


local i=1

forvalues R=0.1(0.1)4{
self_calib_R, alpha(`alpha') tau(`tau') rval(`R')
local lam_SQ_calib=r(lam_SQ_calib)
inc_flag_SQ, lam_SQ(`lam_SQ_calib') rval(`R') tau(`tau')
local prob_inc_flag=r(prob_inc_flag)
replace R=`R' in `i'
replace self_lam=`lam_SQ_calib' in `i'
replace self_inc_flag=`prob_inc_flag' in `i'
local i=`i'+1
}

*Since increasing lambda increases the incorrect flagging rate, the previous
*work is based on the minimum value of R.  Create a column with that
gen min_lam=.
gen min_inc_flag=.
gen min_cor_flag=.
summarize self_lam
replace min_lam=r(min)

forvalues i=1/39 {
*Calculate incorrect flagging for min R
local min_lami=min_lam[`i']
local Ri=R[`i']
inc_flag_SQ, lam_SQ(`min_lami') rval(`Ri') tau(`tau')
local prob_inc_flag=r(prob_inc_flag)
replace min_inc_flag=`prob_inc_flag' in `i'
}

forvalues i=1/39 {
*Calculate correct flagging rates for both
local self_lami=self_lam[`i']
local min_lami=min_lam[`i']
local Ri=R[`i']
cor_flag_SQ, lam_SQ(`self_lami') rval(`Ri') tau(`tau')
local prob_cor_flag=r(prob_cor_flag)
replace self_cor_flag=`prob_cor_flag' in `i'
cor_flag_SQ, lam_SQ(`min_lami') rval(`Ri') tau(`tau')
local prob_cor_flag=r(prob_cor_flag)
replace min_cor_flag=`prob_cor_flag' in `i'
}

keep if R<.
keep R self_lam min_lam self_inc_flag min_inc_flag self_cor_flag min_cor_flag
save "`rnname'.dta", replace

list 

#delimit ;
twoway 
(connected self_inc_flag R, msymbol(X) mcolor(green) lcolor(green))
(connected min_inc_flag R, msymbol(X) mcolor(black) lcolor(black))
(connected self_cor_flag R, msymbol(circle) mcolor(green) lcolor(green))
(connected min_cor_flag R, msymbol(circle) mcolor(black) lcolor(black)),
title("Flagging rates for z_SQ, self-calibrate to R or overall")
ytitle("Flagging rate")
xtitle(" " "R=({&sigma}{sub:{&epsilon}}{sup:2}/n)/{&sigma}{sub:u}{sup:2}", margin(zero))
subtitle("For alpha of `alpha' and tau of `tau'")
legend(order(
1 "Specific R, incorrect" 2 "For all R, incorrect"  3 "Specific R, correct" 4 "For all R, correct") rows(2))
;
#delimit cr;
graph export "Compare self calibrate to specific R vs overall for alpha of `alpha' and tau of `tau'.pdf", as(pdf) replace


