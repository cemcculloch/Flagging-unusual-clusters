clear all
set more off

*Change line below to directory with Excel file
cd "C:\Users\CMcCulloch\Documents\Grants\BLUP R01 application\Flagging methods\CABG data"

import excel "CABG-Surgeon-Summary-Table-2015-2016-Final.xlsx", sheet("CCORP Surgerons Outcome Summary") firstrow clear

table Hospital
keep if Hospital=="Surgeon Overall"

codebook
table E

*Clean up data
drop F
rename Oper Cases_deaths
gen Risk_adj_rate=real(D)
drop D
rename E Rating
codebook

drop if Rating=="NA"

tabulate Rating


