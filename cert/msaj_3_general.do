* Version 0.1, 11/02/2020

*****************************
*  			EBMT  			*
*****************************

cap program drop ebmt
program define ebmt
	use "./data/ebmt", replace
	mat tmat = (.,1,2,.,3,4\.,.,.,5,6,7\.,.,.,8,9,10\.,.,.,.,11,12\.,.,.,.,.,.\.,.,.,.,.,.)
end

cap program drop ebmt_cov
program define ebmt_cov
	use "./data/ebmt", replace
	mat tmat = (.,1,2,.,3,4\.,.,.,5,6,7\.,.,.,8,9,10\.,.,.,.,11,12\.,.,.,.,.,.\.,.,.,.,.,.)
	qui gen cov = 10
	qui replace cov = -3 if id > 500
	qui replace cov = 0 if id > 1000
	qui replace cov = . if id > 1500
end


** #1 Check output given when it should be

* #1a Check entries given for probabilities, se and ci when _d == 1 
ebmt
msaj, transmatrix(tmat) ci se
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' != . if _d == 1
	assert P_AJ_`i'_se != . if _d == 1
	assert P_AJ_`i'_lci != . if _d == 1
	assert P_AJ_`i'_uci != . if _d == 1
}

* #1b Check entries given when by used (not for missing)
ebmt_cov
msaj, transmatrix(tmat) by(cov) ci se
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' != . if _d == 1 & cov != .
	assert P_AJ_`i'_se != . if _d == 1 & cov != .
	assert P_AJ_`i'_lci != . if _d == 1 & cov != .
	assert P_AJ_`i'_uci != . if _d == 1 & cov != .
	
	assert P_AJ_`i' == . if cov == .
	assert P_AJ_`i'_se == . if cov == .
	assert P_AJ_`i'_lci == . if cov == .
	assert P_AJ_`i'_uci == . if cov == .
}

* #1c Check entries given for all observations for LOS (before last event time)
ebmt
msaj, transmatrix(tmat) los
sum _t if _d == 1
forvalues i=1/6 {
	di "State `i'"
	assert LOS_AJ_`i' != . if _t <= `r(max)'
}	

* #1d Check entries given for all observations for LOS and by
ebmt_cov
msaj, transmatrix(tmat) by(cov) los
foreach i in 10 -3 0 {
	sum _t if _d == 1 & cov == `i'
	forvalues j=1/6 {
		di "State `j'"
		assert LOS_AJ_`j' != . if _t <= `r(max)' & cov == `i'
	}
}	
forvalues i=1/6 {
	di "State `i'"
	assert LOS_AJ_`i' == . if cov == .
}	

* #1e Check probabilities and LOS before ltruncated / after exit are missing
ebmt
msaj , transmatrix(tmat) ltruncated(100) exit(500) los ci se
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' ==. if _t < 100 | _t > 500 
	assert P_AJ_`i'_se == . if  _t < 100 | _t > 500 
	assert P_AJ_`i'_lci == . if _t < 100 | _t > 500 
	assert P_AJ_`i'_uci == . if _t < 100 | _t > 500 
	assert LOS_AJ_`i' ==. if _t < 100 | _t > 500 
}

* #1f Check probabilities and LOS before ltruncated / after exit are missing when by used
ebmt_cov
msaj , transmatrix(tmat) by(cov) ltruncated(100) exit(500) los ci se
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' ==. if (_t < 100 | _t > 500) & cov != . 
	assert P_AJ_`i'_se == . if (_t < 100 | _t > 500) & cov != . 
	assert P_AJ_`i'_lci == . if (_t < 100 | _t > 500) & cov != .
	assert P_AJ_`i'_uci == . if (_t < 100 | _t > 500) & cov != .
	assert LOS_AJ_`i' ==. if (_t < 100 | _t > 500) & cov != .
}



** #2 Common sense checks

* #2a Check probabilities add up to 1
ebmt
msaj , transmatrix(tmat) 
egen double testsum=rowtotal(P_AJ_*)
assert sqrt((testsum-1)^2) < $tol if P_AJ_1 != .

* #2b Check for by option
ebmt_cov
msaj , transmatrix(tmat) by(cov)
egen double testsum=rowtotal(P_AJ_*)
assert sqrt((testsum-1)^2) < $tol if P_AJ_1 != .

* #2c Check all probabilities and CI are within [0,1]
ebmt
msaj, transmatrix(tmat) ci
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' >=0 & P_AJ_`i' <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_lci >=0 & P_AJ_`i'_lci <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_uci >=0 & P_AJ_`i'_uci <=1 if P_AJ_`i' != .
}

* #2d Check for by option
ebmt_cov
msaj, transmatrix(tmat) by(cov) ci
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' >=0 & P_AJ_`i' <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_lci >=0 & P_AJ_`i'_lci <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_uci >=0 & P_AJ_`i'_uci <=1 if P_AJ_`i' != .
}

* #2e Check probabilities and LOS are 0 that can never be reached with from
ebmt
msaj , transmatrix(tmat) from(4) los ci se
sum _t if _d == 1
foreach i in 1 2 3 {
	di "State `i'"
	assert P_AJ_`i' == 0 if _d == 1
	assert P_AJ_`i'_se == 0 if _d == 1
	assert P_AJ_`i'_lci == 0 if _d == 1
	assert P_AJ_`i'_uci == 0 if _d == 1
	assert LOS_AJ_`i' == 0 if _t <= `r(max)'
}

* #2f Check LOS adds up to time
ebmt
msaj , transmatrix(tmat) los
egen double testsum = rowtotal(LOS*)
assert sqrt((testsum-_t)^2) < $tol if LOS_AJ_1! = .

* #2g Check for by option
ebmt_cov
msaj , transmatrix(tmat) by(cov) los
egen double testsum = rowtotal(LOS*)
assert sqrt((testsum-_t)^2) < $tol if LOS_AJ_1! = .

* #2h Check all answers the same if cov = 1 as for all
ebmt 
gen cov = 1
msaj, transmatrix(tmat) by(cov) ci se los
rename (P_AJ_*) (by_P_AJ_*)
rename (LOS_AJ*) (by_LOS_AJ*)
msaj, transmatrix(tmat) ci se los
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' == by_P_AJ_`i'
	assert P_AJ_`i'_se == by_P_AJ_`i'_se
	assert P_AJ_`i'_lci == by_P_AJ_`i'_lci
	assert P_AJ_`i'_uci == by_P_AJ_`i'_uci
	assert LOS_AJ_`i' == by_LOS_AJ_`i'
}
 


** 3 Check answers are the same with all the options compared to without

* #3a Overall
ebmt
msaj, transmatrix(tmat)
rename (P_AJ_*) (pP_AJ_*)
msaj, transmatrix(tmat) ci
rename (P_AJ_*ci) (cP_AJ_*ci)
drop P_AJ_*
msaj,  transmatrix(tmat) se
rename (P_AJ_*se) (sP_AJ_*se)
drop P_AJ_*
msaj, transmatrix(tmat) los
rename (LOS*) (lLOS*)
drop P_AJ_*
msaj, transmatrix(tmat) ci se los
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' == pP_AJ_`i'
	assert P_AJ_`i'_se == sP_AJ_`i'_se
	assert P_AJ_`i'_lci == cP_AJ_`i'_lci
	assert P_AJ_`i'_uci == cP_AJ_`i'_uci
	assert LOS_AJ_`i' == lLOS_AJ_`i'
}
 
* #3b With by option
ebmt_cov
msaj, transmatrix(tmat) by(cov)
rename (P_AJ_*) (pP_AJ_*)
msaj, transmatrix(tmat) by(cov) ci 
rename (P_AJ_*ci) (cP_AJ_*ci)
drop P_AJ_*
msaj,  transmatrix(tmat) by(cov) se
rename (P_AJ_*se) (sP_AJ_*se)
drop P_AJ_*
msaj, transmatrix(tmat) by(cov) los
rename (LOS*) (lLOS*)
drop P_AJ_*
msaj, transmatrix(tmat) by(cov) ci se los
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' == pP_AJ_`i'
	assert P_AJ_`i'_se == sP_AJ_`i'_se
	assert P_AJ_`i'_lci == cP_AJ_`i'_lci
	assert P_AJ_`i'_uci == cP_AJ_`i'_uci
	assert LOS_AJ_`i' == lLOS_AJ_`i'
}
 

 
*****************************
*  			HAI  			*
*****************************

cap program drop hai
program define hai
	use "./data/hai", replace
	mat tmat = (.,1,2,3,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
end


** #2 Common sense checks

* #2a Check probabilities add up to 1
hai
msaj , transmatrix(tmat) 
egen double testsum=rowtotal(P_AJ_*)
assert sqrt((testsum-1)^2) < $tol if P_AJ_1 != .

* #2b Check all probabilities and CI are within [0,1]
hai
msaj, transmatrix(tmat) ci
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' >=0 & P_AJ_`i' <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_lci >=0 & P_AJ_`i'_lci <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_uci >=0 & P_AJ_`i'_uci <=1 if P_AJ_`i' != .
}

* #2c Check probabilities and LOS are 0 that can never be reached with from
hai
msaj , transmatrix(tmat) from(2) los ci se
sum _t if _d == 1
foreach i in 1 3 4 {
	di "State `i'"
	assert P_AJ_`i' == 0 if _d == 1
	assert P_AJ_`i'_se == 0 if _d == 1
	assert P_AJ_`i'_lci == 0 if _d == 1
	assert P_AJ_`i'_uci == 0 if _d == 1
	assert LOS_AJ_`i' == 0 if _t <= `r(max)'
}

* #2d Check LOS adds up to time
hai
msaj , transmatrix(tmat) los
egen double testsum = rowtotal(LOS*)
assert sqrt((testsum-_t)^2) < $tol if LOS_AJ_1! = . 


** 3 Check answers are the same with all the options compared to without

* #3a Overall
hai
msaj, transmatrix(tmat)
rename (P_AJ_*) (pP_AJ_*)
msaj, transmatrix(tmat) ci
rename (P_AJ_*ci) (cP_AJ_*ci)
drop P_AJ_*
msaj,  transmatrix(tmat) se
rename (P_AJ_*se) (sP_AJ_*se)
drop P_AJ_*
msaj, transmatrix(tmat) los
rename (LOS*) (lLOS*)
drop P_AJ_*
msaj, transmatrix(tmat) ci se los
forvalues i=1/6 {
	di "State `i'"
	assert P_AJ_`i' == pP_AJ_`i'
	assert P_AJ_`i'_se == sP_AJ_`i'_se
	assert P_AJ_`i'_lci == cP_AJ_`i'_lci
	assert P_AJ_`i'_uci == cP_AJ_`i'_uci
	assert LOS_AJ_`i' == lLOS_AJ_`i'
}
 


*****************************
*  		EBMT - CR  			*
*****************************

** #2 Common sense checks

* #2a Check probabilities add up to 1
use "./data/ebmt_cr", replace
msaj , cr
egen double testsum=rowtotal(P_AJ_*)
assert sqrt((testsum-1)^2) < $tol if P_AJ_1 != .

* #2b Check all probabilities and CI are within [0,1]
use "./data/ebmt_cr", replace
msaj, cr ci
forvalues i=1/5 {
	di "State `i'"
	assert P_AJ_`i' >=0 & P_AJ_`i' <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_lci >=0 & P_AJ_`i'_lci <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_uci >=0 & P_AJ_`i'_uci <=1 if P_AJ_`i' != .
}

* #2d Check LOS adds up to time
use "./data/ebmt_cr", replace
msaj , cr los
egen double testsum = rowtotal(LOS*)
assert sqrt((testsum-_t)^2) < $tol if LOS_AJ_1 != . 


** 3 Check answers are the same with all the options compared to without

* #3a Overall
use "./data/ebmt_cr", replace
msaj, cr
rename (P_AJ_*) (pP_AJ_*)
msaj, cr ci
rename (P_AJ_*ci) (cP_AJ_*ci)
drop P_AJ_*
msaj,  cr se
rename (P_AJ_*se) (sP_AJ_*se)
drop P_AJ_*
msaj, cr los
rename (LOS*) (lLOS*)
drop P_AJ_*
msaj, cr ci se los
forvalues i=1/5 {
	di "State `i'"
	assert P_AJ_`i' == pP_AJ_`i'
	assert P_AJ_`i'_se == sP_AJ_`i'_se
	assert P_AJ_`i'_lci == cP_AJ_`i'_lci
	assert P_AJ_`i'_uci == cP_AJ_`i'_uci
	assert LOS_AJ_`i' == lLOS_AJ_`i'
}



*****************************
*  			Simbi 			*
*****************************

cap program drop simbi
program define simbi
	use "./data/simbi", replace
	mat tmat = (.,1,2\3,.,4\.,.,.)
end


** #2 Common sense checks

* #2a Check probabilities add up to 1
simbi
msaj , transmatrix(tmat) 
egen double testsum=rowtotal(P_AJ_*)
assert sqrt((testsum-1)^2) < $tol if P_AJ_1 != .

* #2b Check all probabilities and CI are within [0,1]
simbi
msaj, transmatrix(tmat) ci
forvalues i=1/3 {
	di "State `i'"
	assert P_AJ_`i' >=0 & P_AJ_`i' <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_lci >=0 & P_AJ_`i'_lci <=1 if P_AJ_`i' != .
	assert P_AJ_`i'_uci >=0 & P_AJ_`i'_uci <=1 if P_AJ_`i' != .
}

* #2d Check LOS adds up to time
simbi
msaj , transmatrix(tmat) los
egen double testsum = rowtotal(LOS*)
assert sqrt((testsum-_t)^2) < $tol if LOS_AJ_1! = . 


** 3 Check answers are the same with all the options compared to without

* #3a Overall
simbi
msaj, transmatrix(tmat)
rename (P_AJ_*) (pP_AJ_*)
msaj, transmatrix(tmat) ci
rename (P_AJ_*ci) (cP_AJ_*ci)
drop P_AJ_*
msaj,  transmatrix(tmat) se
rename (P_AJ_*se) (sP_AJ_*se)
drop P_AJ_*
msaj, transmatrix(tmat) los
rename (LOS*) (lLOS*)
drop P_AJ_*
msaj, transmatrix(tmat) ci se los
forvalues i=1/3 {
	di "State `i'"
	assert P_AJ_`i' == pP_AJ_`i'
	assert P_AJ_`i'_se == sP_AJ_`i'_se
	assert P_AJ_`i'_lci == cP_AJ_`i'_lci
	assert P_AJ_`i'_uci == cP_AJ_`i'_uci
	assert LOS_AJ_`i' == lLOS_AJ_`i'
}
