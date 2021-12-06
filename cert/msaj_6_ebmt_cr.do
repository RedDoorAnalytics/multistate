* Version 0.1, 12/02/2020


*****************************
*  			Itself	 		*
*****************************

* #1 Check by option 

* #1a LOS option
use "./data/ebmt_cr", replace
qui gen cov = mod(id,10)
msaj , cr by(cov) los
rename (LOS_AJ_*) (l*)
forvalues j=0/9 {
	di "By value `j'"
	cap drop *_AJ_*
	msaj if cov == `j', cr los
	forvalues i=1/5 {
		di "State `i'"
		assert l`i' == LOS_AJ_`i' if LOS_AJ_`i' != .
	}
}

* #1b All options
use "./data/ebmt_cr", replace
qui gen cov = mod(id,10)
msaj , cr by(cov) ltruncated(25) ci se los
rename (P_AJ_*) (p*)
rename (LOS_AJ_*) (l*)
forvalues j=0/9 {
	di "By value `j'"
	cap drop *_AJ_*
	msaj if cov == `j', cr ltruncated(25) ci se los
	forvalues i=1/5 {
		di "State `i'"
		assert p`i' == P_AJ_`i' if P_AJ_`i' != .
		assert p`i'_se == P_AJ_`i'_se if P_AJ_`i'_se != .
		assert p`i'_lci == P_AJ_`i'_lci if P_AJ_`i'_lci != .
		assert p`i'_uci == P_AJ_`i'_uci if P_AJ_`i'_uci != .
		assert l`i' == LOS_AJ_`i' if LOS_AJ_`i' != .
	}
}



*****************************
*  			mstate 			*
*****************************

** #2 Probabilities & CI

* #2a Entry time 0, from state 1
use "./data/ebmt_cr", replace
msaj, cr se
merge m:1 _t using "./data/EBMT_CR_prob_e0.dta"
drop if _t == 0 & _merge == 2													// Not in original data
forvalues i=1/5 {
	di "State `i'"
	assert sqrt((P_AJ_`i' - p1`i')^2) < $tol if _d == 1
	assert sqrt((P_AJ_`i'_se - se1`i')^2) < $tol if _d == 1
}


* #2c Additional entry times
cap program drop ebmt_cr_prob
program define ebmt_cr_prob
	syntax , [ num(integer 0) ltruncated(real 0) ]
			
	* From state 1
	use "./data/ebmt_cr", replace
	msaj, cr ltruncated(`ltruncated') se
	merge m:1 _t using "./data/EBMT_CR_prob_e`num'.dta"
	assert P_AJ_1 == . if _merge == 1															
	drop if _merge == 2 & _t == `ltruncated'													// Not in original data	
	forvalues i=1/5 {
		di "State `i'"
		assert sqrt((P_AJ_`i' - p1`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
		assert sqrt((P_AJ_`i'_se - se1`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
	} 
end

ebmt_cr_prob, num(1) ltruncated(0.03)
ebmt_cr_prob, num(2) ltruncated(0.05)
ebmt_cr_prob, num(3) ltruncated(1)
ebmt_cr_prob, num(4) ltruncated(25)
ebmt_cr_prob, num(5) ltruncated(29.5)
ebmt_cr_prob, num(6) ltruncated(121)
ebmt_cr_prob, num(7) ltruncated(397)
ebmt_cr_prob, num(8) ltruncated(1313)


** #3 LOS

* #3a Entry time 0, from state 1
use "./data/ebmt_cr", replace
msaj, cr los
merge m:1 _t using "./data/EBMT_CR_los_e0.dta", keep(3) nogen
sum _t if _d == 1
forvalues i=1/5 {
	di "State `i'"
	assert sqrt((LOS_AJ_`i' - l1`i')^2) < $tol  if _t <= `r(max)'
}


* #3c Additional entry times
cap program drop ebmt_cr_los
program define ebmt_cr_los
	syntax , [ num(integer 0) ltruncated(real 0) ]
			
	* From state 1
	use "./data/ebmt_cr", replace
	msaj, cr ltruncated(`ltruncated') los
	merge m:1 _t using "./data/EBMT_CR_los_e`num'.dta", keep(3)	
	sum _t if _d == 1
	forvalues i=1/5 {
		di "State `i'"
		assert sqrt((LOS_AJ_`i' - l1`i')^2) < $tol if _t <= `r(max)'
	}
end

ebmt_cr_los, num(1) ltruncated(0.03)
ebmt_cr_los, num(2) ltruncated(0.05)
ebmt_cr_los, num(3) ltruncated(1)
ebmt_cr_los, num(4) ltruncated(25)
ebmt_cr_los, num(5) ltruncated(29.5)
ebmt_cr_los, num(6) ltruncated(121)
