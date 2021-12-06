* Version 0.1, 12/02/2020

cap program drop hai
program define hai
	use "./data/hai", replace
	mat tmat = (.,1,2,3,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
end

cap program drop hai_cov
program define hai_cov
	use "./data/hai", replace
	mat tmat = (.,1,2,3,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
	qui gen cov = 0
	qui replace cov = 1 if id > 200
	qui replace cov = -4 if id > 400
end



*****************************
*  			Itself	 		*
*****************************

* #1 Check by option 

* #1a LOS option
hai_cov
msaj , transmatrix(tmat) by(cov) los
rename (LOS_AJ_*) (l*)
foreach j in -4 0 1 {
	di "By value `j'"
	cap drop *_AJ_*
	msaj if cov == `j', transmatrix(tmat) los
	forvalues i=1/6 {
		di "State `i'"
		assert l`i' == LOS_AJ_`i' if LOS_AJ_`i' != .
	}
}

* #1b All options
hai_cov
msaj , transmatrix(tmat) by(cov) from(2) ltruncated(4) ci se los
rename (P_AJ_*) (p*)
rename (LOS_AJ_*) (l*)
foreach j in -4 0 1 {
	di "By value `j'"
	cap drop *_AJ_*
	msaj if cov == `j', transmatrix(tmat) from(2) ltruncated(4)  ci se los
	forvalues i=1/6 {
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
hai
msaj, transmatrix(tmat) se
merge m:1 _t using "./data/HAI_prob_e0.dta"
drop if _t == 0 & _merge == 2													// Not in original data
forvalues i=1/6 {
	di "State `i'"
	assert sqrt((P_AJ_`i' - p1`i')^2) < $tol if _d == 1
	assert sqrt((P_AJ_`i'_se - se1`i')^2) < $tol if _d == 1
}

* #2b Entry time 0, from state 2
hai
msaj, transmatrix(tmat) from(2) se
merge m:1 _t using "./data/HAI_prob_e0.dta"
drop if _t == 0 & _merge == 2													// Not in original data	
foreach i in 2 5 6 {
	di "State `i'"
	assert sqrt((P_AJ_`i' - p2`i')^2) < $tol if _d == 1
	assert sqrt((P_AJ_`i'_se - se2`i')^2) < $tol if _d == 1
}


* #2c Additional entry times
cap program drop hai_prob
program define hai_prob
	syntax , [ num(integer 0) ltruncated(real 0) ]
			
	* From state 1
	hai
	msaj, transmatrix(tmat) ltruncated(`ltruncated') se
	merge m:1 _t using "./data/HAI_prob_e`num'.dta"
	assert P_AJ_1 == . if _merge == 1															
	drop if _merge == 2 & _t == `ltruncated'													// Not in original data	
	forvalues i=1/6 {
		di "State `i' - From 1"
		assert sqrt((P_AJ_`i' - p1`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
		assert sqrt((P_AJ_`i'_se - se1`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
	} 

	* From state 2			
	hai
	msaj, transmatrix(tmat) ltruncated(`ltruncated') from(2) se
	merge m:1 _t using "./data/HAI_prob_e`num'.dta"
	assert P_AJ_4 == . if _merge == 1															
	drop if _merge == 2 & _t == `ltruncated'													// Not in original data	
	foreach i in 2 5 6 {
		di "State `i' - From 2"
		assert sqrt((P_AJ_`i' - p2`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
		assert sqrt((P_AJ_`i'_se - se2`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
	}
end

hai_prob, num(1) ltruncated(3)
hai_prob, num(2) ltruncated(4)
hai_prob, num(3) ltruncated(9)
hai_prob, num(4) ltruncated(42)
hai_prob, num(5) ltruncated(70)



** #3 LOS

* #3a Entry time 0, from state 1
hai
msaj, transmatrix(tmat) los 
merge m:1 _t using "./data/HAI_los_e0.dta", keep(3) nogen
forvalues i=1/6 {
	di "State `i'"
	assert sqrt((LOS_AJ_`i' - l1`i')^2) < $tol 
}

* #3b Entry time 0, from state 2 
hai
msaj, transmatrix(tmat) from(2) los
merge m:1 _t using "./data/HAI_los_e0.dta", keep(3) nogen
foreach i in 2 5 6 {
	di "State `i'"
	assert sqrt((LOS_AJ_`i' - l2`i')^2) < $tol 
}


* #3c Additional entry times
cap program drop hai_los
program define hai_los
	syntax , [ num(integer 0) ltruncated(real 0) ]
			
	* From state 1
	hai
	msaj, transmatrix(tmat) ltruncated(`ltruncated') los
	merge m:1 _t using "./data/HAI_los_e`num'.dta", keep(3)	
	forvalues i=1/6 {
		di "State `i' - From 1"
		assert sqrt((LOS_AJ_`i' - l1`i')^2) < $tol
	}

	* From state 4			
	hai
	msaj, transmatrix(tmat) ltruncated(`ltruncated') from(2) los
	merge m:1 _t using "./data/HAI_los_e`num'.dta", keep(3)	
	forvalues i=1/6 {
		di "State `i' - From 2"
		assert sqrt((LOS_AJ_`i' - l2`i')^2) < $tol
	}
end

hai_los, num(1) ltruncated(3)
hai_los, num(2) ltruncated(4)
hai_los, num(3) ltruncated(9)
hai_los, num(4) ltruncated(42)
hai_los, num(5) ltruncated(70)
