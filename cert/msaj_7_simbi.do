* Version 0.1, 12/02/2020

cap program drop simbi
program define simbi
	use "./data/simbi", replace
	mat tmat = (.,1,2\3,.,4\.,.,.)
end

cap program drop simbi_cov
program define simbi_cov
	use "./data/simbi", replace
	mat tmat = (.,1,2\3,.,4\.,.,.)
	qui gen cov = mod(id,6)*2 - 3
end


*****************************
*  			Itself	 		*
*****************************

* #1 Check by option

* #1a LOS option
simbi_cov
msaj , transmatrix(tmat) by(cov) los
rename (LOS_AJ_*) (l*)
foreach j in -3 -1 1 3 5 7 {
	di "By value `j'"
	cap drop *_AJ_*
	msaj if cov == `j', transmatrix(tmat) los
	forvalues i=1/3 {
		di "State `i'"
		assert l`i' == LOS_AJ_`i' if LOS_AJ_`i' != .
	}
}

* #1b All options
simbi_cov
msaj , transmatrix(tmat) by(cov) from(2) ltruncated(2) ci se los
rename (P_AJ_*) (p*)
rename (LOS_AJ_*) (l*)
foreach j in -3 -1 1 3 5 7 {
	di "By value `j'"
	cap drop *_AJ_*
	msaj if cov == `j', transmatrix(tmat) from(2) ltruncated(2) ci se los
	forvalues i=1/3 {
		di "State `i'"
		assert p`i' == P_AJ_`i' if P_AJ_`i' != .
		assert p`i'_se == P_AJ_`i'_se if P_AJ_`i'_se != .
		assert p`i'_lci == P_AJ_`i'_lci if P_AJ_`i'_lci != .
		assert p`i'_uci == P_AJ_`i'_uci if P_AJ_`i'_uci != .
		assert l`i' == LOS_AJ_`i' if LOS_AJ_`i' != .
	}
}



*****************************
*  			mstate	 		*
*****************************

** #2 Probabilities & CI

* #2a Entry time 0, from state 1
simbi
msaj, transmatrix(tmat) se
merge m:1 _t using "./data/SIMBI_prob_e0.dta"
drop if _t == 0 & _merge == 2													// Not in original data
forvalues i=1/3 {
	di "State `i'"
	assert sqrt((P_AJ_`i' - p1`i')^2) < $tol if _d == 1
	assert sqrt((P_AJ_`i'_se - se1`i')^2) < $tol if _d == 1
}

* #2b Entry time 0, from state 2
simbi
msaj, transmatrix(tmat) from(2) se
merge m:1 _t using "./data/SIMBI_prob_e0.dta"
drop if _t == 0 & _merge == 2													// Not in original data	
forvalues i=1/3 {
	di "State `i'"
	assert sqrt((P_AJ_`i' - p2`i')^2) < $tol if _d == 1
	assert sqrt((P_AJ_`i'_se - se2`i')^2) < $tol if _d == 1
}


* #2c Additional entry times
cap program drop simbi_prob
program define simbi_prob
	syntax , [ num(integer 0) ltruncated(real 0) ]
			
	* From state 1
	simbi
	msaj, transmatrix(tmat) ltruncated(`ltruncated') se
	merge m:1 _t using "./data/SIMBI_prob_e`num'.dta"
	assert P_AJ_1 == . if _merge == 1															
	drop if _merge == 2 & _t == `ltruncated'													// Not in original data	
	forvalues i=1/3 {
		di "State `i' - From 1"
		assert sqrt((P_AJ_`i' - p1`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
		assert sqrt((P_AJ_`i'_se - se1`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
	} 

	* From state 2			
	simbi
	msaj, transmatrix(tmat) ltruncated(`ltruncated') from(2) se
	merge m:1 _t using "./data/SIMBI_prob_e`num'.dta"
	assert P_AJ_2 == . if _merge == 1															
	drop if _merge == 2 & _t == `ltruncated'													// Not in original data	
	forvalues i =1/3 {
		di "State `i' - From 2"
		assert sqrt((P_AJ_`i' - p2`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
		assert sqrt((P_AJ_`i'_se - se2`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
	}
end

simbi_prob, num(1) ltruncated(500)
simbi_prob, num(2) ltruncated(1000)
simbi_prob, num(3) ltruncated(2000)
simbi_prob, num(4) ltruncated(3000)
simbi_prob, num(5) ltruncated(4000)
simbi_prob, num(6) ltruncated(5000)



** #3 LOS

* #3a Entry time 0, from state 1
simbi
msaj, transmatrix(tmat) los
merge m:1 _t using "./data/SIMBI_los_e0.dta", keep(3) nogen
sum _t if _d == 1
forvalues i=1/3 {
	di "State `i'"
	assert sqrt((LOS_AJ_`i' - l1`i')^2) < $tol  if _t <= `r(max)'
}

* #3b Entry time 0, from state 2 
simbi
msaj, transmatrix(tmat) from(2) los
merge m:1 _t using  "./data/SIMBI_los_e0.dta", keep(3) nogen
sum _t if _d == 1
forvalues i=1/3 {
	di "State `i'"
	assert sqrt((LOS_AJ_`i' - l2`i')^2) < $tol if _t <= `r(max)'
}


* #3c Additional entry times
cap program drop simbi_los
program define simbi_los
	syntax , [ num(integer 0) ltruncated(real 0) ]
			
	* From state 1
	simbi
	msaj, transmatrix(tmat) ltruncated(`ltruncated') los
	merge m:1 _t using "./data/SIMBI_los_e`num'.dta", keep(3)	
	sum _t if _d == 1
	forvalues i=1/3 {
		di "State `i' - From 1"
		assert sqrt((LOS_AJ_`i' - l1`i')^2) < $tol if _t <= `r(max)'
	}

	* From state 4			
	simbi
	msaj, transmatrix(tmat) ltruncated(`ltruncated') from(2) los
	merge m:1 _t using "./data/SIMBI_los_e`num'.dta", keep(3)	
	sum _t if _d == 1
	forvalues i =1/3  {
		di "State `i' - From 2"
		assert sqrt((LOS_AJ_`i' - l2`i')^2) < $tol if _t <= `r(max)'
	}
end

simbi_los, num(1) ltruncated(500)
simbi_los, num(2) ltruncated(1000)
simbi_los, num(3) ltruncated(2000)
simbi_los, num(4) ltruncated(3000)
simbi_los, num(5) ltruncated(4000)
simbi_los, num(6) ltruncated(5000)



*****************************
*  			etm		 		*
*****************************

** #2 Probabilities & CI

* #2a Entry time 0, from state 1
simbi
msaj, transmatrix(tmat) se
merge m:1 _t using "./data/SIMBI_etm_prob_e0.dta"
drop if _d == 0 & _merge == 1													
forvalues i=1/3 {
	di "State `i'"
	assert sqrt((P_AJ_`i' - p1`i')^2) < $tol if _d == 1
	assert sqrt((P_AJ_`i'_se - se1`i')^2) < $tol if _d == 1
}

* #2b Entry time 0, from state 2
simbi
msaj, transmatrix(tmat) from(2) se
merge m:1 _t using "./data/SIMBI_etm_prob_e0.dta"
drop if _d == 0 & _merge == 1													
forvalues i=1/3 {
	di "State `i'"
	assert sqrt((P_AJ_`i' - p2`i')^2) < $tol if _d == 1
	assert sqrt((P_AJ_`i'_se - se2`i')^2) < $tol if _d == 1
}


* #2c Additional entry times
cap program drop simbi_etm_prob
program define simbi_etm_prob
	syntax , [ num(integer 0) ltruncated(real 0) ]
			
	* From state 1
	simbi
	msaj, transmatrix(tmat) ltruncated(`ltruncated') se
	merge m:1 _t using "./data/SIMBI_etm_prob_e`num'.dta"
	assert P_AJ_1 == . if _merge == 1 & _t > `ltruncated'															
	drop if  (_t <= `ltruncated') | (_d == 0 & _merge == 1)	
	assert _merge == 3
	forvalues i=1/3 {
		di "State `i' - From 1"
		assert sqrt((P_AJ_`i' - p1`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
		assert sqrt((P_AJ_`i'_se - se1`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
	} 

	* From state 2			
	simbi
	msaj, transmatrix(tmat) ltruncated(`ltruncated') from(2) se
	merge m:1 _t using "./data/SIMBI_etm_prob_e`num'.dta"
	assert P_AJ_2 == . if _merge == 1 & _t > `ltruncated'															
	drop if  (_t <= `ltruncated') | (_d == 0 & _merge == 1)	
	assert _merge == 3		 	
	forvalues i =1/3 {
		di "State `i' - From 2"
		assert sqrt((P_AJ_`i' - p2`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
		assert sqrt((P_AJ_`i'_se - se2`i')^2) < $tol if _d == 1 & _t >=`ltruncated'
	}
end

simbi_etm_prob, num(1) ltruncated(500)
simbi_etm_prob, num(2) ltruncated(1000)
simbi_etm_prob, num(3) ltruncated(2000)
simbi_etm_prob, num(4) ltruncated(3000)
simbi_etm_prob, num(5) ltruncated(4000)
simbi_etm_prob, num(6) ltruncated(5000)
