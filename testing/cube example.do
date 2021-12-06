// Load data

clear
adopath ++ "Z:\stdev\multistate"
cd "Z:\stdev\multistate"


import delimited using  "los.data.csv"
drop v1
foreach var in j01 j02 j03 j12 j13  {
	replace `var' = "" if `var' == "Inf"
}
destring j*, replace
rename j12 j14
rename j13 j15

gen max0 = max(j01, j02, j03)
forvalues i = 1/3 {
	gen status`i' = j0`i' == max0
	replace j0`i' = max0
}

gen max1 = max(j14, j15)
forvalues i = 4/5 {
	gen status`i' = j1`i' == max1
	replace j1`i' = max1
}

// transition matrix for extended illness death model
matrix transmat = (	., 1, 2, 3, ., . \ ///
					., ., ., ., 4, 5 \ ///)
					., ., ., ., ., . \ ///)
					., ., ., ., ., . \ ///)
					., ., ., ., ., . \ ///)
					., ., ., ., ., .)

// msset the data
msset, id(admid) states(status1 status2 status3 status4 status5) times(j01 j02 j03 j14 j15) transmat(transmat) 					
			
stset _stop, enter(_start) failure(_status = 1) exit(time 40)

// timevar for predictions
range timevar 0 40 100

// Fit exponential model and stpm2 model with 3 df for each transition
forvalues i = 1/5 {
	streg if _trans == `i', dist(exp) nolog noshow
	estimates store m_exp`i'
	stpm2 if _trans==`i', scale(hazard) df(3)
	estimates store m_rp`i'
}

// compare hazard functions for each transition
forvalues i = 1/5 {
	estimates restore m_rp`i'
	capture drop tmph
	predict tmph, hazard timevar(timevar) 
	//estimate restore m`i'
	twoway (line tmph timevar),  yline(`=exp(_b[_cons])') name(haz`i', replace) ///
		 ylab(0(0.025)0.20, angle(h) format(%4.3f)) ///
		 xtitle("Time since hospital admission") ///
		 ytitle("hazard rate")
}
graph combine haz1 haz2 haz3 haz4 haz5, nocopies

// predictions for exponential models
predictms, transmatrix(transmat) models(m_exp1 m_exp2 m_exp3 m_exp4 m_exp5) n(100000) timevar(timevar) ///
		gen(p_exp)

// Aalen-Johansson estimates of transition probabilities		
msAJ, transmat(transmat) 

// Compare AJ and exponential transition probabilities
forvalues i = 1/5 {
	twoway (line p_exp_1_`i' timevar) ///
			(line P`i' t_P, connect(stairstep)) ///
			,  name(transprob_exp`i', replace) ///
		 ylab(0(0.2)1, angle(h) format(%3.1f)) ///
		 xtitle("Time since hospital admission") ///
		 ytitle("Probability") ///
		 legend(order(1 "Exponential" 2 "AJ") ring(0) pos(1) cols(1))
}
graph combine transprob_exp1 transprob_exp2 transprob_exp3 transprob_exp4 transprob_exp5, nocopies

// predictions for stpm2 models
predictms, transmatrix(transmat) models(m_rp1 m_rp2 m_rp3 m_rp4 m_rp5) n(100000) timevar(timevar) ///
		gen(p_rp)

// Compare AJ and both parametric transition probabilities
forvalues i = 1/5 {
	twoway (line p_exp_1_`i' timevar) ///
			(line p_rp_1_`i' timevar) ///
			(line P`i' t_P, connect(stairstep)) ///
			,  name(transprob_rp`i', replace) ///
		 ylab(0(0.2)1, angle(h) format(%3.1f)) ///
		 xtitle("Time since hospital admission") ///
		 ytitle("Probability") ///
		 legend(order(1 "Exponential" 2 "RP" 3 "AJ") ring(0) pos(1) cols(1))
}		
graph combine transprob_rp1 transprob_rp2 transprob_rp3 transprob_rp4 transprob_rp5, nocopies


// *******************************************************
// NEED to pass functions to predictms to do this with CIS
// *******************************************************

// Mortality Risk
// P(D(t) = 1|E(t)=0))
gen mortrisk0_exp = p_exp_1_4 / (p_exp_1_1 + p_exp_1_3 + p_exp_1_4)
gen mortrisk0_np = P4 / (P1 + P3 + P4)

twoway (line mortrisk0_exp timevar) ///
		(line mortrisk0_np t_P, connect(stairstep)) ///
		,

gen mortrisk0_rp = p_rp_1_4 / (p_rp_1_1 + p_rp_1_3 + p_rp_1_4)
twoway (line mortrisk0_exp timevar) ///
		(line mortrisk0_rp timevar) ///
		(line mortrisk0_np t_P, connect(stairstep)) ///
		,
// P(D(t) = 1|E(t)=1))
gen mortrisk1_exp = p_exp_1_6 / (p_exp_1_2 + p_exp_1_5 + p_exp_1_6)
gen mortrisk1_np = P6 / (P2 + P5 + P6)
gen mortrisk1_rp = p_rp_1_6 / (p_rp_1_2 + p_rp_1_5 + p_rp_1_6)
twoway (line mortrisk1_exp timevar) ///
		(line mortrisk1_rp timevar) ///
		(line mortrisk1_np t_P, connect(stairstep)) ///
		,
// Overall Mortality
gen mortall_exp = p_exp_1_4 + p_exp_1_6
gen mortall_rp = p_rp_1_4 + p_rp_1_6
gen mortall_np = P4 + P6
twoway (line mortall_exp timevar) ///
		(line mortall_rp timevar) ///
		(line mortall_np t_P, connect(stairstep)) ///
		,
		
// attributable mortality
gen AM_exp = mortrisk1_exp - mortrisk0_exp		
gen AM_rp = mortrisk1_rp - mortrisk0_rp
gen AM_np = mortrisk1_np - mortrisk0_np

twoway	(line AM_exp AM_rp timevar) ///
		(line AM_np t_P, connect(stairstep)) ///
		,
		
// PAF	
gen PAF_exp = (mortall_exp - mortrisk0_exp) / mortall_exp
gen PAF_rp = (mortall_rp - mortrisk0_rp) / mortall_rp
gen PAF_np = (mortall_np - mortrisk0_np) / mortall_np

twoway	(line PAF_exp PAF_rp timevar) ///
		(line PAF_np t_P, connect(stairstep)) ///
		,



		