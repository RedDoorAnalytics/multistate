//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/multistate/multistate"
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./stms"
adopath ++ "./msboxes"
adopath ++ "./msaj"
clear all

do ./build/buildmlib.do
mata mata clear

use "./data/vonCube",clear

mat tmat = (.,1,2,3,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
mat list tmat

msset, id(admid) states(status2 status3 status4 status5 status6) times(t12 t13 t14 t25 t26) transmat(tmat)
mat list r(freqmatrix)

// msboxes, id(admid) transmat(tmat) xvalues(0.25 0.75 0.13 0.37 0.63 0.87) yvalues(0.8 0.8 0.2 0.2 0.2 0.2) ///
// 	statenames("Admission" "HAI" "Discharge" "Death" "Discharge" "Death")

stset _stop, failure(_status) enter(_start)



****************************
* Choose tranisition models*
****************************

forvalues i = 1/5 {
	display "Transition `i'"
	
	* Exponential
	streg if _trans==`i', dist(exponential)
	estimates store m`i'_exp
	
	* Weibull
	streg if _trans==`i', dist(weibull)
	estimates store m`i'_weib
	
	* Gompertz
	streg if _trans==`i', dist(gompertz)
	estimates store m`i'_gom
	
	* Log logistic
	streg if _trans==`i', dist(loglogistic)
	estimates store m`i'_logl
	
	* Log normal
	streg if _trans==`i', dist(lognormal)
	estimates store m`i'_logn
	
	* Generalised gamma
	streg if _trans==`i', dist(ggamma)
	estimates store m`i'_ggam
	
	* Royston Parmar models
	forvalues j=2/5 {
	
		if `i'==2 & `j'==2 {
			stpm2 if _trans==`i', scale(hazard) df(`j') lininit
			estimates store m`i'_rp`j'
		}
	
		else {
			stpm2 if _trans==`i', scale(hazard) df(`j')
			estimates store m`i'_rp`j'
		}
	}
	
}

/*
* Try extra spline positions

stpm2 if _trans==2, scale(hazard) df(6)
estimates store m2_rp6

forvalues i = 3/5 {
	display "Transition `i'"
	forvalues j=6/10 {
		stpm2 if _trans==`i', scale(hazard) df(`j')
		estimates store m`i'_rp`j'
	}
}


* Try using strcs - RCS on the log hazard scale, needs Stata v15.1
* strcs if _trans==1, df(2)
*/
	
* Display the estimates
forvalues i = 1/5 {
	di "Transition `i'"
	estimates stat m`i'*
}

estimates clear



****************************
*  		Hazard rates	   *
****************************

** Graph of best fitting flexible transitions and CIS

* Create the time variable
range timevar80 0 80 50


* Create estimates from best fitting flexible transitions

stpm2 if _trans==1, scale(hazard) df(4)
estimates store m_f1
predict fhaz1, ci hazard timevar(timevar80)

streg if _trans==2, dist(ggamma)
estimates store m_f2
predict fhaz2,  hazard

streg if _trans==3, dist(ggamma)
estimates store m_f3
predict fhaz3, hazard

streg if _trans==4, dist(lnormal)
estimates store m_f4
predict fhaz4, hazard

streg if _trans==5, dist(ggamma)
estimates store m_f5
predict fhaz5, hazard


** Graph comparing hazards

* Create estimates for the exponential model
forvalues i = 1/5 {
	streg if _trans==`i', dist(e)
	estimates store m_e`i'
}

* Create estimates from RP model with 3 df for all
forvalues i = 1/5 {
	stpm2 if _trans==`i', scale(hazard) df(3)
	estimates store m_rp`i'
}

// * Produce the graph using the epan2 smoother 
// estimates restore m_e1
// sts graph if _trans==1, h kernel(epan2) ///
// 	addplot(line fhaz1 timevar80, ///
// 		legend(order(1 "Empirical" 2 "Flexible") )) ///
// 	yline(`=exp(_b[_cons])', lpattern(dash) lcolor(black)) ///	
// 	xlab(0(20)80) ///
//     ylab( , angle(h) format(%4.2f)) ///
// 	xtitle("Time since hospital admission") ///
//     ytitle("Hazard rate") ///
// 	title("Admission to HAI") ///
// 	name(haz1, replace)
//
// local title2 "Admission to Discharge"
// local title3 "Admission to Death"
// local title4 "HAI to Discharge"
// local title5 "HAI to Death"
// forvalues i=2/5 {
// 	estimates restore m_e`i'
// 	sts graph if _trans==`i' , h kernel(epan2) ///
// 		addplot(line fhaz`i' _t, sort ///
// 			legend(order(1 "Empirical" 2 "Flexible") )) ///
// 		yline(`=exp(_b[_cons])', lpattern(dash) lcolor(black)) ///	
// 		xlab(0(20)80) ///
// 		ylab( , angle(h) format(%4.2f)) ///
// 		xtitle("Time since hospital admission") ///
// 		ytitle("Hazard rate") ///
// 		title(`title`i'') ///
// 		name(haz`i', replace)
// }
//
// // grc1leg haz1 haz2 haz3 haz4 haz5 , legendfrom(haz1) xcommon ycommon rows(2)
// graph export "../Figures/Hazards.png", replace



****************************
* Transition Probabilities *
****************************

* Exponential model
predictms, transm(tmat) timevar(timevar80) models(m_e1 m_e2 m_e3 m_e4 m_e5)
rename (_prob_at1_1_*) (eprob*)

* Royston-Parmar model
predictms, transm(tmat) timevar(timevar80) models(m_rp1 m_rp2 m_rp3 m_rp4 m_rp5) //ci m(100)
rename (_prob_at1_1_*) (rpprob*)

* Best fitting flexible
predictms, transm(tmat) timevar(timevar80) models(m_f1 m_f2 m_f3 m_f4 m_f5) //ci m(100)
rename (_prob_at1_1_*) (fprob*)


* Aalen-Johansen estimates (non-para)
msaj, transm(tmat) id(admid)
rename (P_AJ_*) (ajprob*)

local title1 "Admission"
local title2 "HAI"
local title3 "Discharge no HAI"
local title4 "Death no HAI"
local title5 "Discharge after HAI"
local title6 "Death after HAI"

forvalues i=1/6 {
	twoway (line  eprob`i' timevar80, sort ) ///
	(line  rpprob`i' timevar80, sort ) ///
	(line  fprob`i' timevar80, sort ) ///
	(line  ajprob`i' _t, sort connect(stepstair) ) , ///
	name(prob`i', replace) ///
	ylab(0(0.2)1, angle(h) format(%3.1f)) ///
	xtitle("Time since hospital admission") ///
	ytitle("Probability") ///
	title(`title`i'') ///
	legend(order(1 "Exponential" 2 "RP" 3 "Best" 4 "AJ") rows(1))
}

grc1leg prob1 prob3 prob4 prob2 prob5 prob6, xcommon  ycommon rows(2) legendfrom(prob1)
// graph export "./Figures/TransProbs.png", replace


