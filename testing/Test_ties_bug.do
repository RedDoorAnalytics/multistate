clear all 
adopath ++ "C:\multistate\multistate\msAJ"
cd "R:\CanStat\si144\Presentations\Stata UK Meeting 2018\stata_do"
use "http://www.stata-journal.com/software/sj4-2/st0059/prostatecancer", clear
drop touse

gen agegrp = 0 if age < 75
replace agegrp = 1 if age >= 75 	

//gen noise = runiform()/100
//replace time = time + noise
set seed 1323
keep if agegrp==1 & treatment==1
**keep if runiform()<0.2
replace time = 12 if _n==2
stset time, f(status==1) id(id) exit(time 60) scale(12)
capture drop CIF
stcompet CIF = ci if agegrp == 1 & treatment == 1, compet1(2) compet2(3) 
/* 
tw (line CIF _t if status==1, sort connect(stepstair)), plotregion(margin(zero)) ytitle("Probability of death") ///
xtitle("Years since diagnosis") ylab(0(0.2)1, format(%9.2f)) title("Aalen-Johansen Estimator") ///
subtitle("Patients aged over 75 years old and on treatment", size(small)) name(g1, replace) legend(on order(1 "stcompet") cols(1) ring(0) pos(11) size(small)) 
 */
//gen oldt = _t
gen CIFflag = status==1 & agegrp==1 & treatment==1

gen cancer = status == 1
gen cvd = status == 2
gen other = status == 3
msset, id(id) states(cancer cvd other) times(time time time) cr
stset _stop, enter(_start) f(_status == 1) exit(time 60) scale(12)

**set seed 781236
capture drop P_AJ*
msaj if agegrp == 1 & treatment == 1, cr id(id)
  
sort P_AJ_2 _t
list P_AJ_2 _t if P_AJ_2 != .

**set seed 319387
capture drop P_AJ*
msaj if agegrp == 1 & treatment == 1, cr  id(id) 
sort P_AJ_2 _t
list P_AJ_2 _t if P_AJ_2 != .


sort P_AJ_2 _t
list P_AJ_2 _t

capture drop first

bysort P_AJ_2 (_t): gen first = _n==1
 
	tw (line CIF _t if CIFflag, sort connect(stairstep)) ///
	(line P_AJ_2 P_AJ_3 P_AJ_4 _t if agegrp == 1 & treatment == 1 /*& first*/, sort connect(stairstep..) lpattern(dash..)) ///
	, plotregion(margin(zero)) ytitle("Probability of death") ///
	xtitle("Years since diagnosis") ylab(0(0.2)1, format(%9.2f)) title("Aalen-Johansen Estimator") ///
	subtitle("Patients aged over 75 years old and on treatment", size(small)) name(g2, replace) legend(order(1 "stcompet" 2 "msaj") cols(1) ring(0) pos(11) size(small)) 

sort _d _t id
list id _trans _t _d P_AJ*
