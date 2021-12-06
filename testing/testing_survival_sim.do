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

pr define sim1,rclass

	clear
	set obs 1000
	survsim stime died, dist(weibull) lambda(0.1) gamma(1.2) maxt(10)
	
	stset stime, failure(died)

	cap range temptime 0 5 100

	stmerlin , dist(weib)
	est store m1
	
	predictms, singleevent probability models(m1) timevar(temptime) aj ci 
	return scalar b = _prob_at1_1_1[100]
	return scalar b_lci = _prob_at1_1_1_lci[100]
	return scalar b_uci = _prob_at1_1_1_uci[100]
	
	predictms, singleevent probability models(m1) timevar(temptime) ci n(1000) m(200)
	return scalar bsim = _prob_at1_1_1[100]
	return scalar bsim_lci = _prob_at1_1_1_lci[100]
	return scalar bsim_uci = _prob_at1_1_1_uci[100]
	
	predictms, singleevent probability models(m1) timevar(temptime) ci n(1000) m(200) percentile
	return scalar bsim2 = _prob_at1_1_1[100]
	return scalar bsim2_lci = _prob_at1_1_1_lci[100]
	return scalar bsim2_uci = _prob_at1_1_1_uci[100]
	
	mat tmat = (.,1\.,.)
	gen _trans=1
	gen _from = 1
	gen _to = 2
	msaj, transmat(tmat) ci
	
	keep if _d==1 & _t<5
	sort _t
	return scalar baj = P_AJ_1[_N]
	return scalar baj_lci = P_AJ_1_lci[_N]
	return scalar baj_uci = P_AJ_1_uci[_N]
end

set seed 245987
simulate b=r(b) b_lci=r(b_lci) b_uci=r(b_uci) ///
bsim=r(bsim) bsim_lci=r(bsim_lci) bsim_uci=r(bsim_uci) ///
bsim2=r(bsim) bsim2_lci=r(bsim2_lci) bsim2_uci=r(bsim2_uci) ///
 baj=r(baj) baj_lci=r(baj_lci) baj_uci=r(baj_uci) ///
 , reps(200) : sim1

local truth = exp(-0.1*5^1.2)
di `truth'
gen cov = (`truth'>b_lci & `truth'<b_uci)
tab cov
gen covsim = (`truth'>bsim_lci & `truth'<bsim_uci)
tab covsim
gen covsim2 = (`truth'>bsim2_lci & `truth'<bsim2_uci)
tab covsim2
gen covaj = (`truth'>baj_lci & `truth'<baj_uci)
tab covaj
