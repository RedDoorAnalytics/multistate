//local drive Z:\
local drive /Users/Michael/Documents
//local drive c:
cd "`drive'/multistate/multistate"
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./stms"
adopath ++ "./msboxes"
adopath ++ "./msaj"
adopath ++ "./graphms"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 245987

/* q9 Probability of death in a competing risks framework (cause-specific survival) */
// Part B using predictms

use /Users/Michael/Documents/multistate/teaching/wengen/data/colon, clear
drop if stage ==0
gen female = sex==2

// (a) use msset
gen d_cancer = status==1
gen d_other = status==2

matrix transmat = .,1,2 \ .,.,. \.,.,.
matrix list transmat
global statenames alive cancer other
gen surv = surv_mm/12

msset, id(id) states(d_cancer d_other) times(surv surv) ///
  transmat(transmat)

// (b) Fit separate models for the effect of sex
// These are similar to fitting the models simultaneously in Part A
stset _stop, enter(_start) failure(_status=1) exit(time 10)

range temptime 0 10 51
tab agegrp, gen(agegrp)

stmerlin female agegrp2-agegrp4 if _trans==1, dist(rp) df(3) ///
  tvc(agegrp2-agegrp4) dftvc(1)
estimates store cancer_tvc


// (f) Optional extension
// Model other causes using attained age
gen _start_age = age + _start
gen _stop_age = age + _stop

bysort id: gen first = _n==1
rcsgen age, gen(agercs) df(3) if2(first)
global ageknots `r(knots)'

// include splines for age and make time-dependent
stmerlin female agercs* if _trans==1, dist(rp) df(3) ///
  tvc(agercs*) dftvc(1)
estimate store cancer_agercs

stmerlin female agercs* if _trans==2, dist(rp) df(1) ///
  //tvc(agercs*) dftvc(1) 
estimate store other_agercs

// Model with attained age as the time scale for other causes
// Need to re-stset, but be careful if fitting more models!
stset _stop_age, enter(_start_age) failure(_status=1) exit(time _start_age + 10)
stmerlin female if _trans==2, dist(rp) df(3)
estimates store other_agetime

// compare a 40 year old with a 70 year old
rcsgen , scalar(40) gen(a40_) knots(${ageknots})
rcsgen , scalar(70) gen(a70_) knots(${ageknots})

// Note M is made small so runs quickly for teaching purposes
predictms, probability models(cancer_agercs other_agercs) timevar(temptime) transmatrix(transmat) ///
  at1(female 0 agercs1 `=a40_1' agercs2 `=a40_2' agercs3 `=a40_3') ///
  at2(female 0 agercs1 `=a70_1' agercs2 `=a70_2' agercs3 `=a70_3') difference latent
// Change from ugly to nicer names
forvalues i = 1/3 {
  local name = word("${statenames}",`i')
  rename _prob_at1_1_`i'* P_`name'_F_a40a*
  rename _prob_at2_1_`i'* P_`name'_F_a70a*
  rename _diff_prob_at2_1_`i'* Diffa_`name'*
}

predictms, probability models(cancer_agercs other_agetime) timevar(temptime) transmatrix(transmat) ///
  at1(female 0 agercs1 `=a40_1' agercs2 `=a40_2' agercs3 `=a40_3') ///
  at2(female 0 agercs1 `=a70_1' agercs2 `=a70_2' agercs3 `=a70_3') difference tscale2(2) time2(40 70) latent
// Change from ugly to nicer names
forvalues i = 1/3 {
  local name = word("${statenames}",`i')
  rename _prob_at1_1_`i'* P_`name'_F_a40b*
  rename _prob_at2_1_`i'* P_`name'_F_a70b*
  rename _diff_prob_at2_1_`i'* Diffb_`name'*
}


foreach name in $statenames {
  twoway (line P_`name'_F_a40?  temptime, color(red)) ///
         , name(P_`name'_40, replace)
  twoway (line P_`name'_F_a70?  temptime, color(red)) ///
         , name(P_`name'_70, replace)
}	  
