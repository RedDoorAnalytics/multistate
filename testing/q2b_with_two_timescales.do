/* q2 Probability of death in a competing risks framework (cause-specific survival) */
// Part B using predictms

use colon, clear
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
stpm2 female if _trans==1, scale(hazard) df(5)
estimates store cancer
stpm2 female if _trans==2, scale(hazard) df(4)
estimates store other

// (c) Obtain CIFs by sex
// equivalent to using stcompet
msaj , transmat(transmat) by(female)
forvalues i = 1/3 {
  local name = word("${statenames}",`i')
  rename P_AJ_`i' P_`name'_AJ
}

// (d) Predictions using predictms
range temptime 0 10 51
predictms, models(cancer other) timevar(temptime) transmatrix(transmat) ///
  at1(female 0) at2(female 1) 
// give nicer names			
forvalues i = 1/3 {
  local name = word("${statenames}",`i')
  rename _prob_at1_1_`i' P_`name'_model_male
  rename _prob_at2_1_`i' P_`name'_model_female
}

// Compare model based with AJ estimates
foreach name in $statenames {
  twoway  (line P_`name'_AJ _t if female==0, connect(stepstair) sort) ///
    (line P_`name'_model_male temptime) ///
    , name(`name', replace)
}


// (e) Now extend the models to more covariates and time-dependent effects
// more covariates and time-dependent effects 
tab agegrp, gen(agegrp)

stpm2 female agegrp2-agegrp4 if _trans==1, scale(hazard) df(5) ///
  tvc(agegrp2-agegrp4) dftvc(3)
estimates store cancer_tvc

stpm2 female agegrp2-agegrp4 if _trans==2, scale(hazard) df(4)
estimates store other_tvc

// Note M is made small so runs quickly for teaching purposes
predictms, models(cancer_tvc other_tvc) timevar(temptime) transmatrix(transmat) ///
  at1(female 0) at2(female 0 agegrp4 1) difference ci m(20)
// Change from ugly to nicer names
forvalues i = 1/3 {
  local name = word("${statenames}",`i')
  rename _prob_at1_1_`i'* P_`name'_F_young*
  rename _prob_at2_1_`i'* P_`name'_F_old*
  rename _diff_prob_at2_1_`i'* Diff_`name'*
}

// Plot difference in CIFs
foreach name in $statenames {
  twoway (rarea Diff_`name'_lci Diff_`name'_uci temptime, color(red%30)) ///
         (line Diff_`name' temptime, color(red)) ///
         , name(diff_CIF_`name', replace)
}

