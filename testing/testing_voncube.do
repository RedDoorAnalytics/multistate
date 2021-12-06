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

use `drive'/multistate/teaching/wengen/data/VonCube, clear 

// transition matrix for extended illness death model
matrix transmat = (., 1, 2, 3, ., . \ ///
	., ., ., ., 4, 5 \ ///
	., ., ., ., ., . \ ///	
	., ., ., ., ., . \ ///
	., ., ., ., ., . \ ///
	., ., ., ., ., .)

msset, id(admid) states(status2 status3 status4 status5 status6) ///
	times(t12 t13 t14 t25 t26) transmat(transmat) 				

stset _stop, enter(_start) failure(_status = 1) 

// timevar for predictions
range timevar 0 80 100
qui {
// Fit exponential model for each transition
forvalues i = 1/5 {
  display "Transition `i'"
  stmerlin if _trans == `i', dist(exp) nolog
  estimates store m_exp`i'
}

// Fit stpm2 model with 3 df for each transition
forvalues i = 1/5 {
  display "Transition `i'"
  stmerlin if _trans==`i', dist(rp) df(3)
  estimates store m_rp`i'
  predict h_rp`i', hazard timevar(timevar) 
}
}

set seed 923875
predictms, probability transmatrix(transmat)  ///
  models(m_exp1 m_exp2 m_exp3 m_exp4 m_exp5) ///
  n(10000) timevar(timevar) 
su _prob*

set seed 923875
predictms, probability transmatrix(transmat) ///
  models(m_rp1 m_rp2 m_rp3 m_rp4 m_rp5) ///
  timevar(timevar) n(10000)
su _prob*

/*// functions of transition probabilities
mata:
real matrix ProbDead(M)
{
	p4 = ms_user_prob(M,4)
	p6 = ms_user_prob(M,6)
	return(p4 + p6)
}
end


predictms, transmatrix(transmat) ///
	models(m_exp1 m_exp2 m_exp3 m_exp4 m_exp5) ///
	n(10000) timevar(timevar) userf(ProbDead) ci

			
// now do various calculations at the same time.	
mata:
real matrix GenRes(M)
{
// extract transition probabilities
	p1 = ms_user_prob(M,1)
	p2 = ms_user_prob(M,2)
	p3 = ms_user_prob(M,3)
	p4 = ms_user_prob(M,4)
	p5 = ms_user_prob(M,5)
	p6 = ms_user_prob(M,6)
// Do some calculations	
	Pdead = p4 + p6
	Pdead_HAI0 = p4:/(p1 :+ p3 :+ p4)
	Pdead_HAI1 = p6:/(p2 :+ p5 :+ p6)
	AM = Pdead_HAI1 :- Pdead_HAI0
	PAF = (Pdead :- Pdead_HAI0) :/ Pdead
// Return results	
	return(Pdead,Pdead_HAI0,Pdead_HAI1,AM,PAF)
}
end		


predictms, transmatrix(transmat) ///
	models(m_rp1 m_rp2 m_rp3 m_rp4 m_rp5) ///
	n(10000) timevar(timevar) userf(GenRes) ci
*/



		
