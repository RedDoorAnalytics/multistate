

use "./data/multistate_example",clear

msset, id(pid) states(rfi osi) times(rf os) cov(age chemo)

mat tmat = r(transmatrix)
mat list r(freqmatrix)

stset _stop, enter(_start) failure(_status==1) scale(12) 

tab size, gen(sz)

//=====================================================================================================================================//
// Base models, no covariates -> transition probabilities and ci's

//Stacked data
stmerlin _trans2 _trans3, dist(rp) df(3)

predictms, probability transmat(tmat)
predictms, probability transmat(tmat) ci n(1000) m(50) seed(91850943) simulate

stmerlin _trans2 _trans3, dist(weib) 

predictms, probability transmat(tmat)
predictms, probability transmat(tmat) ci

//Separate models

stmerlin if _trans1==1, dist(rp) df(3)
est store m1
stmerlin if _trans2==1, dist(rp) df(3)
est store m2
stmerlin if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) 
predictms , probability transmat(tmat) models(m1 m2 m3) ci 

predictms , probability transmat(tmat) models(m1 m2 m3) aj
predictms , probability transmat(tmat) models(m1 m2 m3) ci aj

stmerlin if _trans1==1, dist(weib)
est store m1
stmerlin if _trans2==1, dist(weib)
est store m2
stmerlin if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3)
predictms , probability transmat(tmat) models(m1 m2 m3) ci 



//=====================================================================================================================================//
// Covariates -> transition probabilities and ci's

//Stacked data

stmerlin age sz2 sz3 _trans2 _trans3, df(3) dist(rp)

predictms, probability transmat(tmat) at1(age 55 sz2 1)
predictms, probability transmat(tmat) at1(age 55 sz2 1) ci 

stmerlin  age sz2 sz3 _trans2 _trans3, dist(weib) 

predictms, probability transmat(tmat) at1(age 55 sz2 1) 
predictms, probability transmat(tmat) at1(age 55 sz2 1) ci 

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) ci 

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) ci 



//=====================================================================================================================================//
// Base models, no covariates -> LoS and ci's

//Separate models

stmerlin if _trans1==1, dist(rp) df(3)
est store m1
stmerlin if _trans2==1, dist(rp) df(3)
est store m2
stmerlin if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) los 
predictms , probability transmat(tmat) models(m1 m2 m3) los ci 

stmerlin if _trans1==1, dist(weib)
est store m1
stmerlin if _trans2==1, dist(weib)
est store m2
stmerlin if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) los 
predictms , probability transmat(tmat) models(m1 m2 m3) los ci 



//=====================================================================================================================================//
// Covariates -> LoS and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los ci 

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los ci 



//=====================================================================================================================================//
// Covariates -> transition probability differences and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ci 

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ci 

//=====================================================================================================================================//
// Covariates -> transition probability ratios and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio ci 

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio ci 


//=====================================================================================================================================//
// Covariates -> LoS differences and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los diff
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ci  diff

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los diff
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ci diff

//=====================================================================================================================================//
// Covariates -> LoS ratios and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio los 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio los ci 

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio los 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio los ci 



//=====================================================================================================================================//
// Different models 

//trans 1
stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==1, df(3) dist(rp)	///
	tvc(sz2 sz3 pr_1) dftvc(1)
est store m1
	
//trans 2
stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==2, dist(weib) 
est store m2

//trans 3
stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==3, df(3) dist(rp)	///
	tvc(pr_1) dftvc(1)
est store m3


//-> trans probs and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) ci 

//-> LoS and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los ci 

//-> diffs in trans probs and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) diff
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ci diff

//-> ratios of trans probs and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio ci 

//-> diffs in LoS and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los diff
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ci diff

//-> ratios of LoS and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ratio 
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ratio ci 




//=====================================================================================================================================//
// User function


mata:
real matrix myf(S)
{
	p1 = ms_user_prob(S,1)
 	p2 = ms_user_prob(S,2)	
 	p3 = ms_user_prob(S,3)	
	return(p1,p2,p3)
}
end

range tvar 0 5 100

predictms , probability at1() at2(age 35) transmat(tmat) models(m1 m2 m3) userf(myf) timevar(tvar) 
predictms , probability at1() at2(age 35) transmat(tmat) models(m1 m2 m3) userf(myf) timevar(tvar) ci 

