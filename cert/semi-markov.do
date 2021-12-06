
use "./data/multistate_example",clear

msset, id(pid) states(rfi osi) times(rf os) cov(age chemo)

mat tmat = r(transmatrix)
mat list r(freqmatrix)

gen _time = _stop - _start

stset _time, failure(_status==1) scale(12) 

tab size, gen(sz)

//=====================================================================================================================================//
// Base models, no covariates -> transition probabilities and ci's

//Stacked data

stmerlin _trans2 _trans3, dist(rp) df(3) 

predictms , probability transmat(tmat) reset
predictms , probability transmat(tmat) ci  reset m(50)

stmerlin _trans2 _trans3, dist(weib) 

predictms , probability transmat(tmat) reset
predictms , probability transmat(tmat) ci reset m(50)

//Separate models

stmerlin if _trans1==1, dist(rp) df(3)
est store m1
stmerlin if _trans2==1, dist(rp) df(3)
est store m2
stmerlin if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) reset
predictms , probability transmat(tmat) models(m1 m2 m3) ci reset m(50)

stmerlin if _trans1==1, dist(weib)
est store m1
stmerlin if _trans2==1, dist(weib)
est store m2
stmerlin if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) reset
predictms , probability transmat(tmat) models(m1 m2 m3) ci reset m(50)



//=====================================================================================================================================//
// Covariates -> transition probabilities and ci's

//Stacked data

stmerlin age sz2 sz3 _trans2 _trans3, dist(rp) df(3)

predictms , probability transmat(tmat) at1(age 55 sz2 1) reset
predictms , probability transmat(tmat) at1(age 55 sz2 1) ci reset m(50)

stmerlin  age sz2 sz3 _trans2 _trans3, dist(weib) 

predictms , probability transmat(tmat) at1(age 55 sz2 1) reset
predictms , probability transmat(tmat) at1(age 55 sz2 1) ci reset m(50)

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) ci reset m(50)

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) ci reset m(50)



//=====================================================================================================================================//
// Base models, no covariates -> LoS and ci's

//Separate models

stmerlin if _trans1==1, dist(rp) df(3)
est store m1
stmerlin if _trans2==1, dist(rp) df(3)
est store m2
stmerlin if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) los reset
predictms , probability transmat(tmat) models(m1 m2 m3) los ci reset m(50)

stmerlin if _trans1==1, dist(weib)
est store m1
stmerlin if _trans2==1, dist(weib)
est store m2
stmerlin if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) los reset
predictms , probability transmat(tmat) models(m1 m2 m3) los ci reset m(50)



//=====================================================================================================================================//
// Covariates -> LoS and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los ci reset m(50)

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los ci reset m(50)



//=====================================================================================================================================//
// Covariates -> transition probability differences and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) reset diff
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ci reset diff m(50)

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) reset diff
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ci reset diff m(50)

//=====================================================================================================================================//
// Covariates -> transition probability ratios and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio ci reset m(50)

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio ci reset m(50)


//=====================================================================================================================================//
// Covariates -> LoS differences and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los reset diff
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ci reset diff m(50)

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ci reset m(50)

//=====================================================================================================================================//
// Covariates -> LoS ratios and ci's

//Separate models

stmerlin age sz2 sz3 if _trans1==1, dist(rp) df(3)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(rp) df(3)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(rp) df(3)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio los reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio los ci reset m(50)

stmerlin age sz2 sz3 if _trans1==1, dist(weib)
est store m1
stmerlin age sz2 sz3 if _trans2==1, dist(weib)
est store m2
stmerlin age sz2 sz3 if _trans3==1, dist(weib)
est store m3

predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio los reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio los ci reset m(50)



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
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) reset simulate
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) ci reset m(50)

//-> LoS and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) los ci reset m(50)

//-> diffs in trans probs and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) reset diff
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ci reset diff m(50)

//-> ratios of trans probs and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) ratio ci reset m(50)

//-> diffs in LoS and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los reset diff
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ci reset diff m(50)

//-> ratios of LoS and ci's
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ratio reset
predictms , probability transmat(tmat) models(m1 m2 m3) at1(age 55 sz2 1) at2(age 55 sz3 1) los ratio ci reset m(50)



