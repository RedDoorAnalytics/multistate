local drive /Users/Michael/Documents
//local drive c:
cd "`drive'/multistate/multistate"
adopath ++ "."
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./stms"
adopath ++ "./msboxes"
adopath ++ "./msaj"
adopath ++ "./graphms"
clear all

tr:do ./build/buildmlib.do
mata mata clear

use "`drive'/multistate/multistate/data/multistate_example",clear
cd

msset, id(pid) states(rfi osi) times(rf os)

mat tmat = r(transmatrix)
mat list r(freqmatrix)

stset _stop, enter(_start) failure(_status==1) scale(12) 
tab size, gen(sz)

mat tmat = (.,1,2\.,.,3\.,.,.)

// //trans 1
// merlin (_t age nodes if _trans==1, family(rp, df(1) failure(_d)))
// est store m1
// //trans 2
// merlin (_t age nodes if _trans==2, family(rp, df(1) failure(_d)))
// est store m2
// //trans 3
// merlin (_t 	age	nodes	if _trans==3								///
// 			, family(rp, df(1) failure(_d) ltruncated(_t0)))	
// est store m3

range time 0 10 10
// set seed 1934
gen test = _n<100
// galahad , 	transmat(tmat)				///
// 			models(m1 m2 m3) 			///
// 			at1(age 45)					///
// 			timevar(time)  				///
// 			prob aj stand standif(test)



// merlin (_t 	age#_trans1 age#_trans2 age#_trans3 			///
// 			nodes#_trans1 nodes#_trans2 nodes#_trans3 		///
// 			_trans2 _trans3									///
// 			_trans2#rcs(_t, log event df(1) orthog) 		///
// 			_trans3#rcs(_t, log event df(1) orthog), 		///
// 		family(rp, df(1) failure(_d) ltruncated(_t0))	timevar(_t))
		
merlin (_t 	age#_trans1 age#_trans2 age#_trans3 ///
			nodes#_trans1 nodes#_trans2 nodes#_trans3			///
			_trans2 _trans3,								///
		family(weib, failure(_d) ltruncated(_t0)) )		

set seed 1934
galahad , 	transmat(tmat)				///
			at1(age 45)					///
			timevar(time)  				///
			n(100) ///
			prob stand standif(test)
		
rename _prob* prob*

forvalues i=1/3 {
	gen age_trans`i' = age * _trans`i'
	gen nodes_trans`i' = nodes * _trans`i'
	
}
stmerlin age_trans1 age_trans2 age_trans3 nodes_trans1 nodes_trans2 nodes_trans3 _trans2 _trans3, dist(weibull)

set seed 1934
galahad , 	transmat(tmat)				///
			at1(age 45)					///
			timevar(time)  				///
			n(100) ///
			prob  stand standif(test)
