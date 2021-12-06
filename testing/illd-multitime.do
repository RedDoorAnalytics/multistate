//=============================================================================================================================//
// cert. script for merlin

//source paths
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/merlin/stmerlin"

//build mlib
clear all
do ./build/buildmlib.do
mata mata clear


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

do ./build/buildmlib.do
mata mata clear

use "`drive'/multistate/multistate/data/multistate_example",clear

//illness-death 
msset, id(pid) states(rfi osi) times(rf os)
mat tmat = r(transmatrix)

//stset on main timescale, time since entry to state 1, primary surgery
stset _stop, enter(_start) failure(_status==1) scale(12) 
tab size, gen(sz)

//trans 1: post-surgery to relapse
merlin (_t hormon if _trans==1, family(rp, df(3) failure(_d)))
est store m1
	
//trans 2: post-surgery to death
merlin (_t hormon if _trans==2, family(rp, df(3) failure(_d)))
est store m2

//trans 3: relapse to death
// - main timescale is still time since surgery, so we have delayed entry
merlin (_t 	hormon													///
			if _trans==3											///
			, family(rp, df(3) failure(_d) ltruncated(_t0)))		///
			, 
est store m3a

//trans 3: relapse to death with multiple timescales
// - main timescale is still time since surgery, so we have delayed entry
// - now including time since relapse as a secondary timescale
merlin (_t 	hormon																///
			rcs(_t, df(1) moffset(_t0))											///
			if _trans==3														///
			, family(rp, df(3) failure(_d) ltruncated(_t0)) timevar(_t))		///
			, 
// mat b2 = e(b)
// mat b2[1,2] = 1
// erepost b=b2
est store m3b
// - for the second timescale, you must use the response time and the ltruncated() time variables in the definition of the rcs() element (you could also use fp() with moffset() option)
// - moffset() adds the negative of _t0 from _t, and then builds the spline.


//trans 3: relapse to death with multiple timescales
// - main timescale is still time since surgery, so we have delayed entry
// - now including time since relapse as a secondary timescale
merlin (_t 	hormon																///
			rcs(_t, df(1) moffset(_t0))											///
			if _trans==3														///
			, family(rp, df(3) failure(_d) ltruncated(_t0)) timevar(_t))		///
			, 
// mat b2 = e(b)
// mat b2[1,2] = 1
// erepost b=b2
est store m3c

range time 0 10 100

predictms , 	transmat(tmat)		///
				probability			///
				models(m1 m2 m3c) 	///
				at1(hormon 1)		///
				timevar(time)  		///
				devcode1(41bsjdh82e198ndu3)
								
rename _prob* prob*

predictms , 	transmat(tmat)		///
				probability			///
				models(m1 m2 m3c) 	///
				at1(hormon 1)		///
				timevar(time)  		///
				simulate			///
				n(1000000)				///
				devcode1(41bsjdh82e198ndu3)

su _prob* prob*

merlin 	(_t hormon if _trans==1, family(rp, df(3) failure(_d)))						///
		(_t hormon if _trans==2, family(rp, df(3) failure(_d)))						///
		(_t hormon 																	///
			rcs(_t, df(1) moffset(_t0)) 											///
			if _trans==3															///
			, family(rp, df(3) failure(_d) ltruncated(_t0)))						///
			, 
			
			
// // mat b2 = e(b)
// // mat b2[1,8] = 1
// // erepost b=b2
			
predict prob1, at(hormon 1) transmatrix(tmat) transprob(1) devcode4(310780) timevar(time) 
predict prob2, at(hormon 1) transmatrix(tmat) transprob(2) devcode4(310780) timevar(time) 
predict prob3, at(hormon 1) transmatrix(tmat) transprob(3) devcode4(310780) timevar(time) 

// // predict los1, transmatrix(tmat) los(1) devcode4(310780) timevar(time) 
// // predict los2, transmatrix(tmat) los(2) devcode4(310780) timevar(time) 
// // predict los3, transmatrix(tmat) los(3) devcode4(310780) timevar(time) 

twoway (line _prob_at1_1_* time) (line prob_at1_1_* time) , name(g1,replace)
twoway (line prob_at1_1_* time) (line prob1 prob2 prob3 time), name(g2,replace)
// twoway (line prob_at1_1_* time) (line prob1 prob2 prob3 time), name(g3,replace)

