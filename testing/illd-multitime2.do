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

tr:do ./build/buildmlib.do
mata mata clear

use "http://fmwww.bc.edu/repec/bocode/m/multistate_example", clear
msset, id(pid) states(rfi osi) times(rf os)
mat tmat = r(transmatrix)

//stset on main timescale, time since entry to state 1, primary surgery
stset _stop, enter(_start) failure(_status==1) scale(12) 
tab size, gen(sz)

//trans 1: post-surgery to relapse
merlin (_t hormon if _trans==1, family(rp, df(1) failure(_d)))
est store m1
	
//trans 2: post-surgery to death
merlin (_t hormon if _trans==2, family(rp, df(1) failure(_d)))
est store m2

//trans 3: relapse to death with clock reset + dynamic time of entry (on time since study start timescale)
// - we need to pretend that clock forward is used, hence the ltruncated() option
// - but we actually reset in the rcs() element, so the main timescale is time since entry
// - ltruncated() needs to be there so it gets updated from previous states
// - can now directly include time of entry (on previous main timescale), which now gets appropriately updated
merlin (_t 	hormon															///
			rcs(_t0, df(1))													///	-time of relapse (on previous timescale)-
			rcs(_t, df(3) orthog log event moffset(_t0))					///	-main timescale, time since entry-
			if _trans==3													///
			, family(logchazard, failure(_d) ltruncated(_t0)) timevar(_t))	///	-must pretend that we use clock forward-
			, diff
est store m3


range time 0 10 100

//multiple timescales
predictms , 	transmat(tmat)		///
				models(m1 m2 m3) 	///
				at1(hormon 1)		///
				timevar(time)  		///
				n(100000) 			///
				seed(1934) 			///
				tsreset(3)			// -must tell galahad that we are resetting in transition 3 model-

ghgraph
