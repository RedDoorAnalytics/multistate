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


// Illness-death example, clock-reset, only one main timescale, time since entry to each state

use "http://fmwww.bc.edu/repec/bocode/m/multistate_example", clear

msset, id(pid) states(rfi osi) times(rf os)

mat tmat = r(transmatrix)
mat list r(freqmatrix)
gen _newt = _stop - _start

stset _newt, failure(_status==1) scale(12) 
tab size, gen(sz)

mat tmat = (.,1,2\.,.,3\.,.,.)

//trans 1
merlin (_t age if _trans==1, family(rp, df(1) failure(_d)))
est store m1
	
//trans 2
merlin (_t age if _trans==2, family(rp, df(1) failure(_d)))
est store m2

//trans 3
merlin (_t 	age	if _trans==3					///
			, family(rp, df(1) failure(_d)))
est store m3

range time 0 10 100

galahad , 	transmat(tmat)		///
			models(m1 m2 m3) 	///
			at1(age 45)			///
			timevar(time)  		///
			reset				///
			los
