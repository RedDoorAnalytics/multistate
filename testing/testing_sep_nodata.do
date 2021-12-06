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

use "./data/multistate_example",clear
msset, id(pid) states(rfi osi) times(rf os)

mat tmat = r(transmatrix)
mat list r(freqmatrix)

stset _stop, enter(_start) failure(_status==1) scale(12) 

tab size, gen(sz)

stmerlin chemo age if _trans1==1, dist(rp) df(3)
est store m1

stmerlin chemo age if _trans2==1, dist(rp) df(3)
est store m2

stmerlin chemo age if _trans3==1, dist(rp) df(3)
est store m3

clear
set obs 100
range tvar 0 10 100

predictms , transmat(tmat) seed(28729)				///
			models(m1 m2 m3) 						///
			probability								///
			at1(age 55) timevar(tvar) n(10000) 		///
			outsample
