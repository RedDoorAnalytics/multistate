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
//!!add ind var for any records changed!!

mat tmat = r(transmatrix)
mat list r(freqmatrix)

stset _stop, enter(_start) failure(_status==1) scale(12) 

tab size, gen(sz)

//trans 1
stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==1, dist(rp) df(3) 

stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==1, dist(rp) df(3) ///
			tvc(sz2 sz3 pr_1) dftvc(1)
est store m1
	
//trans 2
stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==2, dist(weib)
est store m2

//trans 3
stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==3, dist(rp) df(1)

stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==3, dist(rp) df(3) ///
	tvc(pr_1) dftvc(1)
est store m3


predictms , transmatrix(tmat) models(m1 m2 m3)		///
			n(500) at1(hormon 0) at2(hormon 1) 		///
			probability los 

