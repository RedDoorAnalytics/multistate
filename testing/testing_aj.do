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

// qui {
	msset, id(pid) states(rfi osi) times(rf os)
	mat tmat = r(transmatrix)

	stset _stop, enter(_start) failure(_status==1) scale(12) 

	cap range temptime 0 10 2000

	stmerlin age nodes pr_1 if _trans1==1, dist(weib) 
	est store m1
	
	stmerlin age nodes pr_1 if _trans2==1, dist(weib) 
	est store m2

	//stpm2 age nodes pr_1  if _trans3==1, df(3) scale(h) 
	stmerlin age nodes pr_1  if _trans3==1, dist(weib)
	est store m3
// }

timer clear
timer on 1
set seed 245987
predictms, 	transmat(tmat) models(m1 m2 m3) 	///
			at1(age 50 pr_1 3) 					///
			aj 									///
			timevar(temptime) 					///
			ci prob
rename _prob* prob*
predictms, 	transmat(tmat) models(m1 m2 m3) 	///
			at1(age 50 pr_1 3) 					///
			aj 									///
			timevar(temptime) 					///
			ci prob bootstrap
