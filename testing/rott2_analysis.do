adopath ++ "C:\stdev\multistate\version_dev"
cd "C:\stdev\multistate\version_dev"
clear all

buildmlib

use "C:\multistate\Data\rott2",clear

msset, id(pid) states(rfi osi) times(rf os)
//!!add ind var for any records changed!!

mat tmat = r(transmatrix)
mat list r(freqmatrix)

stset _stop, enter(_start) failure(_status==1) scale(12) 

tab size, gen(sz)

//trans 1
stpm2 age sz2 sz3 nodes pr_1 hormon if _trans==1, df(3) scale(h)

stpm2 age sz2 sz3 nodes pr_1 hormon if _trans==1, df(3) scale(h)	///
	tvc(sz2 sz3 pr_1) dftvc(1)
est store m1
	
//trans 2
streg age sz2 sz3 nodes pr_1 hormon if _trans==2, dist(weib) 
est store m2

//trans 3
stpm2 age sz2 sz3 nodes pr_1 hormon if _trans==3, df(3) scale(h)

stpm2 age sz2 sz3 nodes pr_1 hormon if _trans==3, df(3) scale(h)	///
	tvc(pr_1) dftvc(1)
est store m3


//trans probs
pr drop _all
set seed 4578789
cap drop prob*	


predictms, 	transmat(tmat) 						///
			at2(age 54 sz2 1 pr_1 3)				///
			at(age 65 sz2 1 pr_1 4)			///
			model1(m1) model2(m2) model3(m3) //ci m(100) n(1000) los


