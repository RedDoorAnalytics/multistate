//local drive Z:/
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear


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

set seed 245987

use "./data/multistate_example",clear

msset, id(pid) states(rfi osi) times(rf os) cov(age chemo)

mat tmat = r(transmatrix)

stset _stop, enter(_start) failure(_status==1) scale(12) 

stmerlin 	chemo_trans1 chemo_trans2 chemo_trans3 age_trans1 age_trans2 age_trans3 _trans2 _trans3, dist(w) ///
			tvc(_trans2 _trans3) dftvc(1 1)

range tvar 0 10 20

predictms, 	transmat(tmat) 				///
			at1(chemo 1 age 55) 		///
			probability					///
			timevar(tvar)	
qui {
rename _prob* prob*			
stmerlin chemo age if _trans==1, dist(w)
est store m1
stmerlin chemo age if _trans==2, dist(w)	
est store m2
stmerlin chemo age if _trans==3, dist(w) tvc(age) dftvc(1)	
est store m3
}
predictms, 	transmat(tmat) 				///
			models(m1 m2 m3)			///
			at1(chemo 1 age 55) 		///
			probability					///
			timevar(tvar)		
