//local drive Z:\
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

use "./data/multistate_example",clear
set seed 98775
keep if runiform()<0.2
// qui {
msset, id(pid) states(rfi osi) times(rf os)  cr
mat tmat = (.,1,2\.,.,3\.,.,.) 
stset _stop, enter(_start) failure(_status==1) scale(12) 
tab size, gen(sz)

stmerlin hormon age if _trans1==1, dist(cox) 
est store m1

stmerlin hormon age if _trans2==1, dist(cox)
est store m2

stmerlin hormon age if _trans1==1, dist(rcs) df(3)
est store m4

stmerlin hormon age if _trans2==1, dist(rcs) df(3)
est store m5

// }
cap range tvar 0 5 100

set seed 1934
timer clear
timer on 1
// predictms , cr models(m1 m2) 					///
// 			probability							///
// 			at1(age 55) 						///
// 			timevar(tvar) simulate //ltruncated(1)
timer off 1

cap rename _prob* prob*
cap rename _los* los*

timer on 2
predictms , cr models(m4 m5) 					///
			probability							///
			at1(age 55) 						///
			timevar(tvar)						///
			los ///
			simulate latent n(10000) //ltruncated(1)
timer off 2

// merlin (_t hormon age if _trans==1, family(rp, df(3) failure(_d))) ///
// 		(_t hormon age if _trans==2, family(rp, df(3) failure(_d)))
		
// predict cif1, cif at(age 55) zeros
// predict cif2, cif at(age 55) zeros outcome(2)
