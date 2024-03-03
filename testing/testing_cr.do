//local drive Z:\
local drive /Users/Michael/My Drive
cd "`drive'/software/multistate"
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

stmerlin hormon age if _trans1==1, dist(rp) df(3) 
est store m1

stmerlin hormon age if _trans2==1, dist(rp) df(3)
est store m2

// stmerlin hormon age if _trans1==1, dist(rcs) df(3)
// est store m4
//
// stmerlin hormon age if _trans2==1, dist(rcs) df(3)
// est store m5

// }
cap drop tvar
cap range tvar 0 5 100

set seed 1934
timer clear
timer on 1
predictms , cr models(m1 m2) 		///
	probability			///
	at1(age 55) 			///
	timevar(tvar) 			///
	standardise
timer off 1

timer on 2
predictms , cr models(m1 m2) 		///
	probability			///
	at1(age 55) 			///
	timevar(tvar) 			///
	standardise aj
timer off 2

timer list
