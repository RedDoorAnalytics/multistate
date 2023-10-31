
local drive /Users/Michael/My Drive/products/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

do ./build/buildmlib.do
mata mata clear

local drive /Users/Michael/My Drive/products/multistate
cd "`drive'"
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

clear all

use "`drive'/data/multistate_example",clear
set seed 98775
// keep if runiform()<0.2
// qui {


msset, id(pid) states(rfi osi) times(rf os) 
mat tmat = (.,1,2\.,.,3\.,.,.) 
// mat tmat = (.,.,.\1,.,2\3,.,.) 
// mat tmat = (.,1,2,.\.,.,.,3\.,.,.,.\.,.,.,.)
mat list tmat
stset _stop, enter(_start) failure(_status==1) scale(12)
tab size, gen(sz)

// gen _reset = _stop - _start
// stset _reset,failure(_status==1) scale(12) 
// stmerlin hormon age if _trans1==1, dist(cox) 
// est store m1

// stmerlin hormon age if _trans2==1, dist(cox)
// est store m2

// stmerlin hormon age  if _trans3==1, dist(cox)
// est store m3

stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==1, dist(rcs) df(3)     ///
         tvc(sz2 sz3 pr_1) dftvc(1)
// stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==1, dist(rcs) df(3)     ///
//          tvc(sz2 sz3 pr_1) dftvc(1)
// stmerlin hormon age if _trans1==1, dist(rp) df(3)
est store m4

stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==2, dist(weib)
// stmerlin hormon age if _trans2==1, dist(rp) df(3)
est store m5

stmerlin age sz2 sz3 nodes pr_1 hormon if _trans==3, dist(rcs) df(3)     ///
       tvc(pr_1) dftvc(1)
//stmerlin age hormon if _trans3==1, dist(rp) df(3)

est store m6
// }
cap range tvar 0 5 100

set seed 19347

timer clear
timer on 1
predictms , transmat(tmat) models(m5 m4 m6) 	///
			probability 						///
			timevar(tvar)						
		// 			ci from(2 3)
timer off 1

graphms

cap rename _prob* prob*

predictms , transmat(tmat) models(m5 m4 m6) 	///
			probability 						///
			at1(age 55)                      	///
			timevar(tvar) simulate n(100000)
// 			ci from(2 3)
