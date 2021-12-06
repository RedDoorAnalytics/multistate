
//source paths
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/merlin/stmerlin"

//build mlib
clear all
tr:do ./build/buildmlib.do
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

use "./data/multistate_example", clear
set seed 897
// keep if runiform()<0.1
// qui {
msset, id(pid) states(rfi osi) times(rf os) 
mat tmat = (.,1,2\.,.,3\.,.,.) 
stset _stop, enter(_start) failure(_status==1) scale(12)
tab size, gen(sz)

stmerlin hormon age if _trans1==1, dist(cox)
predict bh1, basehazard 
est store m1
stcox hormon age if _trans1==1
predict bh2, basehc

stmerlin hormon age if _trans2==1, dist(cox)
est store m2

stmerlin hormon age if _trans1==1, dist(rp) df(3)
est store m3

stmerlin hormon age if _trans2==1, dist(rp) df(3)
est store m4

// merlin 	(_t hormon age if _trans1==1, family(cox, failure(_d)))	///
// 		(_t hormon age if _trans2==1, family(cox, failure(_d)))

// predict cif1, cif outcome(1) at(age 55) zeros //timevar(tvar)
// predict cif2, cif outcome(2) at(age 55) zeros //timevar(tvar)
		
stcox hormon age if _trans3==1
predict ch1, basehc
stmerlin hormon age  if _trans3==1, dist(cox)
predict ch2, basehazard
est store m5

stmerlin hormon age  if _trans3==1, dist(rp) df(3)
est store m6

cap range tvar 0 16 1000

predictms , transmat(tmat) 			///
			models(m1 m2 m5) 		///
			probability				///
			at1(hormon 1 age 55) 	///
			timevar(tvar) 			///
			n(1000)				///
 			devcode9(98342h2r6fwgi240wbg8gghgffghhgfh)
rename _prob* prob*

predictms , transmat(tmat) 			///
			models(m3 m4 m6) 		///
			probability				///
			at1(hormon 1 age 55) 	///
			timevar(tvar) 			///

twoway (line _prob* tvar) (line prob* tvar), legend(off)
