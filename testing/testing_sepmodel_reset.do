
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

do ./build/buildmlib.do
mata mata clear

local drive /Users/Michael/Documents
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
// keep if runiform()<0.2

// qui {
msset, id(pid) states(rfi osi) times(rf os)

mat tmat = r(transmatrix)
mat list r(freqmatrix)

gen reset = _stop-_start

stset reset, failure(_status==1) scale(12) 

tab size, gen(sz)

stmerlin hormon age  if _trans1==1, dist(cox)
est store m1

stmerlin hormon age  if _trans2==1, dist(cox) 
est store m2

stmerlin hormon age if _trans3==1, dist(cox)
est store m3

stmerlin hormon age if _trans1==1, dist(rp) df(3)
est store m4

stmerlin hormon age if _trans2==1, dist(rp) df(3)
est store m5

stmerlin hormon age  if _trans3==1, dist(rp) df(3)
est store m6

// }

cap range temptime 0 5 100

set seed 1934
timer clear
timer on 1
predictms , transmat(tmat) models(m4 m5 m6) 	///
			probability							///
			at1(age 55) 						///
			timevar(temptime)					///
			ci reset //los visit ci m(20)
timer off 1

rename _prob* prob*

timer on 2
// predictms , transmat(tmat) models(m1 m2 m3) 	///
// 			probability							///
// 			at1(age 55) 						///
// 			timevar(temptime)					///
// 			n(100000) reset //los visit ci m(20)
predictms , transmat(tmat) models(m4 m5 m6) 	///
			probability							///
			at1(age 55) 						///
			timevar(temptime)					///
			n(1000000) reset latent //los visit ci m(20)

timer off 2	
timer list

su _prob*
