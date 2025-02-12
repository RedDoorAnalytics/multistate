
local drive /Users/michael/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software

local drive1 `drive'/merlin
cd "`drive1'"
adopath ++ "`drive1'"
adopath ++ "`drive1'/merlin"
adopath ++ "`drive1'/stmerlin"
clear all

do ./build/buildmlib.do
mata mata clear

local drive2 `drive'/multistate
cd "`drive2'"
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

use "`drive2'/data/multistate_example",clear
set seed 98775

// mat tmat = (.,1,2\.,.,3\.,.,.) 
mat tmat = (.,1,2,.\.,.,.,3\.,.,.,.\.,.,.,.)

msset, id(pid) states(rfi osi osi) times(rf os os) transmat(tmat)


stset _stop, enter(_start) failure(_status==1) scale(12)
tab size, gen(sz)

stmerlin hormon age if _trans1==1, dist(rcs) df(3) 
est store m1

stmerlin hormon age if _trans2==1, dist(rp) df(3) 
est store m2

stmerlin _t0 hormon age  if _trans3==1, dist(rcs) df(3) 
est store m3


cap range tvar 0 5 100

predictms , transmat(tmat) models(m1 m2 m3) 		///
		probability at1(hormon 1 age 44)	///
		timevar(tvar) ci

cap rename _prob* prob*
