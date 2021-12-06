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

mat tmat = (.,1,2 \	///
			3,.,4 \ ///
			.,.,.)

use ./data/sir,clear


//trans1
stset time if from==1, enter(start) f(to==2)
stmerlin, dist(weib)
est store m1

//trans2
stset time if from==1, enter(start) f(to==3)
stmerlin, dist(weib)
est store m2

//trans3
stset time if from==2, enter(start) f(to==1)
stmerlin, dist(weib)
est store m3

//trans4
stset time if from==2, enter(start) f(to==3)
stmerlin, dist(weib)
est store m4

range temptime 0 150 1000

set seed 498
// predictms, transmat(tmat) models(m1 m2 m3 m4) timevar(temptime) n(10000) from(1 2) los //ci
// su _prob* _los* 

predictms, transmat(tmat) models(m1 m2 m3 m4) timevar(temptime) ///
			prob aj from(1 2) los //ci
su _prob* _los* 
