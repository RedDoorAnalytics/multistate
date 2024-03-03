//local drive Z:\
local drive /Users/Michael/My Drive/software
//local drive c:
cd "`drive'/multistate"
adopath ++ "."
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
msset, id(pid) states(rfi osi) times(rf os) cov(age chemo)

keep if _trans1==1

stset _stop, enter(_start) f(_status) scale(12)

mat tmat = (.,1\.,.)

stmerlin chemo age, dist(weibull)

range tvar 1 15 100

predict s1, survival at(chemo 1 age 45) timevar(tvar)
predictms, 	at1(chemo 1 age 45) 	///
			singleevent 			///
			probability 			///
			timevar(tvar) ltruncated(1)

rename _prob* prob*

predictms, 	at1(chemo 1 age 45) 	///
			singleevent 			///
			probability 			///
			timevar(tvar) ltruncated(1) simulate

su _prob*
su prob*
