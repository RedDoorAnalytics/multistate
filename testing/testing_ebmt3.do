//local drive Z:\
local drive /Users/Michael/Documents

//merlin
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/merlin/stmerlin"
clear all
cd "`drive'/merlin"
do ./build/buildmlib.do
mata mata clear

//predictms
adopath ++ "`drive'/multistate/multistate/predictms"
adopath ++ "`drive'/multistate/multistate/msset"
adopath ++ "`drive'/multistate/multistate/graphms"
clear all
cd "`drive'/multistate/multistate"

do ./build/buildmlib.do
mata mata clear

use "/Users/Michael/Documents/multistate/multistate/data/ebmt3",clear

tab age, gen(ag)

msset, id(id) states(prstat rfsstat) times(prtime rfstime)
mat tmat = r(transmatrix)

stset _stop, enter(_start) failure(_status==1) scale(365.24)

cap drop tvar
range tvar 0 10 100

strcs ag2 ag3 if _trans1==1, df(4) nohr 
predict s0, surv zeros timevar(tvar)

stmerlin ag2 ag3 if _trans1==1,  dist(rcs) df(4)
predict s1, surv zeros timevar(tvar)
est store m4

stmerlin ag2 ag3 if _trans1==1,  dist(rp) df(4) diff
predict s2, surv zeros timevar(tvar)





// strcs ag2 ag3 if _trans2==1, df(4) nohr
stmerlin ag2 ag3 if _trans2==1,  dist(rcs) df(4) 
est store m5
predict s2, surv zeros timevar(tvar)

// strcs ag2 ag3 if _trans3==1, df(3) nohr
stmerlin ag2 ag3 if _trans3==1,  dist(rcs) df(3) //diff
est store m6
predict s3, surv zeros timevar(tvar)

timer clear
timer on 1
predictms, singleevent prob timevar(tvar) models(m4) 
graphms

