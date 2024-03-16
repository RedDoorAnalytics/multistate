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

do ./build/buildmlib.do
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
cap drop tvar
cap range tvar 0 5 100
timer clear 
timer on 1
// predict s1, at(age 55) standardise rmst timevar(tvar)
timer off 1

stmerlin hormon age if _trans2==1, dist(rp) df(3)
est store m2



/*
- standardise will require timevar() option
- change to looping over timevar, which gets updated with each time, but 
    replicated across the global dataset touse
*/

timer on 2
predictms , cr models(m1 m2) 		///
	los			///
	at1(age 55) 			///
	timevar(tvar) 			///
	 standardise 
timer off 2

/*
tvar 100 obs
dataset - 1108
timer on second run - 22.41 second 

    Variable |        Obs        Mean    Std. dev.       Min        Max
-------------+---------------------------------------------------------
_prob_at1_~1 |        100    .7726354    .1304271   .5930705          1
_prob_at1_~2 |        100    .2247087    .1287533          0   .4009486
_prob_at1_~3 |        100    .0026559    .0017294          0   .0059809
*/
timer list
su _prob*

