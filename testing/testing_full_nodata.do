local drive Z

cd "`drive':\stdev\multistate"
adopath ++ ".\msset"
adopath ++ ".\predictms"
adopath ++ ".\stms"
clear all

do ./build/buildmlib.do

use ".\data\multistate_example",clear

msset, id(pid) states(rfi osi) times(rf os) cov(age)

mat tmat = r(transmatrix)
mat list r(freqmatrix)

stset _stop, enter(_start) failure(_status==1) scale(12) 

tab size, gen(sz)

stpm2 chemo age_trans1 age_trans2 age_trans3 _trans2 _trans3, scale(h) df(3)

clear

timer clear
timer on 1

set obs 100
range tvar 0.1 10 100
predictms , transmat(tmat) seed(28729)				///
			at(age 55) out timevar(tvar)
			
timer off 1
timer list

su pred_*


