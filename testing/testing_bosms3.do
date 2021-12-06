local drive Z

cd "`drive':\stdev\multistate"
adopath ++ ".\msset"
adopath ++ ".\predictms"
adopath ++ ".\stms"
clear all

do ./build/buildmlib.do

use ".\data\bosms3",clear
//data already msset

stset tstop, enter(tstart) f(status==1)

tab trans, gen(_trans)
stpm2 _trans2 _trans3, scale(hazard) df(5)  
streg _trans2 _trans3, dist(weib) //anc(_trans2 _trans3)

mat tmat = (.,1,2\.,.,3\.,.,.)

timer clear
timer on 1
predictms, transmat(tmat) from(1 2)  n(10000) //ci //timevar(tvar) ci //graph
timer off 1
timer on 2
predictms, transmat(tmat) from(1 2)  n(100000) //ci //timevar(tvar) ci //graph
timer off 2
timer on 3
predictms, transmat(tmat) from(1 2)  n(1000000) //ci //timevar(tvar) ci //graph
timer off 3
timer on 4
predictms, transmat(tmat) from(1 2)  n(10000000) //ci //timevar(tvar) ci //graph
timer off 4
timer list 

