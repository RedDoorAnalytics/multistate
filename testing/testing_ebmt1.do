local drive Z

cd "`drive':\stdev\multistate"
adopath ++ ".\msset"
adopath ++ ".\predictms"
adopath ++ ".\stms"
clear all

do ./build/buildmlib.do

use ".\data\ebmt3",clear

msset, id(id) states(prstat rfsstat) times(prtime rfstime)
mat tmat = r(transmatrix)

stset _stop, enter(_start) failure(_status==1) scale(365.24)

//stpm2  if _trans1==1,  df(5) scale(h) lininit
streg if _trans1==1, dist(exp) nohr
stpm2  if _trans1==1,  df(1) scale(h) 
//strcs if _trans1==1, df(1) 
est store m4

streg  if _trans2==1, dist(exp) nohr
stpm2  if _trans2==1,  df(3) scale(h) 
//strcs if _trans2==1, df(5) 
est store m5

streg if _trans3==1, dist(exp) nohr
stpm2  if _trans3==1,  df(3) scale(h) 
//strcs  if _trans3==1, df(5) 
est store m6

cap drop tvar
range tvar 0 6 61
timer clear
timer on 1
predictms, transmat(tmat) models(m4 m5 m6) n(10000) ci //timevar(tvar) ci //graph
timer off 1
timer list 
