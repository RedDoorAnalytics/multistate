// using msboxes with Rott2
// uses both 3 state and 4 state example

//local drive Z:\
//local drive /Users/Michael/Documents
local drive c:
cd "`drive'/multistate/multistate"
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./stms"
adopath ++ "./msboxes"
adopath ++ "./msaj"
clear all

do ./build/buildmlib.do
mata mata clear

cd "`drive'/multistate/multistate"


clear all
use "./data/multistate_example",clear

msset, id(pid) states(rfi osi) times(rf os)
matrix list r(freqmatrix)
mat tmat = r(transmatrix)

msboxes, transmatrix(tmat) id(pid) ///
  xvalues(0.2 0.7 0.45) ///
  yvalues(0.7 0.7 0.2) ///
  statenames("Surgery" "Recurrence" "Dead") ///
  interactive jsonpath("C:\multistate\multistate\msboxes\Interactive")

stset _stop, enter(_start) failure(_status==1) scale(12) 

tab size, gen(size)
forvalues i = 1/3 {
//	streg size2 size3 if _trans==`i', dist(weibull)
	stpm2 size2 size3 if _trans==`i', scale(h) tvc(size2 size3) df(4) dftvc(3)
	estimates store m`i'
}

range timevar 0 7 500
predictms, interactive jsonpath("C:/multistate/multistate/msboxes/Interactive") ///
	models(m1 m2 m3) transmat(tmat) timevar(timevar) ///
	at1(size2 0 size3 0) ///
	at2(size2 1 size3 0) ///
	at3(size2 0 size3 1) aj	

use "./data/multistate_example",clear
matrix tmat = (.,1,2,. \ ///
  .,.,.,3 \ ///
  .,.,.,. \ ///
  .,.,.,.)
matrix list tmat  
  
msset, id(pid) states(rfi osi osi) times(rf os os) transmatrix(tmat)

msboxes, transmatrix(tmat) id(pid) ///
  xvalues(0.2 0.7 0.2 0.7) ///
  yvalues(0.7 0.7 0.2 0.2) ///
  statenames("Surgery" "Recurrence" "Dead" "Dead") ///
  boxheight(0.2) yrange(0.09 0.81) ysize(3) ///
  interactive jsonpath("C:\multistate\multistate\msboxes\Interactive")

stset _stop, enter(_start) failure(_status==1) scale(12) 

tab size, gen(size)
forvalues i = 1/3 {
	stpm2 size2 size3 if _trans==`i', scale(hazard) df(4) //dist(weibull)
	estimates store m`i'
}

range timevar 0 7 500
predictms, interactive jsonpath("C:/multistate/multistate/msboxes/Interactive") /// 
	models(m1 m2 m3) transmat(tmat) timevar(timevar) ///
	at1(size2 0 size3 0) ///
	at2(size2 1 size3 0) ///
	at3(size2 0 size3 1)  aj 
  
