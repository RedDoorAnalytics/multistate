// using msboxes with Rott2
// uses both 3 state and 4 state example

local drive c:\
//local drive z:\
//local drive /Users/Michael/Documents

cd "`drive'/multistate/multistate"
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./msboxes"

clear all

use "./data/multistate_example",clear

msset, id(pid) states(rfi osi) times(rf os)
matrix list r(freqmatrix)
mat tmat = r(transmatrix)

msboxes, transmatrix(tmat) id(pid) ///
  xvalues(0.2 0.7 0.45) ///
  yvalues(0.7 0.7 0.2) ///
  statenames("Surgery bob" "Recurrence" "Dead") 


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
  boxheight(0.2) yrange(0.09 0.81) ysize(3)

