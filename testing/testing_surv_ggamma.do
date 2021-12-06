//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/multistate/multistate"
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./stms"
adopath ++ "./msboxes"
adopath ++ "./msaj"
clear all

do ./build/buildmlib.do
mata mata clear

webuse brcancer,clear
stset rectime, f(censrec)

streg hormon, dist(ggamma)
est store m1

 predict s1, surv

predictms , models(m1) survival at1(hormon 0) at2(hormon 1)

su _prob*
twoway (scatter _prob_at1_1_1 _time)(scatter _prob_at2_1_1 _time)(scatter s1 _t)
