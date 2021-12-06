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

use "./data/multistate_example",clear

// qui {
msset, id(pid) states(rfi osi) times(rf os)

mat tmat = (.,1,2\.,.,3\.,.,.) 
mat list r(freqmatrix)

stset _stop, enter(_start) failure(_status==1) scale(12) 

tab size, gen(sz)


stmerlin age  if _trans1==1, dist(weib) 
est store m1

stmerlin age  if _trans2==1, dist(weib) 
est store m2

stmerlin age  if _trans3==1, dist(weib) 
est store m3


mata:
real matrix myf(S)
{
	p1 = ms_user_prob(S,1)
// 	p2 = ms_user_prob(S,2)	
// 	p3 = ms_user_prob(S,3)	
	return(p1)
}
end

range tvar 0 5 100

predictms , probability at1() at2(age 35) transmat(tmat) models(m1 m2 m3) userf(myf) timevar(tvar)
predictms , probability at1() at2(age 35) transmat(tmat) models(m1 m2 m3) userf(myf) aj ci  timevar(tvar)
			
