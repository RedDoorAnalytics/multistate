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

tr:do ./build/buildmlib.do
mata mata clear
				 
input illt ills dt ds x1 x2 x3
1 1 5 1 1 6 0.1
1 0 1 1 1 5 0.2
6 1 9 1 1 4 0.1
6 1 7 1 0 3 0
8 0 8 1 0 2 0.5
9 1 12 1 0 1 0
end

gen id = _n
msset , times(illt dt) states(ills ds) cov(x1 x2 x3) id(id)

mat tmat = (.,1,2\.,.,3\.,.,.) 
stset _stop, enter(_start) failure(_status==1)

stmerlin x1 if _trans1==1, dist(cox)
est store m1

stmerlin x2 if _trans2==1, dist(cox)
est store m2

stmerlin x1 if _trans3==1, dist(cox)
est store m3
