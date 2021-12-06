* Version 0.1, 11/02/2020

*****************************
*  			HAI  			*
*****************************

cap program drop hai
program define hai
	use "./data/hai", replace
	mat tmat = (.,1,2,3,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
end

* #1 Check _t is specified
hai
drop _t
cap msaj, transmatrix(tmat)
assert _rc == 198

* #2 Check _d is specified
hai
drop _d
cap msaj, transmatrix(tmat)
assert _rc == 198

* #3 Check _t0 is specified
hai
drop _t0
cap msaj, transmatrix(tmat)
assert _rc == 198

* #4 Check _st is specified
hai
drop _st
cap msaj, transmatrix(tmat)
assert _rc == 198

* #5 Check _trans is specified
hai
drop _trans
cap msaj, transmatrix(tmat)
assert _rc == 198

* #6 Check _to is specified if cr is specified
hai
drop _to
cap msaj, transmatrix(tmat) cr
assert _rc == 198

* #7 Check _to is specified if ci or se is specified
hai
drop _to
cap msaj, transmatrix(tmat) ci
assert _rc == 198
cap msaj, transmatrix(tmat) se
assert _rc == 198
cap msaj, transmatrix(tmat) se ci
assert _rc == 198

* #8 Check _from is specified if ci is specified
hai
drop _from
cap msaj, transmatrix(tmat) ci
assert _rc == 198
cap msaj, transmatrix(tmat) se
assert _rc == 198
cap msaj, transmatrix(tmat) ci se
assert _rc == 198

* #9 Check exactly one of transmatrix or cr is specified
hai 
cap msaj
assert _rc == 198
cap msaj, transmatrix(tmat) cr
assert _rc == 198

* #10 At least two states - tmat option
use "./data/hai", replace
mat tmat = (1)
cap msaj, transmatrix(tmat)
assert _rc == 198

* #11 At least two states - CR option
use "./data/hai", replace
replace _to = 1
cap msaj, cr
assert _rc == 198

* #12 Check from state is a state					
hai
cap msaj, transmatrix(tmat) from(0)
assert _rc == 198
cap msaj, transmatrix(tmat) from(-1)
assert _rc == 198
cap msaj, transmatrix(tmat) from(7)
assert _rc == 198

* #13 Check from state is not absorbing
hai 
cap msaj, transmatrix(tmat) from(3)
assert _rc == 198

* #14 If CR is specified from must be 1
hai
cap msaj, cr from(2)
assert _rc == 198

* #15 Check ltruncated time is positive					
hai
cap msaj, transmatrix(tmat) ltruncated(-1)
assert _rc == 198

* #16 Check exit time is after ltruncated time 
* (also checks exit is +ve and ltruncated before last)
hai
cap msaj, transmatrix(tmat) ltruncated(50) exit(49)
assert _rc == 198
cap msaj, transmatrix(tmat) ltruncated(50) exit(50)
assert _rc == 198
cap msaj, transmatrix(tmat) exit(0)
assert _rc == 198
cap msaj, transmatrix(tmat) ltruncated(82) 
assert _rc == 198
cap msaj, transmatrix(tmat) ltruncated(82.1) 
assert _rc == 198

* #17 Exit time is on or before max time										
hai
cap msaj, transmatrix(tmat) exit(82.1) 
assert _rc == 198

* #18 Check diagonal of tmat are missing
use "./data/hai", replace
mat tmat = (1,2,3,4,.,.\.,.,.,.,5,6\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
cap msaj, transmatrix(tmat)
assert _rc == 198

* #19 Check unique entries in transmatrix
use "./data/hai", replace
mat tmat = (.,1,2,2,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
cap msaj, transmatrix(tmat)
assert _rc == 198

* #20 Check tmat sequentially numbered
use "./data/hai", replace
mat tmat = (.,1,3,2,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
cap msaj, transmatrix(tmat)
assert _rc == 198


*****************************
*  			EBMT  			*
*****************************

* #21 Check it works with EBMT data, los and ltruncated / exit combo 0 but warning
use "./data/ebmt", replace
mat tmat = (.,1,2,.,3,4\.,.,.,5,6,7\.,.,.,8,9,10\.,.,.,.,11,12\.,.,.,.,.,.\.,.,.,.,.,.)
cap msaj, transmatrix(tmat) exit(5800) 
assert _rc == 0
