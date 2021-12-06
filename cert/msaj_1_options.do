* Version 0.1, 11/02/2020

*****************************
*  			EBMT  			*
*****************************

* Use EBMT data (includes censoring)
cap program drop ebmt
program define ebmt
	use "./data/ebmt", replace
	mat tmat = (.,1,2,.,3,4\.,.,.,5,6,7\.,.,.,8,9,10\.,.,.,.,11,12\.,.,.,.,.,.\.,.,.,.,.,.)
end

* EBMT covariate data
cap program drop ebmt_cov
program define ebmt_cov
	use "./data/ebmt", replace
	mat tmat = (.,1,2,.,3,4\.,.,.,5,6,7\.,.,.,8,9,10\.,.,.,.,11,12\.,.,.,.,.,.\.,.,.,.,.,.)
	qui gen cov = 10
	qui replace cov = -3 if id > 500
	qui replace cov = 0 if id > 1000
	qui replace cov = . if id > 1500
end


* #1 Check program runs with tmat
ebmt
msaj, transmatrix(tmat)

* #2 Check if CR works
ebmt
msaj, cr

* #3 Check by option works
ebmt_cov
msaj, transmatrix(tmat) by(cov)

* #4 Check ltruncated works
ebmt
msaj, transmatrix(tmat) ltruncated(0)
ebmt
msaj, transmatrix(tmat) ltruncated(365.25)

* #5 Check exit works
ebmt
msaj, transmatrix(tmat) exit(3652.5)
ebmt
msaj, transmatrix(tmat) exit(1)

* #6 Check from works
ebmt 
msaj, transmatrix(tmat) from(1)		
ebmt 
msaj, transmatrix(tmat) from(4)

* #7 Check CI option works
ebmt
msaj, transmatrix(tmat) ci

* #8 Check the SE option works
ebmt
msaj, transmatrix(tmat) se

* #9 Check los works
ebmt
msaj, transmatrix(tmat) los

* #10 Check all options work together, no cr
ebmt_cov
msaj, transmatrix(tmat) by(cov) ltruncated(365.25) exit(3652.5) from(4) ci se los

* #11 Check all options work together, cr (no from)
ebmt_cov
msaj, cr by(cov) ltruncated(365.25) exit(3652.5) ci se los

* #12 Check _to not needed if CR not specified
ebmt
drop _to
msaj, transmatrix(tmat)	

* #13 Check _to and _from not needed if CI/SE not specified
ebmt
drop _to _from
msaj, transmatrix(tmat)
ebmt
drop _from
msaj, cr		


*****************************
*  			HAI  			*
*****************************

* #9 Check all options work together, no cr
use "./data/hai", replace
mat tmat = (.,1,2,3,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
qui gen cov = 0
qui replace cov = 1 if id > 200
qui replace cov = -4 if id > 400
msaj, transmatrix(tmat) by(cov) ltruncated(3) exit(40) from(2) ci se los

* #10 Check all options work together, cr (no from)
use "./data/hai", replace
mat tmat = (.,1,2,3,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
qui gen cov = 0
qui replace cov = 1 if id > 200
qui replace cov = -4 if id > 400
msaj, cr by(cov) ltruncated(3) exit(40) ci se los


*****************************
*  			EBMT CR			*
*****************************

* #10 Check all options work together, cr (no from)
use "./data/ebmt_cr", replace
qui gen cov = mod(id,10)
msaj, cr by(cov) ltruncated(100) exit(1000) ci se los


*****************************
*  			Simbi 			*
*****************************

* #9 Check all options work together, no cr
use "./data/simbi", replace
mat tmat = (.,1,2\3,.,4\.,.,.)
qui gen cov = mod(id,6)*2 - 3
msaj, transmatrix(tmat) by(cov) ltruncated(2) exit(6) from(2) ci se los

* #10 Check all options work together, cr (no from)
use "./data/simbi", replace
mat tmat = (.,1,2\3,.,4\.,.,.)
qui gen cov = mod(id,6)*2 - 3
msaj, cr by(cov) ltruncated(2) exit(6) ci se los
