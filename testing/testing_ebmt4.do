adopath ++ "C:\multistate\Stata_Commands\stms\version_2_0_0_editing"
cd "C:\multistate\Stata_Commands\stms\version_2_0_0_editing"

clear all

buildmlib
use "C:\multistate\Data\ebmt4",clear

mat tmat = (.,1,2,.,3,4\.,.,.,5,6,7\.,.,.,8,9,10\.,.,.,.,11,12\.,.,.,.,.,.\.,.,.,.,.,.)
mat list tmat

msset, id(id) states(recs aes recaes rels srvs) times(rec ae recae rel srv) transmatrix(tmat)

mat list r(transmatrix)
mat list r(freqmatrix)
mat tmat = r(transmatrix)

stset _stop, enter(_start) f(_status=1) scale(365.25)


//separate models
forvalues i=1/12 {
	local trans `trans' _trans`i'
	streg if _trans`i'==1, dist(weib) diff
	est store m`i'
}

forvalues i=1/12 {
	local transopts `transopts' model`i'(m`i')
}

stms, transmat(tmat) `transopts' graph  obs(10) n(1000)



//one model
local trans
forvalues i=1/11 {
	local trans `trans' _trans`i'
}
streg `trans', dist(weib)

local transopts
forvalues i=1/12 {
	local transopts `transopts' trans`i'(_trans`i' 1)
}

stms, transmat(tmat) `transopts' graph  obs(10) n(1000)



