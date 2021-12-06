
//local drive Z:/
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

do ./build/buildmlib.do
mata mata clear

//local drive Z:\
local drive /Users/Michael/Documents
//local drive c:
cd "`drive'/multistate/multistate"
adopath ++ "."
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./stms"
adopath ++ "./msboxes"
adopath ++ "./msaj"
adopath ++ "./graphms"
clear all

tr:do ./build/buildmlib.do
mata mata clear

clear all

** Data prep

use https://www.mjcrowther.co.uk/data/vonCube, replace
mat tmat = (.,1,2,3,.,.\.,.,.,.,4,5\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.\.,.,.,.,.,.)
msset, id(admid) states(status2 status3 status4 status5 status6) times(t12 t13 t14 t25 t26) transmat(tmat)


** ltruncated with ggamma - just flagging, fine as is for analysis

* ltruncated with time 0 works for e.g. gompertz and merlin
merlin (_stop if _trans == 1, family(lognormal, failure(_status) ltruncated(_start)))

* ltruncated with time 0 does not work for ggamma
// merlin (_stop if _trans == 1, family(ggamma, failure(_status) ltruncated(_start)))

* It works without specifying ltruncated though
// merlin (_stop if _trans == 1, family(ggamma, failure(_status)))

* It does work with ltruncated, when it is not 0 (transition #4)
// merlin (_stop if _trans == 4, family(ggamma, failure(_status) ltruncated(_start)))


** Hazard predictions for ggamma - I can just use streg for these

range timevar0 0 82 165

* It errors with predict and predictnl
merlin (_stop if _trans == 1, family(ggamma, failure(_status)))
predict fhaz, hazard timevar(timevar0)
// predictnl fhaz = predict(hazard)


** Predictms errors with ggamma - would need for analysis (unless using old code)

* Get the models, use rp (which works) for 4 transitions and do 1 with ggamma
foreach i in 1 3 4 5 {
	merlin (_stop if _trans == `i', family(rp, df(4) failure(_status) ltruncated(_start)))
	estimates store m`i'
}

merlin (_stop if _trans == 2, family(loglogistic, failure(_status)))
estimates store m2

* Call predictms
predictms, transm(tmat) timevar(timevar0) models(m1 m2 m3 m4 m5)  ///
	seed(3819407) probability 
