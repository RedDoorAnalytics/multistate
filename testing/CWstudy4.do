//local drive Z:\
local drive /Users/Michael/Documents
clear all

cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/merlin/galahad/galahad"
clear all

do ./build/buildmlib.do

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


foreach mod in excessmod5 expmod deadmod5 dcsdeadmod_stpm24 {
	est use "./data/`mod'"
	est store `mod'
}


// Define transition matrix
mat tmat =  (.,1,2,3,.\.,.,.,.,4\.,.,.,.,5\.,.,.,.,.\.,.,.,.,.)

foreach i in 1985 {

	// Generate relevant year variables
	rcsgen, scalar(`i') knots(1985 1990 1996 2002 2008 2013) gen(y`i')
	forvalues l=1/5 {
		local y`i' `y`i'' yearspl`l' `=y`i'`l''
	}
	di in yellow "`y`i''"
	
	foreach j in 18 {

		// Generate relevant age variables
		rcsgen, scalar(`j') knots(18 25 35 50 66 80) gen(a`j')
		forvalues k=1/5 {
			local a`j' `a`j'' agespl`k' `=a`j'`k''
		}
		di in yellow "`a`j''"	
	
		local tempyear = `i' - `j'
	
		cap drop time
		range time 0 10 2
		
		// Predict PROB's for a male aged 30 at diagnosis, diagnosed in 1995
		set seed 20180604
		predictms, transmatrix(tmat) models(excessmod5 expmod deadmod5 dcsdeadmod_stpm24 dcsdeadmod_stpm24) ///
				out timevar(time) ci m(50) ///
				at1(female 0 `y`i'' `a`j'' tempoffset `tempyear' _age `j' _status 1) /// - specify covariate pattern to make prediction for
				tscale2(2) /// - specify which model that has another time scale
				time2(`j') /// - express how much to add on for the transition with different time scale
				tsreset(4 5) /// - specify models for which we want to re-set the clock
				devcode2(678990) /// - hidden option to allow merlin to be used with predictms
				n(10000) /// - increase number of iterations
				novcv(2)
		
	}
}
