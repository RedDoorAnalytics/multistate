*! version 4.3.1 02sep2021 MJC

/*
History - multistate/predictms
02sep2021: version 4.3.1 - bug fix, graphms ignored extra options; now fixed
                         - bug fix, when rmst was called with simulate it would error out; now fixed
14mar2021: version 4.3.0 - accuracy and efficiency of numerical integration based predictions from an illness-death and extended 
                           illness death model are improved
                         - support for family(hazard) now documented
                         - bug fix for numerical integration based predictions of illness-death and extended illness-death, failing 
                           to pass the entry time appropriately would return a vector of 0s instead of appropriate probabilities; 
                           now fixed
                         - bug fix; if reset or tsreset() was used with an illness-death model, then numerical integration was used 
                           to calculate predictions instead of simulation, giving (obvious) incorrect results - this has been fixed
                         - bug fix when using a relative survival model within predictms, would attempt to call bhazard() when it 
                           should only simulate from the relative survival function; now fixed
                         - warning added when passing a relative survival model to -predictms-, to ensure the expected rate has been 
                           modelled.
09feb2021: version 4.2.1 - bug fix; check on transmatrix() skipped final row, so could erroneously exit telling the user the number 
                               of transitions and models didn't match - now fixed
26jan2012: version 4.2.0 - interactive & freqat() options added to msboxes
16dec2020: version 4.1.0 - family(pwexponential) synced & doc'd
28oct2020: version 4.0.0 - delta method is now the default for ci calculation for methods that use numerical integration or the hybrid AJ estimator for predictions
                                -> this is substantially faster than the bootstrap
                                -> includes contrasts
                                -> bootstrap is still the only method available for large-sample simulation based predictions
                         - bootstrap option added to override the default delta method
                         - latent option now documented (added in 3.0.0 to invoke latent times simulation, overiding the default total hazard + multinomial draws method)
                         - log normal, log logistic and general log cumulative hazard scale models added
                         - memory leak fixed, resulting in speed improvements
11sep2020: version 3.0.3 - msaj option enter() renamed to ltruncated() for consistency with survsim & predictms
08sep2020: version 3.0.2 - rmst missed from standardise calculation; now fixed
                         - ci + visit or rmst error'ed; now fixed
                         - storage improved, so bootstrap cis slightly faster
01sep2020: version 3.0.1 - left over print statements when simulating; now removed
24aug2020: version 3.0.0 - predictms rewrite:
                                 - Stata version 15.1 now required
                                 - all models must be fitted using stmerlin or merlin
                                 - transprob renamed to probability
                                 - survival renamed to singleevent for transmat() shorthand
                                 - ltruncated() introduced, which defines the starting time for all predictions
                                 - rmst added to predict restricted mean survival time, i.e. total time spent in non-absorbing states
                                 - hazard and survival options added to predict transition specific functions
                                 - standard survival, competing risks, illness-death and extended illness-death predictions are now 
                                   calculated using numerical integration as the default method, rather than large-sample simulation
                                 - the default method of simulation to obtain predictions uses Beyersmann et al., with a new -latent- option
                                   to use the previous default of latent times
                                 - each simulated dataset can now be saved
                                 - ci now displays a bootstrap output so it can be monitored
                                 - dist(exp) failed; now fixed
                                 - ci's now calculated on log/atanh/logit scale and back transformed to ensure appropriate bounds
                         - graphms added, with help file
                         - help files substantially improved, new main multistate help file added
*/

/*
History - galahad

MJC 21mar2020: version 1.3.0 - some behind the scenes syncing with survsim
MJC 17mar2020: version 1.2.1 - bug fix; singleevent option failed -> now fixed
                                 - help file example fixes for new transprob option
                                 - bug fix; missing error check, Cox model is not currently supported
MJC 07mar2020: version 1.2.0 - hazard option added to predictms to predict each transition-specific hazard function
                                 - survival changed to singleevent
                                 - survival now predicts transition-specific survival functions
                                 - transprob option now required to request transition probabilities
                                 - bug fix: visit with cis incorrectly returned los predictions instead -> now fixed.
                                 - bug fix: contrasts with visit, new variable names were incorrect -> now fixed
                                 - bug fix: userlink, transformation not applied to point estimate when calculating CI -> now fixed
                                 - bug fix: contrasts and ci, point estimates weren't passed -> now fixed
MJC 30jan2020: version 1.1.0 - bug fix displaying dots for multiple ats()
                                 - aj wasn't synced properly - now working
                                 - synced multiple timescales with all survival models
                                 - added timevar to ret list for use with ghgraph
                                 - ghgraph added, including help file
                                 - bug fix; family(logchazard) wasn't synced - now fixed
MJC 28jan2020: version 1.0.0
*/

/*
History - multistate/predictms
MJC 13aug2019: version 2.3.2 - bug fix: missed an edit in 2.3.1 fix preventing aj to work in predictms (merlin function names); now fixed
                             - error check added for streg, d(ggamma) combined with aj - not supported
MJC 25jul2019: version 2.3.1 - bug fix: edits in merlin function names caused predictms to fail; now fixed
                             - error check added for bhtime with strcs models - not currently supported
MJC 03jun2019: version 2.3.0 - bug fix: stpm2 and strcs models failed, introduced in 2.2.0; now fixed
                             - streg, dist(lognormal) now allowed with clock-forward models
                             - bug fix: streg, dist(lognormal) with reset incorrectly assumed loglogistic; now fixed
MJC 28may2019: version 2.2.0 - bug fix: with novcv() when ci not specified introduced in 2.1.0; now fixed
                             - bug in error check for log normal requiring reset; now fixed
                             - Generalised gamma streg model now supported when used in the models() syntax
                             - bug fix in msaj when number of transitions was not equal to number dimension of transmatrix()
MJC 15apr2019: version 2.1.0 - novcv() option added to use mean vector instead of draws in CI calculations
MJC 19mar2019: version 2.0.1 - help file improved for msaj
MJC 23feb2018: version 2.0.0 - bug fix: los was incorrectly scaled to 1, introduced in 1.2.0. Now fixed.
                                 - normal approximation for CIs is now the default. percentile option replaces normal.
                                 - bug fix: level() was ignored in CI calculations, always assumed 95%. Now fixed.
                                 - gen() removed
                                 - bug fix: in some cases mm_root_vec produced an error when solving the root. Now fixed.
                                 - all predictions now available in one call
                                 - infinite at#()s added (limited to 50 for error check reasons)
                                 - at() changed to at1()
                                 - difference option added
                                 - added atref(#) - default 1
                                 - standardise option added to calculate standardised (population-averaged) predictions
                                 - bug fix: with Stata 15 and streg (apart from dist(exp)), with at2(), the ancilary parameter constant was ignored. Now fixed.
                                 - _time = 0 now allowed with normal ci's
                                 - userfunction() added for user-defined prediction function, subroutines for probs and los
                                 - userlink() added for normal cis of userfunction() - default identity, can also be log or logit
                                 - cr option added to avoid specifying a transition matrix. For use with models() only.
                                 - bug fix: reset now synced with root-finding updates for stpm2 and strcs models
                                 - outsample added
                                 - stms added
                                 - enter() removed, now min(timevar or _t)
                                 - aj added for Aalen-Johansen estimator with Markov models
                                 - reversible transitions now allowed
MJC 17nov2017: version 1.2.1 - bug fix with simulation from strcs models with delayed entry
MJC 13nov2017: version 1.2.0 - tscale2() and time2() added for multiple timescales
                                 - probabilities/los scaled to 1/t
release moved to website
MJC 21nov2016: verison 1.1.1 - bug fix
MJC 16nov2016: version 1.1.0 - some re-writes of source code for massive speed improvements when root-finding required
                                 - trans#() syntax removed, now parsing is done automatically and using covariates expanded using msset
                                 - model#() changed to models()
                                 - improvements to help files
                                 - covariates() added to msset to create transition-specific dummies, which must be used in model fitting (if required)
                                 - default n() now 100,000, or 10,000 with ci
                             - survival added to calculate all predictions on a standard single event model
MJC 09aug2016: version 1.0.1 - gen() option to allow a stub for created variables, default is pred
MJC 12jul2016: version 1.0.0 

Development
stms not allowed yet
MJC 19may2015: version 1.1.0 - error check on transmatrix() added
							 - separate models allowed through model#(name, [at()])
MJC 14may2015: verison 1.0.8 - added error check not allowing AFT weibull or exp
MJC 11may2015: verison 1.0.7 - timevar() added
                             - basic graph option added which creates stacked plots
MJC 11may2015: verison 1.0.6 - fixed bug with weibull and forward approach
                             - Error checks improved                                         
MJC 09may2015: version 1.0.5 - fixed bug in from() which occurred when anything but from(1) was used
                             - fixed bug in forward calculations when enter>0
                             - fixed bug that only calculated predictions to states you could go to from first state
MJC 09may2015: version 1.0.4 - clock-forward approach now the default (simulations incorporate delayed entry), reset option added to use clock-reset
                             - only reset approach allowed with streg, dist(lnormal)
MJC 06may2015: version 1.0.3 - now synced with streg, dist(exp|weib|gompertz|llogistic|lnormal)
                             - normal approximation added for CI calculations
MJC 15apr2015: version 1.0.2 - when no ci's calculated it was using first draw from MVN, this has been fixed to be e(b)
MJC 01apr2015: version 1.0.1 - stpm2 simulation improved by creating and passing struct
                             - odds and normal scales added for stpm2 models
MJC 31mar2015: version 1.0.0
*/

program define predictms, sortpreserve properties(st) rclass
	version 15.1
	syntax 						, 		///
                                                                        ///
                [							///
                        TRANSMatrix(string)				///	-transition matrix-
                        SINGLEEVENT					///	-a single event survival model was fitted-
                        CR						/// -competing risks model was fitted-
                                                                        ///
                        PROBability					/// -transition probabilities-
                        LOS								///	-calculate length of stay in each state-
                        RMST							///	-restricted mean survival time-
                        VISIT							/// -prob. of ever visiting each state within time window-
                        HAZard							/// -transition-specific hazard functions- 
                        SURVival						/// -transition-specific survival functions-
                        USERFunction(string)			///	-name of Mata function-
                        USERLink(string)				///	-link function to apply to Mata function-
                                                                                        ///
                        MODels(string)					///	-est store objects-		
                                                                                        ///
                        FROM(numlist >0 asc int)		/// -starting state for predictions-
                        TIMEvar(varname numeric)		///	-prediction times-
                        LTruncated(numlist >=0 max=1)	///	-time at which observations become at risk-
                        OBS(string)						///	-Number of time points to calculate predictions at between mint() and maxt()-
                        MINT(string)					///	-minimum time to calculate predictions at-
                        MAXT(string)					///	-maximum time to calculate predictions at-
                                                                                        ///
                        SEED(string)					///	-pass to set seed-
                        N(string)						/// -sample size-
                        SIMulate						/// -use simulation instead of numerical integration-
                        LATENT							/// -use latent times simulation-
                        AJ								///	-use AJ estimator-
                                                                                        ///
                        M(numlist >=20 int max=1)		/// -number of parameter samples from MVN-
                        CI								/// -calculate confidence intevals-
                        BOOTstrap						///	-use parametric bootstrap for cis-
                        PERCentile						///	-calculate confidence intervals using percentiles-
                        NOVCV(string)					/// -transitions to skip VCV when calculating CIs-
                        Level(cilevel)					/// -level for CIs-
                                                                                        ///
                        DIFFerence						///	-calculate differences of predictions between at() and at2()-
                        RATIO							///	-calculate ratio of predictions between at() and at2()-
                        ATREFerence(string)				///	-reference at for contrasts-
                                                                                        ///
                        STANDardise						/// -standardise predictions-
                        STANDIF(string)					///	-restricts observations to standardise over-
                                                                                        ///
                        RESET							///	-use clock reset approach in simulations-
                        TSCALE2(string)					/// -transition models on a second timescale-
                        TIME2(string)					/// -time to add to main timescale-
                        TSRESET(string)					/// -transition-specific resets-
                                                                                        ///
                        OUTsample						///	-out of sample predictions-
                        SAVE(string)					///	-save each simulated dataset from point estimates and CI calculations-
                                                                                        ///
                                                                                        ///
                        SURVSIM(string)					/// -not documented-
                        SURVSIMTOUSE(string)			/// -not documented-
                                                                                        ///
                                                                                        ///
                        DEVCODE1(string)				/// -not documented-
                        DEVCODE2(string)				/// -not documented-
                        DEVCODE3(string)				/// -not documented-
                        DEVCODE9(string)				/// -not doc-
                                                                                        ///
                        INTERACTIVE JSONPATH(string) 	///
                                                                                        ///
                        *								/// -infinite ats-
                ]
					
	//================================================================================================================================================//
	// Preliminaries
	
	capture which merlin
	if _rc>0 {
		display in yellow "You need to install the command merlin. This can be installed using,"
		display in yellow ". {stata ssc install merlin}"
		exit  198
	}
	
	local ats `options'
	
	if "`probability'`los'`visit'`hazard'`survival'`userfunction'`rmst'"=="" {
		di as error "At least one prediction type must be specified"
		exit 198
	}
	
	if "`aj'"!="" {
		if "`probability'"=="" {
			local probability probability
		}
		if "`tscale2'"!="" {
			di as error "aj not allowed with tscale2()"
			exit 198
		}
		if "`reset'"!="" {
			di as error "aj not allowed with reset"
			exit 198
		}
		if "`n'"!="" {
			di as error "n() not needed with aj"
			exit 198
		}
		if "`visit'"!="" {
			di as error "visit not allowed with aj"
			exit 198
		}
	}
	
	local K = 1
	if "`standardise'"!="" {
		if "`hazard'"!="" {
			di as error "hazard and standardise not supported"
			exit 198
		}
		if "`outsample'"!="" {
			di as error "Not valid"
			exit 1986
		}
		if "`userfunction'"!="" {
			di as error "userfunction() not allowed with standardise"
			exit 198
		}
		tempvar stdtouse
		if "`standif'"!="" {
			qui gen byte `stdtouse' = (`standif')
		}
		else {
			qui gen byte `stdtouse' = 1
		}
		qui count if `stdtouse'
		local K = r(N)	
	}
	
	if "`outsample'"!="" {
		local out = 1
	}
	else local out = 0
	
	if `out' & "`timevar'"=="" {
		di as error "timevar() must be specified with outsample"
		exit 198
	}
	
	if "`visit'"!="" {
		local simulate simulate
	}
	
	if ("`difference'"!="" | "`ratio'"!="") & "`ats'"=="" {
		di as error "at least at2() needed"
		exit 198
	}
	
	if "`seed'"!="" {
		set seed `seed'
	}
	
	local survsimcall = "`survsim'"!=""
	
	if "`singleevent'"!="" {
		if "`transmatrix'"!="" {
			di in yellow  "transmatrix(`transmatrix') ignored"
		}
		tempname transmatrix
		mat `transmatrix' = (.,1\.,.)
		
		local latent latent					//faster with single event
	}
	
	if "`cr'"!="" {
		if "`models'"=="" {
			di as error "models() required with cr"
			exit 198
		}
		if "`transmatrix'"!="" {
			di in yellow  "transmatrix(`transmatrix') ignored"
		}
		local Ntrans : word count `models'
		tempname transmatrix
		mat `transmatrix' = J(`Ntrans'+1,`Ntrans'+1,.)
		forvalues i=1/`Ntrans' {
			mat `transmatrix'[1,`i'+1] = `i'
		}
	}
	
	cap confirm matrix `transmatrix'
	if _rc>0 & "`singleevent'"=="" & "`cr'"=="" {
		di as error "transmatrix(`transmatrix') not found"
		exit 198
	}
	mata: check_transmatrix()
	
	if ("`tscale2'"!="" & "`time2'"=="" ) | ("`tscale2'"=="" & "`time2'"!="" ) {
		di as error "tscale2() and time2() must both be specified"
		exit 198
	}
	
	if "`tscale2'"!="" & "`reset'"!="" {
		di as error "reset not allowed with tscale2()"
		exit 198
	}
	
	if "`ci'"=="" & "`m'"!="" {
		di as error "Cannot specify m() without ci"
		exit 198
	}

	if "`save'"!="" {
		if "`aj'"!="" {
			di as error "save() not allowed with aj"
			exit 198
		}
		local 0 `save'
		syntax anything(name=savestub) , [REPlace]
		local savereplace `replace'
	}
	
	if "`novcv'"!="" & "`models'"=="" {
		di as error "novcv() requires models()"
		exit 198
	}
	
	if "`models'"!="" {
		//pase models()
		local Nmodels : word count `models'
		if `Ntrans'!=`Nmodels' {
			di as error "Number of estimates objects in model() must be equal to number of transitions" 
			exit 198
		}
		forvalues i=1/`Ntrans' {
			local model`i' : word `i' of `models'
		}
			
		forvalues i=1/`Ntrans' {
			local 0 `model`i''
			syntax name(id="model estimates required") 
			local modelests`i' `namelist'
		}				
	}

	// Checks for interactive options
	if "`jsonpath'" != "" & "`interactive'" == "" {
		di as error "You have used the jsonpath option without using the interactive option."
		exit 198
	}
	if "`interactive'" != "" {
		if "`jsonpath'" != "" {
			mata st_local("direxists",strofreal(direxists("`jsonpath'")))
			if !`direxists' {
				di as error "Folder `jsonpath' does not exist."
				exit 198
			}
			mata st_local("jsonfile",pathjoin("`jsonpath'","msboxes_predictions.json"))   
		}
		else {
			local jsonfile msboxes_predictions.json
		}
		capture confirm file "`jsonfile'"
		if !_rc {
			capture erase "`jsonfile'"
			if _rc {
				display as error "`jsonfile' cannot be deleted'"
			}
		}
	}

	
	//===================================================================================================================//
	
	//core method choices
	
	if "`visit'"!="" {
		local simulate simulate
		local bootstrap bootstrap
	}
	
	//===================================================================================================================//
	//parse ats
	
	//check at#()
	
	local atind = 1
	while "`ats'"!="" {
		
		if `atind'>50 {
			di as error "at#() limit reached, or unrecognised option"
			exit 198
		}
		
		local 0 , `ats'
		syntax , [at`atind'(string) *]
		local ats `options'
		
		if !`out' {
			local varcount : word count `at`atind''
			local count = `varcount'/2			//!!add error check = integer
			tokenize ``at`atind'''
			while "`1'"!="" {
				unab 1: `1'
				cap confirm var `1'
				if _rc {
					di in red "invalid at`atind'(... `1' `2' ...)"
					exit 198
				}
				forvalues i=1/`Ntrans' {
					if "`1'"=="_trans`i'" {
						di as error "Cannot specify _trans# variables in at`atind'()"
						exit 198				
					}
				}
				cap confirm num `2'
				if _rc {
					di in red "invalid at`atind'(... `1' `2' ...)"
					exit 198
				}
				mac shift 2
			}  
		}
		local atind = `atind'+1

	}
	if `atind'>1 {
		local Nats = `atind'-1
	}
	else {
		local Nats = 1
	}

	if `Nats'==1 & ("`difference'"!="" | "`ratio'"!="") {
		di as error "difference/ratio require at least 2 at#() statements"
		exit 198
	}
	
	if "`atreference'"=="" {
		local atref = 1
	}
	else {
		cap confirm integer number `atreference'
		if _rc {
			di as error "atreference() must be an integer"
			exit 198
		}
		if `atreference'<1 {
			di as error "atreference() must be >=1"
			exit 198		
		}
		local atref = `atreference'
	}
	
	if `atref'>`Nats' & `Nats'>1 {
		di as error "atreference(#) must be an at#()"
		exit 198
	}
	
	//==//
			
	if "`from'"=="" local from 1
	
	// prediction time variable
	if "`timevar'"!="" & ("`mint'"!="" | "`maxt'"!="" | "`obs'"!="") {
		di as error "timevar() cannot be specified with mint()/maxt()/obs()"
		exit 198
	}
	
	if "`ltruncated'"=="" {
		local ltruncated = 0
	}
	
	if "`timevar'"=="" {
		if "`maxt'"=="" {
			qui su _t, meanonly
			local maxt = `r(max)'
		}
		if "`mint'"=="" {
			local mint = 0
		}
		
		if "`obs'"=="" & "`aj'"=="" local obs = 20
		else if "`obs'"=="" & "`aj'"!="" local obs = 500
		local timevar _time
		cap drop _time
		cap range2 _time `mint' `maxt' `obs'	
		label var _time "Follow-up time"
		
		//touse variable for predictions etc.
		tempvar touse
		qui gen byte `touse' = _n<= `obs'
	}
	else {
		//touse variable for predictions etc.
		tempvar touse
		qui gen byte `touse' = `timevar'!=.
		qui count if `touse'==1
		local obs = `r(N)'
		qui su `timevar', meanonly
		if `r(max)'<`ltruncated' {
			di as error "max(`timevar') must be > ltruncated()"
			exit 198
		}
		if `r(min)'<`ltruncated' {
			di as error "min(`timevar') must be >= ltruncated()"
			exit 198
		}
	}

	if "`aj'"!="" {
		qui su `timevar' if `touse'
		if `r(min)'!=`ltruncated' {
			di as error "mininum of timevar() must be equal to ltruncated() when using aj"
			exit 198
		}
	}
	
	//get core stuff

	//possible transitions from each state
	local Nstates = colsof(`transmatrix')
	if `Nstates'<2 {
		di as error "Must be at least 2 possible states, including starting state"
		exit 198
	}
	
	forvalues i=1/`Nstates' {
		forvalues j=1/`Nstates' {
			if (`transmatrix'[`i',`j']!=.) {
				local row`i'trans `row`i'trans' `=`transmatrix'[`i',`j']'
				local row`i'next `row`i'next' `j'
			}
		}
	}
	
	//check somewhere to go
	foreach frm in `from' {
		if "`row`frm'next'"=="" {
			di as error "No possible next states from(`frm')"
			exit 198
		}
	}
	
	
	//=====================================================================================================================================================//
	//CORE
		
		// models framework
		
		if "`models'"!="" {
		
			forvalues i=1/`Ntrans' {
				
				//error checks
				cap estimates restore `modelests`i''
					
				if _rc>0 {
					di as error "model estimates `modelests`i'' not found"
					exit 198			
				}
				predictms_modelcheck
				local familys `familys' `r(family)'
				
				//get estimates and variances
				tempname emat`i'
				mat `emat`i'' = e(b)
				
				if "`ci'"!="" {
					tempname  evmat`i'
					mat `evmat`i'' = e(V)
				}
				
				local Nparams`i' = colsof(`emat`i'')
				
			}
		
			//to sync with out of sample predictions - all vars must be present
			if `out' {
				forvalues i=1/`Ntrans' {
					cap estimates restore `modelests`i''
					foreach corevar in `e(allvars)' {
						cap gen `corevar' = 1 in 1
						if !_rc {
							local todrop `todrop' `corevar'
						}
					}
				}
			}			
		
			forval a=1/`Nats' {

				forvalues i=1/`Ntrans' {
					
					//handle outsample and extract standardise variables
					predictms_model, 	trans(`i') 					///
										nparams(`Nparams`i'') 		///
										ntrans(`Ntrans') 			///
										at(`at`a'')					///
										aind(`a')					///
										`standardise'
					
					if "`standardise'"!="" {
						local at`a'stdvars`i' `r(stdvars)'
					}
					
				}
			}
			
		}
		
		//stacked single model framework
		
		else {
			
			predictms_modelcheck
			local familys `r(family)'
		
			//get estimates and variances
			tempname emat
			mat `emat' = e(b)
			
			if "`ci'"!="" {
				tempname  evmat
				mat `evmat' = e(V)
			}
			
			local Nparams = colsof(`emat')
		
			//to sync with out of sample predictions - all vars must be present
			if `out' {
				foreach corevar in `e(allvars)' {
					cap gen `corevar' = 1 in 1
					if !_rc {
						local todrop `todrop' `corevar'
					}
				}
			}			
		
			forval a=1/`Nats' {

				forvalues i=1/`Ntrans' {
				
					//handle outsample and extract standardise variables
					predictms_model, 	trans(`i') 					///
										nparams(`Nparams`i'') 		///
										ntrans(`Ntrans') 			///
										at(`at`a'')					///
										aind(`a')					///
										`standardise'

					local toupdate `toupdate' `r(toupdate)'
					
					if "`standardise'"!="" {
						local at`a'stdvars`i' `r(stdvars)'
					}
					
				}
			}		
		
		}
		
		if "`standardise'"!="" {
			local stdcheck = 0
			forval a=1/`Nats' {
				forvalues i=1/`Ntrans' {
					di as text "Transition `i', at`a'(): Standardising over -> `at`a'stdvars`i''"
					local stdcheck = `stdcheck' + ("`at`a'stdvars`i''"=="")
				}
			}
			if `stdcheck' {
				di as error "No variables for standardising"
				exit 1986
			}
		}
	
		mata: predictms()
		
		if "`interactive'" != "" {
			mata: predictms_writejson2()
		}
	
	//=====================================================================================================================//
	//finish
	
	//tidy up
	cap drop `todrop'	
	
	if `survsimcall' {
		cap drop _time
		exit
	}
	
	// return list
	return matrix transmatrix = `transmatrix', copy
	return local Nstates = `Nstates'
	return local from `from'
	return local timevar `timevar'
	
	//Done

end

// double added to range
program define range2
        version 3.1
        if "`3'"=="" | "`5'"!="" { error 198 }
        confirm new var `1'
        if _N==0 { 
                if "`4'"=="" { error 198 } 
                set obs `4'
                local o "`4'"
        }
        else { 
                if "`4'"!="" { 
                        local o "`4'"
                        if `4' > _N { set obs `4' }
                }
                else    local o=_N
        }
        gen double `1'=(_n-1)/(`o'-1)*((`3')-(`2'))+(`2') in 1/`o'
end

program predictms_modelcheck, rclass

	if "`e(cmd)'"!="merlin" {
		di as error "Only merlin models are supported"
		exit 198	
	}	
	else {
	
		if `e(Nlevels)'>1 {
			di as error "Multilevel models are not supported"
			exit 198
		}
		if "`e(failure1)'"=="" {
			di as error "The merlin model must be a survival/time-to-event model"
			exit 198
		}
		if "`e(family2)'"!="" {
			di as error "Only univariate merlin models are currently supported"
			exit 198
		}
		
		if "`e(bhazard1)'"!="" {
			di as text "warning -> relative survival model detected; only valid when modelling the expected rate"
			di as text "        -> see " in smcl "{helpb exptorcs:help exptorcs}" as text " for more information"
		}
		
	}
	return local family "`e(family1)'"
	
end

mata
void check_transmatrix()
{
	tmat = st_matrix(st_local("transmatrix"))
	tmat_ind = tmat:!=.							//indicator matrix
	
	//Error checks
	if (max(diagonal(tmat_ind))>0) {
		errprintf("All elements on the diagonal of transmatrix() must be coded missing = .\n")
		exit(198)
	}
	row = 1
	rtmat = rows(tmat)
	trans = 1
	while (row<=rtmat) {
		for (i=1;i<=rtmat;i++) {
			if (sum(tmat:==tmat[row,i])>1 & tmat[row,i]!=.) {
				errprintf("Elements of transmatrix() are not unique\n")
				exit(198)
			}
			if (tmat[row,i]!=. & tmat[row,i]!=trans){
				errprintf("Elements of transmatrix() must be sequentially numbered from 1,...,K, where K = number of transitions\n")
				exit(198)
			}		
			if (tmat[row,i]!=.) trans++
		}
		row++
	}
	st_local("Ntrans",strofreal(trans-1))

}

end

program predictms_model, rclass
	syntax [, TRANS(string) NPARAMS(string) Ntrans(real 1) AT(string) STANDardise aind(string)]

	foreach corevar in `e(allvars)' {	
	
		if !strpos("`e(response1)'","`corevar'") {
			local inat = 0
			predictms_atparse, 	corevar(`corevar') 		///
								at(`at') 			 	///
								ntrans(`ntrans') 		///
								i(`trans')
			local inat = r(inat)
			local toupdate `toupdate' `r(toupdate)'
			if !`inat' & "`standardise'"!="" {
				predictms_stdparse `corevar' `trans' `ntrans'
				if r(include) {
					local stdvars `stdvars' `r(stdvar)'
					local toupdate `toupdate' `r(toupdate)'
				}
			}
		}
	}
	
	return local stdvars `stdvars'
	return local toupdate `toupdate'

end

program predictms_atparse, rclass
	syntax, [	corevar(string) 	///
				i(real 0) 			///
				Ntrans(real 1) 		///
				at(string)]
	
	local inat = 0
	tokenize `at'
	while "`1'"!="" {
		unab 1: `1'
		if "`corevar'"=="`1'_trans`i'" {
			local inat = 1
			local toupdate `toupdate' `1'
		}
		else if "`corevar'"=="`1'" {
			local inat = 1
		}
		mac shift 2
	} 
	
	//this makes sure you can't use any of the _trans# vars in at1(), at2() or standardising
	forvalues j=1/`ntrans' {
		if "`corevar'"=="_trans`j'" {
			local inat = 1
		}
	}
	return scalar inat = `inat'
	return local toupdate `toupdate'
end

program predictms_stdparse, rclass
	args corevar i Ntrans
	
	local include = 0
	//see if _trans`i' is in it
	local has_trans = strpos("`corevar'","_trans`i'")
	if `has_trans' {
	
		local strlength = strlen("_trans`i'")
		local strlength2 = strlen("`corevar'")
		local corevarstub = substr("`corevar'",1,`=`strlength2'-`strlength'')
																
		//need to check that each var matches with *_trans`i', others not included in trans specific design matrix
		local extracttrans = substr("`corevar'",`=`strlength2'-`strlength'+1',.)
		if "`extracttrans'"=="_trans`i'" {
			local include = 1
			local toupdate `corevarstub'
		}
										
	}
	else {
		//could still be ?_transj
		local hastransj = 0
		forvalues j=1/`Ntrans' {
			if `i'!=`j' {
				if strpos("`corevar'","_trans`j'") {
					local hastransj = 1
				}
			}
		}
		if !`hastransj' {
			local include = 1
			local corevarstub `corevar'
		}
	}

	return local stdvar `corevarstub'
	return local toupdate `toupdate'
	return scalar include = `include'
	
end

