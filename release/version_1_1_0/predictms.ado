*! version 1.1.0 16nov2016 MJC

/*
History
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
	version 14.1
	syntax 						, 											///
																			///
										[									///
											TRANSMatrix(string)				///	-transition matrix-
											MODels(string)					///	-est store objects-
											GEN(string)						///	-new variable stub-
											RESET							///	-use clock reset approach in simulations-
											FROM(numlist >0 asc int)		/// -starting state for predictions-
											OBS(string)						///	-Number of time points to calculate predictions at between mint() and maxt()-
											MINT(string)					///	-minimum time to calculate predictions at-
											MAXT(string)					///	-maximum time to calculate predictions at-
											TIMEvar(varname numeric)		///	-prediction times-
											enter(string)					///	-time that patients enter model, default 0, for forward predictions-
											exit(string)					/// -time patients exit, for fixed horizon-
											SEED(string)					///	-pass to set seed-
											SURVival						///	-a single event survival model was fitted-
											CR(string)						///
																			///
											N(string)						/// -sample size-
											M(numlist >=20 int max=1)		/// -number of parameter samples from MVN-
											CI								/// -calculate confidence intevals for transprobs-
											NORMAL							///	-calculate confidence intervals using a normal approximation-
											Level(cilevel)					/// -level for normal approx. for CIs-
																			///
											LOS								///	-calculate length of stay in each state-
																			///
											AJ								/// -not documented-
											STD								/// -not documented-
											STDIF(string)					///	-not documented-
											STDIF2(string)					///	-not documented-
											CHECK							///	-not documented-
											OUT								///	-not documented-
											PROB							///	-not documented-
											SURVSIM(string)					/// -not documented-
											SURVSIMTOUSE(string)			/// -not documented-
																			///
											GRAPH							///
											GRAPHOPTS(string)				///
																			///
											AT(string)						/// 
																			///
																			/// // Comparison stats - default is difference
											AT2(string)						/// -second covariate pattern for contrast-
											RATIO							///	-calculate ratio of probabilities between at() and at2()-
																			///
										]
					
	//================================================================================================================================================//
	// Preliminaries
	
	local K = 1
	if "`std'"!="" {
		if "`out'"!="" | "`prob'"!="" {
			exit 1986
		}
		tempvar stdtouse
		if "`stdif'"!="" {
			qui gen byte `stdtouse' = _trans1==1 & `stdif'==1
			qui count if _trans1==1 & `stdif'==1
		}
		else {
			qui gen byte `stdtouse' = _trans1==1
			qui count if _trans1==1
		}
		local K = r(N)	
	}
	
	if "`out'"!="" {
		local out = 1
	}
	else local out = 0
	
	if `out' & "`timevar'"=="" {
		di as error "timevar() must be specified with out"
		exit 1986
	}
	
	if "`n'"=="" & "`ci'"=="" {
		local n = 100000
	}
	else if "`n'"=="" & "`ci'"!="" {
		local n = 10000
	}
	
	if "`prob'"!="" & "`at2'"=="" {
		di as error "at2() needed"
		exit 1986
	}
	
	if "`seed'"!="" {
		set seed `seed'
	}
	
	local survsimcall = "`survsim'"!=""
	
	if "`survival'"!="" {
		if "`transmatrix'"!="" {
			di in yellow  "transmatrix(`transmatrix') ignored"
		}
		tempname transmatrix
		mat `transmatrix' = (.,1\.,.)
	
	}
	
	if "`cr'"!="" {
		if "`transmatrix'"!="" {
			di in yellow  "transmatrix(`transmatrix') ignored"
		}
		tempname transmatrix
		mat `transmatrix' = J(`cr'+1,`cr'+1,.)
		forvalues i=1/`cr' {
			mat `transmatrix'[1,`i'+1] = `i'
		}
	}
	
	if "`graph'"!="" & ("`los'"!="" | "`at2'"!="") {
		di as error "graph not allowed with los or at2()"
		exit 198
	}
	
	cap confirm matrix `transmatrix'
	if _rc>0 & "`survival'"=="" {
		di as error "transmatrix(`transmatrix') not found"
		exit 198
	}
	mata: check_transmatrix()
	
	if "`ci'"=="" & "`m'"!="" {
		di as error "Cannot specify m() without ci"
		exit 198
	}
	
	if "`enter'"!="" & "`exit'"!="" {
		di as error "Cannot specify both enter() and exit()"
		exit 198
	}
	
	if "`e(cmd2)'"=="streg" & ("`e(cmd)'"=="weibull" | "`e(cmd)'"=="ereg") & "`e(frm2)'"== "time" {
		di as error "streg, dist(weib|exp) time, not supported"
		exit 198
	}
	
	cap which lmoremata.mlib
	if _rc>0 {
		di as error "You need to install the moremata library from SSC"
		exit 198
	}	
	
	//parse trans#()/model#() -> predictions calculated at these (essentially at() options), with zeros for everything else
	if "`e(cmd)'"!="stms" {
		local model = "`models'"!=""
		if `model' {
			local Nmodels : word count `models'
			if `Ntrans'!=`Nmodels' {
				di as error "Number of estimates objects in model() must be equal to number of transitions" 
				exit 198
			}
			forvalues i=1/`Ntrans' {
				local model`i' : word `i' of `models'
			}
		}
	}
	else {
		exit 1986
	}
	
	if "`options'"!="" {
		di as error "Unknown option(s): `options'"
		exit 198
	}
		
	//extended parsing of model#()
	if `model' {
		forvalues i=1/`Ntrans' {
			local 0 `model`i''
			syntax name(id="model estimates required") 
			local modelests`i' `namelist'
		}				
	}
	
	local hasat = "`at'"!=""
	local hasat2 = "`at2'"!=""
	local sim = "`aj'"==""
	
	//Error checks on at() and at2()
	if `hasat' & !`out' {
		local varcount : word count `at'
		local count = `varcount'/2			//!!add error check = integer
		tokenize `at'
		while "`1'"!="" {
			unab 1: `1'
			cap confirm var `1'
			if _rc {
				di in red "invalid at(... `1' `2' ...)"
				exit 198
			}
			forvalues i=1/`Ntrans' {
				if "`1'"=="_trans`i'" {
					di as error "Cannot specify _trans# variables in at()"
					exit 198				
				}
			}
			cap confirm num `2'
			if _rc {
				di in red "invalid at(... `1' `2' ...)"
				exit 198
			}
			mac shift 2
		}  
	}
	if `hasat2' & !`out' {
		local varcount : word count `at2'
		local count = `varcount'/2			//!!add error check = integer
		tokenize `at2'
		while "`1'"!="" {
			unab 1: `1'
			cap confirm var `1'
			if _rc {
				di in red "invalid at2(... `1' `2' ...)"
				exit 198
			}
			forvalues i=1/`Ntrans' {
				if "`1'"=="_trans`i'" {
					di as error "Cannot specify _trans# variables in at2()"
					exit 198				
				}
			}
			cap confirm num `2'
			if _rc {
				di in red "invalid at2(... `1' `2' ...)"
				exit 198
			}
			mac shift 2
		}  
	}
	
	if "`ci'"=="" & "`aj'"=="" local m = 1
	
	//default m for sims
	if "`ci'"!="" & "`m'"=="" & "`aj'"=="" local m = 200
			
	if "`from'"=="" local from 1
	
	cap set obs `obs'	
	if "`aj'"=="" cap set obs `m'
	
	// prediction time variable
	if "`timevar'"!="" & ("`mint'"!="" | "`maxt'"!="" | "`obs'"!="") {
		di as error "timevar() cannot be specified with mint()/maxt()/obs()"
		exit 198
	}
	
	if "`timevar'"=="" {
		if "`exit'"=="" {
			if "`maxt'"=="" {
				qui su _t, meanonly
				local maxt = `r(max)'
			}
			if "`mint'"=="" {
				if "`enter'"=="" {
					qui su _t, meanonly
					local mint = `r(min)'
					local enter = 0
				}
				else {
					local mint = `enter'
				}
			}
			else {
				if "`enter'"=="" local enter = 0
			}
			if `maxt'<`enter' {
				di as error "maxt() must be > enter()"
				exit 198
			}
			if `mint'<`enter' {
				di as error "mint() must be >= enter()"
				exit 198
			}
			if `mint'==`enter' & "`ci'"!="" & "`normal'"!="" & "`aj'"=="" {
				di as error "mint() > enter() when confidence intervals with normal approximation are required"
				exit 198
			}			
		}
		else {
			if "`aj'"!="" local enter = 0
			if "`maxt'"=="" {
				local maxt = `exit'
			}
			if "`mint'"=="" {
				qui su _t, meanonly
				local mint = `r(min)'
			}
			if `maxt'>`exit' {
				di as error "maxt() must be <= exit()"
				exit 198
			}
			if `mint'>=`exit' {
				di as error "mint() must be < exit()"
				exit 198
			}
			if `maxt'==`exit' & "`ci'"!="" & "`normal'"!="" & "`aj'"=="" {
				di as error "maxt() < exit() when confidence intervals with normal approximation are required"
				exit 198
			}			
		}
		
		if "`obs'"=="" & "`aj'"=="" local obs = 20
		else if "`obs'"=="" & "`aj'"!="" local obs = 100
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
		if "`enter'"=="" local enter = 0
		qui su `timevar', meanonly
		if "`exit'"=="" {
			if `r(max)'<`enter' {
				di as error "max(`timevar') must be > enter()"
				exit 198
			}
			if `r(min)'<`enter' {
				di as error "min(`timevar') must be >= enter()"
				exit 198
			}
			if `r(min)'==`enter' & "`ci'"!="" & "`normal'"!="" {
				di as error "min(`timevar') > enter() when confidence intervals with normal approximation are required"
				exit 198
			}
		}
		else {
			if "`aj'"!="" local enter = 0
			if `r(max)'>`exit' {
				di as error "max(`timevar') must be <= exit()"
				exit 198
			}
			if `r(min)'>=`exit' {
				di as error "min(`timevar') must be < exit()"
				exit 198
			}
			if `r(max)'==`exit' & "`ci'"!="" & "`normal'"!="" {
				di as error "max(`timevar') < exit() when confidence intervals with normal approximation are required"
				exit 198
			}			
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
	
	//create variables to hold predictions
	if !`survsimcall' {
		if "`gen'"=="" {
			cap drop pred_*
			local stub pred
		}
		else {
			cap drop `gen'_*
			local stub `gen'
		}
		foreach fromstate in `from' {
			forvalues i=1/`Nstates' {
				qui gen double `stub'_`fromstate'_`i' = 0 if _n<= `obs'
				local fromvars`fromstate' `fromvars`fromstate'' `stub'_`fromstate'_`i'
				if "`los'"!="" {
					if "`at2'"=="" {
						label var `stub'_`fromstate'_`i' "Length of stay in state=`i', from state=`fromstate'"
					}
					else {
						if "`ratio'"!="" {
							label var `stub'_`fromstate'_`i' "LoS(at()) / LoS(at2()), in state=`i', from state=`fromstate'"
						}
						else {
							label var `stub'_`fromstate'_`i' "LoS(at()) - LoS(at2()), in state=`i', from state=`fromstate'"
						}
					}
				}
				
				if "`los'"=="" {
					if "`at2'"=="" {
						label var `stub'_`fromstate'_`i' "Transition prob. state=`i', from state=`fromstate'"
					}
					else {
						if "`ratio'"=="" {
							label var `stub'_`fromstate'_`i' "Transition prob(at()) - Transition prob(at2()) state=`i', from state=`fromstate'"
						}
						else {
							label var `stub'_`fromstate'_`i' "Transition prob(at()) / Transition prob(at2()) state=`i', from state=`fromstate'"
						}
					}
				}			
				
				if "`ci'"!="" & "`aj'"=="" {
					if "`normal'"=="" {
						qui gen double `stub'_`fromstate'_`i'_lci = 0 if _n<= `obs'
						qui gen double `stub'_`fromstate'_`i'_uci = 0 if _n<= `obs'
						local fromvarslci`fromstate' `fromvarslci`fromstate'' `stub'_`fromstate'_`i'_lci
						local fromvarsuci`fromstate' `fromvarsuci`fromstate'' `stub'_`fromstate'_`i'_uci
					}
					else {
						tempvar pred_`fromstate'_`i'_se
						qui gen double `pred_`fromstate'_`i'_se' = 0 if _n<= `obs'
						local probvars_`fromstate'_`i' `stub'_`fromstate'_`i' `pred_`fromstate'_`i'_se'
					}
				}
				else if "`ci'"!="" & "`aj'"!="" {
					qui gen double se_`stub'_`fromstate'_`i' = 0 if _n<= `obs'
					local se_fromvars`fromstate' `se_fromvars`fromstate'' se_`stub'_`fromstate'_`i'
				}
			}	
		}
	}
	
	if "`ci'"!="" & "`aj'"=="" {
		tempvar mvnind
		gen byte `mvnind' = _n<=`m'
	}
	
	
	//=====================================================================================================================================================//
	//CORE
	
	//======================================================//
	// stacked model 
	if !`model' {
		
		predictms_modelcheck			
				
		if `sim' {
		
			//coefficients
			tempname emat evmat
			mat `emat' = e(b)
			mat `evmat' = e(V)
			
			if "`e(cmd)'"=="stpm2" {
				mat `emat' = `emat'[1,"xb:"]
				mat `evmat' =`evmat'["xb:","xb:"]
			}
			
			forvalues i=1/`Ntrans' {
				tempname emat`i'
				mat `emat`i'' = `emat'
				local cmds `cmds' `e(cmd)'				//to sync
			}
			local Nparams = colsof(`emat')
			if "`ci'"!="" {
				forvalues i=1/`Nparams' {
					tempvar draws`i'
					local drawvars `drawvars' `draws`i''
				}
				qui drawnorm `drawvars', means(`emat') cov(`evmat')	//!! this is simulating more than I need when _N or `obs' > m
				forvalues i=1/`Ntrans' {
					local drawvars`i' `drawvars'
				}
			}	
		
		
			if "`e(cmd2)'"=="streg" {
				
				if "`e(cmd)'"=="ereg" {
					local cmdline `e(cmdline)'
					gettoken cmd 0 : cmdline
					syntax varlist, [NOCONstant *]
					if "`noconstant'"=="" local addcons _cons
					local corevars `varlist' `addcons'
				}
				else {
					local cmdline `e(cmdline)'
					gettoken cmd 0 : cmdline
					syntax varlist, [NOCONstant ANCillary(varlist) *]
					if "`noconstant'"=="" local addcons _cons
					local corevars `varlist' `addcons' `ancillary' _cons
					//indices for lambda and gamma
					tempname indices			
					mat `indices' = J(2,2,1)
					local colindex = `: word count `varlist'' + ("`noconstant'"=="")
					mat `indices'[2,1] = `colindex'
					mat `indices'[1,2] = `colindex' + 1
					mat `indices'[2,2] = `: word count `corevars''
					forvalues i=1/`Ntrans' {
						tempname indices`i'										//extension to index dm and b
						mat `indices`i'' = `indices'
					}
				}
				
				//now loop over vars and update DM and indices
				//can match variables in at() with varlist
				forvalues i=1/`Ntrans' {
										
					//design matrix for each transition
					tempname dm`i'
					mat `dm`i'' = J(1,`Nparams',0)
					if "`at2'"!="" {
						tempname at2dm`i'
						mat `at2dm`i'' = `dm`i''
					}
				
					local colindex = 1
					foreach corevar in `corevars' {
					
						if "`corevar'"=="_cons" {
							mat `dm`i''[1,`colindex'] = 1
							if "`at2'"!="" mat `at2dm`i''[1,`colindex'] = 1
						}
						else {
							local inat = 0
							predictms_atparse, 	corevar(`corevar') colindex(`colindex') dmmat(`dm`i'') 		///
										i(`i') ntrans(`Ntrans') at(`at') out(`out')
							local inat = r(inat)
							local todrop `todrop' `r(todrop)'
							local inat2 = 0
							if `hasat2' {
								predictms_atparse, 	corevar(`corevar') colindex(`colindex') dmmat(`at2dm`i'') 	///
											i(`i') ntrans(`Ntrans') at(`at2') out(`out')
								local inat2 = r(inat)
								local todrop `todrop' `r(todrop)'
							}
							
							if !`inat' & "`std'"!="" {
								predictms_stdparse `corevar' `i' `Ntrans'
								if r(include) {
									local stdvars`i' `stdvars`i'' `r(stdvar)'
									local stdvarsindex`i' `stdvarsindex`i'' `colindex'
								}
							}
							if `hasat2' & !`inat2' & "`std'"!="" {
								predictms_stdparse `corevar' `i' `Ntrans'
								if r(include) {
									local at2stdvars`i' `sat2tdvars`i'' `r(stdvar)'
									local at2stdvarsindex`i' `at2stdvarsindex`i'' `colindex'
								}
							}
						}
						local colindex = `colindex' + 1
						
					}
				
				}
			
			}
			else if "`e(cmd)'"=="stpm2" | "`e(cmd)'"=="strcs" {
			
				//DM is only for varlist, tvc splines and base splines are handled separately
				local corevars `e(varlist)'
				local Ncovs : word count `corevars'
				local nocons `e(noconstant)'
				local orthog `e(orthog)'
				
				//design matrix for each transition
				if `Ncovs' > 0 {

					//now loop over trans# and update DM and indices
					//can match variables in trans#() with varlist and ancillary
					forvalues i=1/`Ntrans' {
						tempname dm`i'
						mat `dm`i'' = J(1,`Ncovs',0)
						if "`at2'"!="" {
							tempname at2dm`i'
							mat `at2dm`i'' = `dm`i''
						}
						local colindex = 1
						foreach corevar in `corevars' {
						
							local inat = 0
							predictms_atparse, corevar(`corevar') colindex(`colindex') dmmat(`dm`i'') i(`i') ntrans(`Ntrans') at(`at') out(`out')
							local inat = r(inat)
							local todrop `todrop' `r(todrop)'
							local inat2 = 0
							if `hasat2' {
								predictms_atparse, corevar(`corevar') colindex(`colindex') dmmat(`at2dm`i'') i(`i') ntrans(`Ntrans') at(`at2') out(`out')
								local inat2 = r(inat)
								local todrop `todrop' `r(todrop)'
							}
							
							//standardising	-> kept out here as could standardise without at() and at2()
							if !`inat' & "`std'"!="" {
								predictms_stdparse `corevar' `i' `Ntrans'
								if r(include) {
									local stdvars`i' `stdvars`i'' `r(stdvar)'
									local stdvarsindex`i' `stdvarsindex`i'' `colindex'
								}
							}
							if `hasat2' & !`inat2' & "`std'"!="" {
								predictms_stdparse `corevar' `i' `Ntrans'
								if r(include) {
									local at2stdvars`i' `at2stdvars`i'' `r(stdvar)'
									local at2stdvarsindex`i' `at2stdvarsindex`i'' `colindex'
								}
							}
							local colindex = `colindex' + 1
						}
					}
				}
				
				local rcsbaseoff `e(rcsbaseoff)'										//always empty for strcs
				if "`rcsbaseoff'"=="" {
					local Nsplines : word count `e(rcsterms_base)'
					if "`e(cmd)'"=="stpm2" {
						local ln_bknots `e(ln_bhknots)'										//all log baseline knots including boundary knots
						if "`ln_bknots'"=="" {	//this is empty when df(1)
							local ln_bknots `=log(`: word 1 of `e(boundary_knots)'')' `=log(`: word 2 of `e(boundary_knots)'')'
						}
					}
					else {
						local ln_bknots `e(bhknots)'										//all log baseline knots including boundary knots
						if "`ln_bknots'"=="" {	//this is empty when df(1)
							local ln_bknots -5 10 //fudge - these values are not used
						}
					}
					if "`orthog'"!="" {
						tempname rmat
						mat `rmat' = e(R_bh)
						local rmatopt rmatrix(`rmat')
					}			
				}
				
				//to sync with separate model sim
				forvalues i=1/`Ntrans' {
					local rcsbaseoff`i' `rcsbaseoff'
					local scale`i' `e(scale)'
					local orthog`i' `e(orthog)'
					local nocons`i' `e(noconstant)'
					if "`rcsbaseoff`i''"=="" {
						local ln_bknots`i' `ln_bknots'
						if "`orthog`i''"!="" {
							tempname rmat`i'
							mat `rmat`i'' = `rmat'
						}
					}
				}
						
				local tvc `e(tvc)'
				forvalues i=1/`Ntrans' {
					local tvc`i' `tvc'
				}
				if "`tvc'"!="" {
					local i = 1
					foreach tvcvar in `tvc' {
						if "`e(cmd)'"=="stpm2" {
							local boundary_knots_`i' `e(boundary_knots_`tvcvar')'
							local ln_tvcknots_`i' `e(ln_tvcknots_`tvcvar')'
							if "`ln_tvcknots_`i''"=="" {
								local ln_tvcknots_`i' `=log(`: word 1 of `boundary_knots_`i''')' `=log(`: word 2 of `boundary_knots_`i''')'
							}		
						}
						else {
							local boundary_knots_`i' `e(boundary_knots_`tvcvar')'
							local ln_tvcknots_`i' `e(tvcknots_`tvcvar')'
							if "`ln_tvcknots_`i''"=="" {
								local ln_tvcknots_`i' -5 10	//fudge - these values are not used
							}	
						
						}
						forvalues j=1/`Ntrans' {
							local ln_tvcknots`j'_`i' `ln_tvcknots_`i''
						}
						if "`orthog'"!="" {
							tempname R_`i'
							mat `R_`i'' = e(R_`tvcvar')
							forvalues j=1/`Ntrans' {
								tempname R`j'_`i'
								mat `R`j'_`i'' = `R_`i''
							}
						}
						local i = `i' + 1
					}
					local Ntvcvars = `i' - 1
					forvalues i=1/`Ntrans' {
						local Ntvcvars`i' `Ntvcvars'
					}
					//tvc DM
					forvalues i=1/`Ntrans' {
						tempname dmtvc`i'
						mat `dmtvc`i'' = J(1,`Ntvcvars`i'',0)
						if `hasat2' {
							tempname at2dmtvc`i'
							mat `at2dmtvc`i'' = `dmtvc`i''
						}
						local colindex = 1
						foreach corevar in `tvc`i'' {
						
							local inat = 0
							predictms_atparse, corevar(`corevar') colindex(`colindex') dmmat(`dmtvc`i'') i(`i') ntrans(`Ntrans') at(`at') out(`out')
							local inat = r(inat)
							local todrop `todrop' `r(todrop)'
							local inat2 = 0
							if `hasat2' {
								predictms_atparse, corevar(`corevar') colindex(`colindex') dmmat(`at2dmtvc`i'') i(`i') ntrans(`Ntrans') at(`at2') out(`out')
								local inat2 = r(inat)
								local todrop `todrop' `r(todrop)'
							}
							
							//standardising
							if "`std'"!="" {
								if !`inat' {
									predictms_stdparse `corevar' `i' `Ntrans'
									if r(include) {
										local tvcstdvars`i' `tvcstdvars`i'' `r(stdvar)'
										local tvcstdvarsindex`i' `tvcstdvarsindex`i'' `colindex'
									}
								}
								if `hasat2' & !`inat2' {
									predictms_stdparse `corevar' `i' `Ntrans'
									if r(include) {
										local tvcat2stdvars`i' `tvcat2stdvars`i'' `r(stdvar)'
										local tvcat2stdvarsindex`i' `tvcat2stdvarsindex`i'' `colindex'
									}
								}							
							}
							local colindex = `colindex' + 1
						}
					}
					
				}
			
			}
			else if ""=="stms" {
			
				if "`ci'"!="" local count = 1
			
				forvalues i=1/`Ntrans' {
			
					local stmodel : word `i' of `e(models)'		
					
					predictms_stms_model `i' `stmodel' `emat' `at'
				
					tempname dm`i' emat`i' 
					mat `dm`i'' = r(dm)
					
					local Nparams = r(Nparams)
					
					if "`ci'"=="" {
						mat `emat`i'' = r(b)
					}
					else {
						//get drawvar names for each transition model 
						forvalues j=1/`Nparams' {
							local drawvars`i' `drawvars`i'' `draws`count''
							local count = `count' + 1
						}
						//need to skip dxb and s0xb
						if "`stmodel'"=="stpm2" {
							local count = `count' + `Nparams' - `: word count `e(varlist`i')'' - 1
							if `: word `i' of `e(delentry)'' {
								local count = `count' + `Nparams'
							}
						}
					
					}
				
					if "`stmodel'"!="stpm2" & "`stmodel'"!="ereg"{
						tempname indices`i'
						mat `indices`i'' = r(indices)
					}
					if "`stmodel'"=="stpm2" {
						if "`e(orthog`i')'"=="orthog" {
							tempname rmat`i'
							mat `rmat`i'' = r(rmat)
							if "`e(tvc`i')'"!="" {
								forvalues j=1/`Ntvcvars`i'' {
									tempname R`i'_`j'
									mat `R`i'_`j'' = r(R_`j')
								}
							}
						}
						if "`e(tvc`i')'"!="" {
							tempname dmtvc`i'
							mat `dmtvc`i'' = r(dmtvc)
						}
					}
					
					local cmds `cmds' `cmdname'
					
					//second DM for at2()
					if `hasat2' {
						predictms_stms_model_at2 `i' `stmodel' `Nparams' `at2'
						tempname at2dm`i'
						mat `at2dm`i'' = r(dm)
						if "`stmodel'"=="stpm2" & "`e(tvc`i')'"!="" {
							tempname at2dmtvc`i'
							mat `at2dmtvc`i'' = r(dmtvc)
						}
					}
					
				}
			
			}
			
		}
		
		// stacked + AJ
		else {
		
			if "`e(cmd2)'"=="streg" {
				cap set obs `=`obs'+1'
				tempvar ajtimevar ajtouse ajposttouse
				qui gen `ajtimevar' = `timevar'
				qui replace `ajtimevar' = `enter' if _n==`=`obs'+1'
				qui gen byte `ajtouse' = _n<=`=`obs'+1'
				qui gen byte `ajposttouse' = _n<=`obs'
				
				forvalues i=1/`Ntrans' {
					tempvar ch`i'
					if "`ci'"!="" {
						tempvar varch`i' 
						local varopt variance(`varch`i'')				
						local varchvars `varchvars' `varch`i''
					}
					local trans`i' _trans`i' 1
					qui streg_predict_ch `ch`i'' , touse2(`ajtouse') at(`trans`i'') timevar(`ajtimevar') zeros `varopt'
					qui replace `ch`i'' = 0 if `ajtouse' & `ch`i''==.
					if "`ci'"!="" qui replace `varch`i'' = 0 if `ajtouse' & `varch`i''==.
					local chvars `chvars' `ch`i''
				}	
			}
			else if "`e(cmd)'"=="stpm2" {
			
				cap set obs `=`obs'+1'
				tempvar ajtimevar ajtouse ajposttouse
				qui gen `ajtimevar' = `timevar'
				qui replace `ajtimevar' = `enter' if _n==`=`obs'+1'
				qui gen byte `ajtouse' = _n<=`=`obs'+1'
				qui gen byte `ajposttouse' = _n<=`obs'
				
				//AJ
				forvalues i=1/`Ntrans' {
					tempvar ch`i'
					if "`ci'"!="" {
						tempvar varch`i' 
						qui predictnl double `ch`i'' = predict(cumhazard at(`trans`i'') timevar(`ajtimevar') zeros) if `ajtouse', variance(`varch`i'')
						local varchvars `varchvars' `varch`i''
					}
					else {
						qui predict `ch`i'' if `ajtouse' , cumhazard at(`trans`i'') timevar(`ajtimevar') zeros
					}
					local chvars `chvars' `ch`i''
				}	
			}
			else if "`e(cmd)'"=="cox" {
			
				
				cap set obs `=`obs'+1'
				tempvar ajtimevar ajtouse ajposttouse
				qui gen `ajtimevar' = `timevar'
				qui replace `ajtimevar' = `enter' if _n==`=`obs'+1'
				qui gen byte `ajtouse' = _n<=`=`obs'+1'
				qui gen byte `ajposttouse' = _n<=`obs'

				//Get predicted baseline cumulative hazard function
				tempvar basech
				predict `basech', basech
				
				//sample and events
				tempvar eventind
				gen byte `eventind' = e(sample) & _d==1
				
				//indicator for each trans and events to read in basech
				forvalues i= 1/`Ntrans' {
					 tempvar chind`i'
					 gen byte `chind`i'' = (e(sample) & `e(strata)'==`i' & _d==1)
				}
				
				tempname emat
				mat `emat' = e(b)
				
				//linear predictors
				//overall design matrix for each transition, stacked
				if `Nparams'>0 {
					tempname dm indices
					mat `dm' = J(`Ntrans',`Nparams',0)
					if "`at2'"!="" {
						tempname at2dm
						mat `at2dm' = `dm'
					}
				}
				local cmdline `e(cmdline)'
				gettoken cmd 0 : cmdline
				syntax [varlist(default=empty)] [if] [in], [*]
				
				//now loop over vars and update DM and indices
				//can match variables in at() with varlist
				forvalues i=1/`Ntrans' {
					local colindex = 1
					foreach corevar in `varlist' {
						local inat = 0
						if "`at'"!="" {
							tokenize `at'
							while "`1'"!="" {
								unab 1: `1'
								if "`corevar'"=="`1'_trans`i'" | "`corevar'"=="`1'" {
									mat `dm'[`i',`colindex'] = `2'
									local inat = 1
								}
								mac shift 2
							} 
						}
						local inat2 = 0
						if "`at2'"!="" {
							tokenize `at2'
							while "`1'"!="" {
								unab 1: `1'
								if "`corevar'"=="`1'_trans`i'" | "`corevar'"=="`1'" {
									mat `at2dm'[`i',`colindex'] = `2'
									local inat2 = 1
								}
								mac shift 2
							} 
						}
						
						if "`corevar'"=="_trans`i'" {
							mat `dm'[`i',`colindex'] = 1
							if "`at2'"!="" mat `at2dm'[`i',`colindex'] = 1
						}
						forvalues j=1/`Ntrans' {
							if "`corevar'"=="_trans`j'" {
								local inat = 1
								if "`at2'"!="" local inat2 = 1
							}
						}
						
						if `hasat' & !`inat' & "`std'"!="" {
							predictms_stdparse `corevar' `i' `Ntrans'
							if r(include) {
								local stdvars`i' `stdvars`i'' `r(stdvar)'
								local stdvarsindex`i' `stdvarsindex`i'' `colindex'
							}
						}
						if `hasat2' & !`inat2' & "`std'"!="" {
							predictms_stdparse `corevar' `i' `Ntrans'
							if r(include) {
								local at2stdvars`i' `sat2tdvars`i'' `r(stdvar)'
								local at2stdvarsindex`i' `at2stdvarsindex`i'' `colindex'
							}
						}
						
						local colindex = `colindex' + 1
					}

				}
			}
		}

	}	
	
	//======================================================//
	// models framework
	else {
		
		if `sim' {
		
			forvalues i=1/`Ntrans' {
			
				capture estimates restore `modelests`i''								//could also have use
				if _rc {
					di as error "model estimates `modelests`i'' not found"
					exit 198			
				}
				predictms_modelcheck			
				
				tempname emat`i' evmat`i'
				mat `emat`i'' = e(b)
				mat `evmat`i'' = e(V)
				
				if "`e(cmd)'"=="stpm2" {
					mat `emat`i'' = `emat`i''[1,"xb:"]
					mat `evmat`i'' =`evmat`i''["xb:","xb:"]
				}
				
				local Nparams`i' = colsof(`emat`i'')
				if "`ci'"!="" {
					forvalues j=1/`Nparams`i'' {
						tempvar draws`i'`j'
						local drawvars`i' `drawvars`i'' `draws`i'`j''
					}
					qui drawnorm `drawvars`i'', means(`emat`i'') cov(`evmat`i'')	//!! this is simulating more than I need when _N or `obs' > m -> move to Mata
				}	
				
				predictms_model, trans(`i') nparams(`Nparams`i'') ntrans(`Ntrans') at(`at') `std' out(`out')
				local todrop `todrop' `r(todrop)'
				
				if "`std'"!="" {
					local stdvars`i' `r(stdvars)'
					local stdvarsindex`i' `r(stdvarsindex)'
				}
				
				tempname dm`i'
				mat `dm`i'' = r(dm)

				if "`e(cmd)'"!="stpm2" & "`e(cmd)'"!="ereg" & "`e(cmd)'"!="strcs" {
					tempname indices`i'
					mat `indices`i'' = r(indices)
				}
				if "`e(cmd)'"=="stpm2" | "`e(cmd)'"=="strcs" {
					if "`e(orthog)'"=="orthog" {
						tempname rmat`i'
						mat `rmat`i'' = r(rmat)
						if "`e(tvc)'"!="" {
							forvalues j=1/`Ntvcvars`i'' {
								tempname R`i'_`j'
								mat `R`i'_`j'' = r(R_`j')
							}
						}
					}
					if "`e(tvc)'"!="" {
						tempname dmtvc`i'
						mat `dmtvc`i'' = r(dmtvc)
						if "`std'"!="" {
							local tvcstdvars`i' `r(tvcstdvars)'
							local tvcstdvarsindex`i' `r(tvcstdvarsindex)'
						}
					}
				}
				local cmds `cmds' `cmdname'
				
				//second DM for at2()
				if "`at2'"!="" {
				
					predictms_model_at2, model(`e(cmd)') trans(`i') ntrans(`Ntrans') nparams(`Nparams`i'') at(`at2') `std' out(`out')

					if "`std'"!="" {
						local at2stdvars`i' `r(stdvars)'
						local at2stdvarsindex`i' `r(stdvarsindex)'
					}
					
					tempname at2dm`i'
					mat `at2dm`i'' = r(dm)
					if ("`e(cmd)'"=="stpm2" | "`e(cmd)'"=="strcs") & "`e(tvc)'"!="" {
						tempname at2dmtvc`i'
						mat `at2dmtvc`i'' = r(dmtvc)
						if "`std'"!="" {
							local tvcat2stdvars`i' `r(tvcstdvars)'
							local tvcat2stdvarsindex`i' `r(tvcstdvarsindex)'
						}
					}
					
				}
				
			}
		
		}	
		else {
			
			cap set obs `=`obs'+1'
			tempvar ajtimevar ajtouse ajposttouse
			qui gen `ajtimevar' = `timevar'
			qui replace `ajtimevar' = `enter' if _n==`=`obs'+1'
			qui gen byte `ajtouse' = _n<=`=`obs'+1'
			qui gen byte `ajposttouse' = _n<=`obs'
			
			forvalues i=1/`Ntrans' {
				
				capture estimates restore `modelests`i''								//could also have use
				if _rc {
					di as error "model`i'() estimates not found"
					exit 198			
				}
				tempvar ch`i'
				if "`e(cmd2)'"=="streg" {
					qui streg_predict_ch `ch`i'' , touse2(`ajtouse') at(`at') timevar(`ajtimevar') zeros
				}
				else if "`e(cmd)'"=="stpm2" {
					qui predict `ch`i'' if `ajtouse' , cumhazard at(`at') timevar(`ajtimevar') zeros
				}
				else if "`e(cmd)'"=="stgenreg" {
					qui predict `ch`i'' if `ajtouse', cumhazard at(`at') timevar(`ajtimevar') zeros
					qui replace `ch`i'' = 0 if `ajtouse' & `ch`i''==.
				
				}
				local chvars `chvars' `ch`i''
			
			}
		}
	
	}
	
	if "`std'"!="" {
		local stdcheck = 0
		forvalues i=1/`Ntrans' {
			if "`check'"!="" di in yellow "Transition `i': Standardising over -> `stdvars`i''"
			local stdcheck = `stdcheck' + ("`stdvars`i''"=="")
		}
		if `stdcheck' {
			di as error "No variables for standardising"
			exit 1986
		}
	}
	
	cap drop `todrop'	
	if "`aj'"=="" & "`e(cmd)'"!="cox" {
		mata: predictms()				
	}
	if `survsimcall' {
		cap drop _time
		exit
	}
	
	//=====================================================================================================================================================//
	//normal approximation stuff
	//transform back to prob scales
	//if prediction = 1, need to replace uci and lci with 1's, as logit calc. will give missings
	
	quietly {
		
		if "`ci'"!="" & "`normal'"!="" & "`aj'"=="" {	
			if "`los'"!="" {
				local func exp
				if "`at2'"!="" & "`ratio'"=="" local func
			}
			else {
				local func invlogit
				if "`at2'"!="" {
					if "`ratio'"!="" local func exp
					else local func tanh
				}
			}
		
			local siglev = abs(invnormal((100-`level')/200))
			foreach fromstate in `from' {

				forvalues i=1/`Nstates' {
						
					gen double `stub'_`fromstate'_`i'_lci = `stub'_`fromstate'_`i' - `siglev'* sqrt(`pred_`fromstate'_`i'_se') if _n<=`obs'
					gen double `stub'_`fromstate'_`i'_uci = `stub'_`fromstate'_`i' + `siglev'* sqrt(`pred_`fromstate'_`i'_se') if _n<=`obs'
					
					replace `stub'_`fromstate'_`i' = `func'(`stub'_`fromstate'_`i') if _n<=`obs'
					replace `stub'_`fromstate'_`i'_lci = `func'(`stub'_`fromstate'_`i'_lci) if _n<=`obs'
					replace `stub'_`fromstate'_`i'_uci = `func'(`stub'_`fromstate'_`i'_uci) if _n<=`obs'
					
				}	
			}
		}
		
		if "`ci'"!="" & "`aj'"!="" {	
	
			local siglev = abs(invnormal((100-`level')/200))
			foreach fromstate in `from' {

				forvalues i=1/`Nstates' {
					gen double `stub'_`fromstate'_`i'_lci = `stub'_`fromstate'_`i' - `siglev'* se_`stub'_`fromstate'_`i' if _n<=`obs'
					gen double `stub'_`fromstate'_`i'_uci = `stub'_`fromstate'_`i' + `siglev'* se_`stub'_`fromstate'_`i' if _n<=`obs'
				}	
			}
		}
		
	}
	
	if "`graph'"!="" {
		if "`e(cmd)'"=="cox" local coxgraph cox
		msgraph, from(`from') nstates(`Nstates') timevar(`timevar') enter(`enter') gen(`gen') `graphopts' `coxgraph'
	}
	
	// return list
	return matrix transmatrix = `transmatrix', copy
	return local Nstates = `Nstates'
	return local from `from'
	return local stub `stub'
	
	cap mata: rmexternal(st_local("predictms_struct"))
	//Done

end

//for stgenreg
program updatebmat, eclass
	syntax, N(int) j(int) draws(string)
	
	tempname bmat
	matrix `bmat' = e(b)
	
	forvalues i=1/`n' {
		local var`i': word `i' of `draws'
	}
	
	forvalues i=1/`n' {
		matrix `bmat'[1,`i'] = `var`i''[`j']
	}
	ereturn repost b = `bmat'
end

* double added to range
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

program predictms_modelcheck

	if "`e(cmd)'"!="stms" & "`e(cmd)'"!="stpm2" & "`e(cmd)'"!="stgenreg" & "`e(cmd2)'"!="streg" & "`e(cmd)'"!="cox" & "`e(cmd)'"!="strcs" {
		di as error "Last estimates not found"
		exit 198	
	}	
	
	if "`e(cmd2)'"=="streg" & "`e(cmd)'"=="gamma" {
		di as error "Generalised gamma not supported"
		exit 198
	}

	if "`e(cmd2)'"=="streg" & "`e(cmd)'"=="lnormal" & "`reset'"=="" {
		di as error "Clock forward approach not allowed with streg, dist(lnormal), you must use the clock reset approach"
		exit 198
	}	
	
	if "`e(cmd2)'"=="streg" & ("`e(cmd)'"=="weibull" | "`e(cmd)'"=="ereg") & "`e(frm2)'"== "time" {
		di as error "streg, dist(weib|exp) time, not supported"
		exit 198
	}
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
	if (max(lowertriangle(tmat_ind))>0) {
		errprintf("All elements of the lower triangle of transmatrix() must be coded missing = .\n")
		exit(198)
	}
	row = 1
	rtmat = rows(tmat)
	trans = 1
	while (row<rtmat) {
		for (i=row+1;i<=rtmat;i++) {
			if (sum(tmat:==tmat[row,i])>1 & tmat[row,i]!=.) {
				errprintf("Elements in the upper triangle of transmatrix() are not unique\n")
				exit(198)
			}
			if (tmat[row,i]!=. & tmat[row,i]!=trans){
				errprintf("Elements in the upper triangle of transmatrix() must be sequentially numbered from 1,...,K, where K = number of transitions\n")
				exit(198)
			}		
			if (tmat[row,i]!=.) trans++
		}
		row++
	}
	st_local("Ntrans",strofreal(trans-1))

}
end

