*! Version 1.0.1 11sep2020 MJC

/* 
History
MJC 11sep2020: version 1.0.1 - enter() renamed to ltruncated() for consistency with survsim & predictms
MH 02mar2020: version 1.0.0  - restructured so calculations are all performed in mata (improved efficiency)
							 - bug fix: confidence intervals caused a conformability error and kronecker delta function incorrectly defined; now fixed
							 - removal of gen option
							 - id option no longer required due to restructure
							 - enter time option added
							 - exit time option added
							 - from state option added
							 - los option added
							 - se option added
							 - program now works for bi-directional models
							 - additional error checks							 
PCL 30apr2019: version 0.6.0 - 
*/


program define msaj
	version 14.2
	syntax [if] [in] , 	[		TRANSMatrix(name)	///
								BY(varname) 		///
								LTruncated(real 0)	///
								EXIT(real -99)		///
								FROM(integer 1)		///
								CR					///
								LOS					///
								CI					///
								SE					///
				]

	marksample touse
	if "`by'" != "" {
		qui replace `touse' = 0 if `by' == . 
		qui levelsof `by' if `touse', local(bylevels)
	}
		
		
	// Error checks

	cap confirm variable _t _d _st _t0, exact
	if _rc {
		di as error "Data must be stset (at least one variable of _st, _t, _t0 or _d missing)"
		exit 198
	}
	
	cap confirm variable _trans, exact
	if _rc {
		di as error "Data must be msset (_trans variable missing)"
		exit 198
	}
	
	cap confirm variable _to, exact
	if _rc & "`cr'"!="" {
		di as error "Data must be msset (_to variable missing, needed for CR option)"
		exit 198
	}
	
	cap confirm variable _to _from, exact
	if ("`ci'" != "" | "`se'" != "") & _rc {
		di as error "_from and _to variables must be specified for CI or SE option"
		exit 198
	}	
	
	if ("`transmatrix'" == ""  & "`cr'" == "") | ("`transmatrix'" != ""  & "`cr'" != "") {
		di as error "You must specify either the transition matrix or use the cr option"
		exit 198
	}
	
	if "`transmatrix'" != "" {
		local Nstates = colsof(`transmatrix')
		
		if `Nstates' < 2 {
			di as error "Must be at least 2 possible states, including starting state"
			exit 198
		}	
		
		if `from' < 1 | `from' > `Nstates' {
			di as error "From state must be between 1 and the maximum number of states"
			exit 198
		}
	}
	
	if "`cr'"!="" & `from'!=1 {
		di as error "If CR is specified then from must be 1"
		exit 198
	}
	
	if `ltruncated' < 0 {
		di as error "ltruncated() time must be positive"
		exit 198
	}
	
		
	// Set exit time if not given and check it is after ltruncated time

	summ _t if `touse', meanonly
	local max_t = `r(max)'
	summ _t if `touse' & _d ==1, meanonly
	local max_event = `r(max)'
	
	if `exit' == -99 local exit `max_event'

	if `exit' <= `ltruncated' {
		di as error "ltruncated() time (default 0) is greater than or equal to exit time (default max any event time)"
		exit 198
	}
	
	if `exit' > `max_event' & `exit' <= `max_t' {
		di as error "Warning: Exit time is greater than last event time"
	}	
	
	if `exit' > `max_t' {
		di as error "Exit time is greater than the maximum time"
		exit 198
	}	


	// Build the transition matrix if CR is specified

	if "`cr'" != "" {
		
		summ _to if `touse', meanonly
		local Nstates `r(max)'
		
		if `Nstates' < 2 {
			di as error "Must be at least 2 possible states, including starting state (max _to is 1)"
			exit 198
		}	
		
		tempname transmatrix
		matrix `transmatrix' = J(`Nstates',`Nstates',.)
		forvalues i = 2/`Nstates' {
			local tmptrans = `i' - 1
			matrix `transmatrix'[1,`i'] = `tmptrans'
		}			
	}
	

	// Prepare the output variables

	forvalues i = 1/`Nstates' {
		local newvars `newvars' P_AJ_`i'
		if "`ci'" != "" {
			local newvars `newvars' P_AJ_`i'_lci P_AJ_`i'_uci
		}
		if "`se'" != "" {
			local newvars_se `newvars_se' P_AJ_`i'_se
		}
		if "`los'" != "" {
			local newvars_LOS `newvars_LOS' LOS_AJ_`i'
		}
	}

	local Nnewvars = wordcount("`newvars'")
	forvalues i = 1/`Nnewvars' {
		local tmp = word("`newvars'",`i')
		qui gen double `tmp'=.
	}
	
	if "`se'" != "" {
		local Nnewvars_se = wordcount("`newvars_se'")
		forvalues i = 1/`Nnewvars_se' {
			local tmp = word("`newvars_se'",`i')
			qui gen double `tmp'=.
		}
	}
	
	if "`los'" != "" {
		local Nnewvars_LOS = wordcount("`newvars_LOS'")
		forvalues i = 1/`Nnewvars_LOS' {
			local tmp = word("`newvars_LOS'",`i')
			qui gen double `tmp'=.
		}
	}
	

	// Get the hazards

	tempvar Nrisk Nevents

	sts gen `Nrisk'=n if `touse', by(_trans `by')
	sts gen `Nevents'=d if `touse', by(_trans `by')

	
	// Call mata to apply AJ equations
	mata AJ()

end				
