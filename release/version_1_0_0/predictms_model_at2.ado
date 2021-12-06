
program predictms_model_at2, rclass
	gettoken model 0 : 0
	gettoken Nparams at : 0
		
	tempname dm
	
	local varlist `e(varlist)'

	if "`model'"!="stpm2" {
		
		mat `dm' = J(1,`Nparams',0)
		
		local cmdline `e(cmdline)'
		gettoken cmd 0 : cmdline
		syntax [varlist(default=empty)] [if] [in], [NOCONstant ANCillary(varlist) ANC2(varlist) *]		//!!change to substr check that nocons is in it
		
		local colindex = 1
		foreach corevar in `varlist' {
			tokenize `at'
			while "`1'"!="" {
				unab 1: `1'
				if "`corevar'"=="`1'" {
					mat `dm'[1,`colindex'] = `2'
				}
				mac shift 2
			} 
			local colindex = `colindex' + 1
		}
		if "`e(noconstant)'"=="" {
			mat `dm'[1,`colindex'] = 1
			local colindex = `colindex' + 1
		}
		
		if "`model'"!="exp" {
			local Nmleqns = 2
			foreach corevar in `ancillary' {
				tokenize `at'
				while "`1'"!="" {
					unab 1: `1'
					if "`corevar'"=="`1'" {
						mat `dm'[1,`colindex'] = `2'
					}
					mac shift 2
				} 
				local colindex = `colindex' + 1
			}
			mat `dm'[1,`colindex'] = 1
			local colindex = `colindex' + 1
				
			if "`model'"=="gamma" {
				local Nmleqns = 3
				foreach corevar in `anc2' {
					tokenize `at'
					while "`1'"!="" {
						unab 1: `1'
						if "`corevar'"=="`1'" {
							mat `dm'[1,`colindex'] = `2'
						}
						mac shift 2
					} 
					local colindex = `colindex' + 1
				}
				mat `dm'[1,`colindex'] = 1
			}
			
		}
		
		return matrix dm = `dm'		
	}
	else if "`model'"=="stpm2" {
	
		//DM is only for varlist, tvc splines and base splines are handled separately
		local Ncovs : word count `varlist'

		if `Ncovs' > 0 {
			tempname dm 
			mat `dm' = J(1,`Ncovs',0)
			//now update DM and indices
			//can match variables in trans#() with varlist and ancillary
			local colindex = 1
			foreach corevar in `varlist' {
				tokenize `at'
				while "`1'"!="" {
					unab 1: `1'
					if "`corevar'"=="`1'" {
						mat `dm'[1,`colindex'] = `2'
					}
					mac shift 2
				} 
				local colindex = `colindex' + 1
			}
			return matrix dm = `dm'
		}
				
		if "`e(tvc)'"!="" {

			local Ntvcvars : word count `e(tvc)'
			tempname dmtvc
			mat `dmtvc' = J(1,`Ntvcvars',0)
			local colindex = 1
			foreach corevar in `e(tvc)' {
				tokenize `at'
				while "`1'"!="" {
					unab 1: `1'
					if "`corevar'"=="`1'" {
						mat `dmtvc'[1,`colindex'] = `2'
					}
					mac shift 2
				} 
				local colindex = `colindex' + 1
			}
			
			return matrix dmtvc = `dmtvc'
		}
	
	}



end

