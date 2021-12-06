
program predictms_model, rclass
	gettoken trans 0 : 0
	gettoken Nparams at : 0

	local cmdname "`e(cmd)'"
	if "`e(cmd2)'"=="streg" {
		
		//design matrix for transition
		tempname dm indices
		mat `dm' = J(1,`Nparams',0)
	
		if "`e(cmd)'"=="ereg" {						//!! don't need index for exp
			
			local cmdline `e(cmdline)'
			gettoken cmd 0 : cmdline
			syntax [varlist(default=empty)] [if] [in], [NOCONstant *]

			//now update DM and indices
			//can match variables in trans#() with varlist
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
			if "`noconstant'"=="" {
				mat `dm'[1,`colindex'] = 1
			}
	
		}
		else {
		
			local Nmleqns = 2
			local cmdline `e(cmdline)'
			gettoken cmd 0 : cmdline
			syntax [varlist(default=empty)] [if] [in], [NOCONstant ANCillary(varlist) *]		//!!change to substr check that nocons is in it
			local corevars1 `varlist'		
			local corevars2 `ancillary'
			
			//now update DM and indices
			//can match variables in trans#() with varlist and ancillary
			local colindex = 1
			foreach corevar in `corevars1' {
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
			if "`noconstant'"=="" {
				mat `dm'[1,`colindex'] = 1
				local colindex = `colindex' + 1
			}
			foreach corevar in `corevars2' {
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
			
			//indices
			mat `indices' = J(2,`Nmleqns',1)
			local colindex = 1
			foreach corevar in `corevars1' {
				local colindex = `colindex' + 1
			}
			if "`noconstant'"=="" {
				local colindex = `colindex' + 1
			}
			mat `indices'[2,1] = `colindex' - 1
			mat `indices'[1,2] = `colindex'
			mat `indices'[2,2] = `Nparams'
			
			return matrix indices = `indices'
		}
		
		return matrix dm = `dm'
			
	}
	else if "`e(cmd)'"=="stpm2" {
	
		//DM is only for varlist, tvc splines and base splines are handled separately
		local corevars `e(varlist)'
		local Ncovs : word count `corevars'
		c_local Ncovs`trans' `Ncovs'
		c_local nocons`trans' `e(noconstant)'
		c_local orthog`trans' `e(orthog)'
		c_local scale`trans' `e(scale)'
		
		//overall design matrix for each transition, stacked
		if `Ncovs' > 0 {
			tempname dm 
			mat `dm' = J(1,`Ncovs',0)
			//now update DM and indices
			//can match variables in trans#() with varlist and ancillary
			local colindex = 1
			foreach corevar in `corevars' {
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
		
		c_local rcsbaseoff`trans' `e(rcsbaseoff)'
		if "`e(rcsbaseoff)'"=="" {
			local Nsplines : word count `e(rcsterms_base)'
			c_local ln_bknots`trans' `e(ln_bhknots)'										//all log baseline knots including boundary knots
			if "`e(ln_bhknots)'"=="" {	//this is empty when df(1)
				c_local ln_bknots`trans' `=log(`: word 1 of `e(boundary_knots)'')' `=log(`: word 2 of `e(boundary_knots)'')'
			}
			if "`e(orthog)'"=="orthog" {
				tempname rmat
				matrix `rmat' = e(R_bh)
				return matrix rmat = `rmat'
			}			
		}
				
		c_local tvc`trans' `e(tvc)'
		if "`e(tvc)'"!="" {
			local i = 1
			foreach tvcvar in `e(tvc)' {
				local boundary_knots_`i' `e(boundary_knots_`tvcvar')'
				local ln_tvcknots`trans'_`i' `e(ln_tvcknots_`tvcvar')'
				if "`ln_tvcknots`trans'_`i''"=="" {
					local ln_tvcknots`trans'_`i' `=log(`: word 1 of `boundary_knots_`i''')' `=log(`: word 2 of `boundary_knots_`i''')'
				}				
				if "`e(orthog)'"=="orthog" {
					tempname R_`i'
					mat `R_`i'' = e(R_`tvcvar')
					return matrix R_`i' = `R_`i''
				}
				c_local ln_tvcknots`trans'_`i' `ln_tvcknots`trans'_`i''
				local i = `i' + 1
			}
			local Ntvcvars = `i' - 1
			c_local Ntvcvars`trans' = `Ntvcvars'

			//tvc DM
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

	c_local cmdname `cmdname'


end

