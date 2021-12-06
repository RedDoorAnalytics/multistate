*! version 1.0.1 03feb2020 MJC

/*
Graph utility for predictms to plot stacked transition probabilities
*/


program graphms
	syntax  ,	[								///
					FROM(numlist >0 asc int)	/// -starting state for predictions-
					AT(numlist >0 asc int)		/// -at# patterns to plot (separately)-
					Nstates(numlist max=1 int)	///	-# of potential states to be in-
					TIMEvar(string)				///	-timevar-
					*							///	-graph opts-
				]

	if "`from'"=="" {
		local from `r(from)'
		if "`from'"=="" {
			di as error "from() required"
			exit 198
		}
	}
	if "`at'"=="" {
		local at = 1
	}
	
	if "`nstates'"=="" { 
		local nstates `r(Nstates)'
		if "`nstates'"=="" {
			di as error "nstates() required"
			exit 198
		}
	}
	
	local stub _prob_at
	if "`timevar'"=="" {
		local timevar `r(timevar)'
		if "`timevar'"=="" {
			di as error "timevar() required"
			exit 198
		}
	}
	
	qui su `timevar'
	local mint = `r(min)'
	
	cap numlist "`nstates'/1"
	local legorder `r(numlist)'
	
	foreach a in `at' {
		foreach f in `from' {
			tempvar plotvars1_`a'_`f'
			qui gen double `plotvars1_`a'_`f'' = `stub'`a'_`f'_1
			label var `plotvars1_`a'_`f'' "P(Y(t) = 1 | Y(`mint') = `f')"
			
			local plots_`a'_`f' `plots_`a'_`f'' (area `plotvars1_`a'_`f'' `timevar', base(0)) 
			forvalues i=2/`nstates' {
				tempvar plotvars`i'_`a'_`f'
				qui gen double `plotvars`i'_`a'_`f'' = `plotvars`=`i'-1'_`a'_`f'' + `stub'`a'_`f'_`i'
				label var `plotvars`i'_`a'_`f'' "P(Y(t) = `i' | Y(`mint') = `f')"
				local plots_`a'_`f' (area `plotvars`i'_`a'_`f'' `timevar', base(0)) `plots_`a'_`f'' 
			}                   
			
			twoway `plots_`a'_`f'', 								///
					ylabel(0(0.1)1, angle(h) format(%2.1f)) 		///
					ytitle("Transition Probability") 				///
					xtitle("Follow-up time (t)")					///
					title("P(Y(t) = b | Y(`mint') = `f', at`a')")	///
					plotregion(m(zero))								///
					legend(order(`legorder')) 						///
					name(gh_`a'_`f', replace)						///
					`options'										//
		}
	}
                        
end


