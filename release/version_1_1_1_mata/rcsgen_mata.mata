version 12.1

mata:

real matrix rcsgen_core(	real colvector variable,	///
							real rowvector knots, 		///
							real scalar deriv,|			///
							real matrix rmatrix			///
						)
{
	real scalar Nobs, Nknots, kmin, kmax, interior, Nparams
	real matrix splines, knots2

	//======================================================================================================================================//
	// Extract knot locations

		Nobs 	= rows(variable)
		Nknots 	= cols(knots)
		kmin 	= knots[1,1]
		kmax 	= knots[1,Nknots]
	
		if (Nknots==2) interior = 0
		else interior = Nknots - 2
		Nparams = interior + 1
		
		splines = J(Nobs,Nparams,.)
		
	//======================================================================================================================================//
	// Calculate splines

		if (Nparams>1) {
			lambda = J(Nobs,1,(kmax:-knots[,2..Nparams]):/(kmax:-kmin))
			knots2 = J(Nobs,1,knots[,2..Nparams])
		}

		if (deriv==0) {
			splines[,1] = variable
			if (Nparams>1) {
				splines[,2..Nparams] = (variable:-knots2):^3 :* (variable:>knots2) :- lambda:*((variable:-kmin):^3):*(variable:>kmin) :- (1:-lambda):*((variable:-kmax):^3):*(variable:>kmax) 
			}
		}
		else if (deriv==1) {
			splines[,1] = J(Nobs,1,1)
			if (Nparams>1) {
				splines[,2..Nparams] = 3:*(variable:-knots2):^2 :* (variable:>knots2) :- lambda:*(3:*(variable:-kmin):^2):*(variable:>kmin) :- (1:-lambda):*(3:*(variable:-kmax):^2):*(variable:>kmax) 	
			}
		}
		else if (deriv==2) {
			splines[,1] = J(Nobs,1,0)
			if (Nparams>1) {
				splines[,2..Nparams] = 6:*(variable:-knots2) :* (variable:>knots2) :- lambda:*(6:*(variable:-kmin)):*(variable:>kmin) :- (1:-lambda):*(6:*(variable:-kmax)):*(variable:>kmax) 	
			}
		}
		else if (deriv==3) {
			splines[,1] = J(Nobs,1,0)
			if (Nparams>1) {
				splines[,2..Nparams] = 6:*(variable:>knots2) :- lambda:*6:*(variable:>kmin) :- (1:-lambda):*6:*(variable:>kmax)
			}
		}
		
		//orthog
		if (args()==4) {
			real matrix rmat
			rmat = luinv(rmatrix)
			if (deriv==0) splines = (splines,J(Nobs,1,1)) * rmat[,1..Nparams]
			else splines = splines * rmat[1..Nparams,1..Nparams]
		}
		return(splines)
}
end
