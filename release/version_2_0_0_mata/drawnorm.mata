



version 14.1

local SS 	string scalar
local RS	real scalar
local RM	real matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector

mata:

`RM' drawnormm( 					///
					`RC' mvec,		///
					`RM' V,			///
					`RS' n			///
					)	
{
	`RS' nvars
	`RM' cholV, z, res
	
	nvars = rows(mvec)
	if (nvars!=rows(V) & nvars!=cols(V)) {
		errprintf("Dimensions not valid\n")
		exit(198)
	}
	
	cholV = cholesky(V)
	z = rnormal(nvars,n,0,1)

	/* Draw via x = mu + A*z, repeat this m times */
	res = J(nvars,n,.)
	for (i=1; i<=n; i++) res[.,i] = mvec + A*z[., i]
	
	return(res')
}

end
