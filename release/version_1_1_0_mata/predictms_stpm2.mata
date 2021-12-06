*! version 1.0.0 19may2015 MJC

/*
Mata functions for predictms
-> predictms_stpm2 is a structure to hold all info about an stpm2 model (full or separate)
-> predictms_stpm2_sepmodel_setup fills a struct for a particular transition
-> stpm2_sim_sep function to simulate from an stpm2 model, for a particular transition
*/

/*
History
19may 2015 version 1.0.0
*/

version 12.1

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector

mata:

struct predictms_stpm2
{
	`RS' hasvars, hasrcsbase, hasorthog, hastvcs, hascons, Ntvcs, hazard, odds, normal, at2
	real rowvector ln_bknots
	`NM' X, rmat, Xtvcs, Xat2, Xtvcsat2
	`RS' hasXtvcs, hasXtvcsat2
	`TR' tvcs_ln_knots, tvcs_rmats
}

pointer(struct predictms_stpm2 scalar) scalar predictms_stpm2_setup(`SS' trans)
{
	struct predictms_stpm2 scalar S
	pointer scalar p

	stata("tempname predictms_struct"+trans)
	rmexternal(st_local("predictms_struct"+trans))
	p = crexternal(st_local("predictms_struct"+trans))

	S.at2 = st_local("at2")!=""
	S.hasvars = strtoreal(st_local("Ncovs"+trans))
	K = strtoreal(st_local("K"))
	if (S.hasvars) {
		S.X = st_matrix(st_local("dm"+trans))
		if (K>1) {
			S.X = J(K,1,S.X)
			S.X[,strtoreal(tokens(st_local("stdvarsindex"+trans)))] = st_data(.,tokens(st_local("stdvars"+trans)),st_local("stdtouse"))
		}

		if (S.at2) {
			S.Xat2 = st_matrix(st_local("at2dm"+trans))
			if (K>1) {
				S.Xat2 = J(K,1,S.Xat2)
				S.Xat2[,strtoreal(tokens(st_local("at2stdvarsindex"+trans)))] = st_data(.,tokens(st_local("at2stdvars"+trans)),st_local("stdtouse"))
			}	
		}
	}

	S.hasrcsbase = st_local("rcsbaseoff"+trans)==""
	S.hasorthog = st_local("orthog"+trans)!=""
	if (S.hasrcsbase) {
		S.ln_bknots = strtoreal(tokens(st_local("ln_bknots"+trans)))
		if (S.hasorthog) S.rmat = st_matrix(st_local("rmat"+trans))
	}
	S.hastvcs = st_local("tvc"+trans)!=""
	S.hascons = st_local("nocons"+trans)==""

	if (S.hastvcs) {
		S.Ntvcs = strtoreal(st_local("Ntvcvars"+trans))
		S.Xtvcs = st_matrix(st_local("dmtvc"+trans))	
		if (K>1) {
			S.hasXtvcs = st_local("tvcstdvarsindex"+trans)!=""
			S.Xtvcs = J(K,1,S.Xtvcs)
			if (S.hasXtvcs) S.Xtvcs[,strtoreal(tokens(st_local("tvcstdvarsindex"+trans)))] = st_data(.,tokens(st_local("tvcstdvars"+trans)),st_local("stdtouse"))
		}
		if (S.at2) {
			S.Xtvcsat2 = st_matrix(st_local("at2dmtvc"+trans))
			if (K>1) {
				S.hasXtvcsat2 = st_local("tvcat2stdvarsindex"+trans)!=""
				S.Xtvcsat2 = J(K,1,S.Xtvcsat2)
				if (S.hasXtvcsat2) S.Xtvcsat2[,strtoreal(tokens(st_local("tvcat2stdvarsindex"+trans)))] = st_data(.,tokens(st_local("tvcat2stdvars"+trans)),st_local("stdtouse"))
			}
		}
		S.tvcs_ln_knots = asarray_create("real",1)
		for (i=1;i<=S.Ntvcs;i++) {
			asarray(S.tvcs_ln_knots,i,strtoreal(tokens(st_local("ln_tvcknots"+trans+"_"+strofreal(i)))))
		}
		if (S.hasorthog) {
			S.tvcs_rmats = asarray_create("real",1)
			for (i=1;i<=S.Ntvcs;i++) {
				asarray(S.tvcs_rmats,i,st_matrix(st_local("R"+trans+"_"+strofreal(i))))
			}		
		}
	}

	S.hazard = st_local("scale"+trans)=="hazard"
	S.odds = st_local("scale"+trans)=="odds"
	S.normal = st_local("scale"+trans)=="normal"
	
	swap((*p), S)
	return(p)
}

function stpm2_sim(	
					`RC' x,												/// -simulated survival time-
					`RC' U,												///	-random draw from U(0,1)- 
					`RC' bmat,											///	-coefficient vector-
					`RS' dmind,											/// -at or at2 = 1 or 5-
					`RS' stdind,										///	-index of DM for standardisation-
					pointer(struct predictms_stpm2 scalar) scalar p, |	///
					`RC' t0												///	-entry time-
					)
{
	struct predictms_stpm2 S
	`RR' DM
	`RC' ret1

	S = *p
	nobs = rows(x)
	DM = J(nobs,0,.)
	if (S.hasvars) {
		if (dmind==5) DM = DM,J(nobs,1,S.Xat2[stdind,])
		else DM = DM,J(nobs,1,S.X[stdind,])
	}

	if (S.hasrcsbase) {
		if (S.hasorthog) DM = DM,rcsgen_core(log(x),S.ln_bknots,0,S.rmat)
		else DM = DM,rcsgen_core(log(x),S.ln_bknots,0)
	}

	if (S.hastvcs) {
		if (dmind==5) Xtvctemp = J(nobs,1,S.Xtvcsat2[stdind,])
		else Xtvctemp = J(nobs,1,S.Xtvcs[stdind,])
		if (S.hasorthog) {
			for (i=1;i<=S.Ntvcs;i++) DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
		}
		else {
			for (i=1;i<=S.Ntvcs;i++) DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0))
		}
	}		
	if (S.hascons) DM = DM,J(nobs,1,1)

	if (S.hazard) 		ret1 = exp(DM * bmat) :+ log(U)
	else if (S.odds) 	ret1 = 1:/(exp(DM * bmat):+1) :- U
	else if (S.normal) 	ret1 = (1:-normal(DM * bmat)) :- U

	if (args()>6) {

		`RC' index

		index = 1::rows(x)
		index = select(index,t0:>0)
		nobs2 = rows(index)

		if (nobs2) {
			`RM' DM0
			DM0 = J(nobs2,0,.)
			if (S.hasvars) {
				if (dmind==5) DM0 = DM0,J(nobs2,1,S.Xat2[stdind,])
				else DM0 = DM0,J(nobs2,1,S.X[stdind,])
			}
	
			if (S.hasrcsbase) {
				if (S.hasorthog) DM0 = DM0,rcsgen_core(log(t0[index,]),S.ln_bknots,0,S.rmat)
				else DM0 = DM0,rcsgen_core(log(t0[index,]),S.ln_bknots,0)
			}
			if (S.hastvcs) {
				if (dmind==5) Xtvctemp = J(nobs2,1,S.Xtvcsat2[stdind,])
				else Xtvctemp = J(nobs2,1,S.Xtvcs[stdind,])
				if (S.hasorthog) {
					for (i=1;i<=S.Ntvcs;i++) DM0 = DM0,(Xtvctemp[,i] :* rcsgen_core(log(t0[index,]),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
				}
				else {
					for (i=1;i<=S.Ntvcs;i++) DM0 = DM0,(Xtvctemp[,i] :* rcsgen_core(log(t0[index,]),asarray(S.tvcs_ln_knots,i),0))
				}
			}		
			if (S.hascons) DM0 = DM0,J(nobs2,1,1)

			if (S.hazard) 		ret1[index,] = ret1[index,] :- exp(DM0 * bmat)
			else if (S.odds) 	ret1[index,] = ret1[index,] :+ (1:-U[index,]) :/ (exp(DM0 * bmat):+1)
			else if (S.normal) 	ret1[index,] = ret1[index,] :+ (1:-U[index,]) :* (1:-normal(DM0 * bmat))
		}
	}
	return(ret1)
}

/*
function stpm2_sim_sep_vec_shell(	
							`RC' x,												/// -simulated survival time-
							`RC' U,												///	-random draw from U(0,1)- 
							`RM' bmat,											///	-coefficient vector-
							`RS' dmind,											/// -at or at2 = 1 or 5-
							`RS' stdind,										///	-index of DM for standardisation-
							pointer(struct predictms_stpm2 scalar) colvector p, |	///
							`RC' t0												///	-entry time-
						)
{

	ind = rows(p)
	ret1 = log(U)
	for (i=1;i<=ind;i++) ret1 = ret1 :+ stpm2_sim_sep_vec(x,bmat[,i],dmind,stdind,p[i],t0)
	return(ret1)

}
*/

end

