*! version 1.0.0 19may2015 MJC

/*
Mata functions for predictms
-> predictms_stpm2 is a structure to hold all info about an stpm2 model (full or separate)
-> predictms_stpm2_fullmodel_setup fills the above struct when using a full model
-> predictms_stpm2_sepmodel_setup fills a struct for a particular transition
-> stpm2_sim function to simulate from an stpm2 model, transition DMs stacked
-> stpm2_sim_sep function to simulate from an stpm2 model, for a particular transition	-> should combine with above

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
	`RS' hasvars, hasrcsbase, hasorthog, hastvcs, hascons, Ntvcs, hazard, odds, normal, reset, at2
	real rowvector ln_bknots, Xat2, Xtvcsat2
	`NM' X, rmat, Xtvcs
	`TR' tvcs_ln_knots, tvcs_rmats
}

pointer(struct predictms_stpm2 scalar) scalar predictms_stpm2_fullmodel_setup()
{
	struct predictms_stpm2 scalar S
	pointer scalar p

	stata("tempname predictms_struct")
	rmexternal(st_local("predictms_struct"))
	p = crexternal(st_local("predictms_struct"))

	S.reset = st_local("reset")!=""
	S.hasvars = st_local("dm")!=""
	if (S.hasvars) S.X = st_matrix(st_local("dm"))
	
	S.hasrcsbase = st_local("rcsbaseoff")==""
	S.hasorthog = st_local("orthog")!=""
	if (S.hasrcsbase) {
		S.ln_bknots = strtoreal(tokens(st_local("ln_bknots")))
		if (S.hasorthog) S.rmat = st_matrix(st_local("rmat"))
	}
	
	S.hastvcs = st_local("tvc")!=""
	S.hascons = st_local("nocons")==""
	
	if (S.hastvcs) {
		S.Ntvcs = strtoreal(st_local("Ntvcvars"))
		S.Xtvcs = st_matrix(st_local("dmtvc"))		
		S.tvcs_ln_knots = asarray_create("real",1)
		for (i=1;i<=S.Ntvcs;i++) {
			asarray(S.tvcs_ln_knots,i,strtoreal(tokens(st_local("ln_tvcknots_"+strofreal(i)))))
		}
		if (S.hasorthog) {
			S.tvcs_rmats = asarray_create("real",1)
			for (i=1;i<=S.Ntvcs;i++) {
				asarray(S.tvcs_rmats,i,st_matrix(st_local("R_"+strofreal(i))))
			}		
		}
	}
	
	S.hazard = st_global("e(scale)")=="hazard"
	S.odds = st_global("e(scale)")=="odds"
	S.normal = st_global("e(scale)")=="normal"
	
	swap((*p), S)
	return(p)
}


pointer(struct predictms_stpm2 scalar) scalar predictms_stpm2_sepmodel_setup(`SS' trans)
{
	struct predictms_stpm2 scalar S
	pointer scalar p

	stata("tempname predictms_struct"+trans)
	rmexternal(st_local("predictms_struct"+trans))
	p = crexternal(st_local("predictms_struct"+trans))

	S.at2 = st_local("at2")!=""
	S.reset = st_local("reset")!=""
	S.hasvars = strtoreal(st_local("Ncovs"+trans))
	if (S.hasvars) {
		S.X = st_matrix(st_local("dm"+trans))
		if (S.at2) S.Xat2 = st_matrix(st_local("at2dm"+trans))
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
		if (S.at2) S.Xtvcsat2 = st_matrix(st_local("at2dmtvc"+trans))
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
						`RS' x,											/// -simulated survival time-
						`RS' U,											///	-random draw from U(0,1)- 
						`RC' bmat,										///	-coefficient vector-
						`RS' transid,									///	-indicate which transition for DMs-
						pointer(struct predictms_stpm2 scalar) scalar p, |	///
						`RS' t0											///	-entry time-
					)
{
	struct predictms_stpm2 S
	`RR' DM
	
	S = *p
	DM = J(1,0,.)
	if (S.hasvars) DM = DM,S.X[transid,]

	if (args()==6 & t0>0) {
		`RR' DM0
		DM0 = DM
		if (S.hasrcsbase) {
			if (S.hasorthog) {
				DM = DM,rcsgen_core(log(x),S.ln_bknots,0,S.rmat)
				DM0 = DM0,rcsgen_core(log(t0),S.ln_bknots,0,S.rmat)
			}
			else {
				DM = DM,rcsgen_core(log(x),S.ln_bknots,0)
				DM0 = DM0,rcsgen_core(log(t0),S.ln_bknots,0)
			}
		}
		if (S.hastvcs) {
			if (S.hasorthog) {
				for (i=1;i<=S.Ntvcs;i++) {
					DM = DM,(S.Xtvcs[transid,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
					DM0 = DM0,(S.Xtvcs[transid,i] :* rcsgen_core(log(t0),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
				}
			}
			else {
				for (i=1;i<=S.Ntvcs;i++) {
					DM = DM,(S.Xtvcs[transid,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0))
					DM0 = DM0,(S.Xtvcs[transid,i] :* rcsgen_core(log(t0),asarray(S.tvcs_ln_knots,i),0))
				}
			}
		}		
		if (S.hascons) {
			DM = DM,1
			DM0 = DM0,1
		}

		if (S.hazard) return(exp(-exp(DM * bmat)) :- U :* exp(-exp(DM0 * bmat)))
		else if (S.odds) return(1:/(exp(DM * bmat):+1) :- U :/ (exp(DM0 * bmat):+1))
		else if (S.normal) return((1:-normal(DM * bmat)) :- U :* (1:-normal(DM0 * bmat)))
	}
	else {

		if (S.hasrcsbase) {
			if (S.hasorthog) DM = DM,rcsgen_core(log(x),S.ln_bknots,0,S.rmat)
			else DM = DM,rcsgen_core(log(x),S.ln_bknots,0)
		}
		if (S.hastvcs) {
			if (S.hasorthog) {
				for (i=1;i<=S.Ntvcs;i++) {
					DM = DM,(S.Xtvcs[transid,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
				}
			}
			else {
				for (i=1;i<=S.Ntvcs;i++) {
					DM = DM,(S.Xtvcs[transid,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0))
				}
			}
		}		
		if (S.hascons) DM = DM,1
		if (S.hazard) return(exp(-exp(DM * bmat)) :- U)
		else if (S.odds) return(1:/(exp(DM * bmat):+1) :- U)
		else if (S.normal) return((1:-normal(DM * bmat)) :- U)
	}
}

function stpm2_sim_sep(	
							`RS' x,											/// -simulated survival time-
							`RS' U,											///	-random draw from U(0,1)- 
							`RC' bmat,										///	-coefficient vector-
							`RS' dmind,										/// -- at or at2 = 1 or 5
							pointer(struct predictms_stpm2 scalar) scalar p, |	///
							`RS' t0											///	-entry time-
						)
{
	struct predictms_stpm2 S
	`RR' DM
	S = *p
	DM = J(1,0,.)

	if (S.hasvars) {
		if (dmind==5) DM = DM,S.Xat2
		else DM = DM,S.X
	}
	if (args()==6 & t0>0) {
		`RR' DM0
		DM0 = DM
		if (S.hasrcsbase) {
			if (S.hasorthog) {
				DM = DM,rcsgen_core(log(x),S.ln_bknots,0,S.rmat)
				DM0 = DM0,rcsgen_core(log(t0),S.ln_bknots,0,S.rmat)
			}
			else {
				DM = DM,rcsgen_core(log(x),S.ln_bknots,0)
				DM0 = DM0,rcsgen_core(log(t0),S.ln_bknots,0)
			}
		}

		if (S.hastvcs) {
			if (dmind==5) Xtvctemp = S.Xtvcsat2
			else Xtvctemp = S.Xtvcs
			if (S.hasorthog) {
				for (i=1;i<=S.Ntvcs;i++) {
					DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
					DM0 = DM0,(Xtvctemp[1,i] :* rcsgen_core(log(t0),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
				}
			}
			else {
				for (i=1;i<=S.Ntvcs;i++) {
					DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0))
					DM0 = DM0,(Xtvctemp[1,i] :* rcsgen_core(log(t0),asarray(S.tvcs_ln_knots,i),0))
				}
			}
		}		
		if (S.hascons) {
			DM = DM,1
			DM0 = DM0,1
		}

		if (S.hazard) return(exp(-exp(DM * bmat)) :- U :* exp(-exp(DM0 * bmat)))
		else if (S.odds) return(1:/(exp(DM * bmat):+1) :- U :/ (exp(DM0 * bmat):+1))
		else if (S.normal) return((1:-normal(DM * bmat)) :- U :* (1:-normal(DM0 * bmat)))
	}
	else {
		if (S.hasrcsbase) {
			if (S.hasorthog) DM = DM,rcsgen_core(log(x),S.ln_bknots,0,S.rmat)
			else DM = DM,rcsgen_core(log(x),S.ln_bknots,0)
		}
		
		if (S.hastvcs) {
			if (dmind==5) Xtvctemp = S.Xtvcsat2
			else Xtvctemp = S.Xtvcs		
			if (S.hasorthog) {
				for (i=1;i<=S.Ntvcs;i++) DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
			}
			else {
				for (i=1;i<=S.Ntvcs;i++) DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(x),asarray(S.tvcs_ln_knots,i),0))
			}
		}		

		if (S.hascons) DM = DM,1
		if (S.hazard) return(exp(-exp(DM * bmat)) :- U)
		else if (S.odds) return(1:/(exp(DM * bmat):+1) :- U)
		else if (S.normal) return((1:-normal(DM * bmat)) :- U)
	}
}


end

