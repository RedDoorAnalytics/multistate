*! version 1.0.0 19may2015 MJC

/*
Mata functions for predictms with strcs
-> predictms_stpm2 is a structure to hold all info about an stpm2 model (full or separate)
-> predictms_stpm2_sepmodel_setup fills a struct for a particular transition
-> stpm2_sim_sep function to simulate from an stpm2 model, for a particular transition
*/

/*
History
?????2016 version 1.0.0
*/
version 12.1

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix
local RM	real matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector

mata:

struct predictms_strcs
{
	`RS' hasvars, hasrcsbase, hasorthog, hastvcs, hascons, Ntvcs, Ngk
	`RR' ln_bknots
	`RM' X, rmat, Xat2, Xtvcsat2, qd
	`RS' hasXtvcs, hasXtvcsat2
	`TR' tvcs_ln_knots, tvcs_rmats, Xtvcs
	`RS' hasbhazard
}

pointer(struct predictms_strcs scalar) scalar predictms_strcs_setup(`SS' trans)
{
	struct predictms_strcs scalar S
	pointer scalar p

	stata("tempname predictms_struct"+trans)
	rmexternal(st_local("predictms_struct"+trans))
	p = crexternal(st_local("predictms_struct"+trans))

	S.Ngk = st_numscalar("e(nodes)")
	S.qd = predictms_gq(S.Ngk)
	
	
	Nats = strtoreal(st_local("Nats"))
	S.hasvars = strtoreal(st_local("Ncovs"+trans))
	K = strtoreal(st_local("K"))
	
	S.hasrcsbase = st_local("rcsbaseoff"+trans)==""
	S.hasorthog = st_local("orthog"+trans)!=""
	if (S.hasrcsbase) {
		S.ln_bknots = strtoreal(tokens(st_local("ln_bknots"+trans)))
		if (S.hasorthog) S.rmat = st_matrix(st_local("rmat"+trans))
	}
	S.hastvcs = st_local("tvc"+trans)!=""
	S.hascons = st_local("nocons"+trans)==""

	S.hasbhazard = st_local("bhazard"+trans)!=""
	if (S.hasbhazard) {
		
		
		
		
	}
	
	if (S.hastvcs) {
		S.Ntvcs = strtoreal(st_local("Ntvcvars"+trans))
		S.Xtvcs = asarray_create("real",1)
		for (at=1;at<=Nats;at++) {
			tvcdm = st_matrix(st_local("at"+strofreal(at)+"dmtvc"+trans))	
			if (K>1) {
				S.hasXtvcs = st_local("at"+strofreal(at)+"tvcstdvarsindex"+trans)!=""
				tvcdm = J(K,1,tvcdm)
				if (S.hasXtvcs) tvcdm[,strtoreal(tokens(st_local("at"+strofreal(at)+"tvcstdvarsindex"+trans)))] = st_data(.,tokens(st_local("at"+strofreal(at)+"tvcstdvars"+trans)),st_local("stdtouse"))
			}
			asarray(S.Xtvcs,at,tvcdm)
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
	
	swap((*p),S)
	return(p)
}

function strcs_sim(	
					`RC' x,												/// -simulated survival time-
					`RC' U,												///	-random draw from U(0,1)- 
					`RC' bmat,											///	-coefficient vector-
					`RR' dm,											///	-design matrix-
					`RS' at, `RS' std,									///
					pointer(struct predictms_strcs scalar) scalar p, |	///
					`RC' t0												///	-entry time-
					)
{
	struct predictms_strcs S
	`RR' DM
	`RC' ret1

	S 			= *p
	nobs 		= rows(x)
	logh_q 		= J(nobs,S.Ngk,0)
	nodelent 	= args()==7
	
	if (nodelent) 	nodes = x:/2 :* J(nobs,1,S.qd[,1]') :+ x:/2
	else  			nodes = (x:-t0):/2 :* J(nobs,1,S.qd[,1]') :+ (x:+t0):/2
	
	for (q=1;q<=S.Ngk;q++) {
		DM = J(nobs,0,.)
		
		if (S.hasvars) DM = DM,J(nobs,1,dm)

		if (S.hasrcsbase) {
			if (S.hasorthog) DM = DM,rcsgen_core(log(nodes[,q]),S.ln_bknots,0,S.rmat)
			else DM = DM,rcsgen_core(log(nodes[,q]),S.ln_bknots,0)
		}

		if (S.hastvcs) {
			Xtvctemp = J(nobs,1,asarray(S.Xtvcs,at)[std,])
			if (S.hasorthog) {
				for (i=1;i<=S.Ntvcs;i++) DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(nodes[,q]),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
			}
			else {
				for (i=1;i<=S.Ntvcs;i++) DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(nodes[,q]),asarray(S.tvcs_ln_knots,i),0))
			}
		}		
		if (S.hascons) DM = DM,J(nobs,1,1)
		logh_q[,q] = DM * bmat
		
		//if (S.hasbhazard) logh_q[,q] = logh_q[,q] //:+ predictms_bhazard(S,)
	}

	if (nodelent) 	ret1 = x:/2 :* (exp(logh_q) * S.qd[,2]) :+ log(U)
	else  			ret1 = (x:-t0):/2 :* (exp(logh_q) * S.qd[,2]) :+ log(U)
	return(ret1)
}

`RC' strcs_ch(	
				`RC' x,												/// -simulated survival time-
				`RC' bmat,											///	-coefficient vector
				`RR' dm,											///	-design matrix-
				`RS' at, `RS' std,									///
				pointer(struct predictms_strcs scalar) scalar p	///
				)
{
	struct predictms_strcs S
	`RR' DM

	S 			= *p
	nobs 		= rows(x)
	logh_q 		= J(nobs,S.Ngk,0)
	
	nodes 		= x:/2 :* J(nobs,1,S.qd[,1]') :+ x:/2

	for (q=1;q<=S.Ngk;q++) {
		DM = J(nobs,0,.)
		
		if (S.hasvars) DM = DM,J(nobs,1,dm)

		if (S.hasrcsbase) {
			if (S.hasorthog) DM = DM,rcsgen_core(log(nodes[,q]),S.ln_bknots,0,S.rmat)
			else DM = DM,rcsgen_core(log(nodes[,q]),S.ln_bknots,0)
		}

		if (S.hastvcs) {
			Xtvctemp = J(nobs,1,asarray(S.Xtvcs,at)[std,])
			if (S.hasorthog) {
				for (i=1;i<=S.Ntvcs;i++) DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(nodes[,q]),asarray(S.tvcs_ln_knots,i),0,asarray(S.tvcs_rmats,i)))
			}
			else {
				for (i=1;i<=S.Ntvcs;i++) DM = DM,(Xtvctemp[1,i] :* rcsgen_core(log(nodes[,q]),asarray(S.tvcs_ln_knots,i),0))
			}
		}		
		if (S.hascons) DM = DM,J(nobs,1,1)
		logh_q[,q] = DM * bmat
	}
	res1 = x:/2 :* (exp(logh_q) * S.qd[,2])
	//correct for time = 0
	time0 = selectindex(x:==0)
	if (rows(time0) & cols(time0)) res1[time0] = J(rows(time0),1,0)
	return(res1)
}

`RM' predictms_gq(`RS' n)
{
	i = range(1,n,1)'
	i1 = range(1,n-1,1)'
	alpha = strtoreal(st_local("alpha"))
	beta = strtoreal(st_local("beta"))

	muzero = 2
	a = J(1,n,0)
	b = i1:/sqrt(4 :* i1:^2 :- 1)

	A= diag(a)
	for(j=1;j<=n-1;j++){
			A[j,j+1] = b[j]
			A[j+1,j] = b[j]
	}
	symeigensystem(A,vec,nodes)
	weights = (vec[1,]:^2:*muzero)'
	weights = weights[order(nodes',1)]
	nodes = nodes'[order(nodes',1)']
	return(nodes,weights)
}

end
