*! version 1.0.0 ?????2015 MJC

version 12.1

local ss 	string scalar
local RS	real scalar
local NM	numeric matrix
local RM 	real matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local PS	pointer(struct predictms_struct scalar) scalar
local SS	struct predictms_struct scalar
local gml 	struct merlin_struct scalar
local cgml 	struct merlin_struct colvector
local Pcm	pointer(struct merlin_struct scalar) colvector

mata:

void predictms_analytic_cr(`SS' S, `Pcm' Pmerlin, `RS' Nobs)
{
	t = S.predtime
	if (S.getprobs) {
		pred = predictms_analytic_cr_p(S,Pmerlin,Nobs,t)
		if (S.standardise) 	S.pt = S.pt :+ pred
		else 			S.pt = pred
	}

	if (S.getlos | S.getrmst) {
		Nqp 	= S.chips
		gq 	= predictms_gq(Nqp)
		qp	= (t:-S.enter) :/ 2 :* J(Nobs,1,gq[,1]') ///
				:+ (t:+S.enter) :/2
		
		pred	= J(Nobs,S.Nstates,0)
		for (q=1; q<=Nqp; q++) {
			pred = pred :+ 					///
			predictms_analytic_cr_p(S,Pmerlin,Nobs,qp[,q]) 	///
				:* gq[q,2]
		}
		pred = pred :* (t:-S.enter) :/ 2
		
		if (S.standardise) 	S.los = S.los :+ pred
		else 			S.los = pred
	}
	
	if (S.getrmst) {
		if (S.standardise) 	S.rmst = S.rmst :+ pred[,1]
		else 			S.rmst = pred[,1]
	}

	if (S.hasuser) predictms_calc_user(S)
}

`RM' predictms_analytic_cr_p(`SS' S, `Pcm' Pmerlin, `RS' Nobs, `RC' t)
{
	`gml' gml
	Nqp 	= S.chips
	gq 	= predictms_gq(Nqp)
	pred 	= J(Nobs,S.Nstates,0)
	qp	= (t:-S.enter) :/ 2 :* J(Nobs,1,gq[,1]') :+ (t:+S.enter) :/2

	//get total S at qps
	totalS 	= J(Nobs,Nqp,0)
	for (s=1;s<=S.Ntrans;s++) {
		gml = *Pmerlin[s]
		for (q=1; q<=Nqp; q++) {						
			ch = (*gml.Pch[1])(gml,qp[,q])
			_editmissing(ch,0)
			totalS[,q] = totalS[,q] :+ ch
		}
	}
	totalS = exp(-totalS)
	
	nextstates = asarray(S.posnextstates,1)
	for (s=1;s<=S.Ntrans;s++) {	
		gml 	= *Pmerlin[s]
		cifs	= J(Nobs,Nqp,.)
		for (q=1; q<=Nqp; q++) {
			cifs[,q] = exp((*gml.Plogh[1])(gml,qp[,q])) :* totalS[,q]
		}
		pred[,nextstates[s]] =  (cifs * gq[,2]) :* (t:-S.enter):/2
		
		if (missing(pred)==Nobs) merlin_error("Error in numerical integration; try the -simulate- option")
		_editmissing(pred,0)
	}
	
	if (S.enter) {
		for (s=1;s<=S.Ntrans;s++) {
			gml = *Pmerlin[s]
			pred[,2..S.Nstates] = pred[,2..S.Nstates] 	///
				:/ exp(-(*gml.Pch[1])(gml,J(Nobs,1,S.enter)))
		}
	}
	
	pred[,1] = 1:-rowsum(pred[,2..S.Nstates])
	return(pred)
}

void predictms_analytic_cr_stand(`SS' S, `Pcm' Pmerlin, `RS' Nobs)
{
	`RC' 	t
	`RS'	Nt
	`gml' 	gml
	`cgml' 	cgml
	
	t = S.predtime
	Nt = rows(t)
	cgml = J(S.Ntrans,1,gml)
	for (c=1;c<=S.Ntrans;c++) cgml[c] = *Pmerlin[c]
	
	if (S.getprobs) {
		for (i=1;i<=Nt;i++) {
			S.pt[i,] = predictms_analytic_cr_p_stand(S,cgml,Nobs,t[i])
		}
	}
	
	if (S.getlos | S.getrmst) {
		Nqp = S.chips
		gq  = predictms_gq(Nqp)
		qp  = (t:-S.enter) :/ 2 :* J(Nt,1,gq[,1]') :+ (t:+S.enter) :/2
		for (i=1;i<=Nt;i++) {
			res = J(1,S.Nstates,0)
			for (j=1;j<=Nqp;j++) {
				res = res :+ predictms_analytic_cr_p_stand(S,
						cgml,Nobs,qp[i,j])	///
					:* gq[j,2] 
			}
			S.los[i,] = res :* (t[i]:-S.enter):/2
		}	
	}	
	
	if (S.getrmst) S.rmst = S.los[,1]
	if (S.hasuser) predictms_calc_user(S)
}

`RM' predictms_analytic_cr_p_stand(`SS' S, `cgml' cgml, `RS' Nobs, `RS' t)
{
	Nqp 	= S.chips
	gq 	= predictms_gq(Nqp)
	pred 	= J(1,S.Nstates,.)
	
	qp	= (t:-S.enter) :/ 2 :* gq[,1]' :+ (t:+S.enter) :/2
	qp	= J(S.K,1,qp)
	
	//get total S at qps
	totalS 	= J(S.K,Nqp,0)
	
	for (s=1;s<=S.Ntrans;s++) {
		for (q=1; q<=Nqp; q++) {						
			ch = (*cgml[s].Pch[1])(cgml[s],qp[,q])
			_editmissing(ch,0)
			totalS[,q] = totalS[,q] :+ ch
		}
	}
	totalS = exp(-totalS)
	nextstates = asarray(S.posnextstates,1)
	for (s=1;s<=S.Ntrans;s++) {	
		cifs	= J(S.K,Nqp,.)
		for (q=1; q<=Nqp; q++) {
			cifs[,q] = exp((*cgml[s].Plogh[1])(cgml[s],qp[,q])) :* totalS[,q]
		}
		pred[,nextstates[s]] =  mean((cifs * gq[,2]) :* (t:-S.enter):/2)
		
		if (missing(pred)==S.K) merlin_error("Error in numerical integration; try the -simulate- option")
		_editmissing(pred,0)
	}

	if (S.enter) {
		for (s=1;s<=S.Ntrans;s++) {
			pred[,2..S.Nstates] = pred[,2..S.Nstates] 	///
				:/ mean(exp(-(*cgml[s].Pch[1])(cgml[s],J(S.K,1,S.enter))))
		}
	}
	
	pred[,1] = 1:-rowsum(pred[,2..S.Nstates])
	return(pred)
}

end
