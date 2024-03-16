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
local Pcm	pointer(struct merlin_struct scalar) colvector

mata:

void predictms_analytic_surv(`SS' S, `Pcm' Pmerlin, `RS' Nobs)
{
	struct merlin_struct gml
	
	t 	= S.predtime
	gml	= *Pmerlin[1]
	
	if (S.getprobs) {
		pred 		= J(Nobs,2,.)
		ch			= (*gml.Pch[1])(gml,t)
		_editmissing(ch,0)
		pred[,1] 	= exp(-ch)
		if (S.enter) {
			pred[,1] = pred[,1] :/ exp(-(*gml.Pch[1])(gml,J(Nobs,1,S.enter)))
		}
		pred[,2]	= 1:-pred[,1]

		if (S.standardise) 	S.pt = S.pt :+ pred
		else 			S.pt = pred
	}
	
	if (S.getlos | S.getrmst) {
		
		pred	= J(Nobs,2,.)
		Nqp 	= S.chips
		gq 	= predictms_gq(Nqp)
		qp	= (t:-S.enter) :/ 2 :* J(Nobs,1,gq[,1]') :+ (t:+S.enter) :/2
		
		rmst	= J(Nobs,Nqp,.)
		for (q=1; q<=Nqp; q++) {
			rmst[,q] = exp(-(*gml.Pch[1])(gml,qp[,q]))
		}
		pred[,1] = rmst * gq[,2] :* (t:-S.enter) :/ 2
		if (S.enter) pred[,1]  = pred[,1] :/ exp(-(*gml.Pch[1])(gml,J(Nobs,1,S.enter)))
		pred[,2] = t :- S.enter :- pred[,1]
		
		if (S.standardise) 	S.los = S.los :+ pred
		else 				S.los = pred
	}
	
	if (S.getrmst) {
		if (S.standardise) 	S.rmst = S.rmst :+ pred[,1]
		else 			S.rmst = pred[,1]
	}
	
	if (S.hasuser) predictms_calc_user(S)
		
}

void predictms_analytic_s_stand(`SS' S, `Pcm' Pmerlin)
{
	`gml' gml
	`RC' t
	`RS' Nt
	
	gml = *Pmerlin[1]
	t = S.predtime
	Nt = rows(t)
	
	if (S.getprobs) {	
		for (i=1;i<=Nt;i++) {
			S.pt[i,1] = 	///
			    predictms_analytic_s_p_stand(S,gml,t[i])
		}
		S.pt[,2] = 1 :- S.pt[,1]
	}
	
	if (S.getlos | S.getrmst) {
		Nqp = S.chips
		gq  = predictms_gq(Nqp)
		qp  = (t:-S.enter) :/ 2 :* J(Nt,1,gq[,1]') :+ (t:+S.enter) :/2
		res = J(Nt,Nqp,.)
		for (i=1;i<=Nt;i++) {
			for (j=1;j<=Nqp;j++) {
				res[i,j] = predictms_analytic_s_p_stand(S,
						gml,qp[i,j])
			}
		}
		S.los[,1]  = res * gq[,2] :* (t:-S.enter) :/ 2
		if (S.getlos) S.los[,2] = t :- S.enter :- S.los[,1]
	}
	
	if (S.getrmst) S.rmst = S.los[,1]
	if (S.hasuser) predictms_calc_user(S)		
}

`RS' predictms_analytic_s_p_stand(`SS' S, `gml' gml, `RS' t)
{
	ch	= (*gml.Pch[1])(gml,J(S.K,1,t))
	_editmissing(ch,0)
	pred = mean(exp(-ch))
	if (S.enter) {
		pred = pred :/ mean(exp(-(*gml.Pch[1])(gml,J(S.K,1,S.enter))))
	}
	return(pred)
}

end
