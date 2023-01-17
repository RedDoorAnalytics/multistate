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

void predictms_analytic_survival(`SS' S, `Pcm' Pmerlin, `RS' Nobs)
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
		else 				S.pt = pred
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
		else 				S.rmst = pred[,1]
	}
	
	if (S.hasuser) predictms_calc_user(S)
		
}

end
