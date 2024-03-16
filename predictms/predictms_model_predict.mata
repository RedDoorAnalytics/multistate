
version 14.2

local ss 	string scalar
local RS	real scalar
local NM	numeric matrix
local RM 	real matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local PS	pointer(struct merlin_struct scalar) scalar
local SS	struct predictms_struct scalar
local gml 	struct merlin_struct scalar

mata:

void predictms_model_predict(`SS' S, `RS' from)
{
	`gml' 	gml

	Nobs 		= rows(S.predtime)
	postrans 	= asarray(S.postrans,from)
	Npostrans 	= rows(postrans)

	for (c=1; c<=Npostrans; c++) {
		trans 	= postrans[c]
		b	= predictms_get_b(S,trans)
		gml 	= *predictms_merlin_setup(S,b,Nobs,trans)
		if (S.gethazard) {
			S.hazard[,c] = exp((*gml.Plogh[1])(gml,S.predtime))
		}
		if (S.getsurvival) 	{
			survpred = exp(-(*gml.Pch[1])(gml,S.predtime))
			
			if (S.enter) {
				survpred = survpred :/ 	///
				    exp(-(*gml.Pch[1])(gml,J(Nobs,1,S.enter)))
			}
			
			if (min(S.predtime)==0) {
				index0 	= selectindex(S.predtime:==0)
				nind	= rows(index0)
				survpred[index0] = J(nind,1,1)
			}
			S.survival[,c] 	= S.survival[,c] :+ survpred
		}	
		rmexternal(st_local("GML"+strofreal(trans)))
	}
}

void predictms_model_predict_stand(`SS' S, `RS' from)
{
	`PS' 	gml
	`RC'	t
	`RS' 	Nt, Npostrans
	
	t 		= S.predtime
	Nt 		= rows(t)
	postrans 	= asarray(S.postrans,from)
	Npostrans 	= rows(postrans)
	for (c=1; c<=Npostrans; c++) {
		trans 	= postrans[c]
		b	= predictms_get_b(S,trans)
		Pmerlin	= predictms_merlin_setup_stand(S,b,Nt,trans)
		if (S.getsurvival) 	{
			for (i=1;i<=Nt;i++) {
				S.survival[i,c] = 	///
				    predictms_analytic_s_p_stand(S,Pmerlin,t[i])
			}	
			if (min(S.predtime)==0) {
				index0 	= selectindex(S.predtime:==0)
				nind	= rows(index0)
				S.survival[index0,c] = J(nind,1,1)
			}
		}
		if (S.gethazard) {
			for (i=1;i<=Nt;i++) {
				S.hazard[i,c] = 	///
				    predictms_analytic_h_p_stand(S,Pmerlin,t[i])
			}	
			
		}	
		rmexternal(st_local("GML"+strofreal(trans)))
	}
}

`RS' predictms_analytic_h_p_stand(`SS' S, `Pcm' Pmerlin, `RS' t)
{
	`gml' gml
	gml	= *Pmerlin[1]
	pred 	= J(1,2,.)
	ch	= (*gml.Pch[1])(gml,J(S.K,1,t))
	_editmissing(ch,0)
	pred = mean(exp((*gml.Plogh[1])(gml,J(S.K,1,t)) :- ch))
	if (S.enter) {
		pred = pred :/ mean(exp(-(*gml.Pch[1])(gml,J(S.K,1,S.enter))))
	}
	return(pred)
}

end

