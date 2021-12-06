/*
predictms Mata sourcecode
*/

version 14.2

local RM 	real matrix
local ss 	string scalar
local RS	real scalar
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local SS	struct predictms_struct scalar
local PS	pointer(struct predictms_struct scalar) scalar
local PC	pointer(struct predictms_struct scalar) colvector
local Ps	pointer scalar

mata:

void predictms_post_cis_dm(`SS' S, `RS' fr, `RS' predtype)
{
	`RS' Nvars, critval
	Nvars 	= S.Nstates
	critval = invnormal((1-(100-S.level)/200))

	for (at=1;at<=S.Nats;at++) {
		
		if (predtype==0) {
			pest	= asarray(S.probs,(fr,at,1))
			dbetas 	= asarray(S.probs,(fr,at,2))
			ses 	= se_dm_logit(pest,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			stub 	= "_prob_at"
			indexes	= 1::Nvars
		}
		else if (predtype==1) {
			pest	= asarray(S.loss,(fr,at,1))
			dbetas 	= asarray(S.loss,(fr,at,2))
			ses 	= se_dm_log(pest,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			stub 	= "_los_at"
			indexes	= 1::Nvars
		}
		else if (predtype==2) {
			Nvars	= S.Nuservars
			stub 	= "_user_at"
			indexes	= 1::Nvars
			pest 	= asarray(S.users,(fr,at,1))
			dbetas 	= asarray(S.users,(fr,at,2))
			if 		(S.userlink==1) ses 	= se_dm_logit(pest,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			else if (S.userlink==2) ses 	= se_dm_log(pest,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			else 					ses 	= se_dm_iden(pest,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
		}
		else if (predtype==6) {
			stub 	= "_rmst_at"
			Nvars	= 1
			indexes	= 1::Nvars
			pest 	= asarray(S.rmsts,(fr,at,1))
			dbetas 	= asarray(S.rmsts,(fr,at,2))
			ses 	= se_dm_log(pest,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
		}
		else {
			indexes 	= asarray(S.posnextstates,fr)
			Nvars		= rows(indexes)
			if (predtype==4) {
				pest 	= asarray(S.hazards,(fr,at,1))
				dbetas 	= asarray(S.hazards,(fr,at,2))
				ses 	= se_dm_log(pest,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				stub 	= "_hazard_at"
			}
			else {
				pest 	= asarray(S.survivals,(fr,at,1))
				dbetas 	= asarray(S.survivals,(fr,at,2))
				ses 	= se_dm_logit(pest,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				stub 	= "_survival_at"
			}
		}
		
		stub = stub + strofreal(at) + "_" + strofreal(fr)
		
		for (tr=1;tr<=Nvars;tr++) {
			
			name = stub+"_"+strofreal(indexes[tr])
			
			stata("cap drop "+name+"_*")
			names = name+"_lci",name+"_uci"
			(void) st_addvar("double",names)
			postpreds 	= J(S.obs,2,.)
			
			if 		(predtype==0 | predtype==5 | (predtype==2 & S.userlink==1)) 	{
				
				postpreds[,1] 	= logit(pest[,tr]) :- critval :* ses[,tr]
				postpreds[,2] 	= logit(pest[,tr]) :+ critval :* ses[,tr]
				postpreds 		= invlogit(postpreds)
				
				index1			= selectindex(pest[,tr]:==1)				//1 - boundary
				nobs 			= rows_cols(index1)
				if (nobs[1] & nobs[2]) 	postpreds[index1,] = J(nobs[1],2,1)
				index1			= selectindex(pest[,tr]:==0)				//0 - boundary
				nobs 			= rows_cols(index1)
				if (nobs[1] & nobs[2]) 	postpreds[index1,] = J(nobs[1],2,0)
				
			}
			else if (predtype==1 | predtype==6 | predtype==4 | (predtype==2 & S.userlink==2))	{
				
				postpreds[,1] 	= log(pest[,tr]) :- critval :* ses[,tr]
				postpreds[,2] 	= log(pest[,tr]) :+ critval :* ses[,tr]
				postpreds 		= exp(postpreds)
				
				index1			= selectindex(pest[,tr]:==0)				//0 - boundary
				nobs 			= rows_cols(index1)
				if (nobs[1] & nobs[2]) 	postpreds[index1,] = J(nobs[1],2,0)
				
			}
			else {
				postpreds[,1] 	= pest[,tr] :- critval :* ses[,tr]
				postpreds[,2] 	= pest[,tr] :+ critval :* ses[,tr]
			}
			
			st_store(.,names,S.touse,postpreds)
		}
		
	}
}

void predictms_contrast_post_cis_dm(`SS' S, `RS' fr, `RS' predtype, `RS' diff)
{

	`RS' Nvars, critval
	`ss' cont, stub, name
	`TR' preds
	`RM' base, tpred
	
	if (diff) 	cont = "_diff"
	else 		cont = "_ratio"
	
	Nvars 	= S.Nstates
	indexes = 1::Nvars
	critval = invnormal((1-(100-S.level)/200))
	
	if (predtype==0) {
		refpests 	= asarray(S.probs,(fr,S.atref,1))
		refdbetas 	= asarray(S.probs,(fr,S.atref,2))
		stub 		= cont+"_prob_at"
	}
	else if (predtype==1) {
		refpests	= asarray(S.loss,(fr,S.atref,1))
		refdbetas	= asarray(S.loss,(fr,S.atref,2))
		stub 		= cont+"_los_at"
	}
	else if (predtype==2) {
		refpests 	= asarray(S.users,(fr,S.atref,1))
		refdbetas 	= asarray(S.users,(fr,S.atref,2))
		stub 		= cont+"_user_at"
		Nvars 		= S.Nuservars
	}
	else if (predtype==3) {
		refpests 	= asarray(S.visits,(fr,S.atref,1))
		refdbetas 	= asarray(S.visits,(fr,S.atref,2))
		stub 		= cont+"_visit_at"
	}	
	else if (predtype==6) {
		refpests 	= asarray(S.rmsts,(fr,S.atref,1))
		refdbetas 	= asarray(S.rmsts,(fr,S.atref,2))
		Nvars		= 1
		stub 		= cont+"_rmst_at"
	}	
	else {
		if (predtype==4) {
			refpests 	= asarray(S.hazards,(fr,S.atref,1))
			refdbetas 	= asarray(S.hazards,(fr,S.atref,2))
			stub 		= cont+"_hazard_at"
		}
		else {
			refpests 	= asarray(S.survivals,(fr,S.atref,1))
			redbetas 	= asarray(S.survivals,(fr,S.atref,2))
			stub 		= cont+"_survival_at"
		}
		indexes 	= asarray(S.posnextstates,fr)
		Nvars		= rows(indexes)
	}

	for (at=1;at<=S.Nats;at++) {
		
		if (at!=S.atref) {
		
			if (predtype==0) {
				pests 	= asarray(S.probs,(fr,at,1))
				dbetas 	= asarray(S.probs,(fr,at,2))
				if (diff) 	ses = se_dm_atanh_diff(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				else 		ses = se_dm_log_ratio(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			}
			else if (predtype==1) {
				pests 	= asarray(S.loss,(fr,at,1))
				dbetas 	= asarray(S.loss,(fr,at,2))
				if (diff) 	ses = se_dm_diff(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				else 		ses = se_dm_log_ratio(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			}
			else if (predtype==2) {
				pests 	= asarray(S.users,(fr,at,1))
				dbetas 	= asarray(S.users,(fr,at,2))
				if (diff) 	{
					if 	 (S.userlink==1) 	ses = se_dm_atanh_diff(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
					else 					ses = se_dm_diff(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				}
				else {
					if 	 (S.userlink==1) 	ses = se_dm_log_ratio(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
					else if (S.userlink==2)	ses = se_dm_log_ratio(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
					else					ses = se_dm_ratio(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				}
			}
			else if (predtype==3) {
				pests 	= asarray(S.visits,(fr,at,1))
				dbetas 	= asarray(S.visits,(fr,at,2))
				if (diff) 	ses = se_dm_atanh_diff(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				else 		ses = se_dm_log_ratio(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			}
			else if (predtype==4) {
				pests 	= asarray(S.hazards,(fr,at,1))
				dbetas 	= asarray(S.hazards,(fr,at,2))
				if (diff) 	ses = se_dm_diff(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				else 		ses = se_dm_log_ratio(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			}
			else if (predtype==5) {
				pests 	= asarray(S.survivals,(fr,at,1))
				dbetas 	= asarray(S.survivals,(fr,at,2))
				if (diff) 	ses = se_dm_atanh_diff(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				else 		ses = se_dm_log_ratio(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			}
			else {
				pests 	= asarray(S.rmsts,(fr,at,1))
				dbetas 	= asarray(S.rmsts,(fr,at,2))
				if (diff) 	ses = se_dm_diff(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
				else 		ses = se_dm_log_ratio(refpests,pests,refdbetas,dbetas,Nvars,S.obs,S.Nbs,S.VCV)
			}

			for (tr=1;tr<=Nvars;tr++) {
				
				name = stub+strofreal(at)+"_"+strofreal(fr)+"_"+strofreal(indexes[tr])
				stata("cap drop "+name+"_*")
				names = name+"_lci",name+"_uci"
				(void) st_addvar("double",names)	
				
				postpreds 	= J(S.obs,2,.)
				
				if (diff) 	{
					if (predtype==0 | predtype==3 | predtype==5) {				
						postpreds[,1] 	= atanh(pests[,tr] :- refpests[,tr]) :- critval :* ses[,tr]
						postpreds[,2] 	= atanh(pests[,tr] :- refpests[,tr]) :+ critval :* ses[,tr]
						postpreds 		= tanh(postpreds)
					}
					else {
						postpreds[,1] 	= pests[,tr] :- refpests[,tr] :- critval :* ses[,tr]
						postpreds[,2] 	= pests[,tr] :- refpests[,tr] :+ critval :* ses[,tr]
					}
				}
				else {
					postpreds[,1] 	= log(pests[,tr] :/ refpests[,tr]) :- critval :* ses[,tr]
					postpreds[,2] 	= log(pests[,tr] :/ refpests[,tr]) :+ critval :* ses[,tr]
					postpreds 		= exp(postpreds)
				}

				st_store(.,names,S.touse,postpreds)
			}
		}
	}
	
}

`RM' se_dm_logit(`RM' pest, `TR' A, `RS' Nvars, `RS' Nobs, `RS' Nbs, `RM' VCV)
{
	ses = J(Nobs,Nvars,0)
	for (i=1;i<=Nvars;i++) {
		Gb = J(Nobs,Nbs,.)
		for (b=1;b<=Nbs;b++) Gb[,b] = asarray(A,b)[,i]
		Gb = Gb :* (1:/pest[,i] :+ 1:/(1:-pest[,i]))
		ses[,i] = diagonal(Gb * VCV * Gb')
	}
	return(sqrt(ses))
}			

`RM' se_dm_log(`RM' pest, `TR' A, `RS' Nvars, `RS' Nobs, `RS' Nbs, `RM' VCV)
{
	ses = J(Nobs,Nvars,0)
	for (i=1;i<=Nvars;i++) {
		Gb = J(Nobs,Nbs,.)
		for (b=1;b<=Nbs;b++) Gb[,b] = asarray(A,b)[,i]
		Gb = Gb :/ pest[,i]
		ses[,i] = diagonal(Gb * VCV * Gb')
	}
	return(sqrt(ses))
}

`RM' se_dm_iden(`RM' pest, `TR' A, `RS' Nvars, `RS' Nobs, `RS' Nbs, `RM' VCV)
{
	ses = J(Nobs,Nvars,0)
	for (i=1;i<=Nvars;i++) {
		Gb = J(Nobs,Nbs,.)
		for (b=1;b<=Nbs;b++) Gb[,b] = asarray(A,b)[,i]
		ses[,i] = diagonal(Gb * VCV * Gb')
	}
	return(sqrt(ses))
}

`RM' se_dm_atanh_diff(`RM' pest1, `RM' pest2, `TR' A, `TR' B, `RS' Nvars, `RS' Nobs, `RS' Nbs, `RM' VCV)
{
	ses = J(Nobs,Nvars,0)
	for (i=1;i<=Nvars;i++) {
		Gb =  J(Nobs,Nbs,.)
		for (b=1;b<=Nbs;b++) Gb[,b] = asarray(B,b)[,i] :- asarray(A,b)[,i]
		Gb = Gb :/ (cosh(pest2[,i] :- pest1[,i]):^2)
		ses[,i] = diagonal(Gb * VCV * Gb')
	}
	return(sqrt(ses))
}

`RM' se_dm_log_ratio(`RM' pest1, `RM' pest2, `TR' A, `TR' B, `RS' Nvars, `RS' Nobs, `RS' Nbs, `RM' VCV)
{
	ses = J(Nobs,Nvars,0)
	for (i=1;i<=Nvars;i++) {
		Gb =  J(Nobs,Nbs,.)
		for (b=1;b<=Nbs;b++) Gb[,b] = asarray(B,b)[,i]:/pest2[,i] :- asarray(A,b)[,i]:/pest1[,i]
		ses[,i] = diagonal(Gb * VCV * Gb')
	}
	return(sqrt(ses))
}

`RM' se_dm_diff(`RM' pest1, `RM' pest2, `TR' A, `TR' B, `RS' Nvars, `RS' Nobs, `RS' Nbs, `RM' VCV)
{
	ses = J(Nobs,Nvars,0)
	for (i=1;i<=Nvars;i++) {
		Gb =  J(Nobs,Nbs,.)
		for (b=1;b<=Nbs;b++) Gb[,b] = asarray(B,b)[,i] :- asarray(A,b)[,i]
		ses[,i] = diagonal(Gb * VCV * Gb')
	}
	return(sqrt(ses))
}

`RM' se_dm_ratio(`RM' pest1, `RM' pest2, `TR' A, `TR' B, `RS' Nvars, `RS' Nobs, `RS' Nbs, `RM' VCV)
{
	ses = J(Nobs,Nvars,0)
	for (i=1;i<=Nvars;i++) {
		Gb =  J(Nobs,Nbs,.)
		for (b=1;b<=Nbs;b++) Gb[,b] = asarray(B,b)[,i]  :- asarray(A,b)[,i] :* pest2[,i] :/ pest1[,i]:^2
		ses[,i] = diagonal(Gb * VCV * Gb')
	}
	return(sqrt(ses))
}

end
