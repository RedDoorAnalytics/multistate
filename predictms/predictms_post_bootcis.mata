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

void predictms_post_cis(`SS' S, `RS' fr, `RS' predtype)
{
	`RS' Nvars
	Nvars = S.Nstates
	if (!S.percentiles) {
		`RS' critval
		critval = invnormal((1-(100-S.level)/200))
	}

	for (at=1;at<=S.Nats;at++) {
		
		if (predtype==0) {
			pest	= asarray(S.probs,(fr,at,1))
			preds 	= asarray(S.probs,(fr,at,2))
			stub 	= "_prob_at"+strofreal(at)+"_"+strofreal(fr)
		}
		else if (predtype==1) {
			pest 	= asarray(S.loss,(fr,at,1))
			preds 	= asarray(S.loss,(fr,at,2))
			stub 	= "_los_at"+strofreal(at)+"_"+strofreal(fr)
		}
		else if (predtype==2) {
			pest 	= asarray(S.users,(fr,at,1))
			preds 	= asarray(S.users,(fr,at,2))
			stub 	= "_user_at"+strofreal(at)+"_"+strofreal(fr)
			Nvars	= S.Nuservars
		}
		else if (predtype==3) {
			pest 	= asarray(S.visits,(fr,at,1))
			preds 	= asarray(S.visits,(fr,at,2))
			stub 	= "_visit_at"+strofreal(at)+"_"+strofreal(fr)
		}
		else if (predtype==6) {
			pest 	= asarray(S.rmsts,(fr,at,1))
			preds 	= asarray(S.rmsts,(fr,at,2))
			name 	= "_rmst_at"+strofreal(at)+"_"+strofreal(fr)
			
			stata("cap drop "+name+"_*")
			names = name+"_lci",name+"_uci"
			(void) st_addvar("double",names)	
			
			base 		= asarray(preds,1)'
			basepest 	= pest

			if (S.percentiles) predictms_ci_perc(S,base,name)
			else {
				postpreds 	= J(S.obs,2,.)
				predictms_cis_ident(S,base,basepest,critval,postpreds,0,.)	//rmst
				st_store(.,names,S.touse,postpreds)
			}
			
			continue
		}
		else {
			nextstates 	= asarray(S.posnextstates,fr)
			Nvars		= rows(nextstates)
			if (predtype==4) {
				pest 		= asarray(S.hazards,(fr,at,1))
				preds 		= asarray(S.hazards,(fr,at,2))
				stub 		= "_hazard_at"+strofreal(at)+"_"+strofreal(fr)
			}
			else {
				pest 		= asarray(S.survivals,(fr,at,1))
				preds 		= asarray(S.survivals,(fr,at,2))
				stub 		= "_survival_at"+strofreal(at)+"_"+strofreal(fr)
			}
			for (tr=1;tr<=Nvars;tr++) {
			
				name = stub+"_"+strofreal(nextstates[tr])
				stata("cap drop "+name+"_*")
				names = name+"_lci",name+"_uci"
				(void) st_addvar("double",names)	
				
				base 		= asarray(preds,tr)'
				basepest 	= pest[,tr]

				if (S.percentiles) 			predictms_ci_perc(S,base,name)
				else {
					postpreds 	= J(S.obs,2,.)
					if (predtype==4) 	predictms_cis_log(S,base,basepest,critval,postpreds)
					else 				predictms_cis_logit(S,base,basepest,critval,postpreds)
					st_store(.,names,S.touse,postpreds)
				}
				
			}
			
			continue
		}
		
		for (tr=1;tr<=Nvars;tr++) {
			
			name = stub+"_"+strofreal(tr)
			stata("cap drop "+name+"_*")
			names = name+"_lci",name+"_uci"
			(void) st_addvar("double",names)	
			
			base 		= asarray(preds,tr)'
			basepest 	= pest[,tr]

			if (S.percentiles) 			predictms_ci_perc(S,base,name)
			else {
				postpreds 	= J(S.obs,2,.)
				if (predtype==0) 		predictms_cis_logit(S,base,basepest,critval,postpreds)	//probs
				else if (predtype==1) 	predictms_cis_log(S,base,basepest,critval,postpreds)	//los
				else if (predtype==2)	predictms_cis_user(S,base,basepest,critval,postpreds)	//user
				else 					predictms_cis_logit(S,base,basepest,critval,postpreds)	//probs visit
				st_store(.,names,S.touse,postpreds)
			}
			
		}
	}
}

void predictms_contrast_post_cis(`SS' S, `RS' fr, `RS' predtype, `RS' diff)
{

	`RS' Nvars
	`ss' cont, stub, name
	`TR' preds
	`RM' base, tpred
	
	if (diff) 	cont = "_diff"
	else 		cont = "_ratio"
	
	Nvars = S.Nstates
	
	if (!S.percentiles) {
		`RS' critval
		critval = invnormal((1-(100-S.level)/200))
	}
	
	if (predtype==0) {
		refpests	= asarray(S.probs,(fr,S.atref,1))
		refpreds 	= asarray(S.probs,(fr,S.atref,2))
		stub 		= cont+"_prob_at"
	}
	else if (predtype==1) {
		refpests	= asarray(S.loss,(fr,S.atref,1))
		refpreds	= asarray(S.loss,(fr,S.atref,2))
		stub 		= cont+"_los_at"
	}
	else if (predtype==2) {
		refpests 	= asarray(S.users,(fr,S.atref,1))
		refpreds 	= asarray(S.users,(fr,S.atref,2))
		stub 		= cont+"_user_at"
		Nvars 		= S.Nuservars
	}
	else if (predtype==3) {
		refpests 	= asarray(S.visits,(fr,S.atref,1))
		refpreds 	= asarray(S.visits,(fr,S.atref,2))
		stub 		= cont+"_visit_at"
	}	
	else if (predtype==6) {
		refpests 	= asarray(S.rmsts,(fr,S.atref,1))
		refpreds 	= asarray(S.rmsts,(fr,S.atref,2))
		Nvars		= 1
		stub 		= cont+"_rmst_at"
	}	
	else {
		if (predtype==4) {
			refpests 	= asarray(S.hazards,(fr,S.atref,1))
			refpreds 	= asarray(S.hazards,(fr,S.atref,2))
			stub 		= cont+"_hazard_at"
		}
		else {
			refpests 	= asarray(S.survivals,(fr,S.atref,1))
			refpreds 	= asarray(S.survivals,(fr,S.atref,2))
			stub 		= cont+"_survival_at"
		}
		nextstates 	= asarray(S.posnextstates,fr)
		Nvars		= rows(nextstates)
	}
	
	for (at=1;at<=S.Nats;at++) {
		
		if (at!=S.atref) {
		
			if (predtype==0) 		pests 	= asarray(S.probs,(fr,at,1))
			else if (predtype==1) 	pests 	= asarray(S.loss,(fr,at,1))
			else if (predtype==2)	pests 	= asarray(S.users,(fr,at,1))
			else if (predtype==3)	pests 	= asarray(S.visits,(fr,at,1))
			else if (predtype==4)	pests 	= asarray(S.hazards,(fr,at,1))
			else if (predtype==5)	pests 	= asarray(S.survivals,(fr,at,1))
			else 					pests 	= asarray(S.rmsts,(fr,at,1))
			if (predtype==0) 		preds 	= asarray(S.probs,(fr,at,2))
			else if (predtype==1) 	preds 	= asarray(S.loss,(fr,at,2))
			else if (predtype==2)	preds 	= asarray(S.users,(fr,at,2))
			else if (predtype==3)	preds 	= asarray(S.visits,(fr,at,2))
			else if (predtype==4)	preds 	= asarray(S.hazards,(fr,at,2))
			else if (predtype==5)	preds 	= asarray(S.survivals,(fr,at,2))
			else 					preds 	= asarray(S.rmsts,(fr,at,2))
			
			for (tr=1;tr<=Nvars;tr++) {
				
				if (predtype==4 | predtype==5) {
					name = stub+strofreal(at)+"_"+strofreal(fr)+"_"+strofreal(nextstates[tr])
				}
				else if (predtype==6) {
					name = stub+strofreal(at)+"_"+strofreal(fr)
				}
				else name = stub+strofreal(at)+"_"+strofreal(fr)+"_"+strofreal(tr)
				stata("cap drop "+name+"_*")
				names = name+"_lci",name+"_uci"
				(void) st_addvar("double",names)	
				
				if (diff) 	{
					basepest	= pests[,tr] :- refpests[,tr]
					base 		= (asarray(preds,tr) :- asarray(refpreds,tr))'
				}
				else {
					basepest	= pests[,tr] :/ refpests[,tr]
					base 		= (asarray(preds,tr) :/ asarray(refpreds,tr))'
				}
				
				if (S.percentiles) 	predictms_ci_perc(S,base,name)
				else {
					postpreds 	= J(S.obs,2,.)
					if (predtype==0 | predtype==3 | predtype==5) {									//probs, visit, survival
						if (diff) 	predictms_cis_atanh(S,base,basepest,critval,postpreds)
						else 		predictms_cis_log(S,base,basepest,critval,postpreds)
					}
					else if (predtype==1 | predtype==6) {											//los, rmst
						if (diff) 	predictms_cis_ident(S,base,basepest,critval,postpreds,.,.)	
						else 		predictms_cis_log(S,base,basepest,critval,postpreds)
					}
					else if (predtype==2) {												//user
						if (diff) 	predictms_cis_user_diff(S,base,basepest,critval,postpreds)	
						else 		predictms_cis_user_ratio(S,base,basepest,critval,postpreds)	
					}
					else {																//hazard
						if (diff) 	predictms_cis_ident(S,base,basepest,critval,postpreds,.,.)	
						else 		predictms_cis_log(S,base,basepest,critval,postpreds)	
					}
					st_store(.,names,S.touse,postpreds)
				}
			}
		}
	}
	
}

void predictms_ci_perc(`SS' S, `RM' base, `ss' name)
{
	tpred = J(S.obs,2,.)
	for (k=1;k<=S.obs;k++) {
		temp = sort(base[,k],1)
		tpred[k,1] = temp[S.lci,1]
		tpred[k,2] = temp[S.uci,1]
	}
	st_store(.,name+"_lci",S.touse,tpred[,1])
	st_store(.,name+"_uci",S.touse,tpred[,2])
}

void predictms_cis_ident(`SS' S, `RM' base, `RM' basepest, `RS' critval, `RM' postpreds, `RS' lim1, `RS' lim2)
{
	for (j=1;j<=S.obs;j++) {
		tpred = quadvariance(base[,j])'
		postpreds[j,1] = basepest[j] :- critval :* sqrt(tpred)
		postpreds[j,2] = basepest[j] :+ critval :* sqrt(tpred)
		if (lim1!=.) {
			if (postpreds[j,1]<lim1) postpreds[j,1] = lim1
		}
		if (lim2!=.) {
			if (postpreds[j,2]>lim2) postpreds[j,2] = lim2
		}
	}
}

void predictms_cis_log(`SS' S, `RM' base, `RM' basepest, `RS' critval, `RM' postpreds)
{
	for (j=1;j<=S.obs;j++) {
		if (min(base[,j])==0 ) {
			tpred = quadvariance(base[,j])'
			postpreds[j,1] = basepest[j] :- critval :* sqrt(tpred)
			postpreds[j,2] = basepest[j] :+ critval :* sqrt(tpred)
			//re-bound between 0 and +infty
			if (postpreds[j,1]<0) postpreds[j,1] = 0
		}
		else {
			tpred = quadvariance(log(base[,j]))'
			postpreds[j,1] = log(basepest[j]) :- critval :* sqrt(tpred)
			postpreds[j,2] = log(basepest[j]) :+ critval :* sqrt(tpred)
			postpreds[j,] = exp(postpreds[j,])
		}
	}
}

void predictms_cis_logit(`SS' S, `RM' base, `RM' basepest, `RS' critval, `RM' postpreds)
{
	for (j=1;j<=S.obs;j++) {
		if (min(base[,j])==0  | max(base[,j])==1) {
			tpred = quadvariance(base[,j])'
			postpreds[j,1] = basepest[j] :- critval :* sqrt(tpred)
			postpreds[j,2] = basepest[j] :+ critval :* sqrt(tpred)
			//re-bound between 0 and 1
			if (postpreds[j,1]<0) postpreds[j,1] = 0
			if (postpreds[j,2]>1) postpreds[j,2] = 1
		}
		else {
			tpred = quadvariance(logit(base[,j]))'
			postpreds[j,1] = logit(basepest[j]) :- critval :* sqrt(tpred)
			postpreds[j,2] = logit(basepest[j]) :+ critval :* sqrt(tpred)
			postpreds[j,] = invlogit(postpreds[j,])
		}
	}
}

void predictms_cis_atanh(`SS' S, `RM' base, `RM' basepest, `RS' critval, `RM' postpreds)
{
	for (j=1;j<=S.obs;j++) {
		if (min(base[,j])==-1  | max(base[,j])==1) {
			tpred = quadvariance(base[,j])'
			postpreds[j,1] = basepest[j] :- critval :* sqrt(tpred)
			postpreds[j,2] = basepest[j] :+ critval :* sqrt(tpred)
			//re-bound between -1 and 1
			if (postpreds[j,1]<-1) postpreds[j,1] = -1
			if (postpreds[j,2]>1)  postpreds[j,2] = 1
		}
		else {
			tpred = quadvariance(atanh(base[,j]))'
			postpreds[j,1] = atanh(basepest[j]) :- critval :* sqrt(tpred)
			postpreds[j,2] = atanh(basepest[j]) :+ critval :* sqrt(tpred)
			postpreds[j,] = tanh(postpreds[j,])
		}
	}
}

void predictms_cis_user(`SS' S, `RM' base, `RM' basepest, `RS' critval, `RM' postpreds)
{
	if (S.userlink==0)		predictms_cis_ident(S,base,basepest,critval,postpreds,.,.)
	else if (S.userlink==1) predictms_cis_logit(S,base,basepest,critval,postpreds)
	else if (S.userlink==2) predictms_cis_log(S,base,basepest,critval,postpreds)
}

void predictms_cis_user_diff(`SS' S, `RM' base, `RM' basepest, `RS' critval, `RM' postpreds)
{
	if (S.userlink==0)		predictms_cis_ident(S,base,basepest,critval,postpreds,.,.)
	else if (S.userlink==1) predictms_cis_atanh(S,base,basepest,critval,postpreds)
	else if (S.userlink==2) predictms_cis_ident(S,base,basepest,critval,postpreds)
}

void predictms_cis_user_ratio(`SS' S, `RM' base, `RM' basepest, `RS' critval, `RM' postpreds)
{
	if (S.userlink==0)		predictms_cis_ident(S,base,basepest,critval,postpreds,.,.)
	else if (S.userlink==1) predictms_cis_log(S,base,basepest,critval,postpreds)
	else if (S.userlink==2) predictms_cis_log(S,base,basepest,critval,postpreds)
}

end
