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

void predictms_post_predictions(`SS' S, `RS' fr, `RS' predtype)
{	
								predictms_post(S,fr,predtype)
	if (S.getdiffs)				predictms_contrast_post(S,fr,predtype,1)
	if (S.getratios) 			predictms_contrast_post(S,fr,predtype,0)

	if (S.getcis) {
		if (!S.cimethod) {		//DM
								predictms_post_cis_dm(S,fr,predtype)
			if (S.getdiffs)		predictms_contrast_post_cis_dm(S,fr,predtype,1)
			if (S.getratios) 	predictms_contrast_post_cis_dm(S,fr,predtype,0)
		}
		else {					//Bootstrap
								predictms_post_cis(S,fr,predtype)
			if (S.getdiffs)		predictms_contrast_post_cis(S,fr,predtype,1)
			if (S.getratios) 	predictms_contrast_post_cis(S,fr,predtype,0)
		}
	}
}

void predictms_post(`SS' S, `RS' fr, `RS' predtype)
{
	Nvars = S.Nstates

	for (at=1;at<=S.Nats;at++) {
		
		if (predtype==0) {
			preds 	= asarray(S.probs,(fr,at,1))
			stub 	= "_prob_at"+strofreal(at)+"_"+strofreal(fr)
		}
		else if (predtype==1) {
			preds 	= asarray(S.loss,(fr,at,1))
			stub 	= "_los_at"+strofreal(at)+"_"+strofreal(fr)
		}
		else if (predtype==2) {
			preds 	= asarray(S.users,(fr,at,1))
			stub 	= "_user_at"+strofreal(at)+"_"+strofreal(fr)
			Nvars 	= S.Nuservars
		}
		else if (predtype==3) {
			preds 	= asarray(S.visits,(fr,at,1))
			stub 	= "_visit_at"+strofreal(at)+"_"+strofreal(fr)
		}
		else if (predtype==6) {
			preds 	= asarray(S.rmsts,(fr,at,1))
			name 	= "_rmst_at"+strofreal(at)+"_"+strofreal(fr)
			stata("cap drop "+name+"*")
			(void) st_addvar("double",name)	
			st_store(.,name,S.touse, preds)
			continue
		}
		else {
			nextstates = asarray(S.posnextstates,fr)
			Nvars	= rows(nextstates)
			if (predtype==4) {
				preds 	= asarray(S.hazards,(fr,at,1))
				stub 	= "_hazard_at"+strofreal(at)+"_"+strofreal(fr)
			}
			else {
				preds 	= asarray(S.survivals,(fr,at,1))
				stub 	= "_survival_at"+strofreal(at)+"_"+strofreal(fr)
			}
			stata("cap drop "+stub+"*")
			for (tr=1;tr<=Nvars;tr++) {
				name = stub+"_"+strofreal(nextstates[tr])
				(void) st_addvar("double",name)	
				st_store(.,name,S.touse, preds[,tr])
			}
			continue
		}
		
		stata("cap drop "+stub+"*")
		for (tr=1;tr<=Nvars;tr++) {
			name = stub+"_"+strofreal(tr)
			(void) st_addvar("double",name)	
			st_store(.,name,S.touse, preds[,tr])
		}
	}
	
}

void predictms_contrast_post(`SS' S, `RS' fr, `RS' predtype, `RS' diff)
{
	`RS' Nvars
	`ss' cont, stub, name
	`RC' preds
	
	if (diff) 	cont = "_diff"
	else 		cont = "_ratio"
	
	Nvars = S.Nstates
	if (predtype==0) {
		refpreds 	= asarray(S.probs,(fr,S.atref,1))
		stub = cont+"_prob_at"
	}
	else if (predtype==1) {
		refpreds	= asarray(S.loss,(fr,S.atref,1))
		stub = cont+"_los_at"
	}
	else if (predtype==2) {
		refpreds 	= asarray(S.users,(fr,S.atref,1))
		stub 		= cont+"_user_at"
		Nvars 		= S.Nuservars
	}
	else if (predtype==3) {
		refpreds	= asarray(S.visits,(fr,S.atref,1))
		stub = cont+"_visit_at"
	}
	else if (predtype==6) {
		refpreds	= asarray(S.rmsts,(fr,S.atref,1))
		stub = cont+"_rmst_at"
	}
	else {
		if (predtype==4) {
			refpreds 	= asarray(S.hazards,(fr,S.atref,1))
			stub 		= cont+"_hazard_at"
		}
		else {
			refpreds 	= asarray(S.survivals,(fr,S.atref,1))
			stub 		= cont+"_survival_at"
		}
		nextstates = asarray(S.posnextstates,fr)
		Nvars	= rows(nextstates)
	}
	
	for (at=1;at<=S.Nats;at++) {
		
		if (at!=S.atref) {
		
			name = stub+strofreal(at)+"_"+strofreal(fr)
			stata("cap drop "+name+"*")
			
			if (predtype==0) {
				if (diff) 	preds 	= asarray(S.probs,(fr,at,1)) :- refpreds
				else 		preds 	= asarray(S.probs,(fr,at,1)) :/ refpreds
			}
			else if (predtype==1) {
				if (diff)	preds 	= asarray(S.loss,(fr,at,1))  :- refpreds
				else 		preds 	= asarray(S.loss,(fr,at,1))  :/ refpreds
			}
			else if (predtype==2) {
				if (diff) 	preds 	= asarray(S.users,(fr,at,1)) :- refpreds
				else 		preds 	= asarray(S.users,(fr,at,1)) :/ refpreds
			}
			else if (predtype==3) {
				if (diff) 	preds 	= asarray(S.visits,(fr,at,1)) :- refpreds
				else 		preds 	= asarray(S.visits,(fr,at,1)) :/ refpreds
			}
			else if (predtype==6) {
				if (diff) 	preds 	= asarray(S.rmsts,(fr,at,1)) :- refpreds
				else 		preds 	= asarray(S.rmsts,(fr,at,1)) :/ refpreds
				(void) st_addvar("double",name)	
				st_store(.,name,S.touse, preds)
				continue
			}
			else {
				if (predtype==4) {
					if (diff) 	preds 	= asarray(S.hazards,(fr,at,1)) :- refpreds
					else 		preds 	= asarray(S.hazards,(fr,at,1)) :/ refpreds
				}
				else {
					if (diff) 	preds 	= asarray(S.survivals,(fr,at,1)) :- refpreds
					else 		preds 	= asarray(S.survivals,(fr,at,1)) :/ refpreds
				}
				for (tr=1;tr<=Nvars;tr++) {
					name2 = name+"_"+strofreal(nextstates[tr])
					(void) st_addvar("double",name2)	
					st_store(.,name2,S.touse, preds[,tr])					
				}
				continue
			}
			
			for (tr=1;tr<=Nvars;tr++) {
				
				name2 = name+"_"+strofreal(tr)
				(void) st_addvar("double",name2)	
				st_store(.,name2,S.touse, preds[,tr])
				
			}
		}
	}
	
}

end
