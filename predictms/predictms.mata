
version 14.2

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

mata:

void predictms()
{
	`SS' S
	`RC' from
	`RS' Nstarts

	predictms_setup(S)

	from 	= strtoreal(tokens(st_local("from")))'
	Nstarts = rows(from)

	for (fr=1;fr<=Nstarts;fr++) {

		for (at=1;at<=S.Nats;at++) {
			S.at 	= at
			S.Kind 	= 0
			S.std 	= 0
			
			predictms_core(S,from[fr])
		}

		if (S.getprobs)  	predictms_post_predictions(S,from[fr],0)
		if (S.getlos) 	 	predictms_post_predictions(S,from[fr],1)
		if (S.hasuser) 	 	predictms_post_predictions(S,from[fr],2)
		if (S.getvisit)  	predictms_post_predictions(S,from[fr],3)
		if (S.gethazard) 	predictms_post_predictions(S,from[fr],4)
		if (S.getsurvival) 	predictms_post_predictions(S,from[fr],5)
		if (S.getrmst) 		predictms_post_predictions(S,from[fr],6)

	}
	
}

void predictms_core(`SS' S, `RS' from)
{
	
	//========================================================================================================================//
	//point estimates
	
		ptlosvisit = S.getprobs | S.getlos | S.getvisit | S.getrmst
		
		if (S.getprobs)	 	S.pt[,] 		= J(S.obs,S.Nstates,0)											
		if (S.getlos) 	 	S.los[,] 		= J(S.obs,S.Nstates,0)
		if (S.getrmst) 	 	S.rmst[,] 		= J(S.obs,1,0)
		if (S.getvisit)  	S.visit[,] 		= J(S.obs,S.Nstates,0)
		if (S.gethazard) 	S.hazard	 	= J(S.obs,S.Nnextstates[from],0)
		if (S.getsurvival)	S.survival	 	= J(S.obs,S.Nnextstates[from],0)

		//std loop
		for (std=1;std<=S.K;std++) {
			
			S.std = std
			
			if (ptlosvisit) {
				if 		(S.method==0 | S.method==3)	predictms_sim(S,from)		
				else if (S.method==1)				predictms_aj(S,from)		
				else if (S.method==2)				predictms_analytic(S,from)	
			}
			
			if (S.gethazard | S.getsurvival) {
				predictms_model_predict(S,from)
			}

		}

		if (S.standardise) {
			if (S.getprobs)	 	S.pt 		= S.pt :/ S.K
			if (S.getlos) 	 	S.los 	 	= S.los :/ S.K
			if (S.getrmst) 	 	S.rmst 	 	= S.rmst :/ S.K
			if (S.getvisit)  	S.visit 	= S.visit :/ S.K
			if (S.hasuser) 	 	S.user 		= S.user :/ S.K
			if (S.getsurvival) 	S.survival 	= S.survival :/ S.K
		}
		
		//store
		if (S.getprobs) 	{
			if (S.method==0 | S.method==3) S.pt = S.pt :/ rowsum(S.pt)
			asarray(S.probs,(from,S.at,1),S.pt)
		}
		if (S.getlos) 		asarray(S.loss,(from,S.at,1),S.los)
		if (S.getrmst) 		asarray(S.rmsts,(from,S.at,1),S.rmst)
		if (S.getvisit) 	asarray(S.visits,(from,S.at,1),S.visit)
		if (S.hasuser) 		asarray(S.users,(from,S.at,1),S.user)
		if (S.gethazard) 	asarray(S.hazards,(from,S.at,1),S.hazard)
		if (S.getsurvival) 	asarray(S.survivals,(from,S.at,1),S.survival)
		
	//========================================================================================================================//
	//confidence intervals

	if (S.getcis) {
				
		//to hold predictions 
		if (S.getprobs)		A = asarray_create("real",1)	
		if (S.getlos) 		B = asarray_create("real",1)
		if (S.getrmst) 		C = asarray_create("real",1)
		if (S.hasuser) 	 	D = asarray_create("real",1)
		if (S.getvisit)  	E = asarray_create("real",1)
		if (S.gethazard) 	F = asarray_create("real",1)
		if (S.getsurvival) 	G = asarray_create("real",1)
		
		if (!S.cimethod) {
			
			//get global parameter vector & VCV
			if (S.hasmodels) {
				Best 	= J(0,1,.)
				VCV 	= J(0,0,.)
				Bindex 	= J(S.Ntrans,2,.)
				i1 		=  1
				for (i=1;i<=S.Ntrans;i++) {
					Bi			= asarray(S.transinfo,(i,1))'
					NBi			= rows(Bi)
					Bindex[i,1] = i1
					Bindex[i,2] = i1 + NBi - 1
					Best		= Best\Bi
					if (S.novcv[i]) VCV = blockdiag(VCV,J(NBi,NBi,0))
					else			VCV = blockdiag(VCV,asarray(S.transinfo,(i,2)))
					i1 			= Bindex[i,2] + 1
				}
			}
			else {
				Best 	= asarray(S.transinfo,1)'
				VCV 	= asarray(S.transinfo,2)
			}
			S.VCV 	= VCV
			S.Nbs 	= rows(Best)
			
			printf("\n")
			stata(`"_dots 0 , title("Calc. CIs using the delta method for at"'+strofreal(S.at)+`"()") reps("'+strofreal(S.Nbs)+")")
			
			//now derivative with respect to each parameter
			for (b=1;b<=S.Nbs;b++) {

				newB	= Best
				predictms_init_storage(S,from)
				
				if (abs(Best[b])<1) hstep = abs(Best[b])
				else 				hstep = 1
				hstep = hstep * c("epsdouble") :^ (1/3)
				newB[b] = Best[b] + hstep/2
				
				//left update
				
					if (S.hasmodels) {
						for (i=1;i<=S.Ntrans;i++) {
							if (!S.novcv[i]) asarray(S.transinfo,(i,1),newB[|Bindex[i,1],1\Bindex[i,2],1|]')
						}
					}
					else asarray(S.transinfo,1,newB')
					
					for (std=1;std<=S.K;std++) {
						
						S.std = std
						if (ptlosvisit) {
							if 		(S.method==0 | S.method==3)	predictms_sim(S,from)
							else if (S.method==1)				predictms_aj(S,from)	
							else if (S.method==2)				predictms_analytic(S,from)
						}
						if (S.gethazard | S.getsurvival) 		predictms_model_predict(S,from)
						
					}

					if (S.standardise) {
						if (S.getprobs) 	S.pt 		= S.pt 			:/ S.K
						if (S.getlos) 	 	S.los   	= S.los   		:/ S.K
						if (S.getvisit)  	S.visit 	= S.visit  		:/ S.K
						if (S.hasuser) 	 	S.user  	= S.user   		:/ S.K
						if (S.getsurvival) 	S.survival 	= S.survival 	:/ S.K
					}

					//scale probs to total of 1
					if (S.getprobs & (S.method==0 | S.method==3)) 	S.pt = S.pt :/ rowsum(S.pt)

					//store lhs
					if (S.getprobs) 	pt 			= S.pt //logit(S.pt)
					if (S.getlos) 	 	los   		= S.los
					if (S.getrmst)  	rmst 		= S.rmst
					if (S.getvisit)  	visit 		= S.visit
					if (S.hasuser) 		user 		= S.user
					if (S.gethazard) 	hazard		= S.hazard
					if (S.getsurvival) 	survival	= S.survival

				//right update
				
					newB[b] = Best[b] - hstep/2
					
					if (S.hasmodels) {
						for (i=1;i<=S.Ntrans;i++) {
							if (!S.novcv[i]) asarray(S.transinfo,(i,1),newB[|Bindex[i,1],1\Bindex[i,2],1|]')
						}
					}
					else asarray(S.transinfo,1,newB')

					//reset
					predictms_init_storage(S,from)

					//1 or std loop
					for (std=1;std<=S.K;std++) {
						
						S.std = std
						if (ptlosvisit) {
							if 		(S.method==0 | S.method==3)	predictms_sim(S,from)
							else if (S.method==1)				predictms_aj(S,from)	
							else if (S.method==2)				predictms_analytic(S,from)
						}
						if (S.gethazard | S.getsurvival) 		predictms_model_predict(S,from)
						
					}

					if (S.standardise) {
						if (S.getprobs) 	S.pt 		= S.pt 			:/ S.K
						if (S.getlos) 	 	S.los   	= S.los    		:/ S.K
						if (S.getvisit)  	S.visit 	= S.visit  		:/ S.K
						if (S.hasuser) 	 	S.user  	= S.user   		:/ S.K
						if (S.getsurvival) 	S.survival 	= S.survival 	:/ S.K
					}

					//scale probs to total of 1
					if (S.getprobs & (S.method==0 | S.method==3)) 	S.pt = S.pt :/ rowsum(S.pt)	
					

				//deriv
				if (S.getprobs) 	asarray(A,b,(pt:-S.pt):/hstep)
				if (S.getlos) 	 	asarray(B,b,(los:-S.los):/hstep)
				if (S.getrmst) 	 	asarray(C,b,(rmst:-S.rmst):/hstep)
				if (S.hasuser) 		asarray(D,b,(user:-S.user):/hstep)
				if (S.getvisit)  	asarray(E,b,(visit:-S.visit):/hstep)
				if (S.gethazard) 	asarray(F,b,(hazard:-S.hazard):/hstep)
				if (S.getsurvival) 	asarray(G,b,(survival:-S.survival):/hstep)

				stata("_dots "+strofreal(b)+" 0")
			}
			
			//store derivatives
			if (S.getprobs) 	asarray(S.probs,(from,S.at,2),A) 
			if (S.getlos) 		asarray(S.loss,(from,S.at,2),B)
			if (S.getrmst) 		asarray(S.rmsts,(from,S.at,2),C)
			if (S.hasuser) 		asarray(S.users,(from,S.at,2),D)
			if (S.getvisit) 	asarray(S.visits,(from,S.at,2),E)
			if (S.gethazard) 	asarray(S.hazards,(from,S.at,2),F)
			if (S.getsurvival) 	asarray(S.survivals,(from,S.at,2),G)

		}
		else {
		
			printf("\n")
			stata(`"_dots 0 , title("Calculating CIs via parametric bootstrap for at"'+strofreal(S.at)+`"()") reps("'+strofreal(S.M)+")")
		
			// loop over simulations
			for (k=1;k<=S.M;k++) {

				S.Kind = k
				
				//reset to 0 for each bootstrap sample in case of standardising
				if (S.method==1) 	S.pt = J(S.obs,S.Nstates^2,0)
				else 				predictms_init_storage(S,from)

				//1 or std loop
				for (std=1;std<=S.K;std++) {
					
					S.std = std
					if (ptlosvisit) {
						if 		(S.method==0 | S.method==3)	predictms_sim(S,from)
						else if (S.method==1)				predictms_aj(S,from)	
						else if (S.method==2)				predictms_analytic(S,from)
					}
					if (S.gethazard | S.getsurvival) 		predictms_model_predict(S,from)
					
				}

				if (S.standardise) {
					if (S.getprobs) 	S.pt 		= S.pt 			:/ S.K
					if (S.getlos) 	 	S.los   	= S.los    		:/ S.K
					if (S.getrmst) 	 	S.rmst   	= S.rmst    	:/ S.K
					if (S.hasuser) 	 	S.user  	= S.user   		:/ S.K
					if (S.getvisit)  	S.visit 	= S.visit  		:/ S.K
					if (S.getsurvival) 	S.survival 	= S.survival 	:/ S.K
				}

				//scale probs to total of 1
				if (S.getprobs & (S.method==0 | S.method==3)) 	S.pt = S.pt :/ rowsum(S.pt)

				//Store predictions
				//-> time points x simulations, for each current and next state combination
				if (S.getprobs) 	_store_prob(A,from,S)
				if (S.getlos) 		_store_los(B,from,S)
				if (S.getrmst) 		_store_rmst(C,from,S)
				if (S.hasuser) 		_store_user(D,from,S)		
				if (S.getvisit) 	_store_visit(E,from,S)	
				if (S.gethazard) 	_store_hazard(F,from,S)			
				if (S.getsurvival) 	_store_survival(G,from,S)			
				
				stata("_dots "+strofreal(S.Kind)+" 0")
			}
			
			//store bootstraps
			if (S.getprobs)		asarray(S.probs,(from,S.at,2),A)
			if (S.getlos) 		asarray(S.loss,(from,S.at,2),B)
			if (S.getrmst) 		asarray(S.rmsts,(from,S.at,2),C)
			if (S.hasuser) 	 	asarray(S.users,(from,S.at,2),D)
			if (S.getvisit)  	asarray(S.visits,(from,S.at,2),E)
			if (S.gethazard) 	asarray(S.hazards,(from,S.at,2),F)
			if (S.getsurvival) 	asarray(S.survivals,(from,S.at,2),G)
			
		}
		
	}
	
	//tidy up gml structs
	if (S.method==2 | S.gethazard | S.getsurvival) {
		for (trans=1;trans<=S.Ntrans;trans++) {
			if (strtrim(st_local("GML"+strofreal(trans)))!="") {
				rmexternal(st_local("GML"+strofreal(trans)))
			}
		}
	}
	
}

void _store_prob(`TR' A, `RS' from, `SS' S) 
{
	if (S.Kind==1) {
		for (tr=1;tr<=S.Nstates;tr++) {
			asarray(A,tr,S.pt[,tr])								
		}
	}
	else {
		for (tr=1;tr<=S.Nstates;tr++) {
			asarray(A,tr,(asarray(A,tr),S.pt[,tr]))
		}
	}
}

void _store_los(`TR' B, `RS' from, `SS' S) 
{
	if (S.Kind==1) {
		for (tr=1;tr<=S.Nstates;tr++) {
			asarray(B,tr,S.los[,tr])								
		}
	}
	else {
		for (tr=1;tr<=S.Nstates;tr++) {
			asarray(B,tr,(asarray(B,tr),S.los[,tr]))
		}
	}
}

void _store_rmst(`TR' C, `RS' from, `SS' S) 
{
	if (S.Kind==1) 	asarray(C,1,S.rmst)
	else 			asarray(C,1,(asarray(C,1),S.rmst))
}

void _store_user(`TR' D, `RS' from, `SS' S) 
{
	nvars = cols(S.user)
	if (S.Kind==1) {
		for (tr=1;tr<=nvars;tr++) {
			asarray(D,tr,S.user[,tr])								
		}
	}
	else {
		for (tr=1;tr<=nvars;tr++) {
			asarray(D,tr,(asarray(D,tr),S.user[,tr]))
		}
	}
}

void _store_visit(`TR' E, `RS' from, `SS' S) 
{
	if (S.Kind==1) {
		for (tr=1;tr<=S.Nstates;tr++) {
			asarray(E,tr,S.visit[,tr])								
		}
	}
	else {
		for (tr=1;tr<=S.Nstates;tr++) {
			asarray(E,tr,(asarray(E,tr),S.visit[,tr]))
		}
	}
}

void _store_hazard(`TR' F, `RS' from, `SS' S) 
{
	nvars = cols(S.hazard)
	if (S.Kind==1) {
		for (tr=1;tr<=nvars;tr++) {
			asarray(F,tr,S.hazard[,tr])								
		}
	}
	else {
		for (tr=1;tr<=nvars;tr++) {
			asarray(F,tr,(asarray(F,tr),S.hazard[,tr]))
		}
	}
}

void _store_survival(`TR' G, `RS' from, `SS' S) 
{
	nvars = cols(S.survival)
	if (S.Kind==1) {
		for (tr=1;tr<=nvars;tr++) {
			asarray(G,tr,S.survival[,tr])								
		}
	}
	else {
		for (tr=1;tr<=nvars;tr++) {
			asarray(G,tr,(asarray(G,tr),S.survival[,tr]))
		}
	}
}

end

