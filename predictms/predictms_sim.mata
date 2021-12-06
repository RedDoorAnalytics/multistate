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

mata:

void predictms_sim(`SS' S, `RS' from)
{

	states 	= J(S.N,1,from)				
	stimes 	= J(S.N,1,S.enter)
	ind 	= 1	
	done 	= J(S.N,1,0)
	Ndone	= 0
	
	while (Ndone!=S.N) {

		// initialise new states and times
		states = states,states[,ind]
		stimes = stimes,stimes[,ind]

		for (j=1; j<=S.Nstates; j++) {

			// if there are possible next states
			if (S.Nnextstates[j]) {

				// update index for current state and not done
				index 	= select(S.coreindex,(states[,ind]:==j) :* (done:==0))	
				Nsim 	= rows_cols(index)
				if (!(Nsim[1] & Nsim[2])) continue								

				// simulate survival times and events
				if (!S.method) {
					if (S.iscox)	newdata = predictms_sim_cr_cox(S,Nsim[1],stimes[index,ind],j) 
					else			newdata = predictms_sim_cr(S,Nsim[1],stimes[index,ind],j)
				}
				else 				newdata = predictms_sim_cr_latent(S,Nsim[1],stimes[index,ind],j)
				
				if (S.survsim) 	predictms_survsim_post(newdata[,1])

				stimes[index,ind+1] = newdata[,1]
				states[index,ind+1] = newdata[,2]
				done[index] 		= newdata[,3]

			}
			
		}	
		
		Ndone = sum(done)
		ind++

	} 	
	
	if (S.getprobs) 			predictms_calc_prob(S,stimes,states)
	if (S.getlos | S.getrmst) 	predictms_calc_los(S,stimes,states)
	if (S.getrmst) 				predictms_calc_rmst(S,stimes,states)
	if (S.getvisit) 			predictms_calc_visit(S,stimes,states)
	if (S.hasuser) 				predictms_calc_user(S)

	if (S.tosave) predictms_sim_save(S,stimes,states)
	
}

end

