*! version 1.0.0 ?????2015 MJC

/*
Notes
*/

/*
History
05sep2016 - bug fix in los
*/

version 12.1

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix
local RM 	real matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local PS	pointer(struct predictms_struct scalar) scalar

mata:

`TR' predictms_core(	`PS' p,			///
						`RS' dmind,		///
						`RS' from)
{

	struct predictms_struct scalar S
	`RM' posnextstates
	
	S = *p

	if (S.getcis) {
		transmorphic A
		A = asarray_create("real",1)		//to hold predictions across trans and sims
	}

	// loop over simulations
	for (k=1;k<=S.M;k++) {
	
		pt = J(S.obs,S.Nstates,0)				//to hold predictions
		
		//std loop
		for (std=1;std<=S.K;std++) {
		
			for (t=1;t<=S.Nobstosim;t++) {
				
				//all start in state from at time 0 defined in enter()
				states = J(S.N,1,from)
				stimes = J(S.N,1,S.enter[t])

				ind = 1	//indexes the move
				
				for (j=from;j<=S.Nstates;j++) {
				
					//if there are possible next states
					if (S.Nnextstates[j]>0) {
					
						//update index vector
						index = select(S.coreindex,(states[,ind]:==j))
						Nsim = rows(index)										//number of obs to simulate
						if (!Nsim) continue										//no one in current state
						
						//use previous time and state as basis
						stimes = stimes,stimes[,ind]														
						states = states,states[,ind]

						//simulate survival times
						newt = J(Nsim,S.Nnextstates[j,1],0)
						Umat = runiform(Nsim,S.Nnextstates[j,1])		//need to do this here to keep U draws the same from separate model approach
						posnextstates = asarray(S.posnextstatesj,j)		//possible next states
						
						for (ij=1; ij<=S.Nnextstates[j,1]; ij++) {

							trans = asarray(S.postrans,j)[ij]
							
							//RP
							if (S.modeltrans[trans,6]) {
								test1 = asarray(S.transinfo,(trans,2))
								pstpm2 = asarray(S.transinfo,(trans,4))
								if (S.reset) {
									temp3 = S.maxt:-(stimes[index,ind]):+0.01		//plus a bit takes care of censoring in state calc below
									rc = mm_root_vec(x=J(Nsim,1,.),&stpm2_sim(),mint=smallestdouble(),temp3,0,1000,Umat[,ij],test1[k,]',dmind,std,pstpm2) 	
									newt[,ij] = x 
								}
								else {
									test3 = stimes[index,ind]
									rc = mm_root_vec(x=J(Nsim,1,.),&stpm2_sim(),mint=smallestdouble(),S.maxt+0.01,0,1000,Umat[,ij],test1[k,]',dmind,std,pstpm2,test3)	
									newt[,ij] = x 
								}
								continue							
							}
							//strcs
							if (S.modeltrans[trans,7]) {
								test1 = asarray(S.transinfo,(trans,2))
								pstrcs = asarray(S.transinfo,(trans,4))
								if (S.reset) {
									temp3 = S.maxt:-(stimes[index,ind]):+0.01		//plus a bit takes care of censoring in state calc below
									rc = mm_root_vec(x=J(Nsim,1,.),&strcs_sim(),mint=smallestdouble(),temp3,0,1000,Umat[,ij],test1[k,]',dmind,std,pstrcs) 	
									newt[,ij] = x 
								}
								else {
									test3 = stimes[index,ind]
									rc = mm_root_vec(x=J(Nsim,1,.),&strcs_sim(),mint=smallestdouble(),S.maxt+0.01,0,1000,Umat[,ij],test1[k,]',dmind,std,pstrcs,test3)	
									newt[,ij] = x 
								}
								continue							
							}
							//exponential
							if (S.modeltrans[trans,1]) {
								lambda = exp(asarray(S.transinfo,(trans,dmind))[std,] * asarray(S.transinfo,(trans,2))[k,]')	
								if (S.reset) newt[,ij] = log(Umat[,ij]) :/ (-lambda)
								else newt[,ij] = log( exp(-lambda :* stimes[index,ind]) :* Umat[,ij]) :/ (-lambda)
								continue
							}
							//weibull
							if (S.modeltrans[trans,2]) {
								indices = asarray(S.transinfo,(trans,3))
								lambda = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,1]::indices[2,1]][std,] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,1]|])
								gamm = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,2]::indices[2,2]][std,] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,2]|])
								if (S.reset) newt[,ij] = (log(Umat[,ij]) :/ (-lambda)):^(1:/gamm)
								else newt[,ij] = ((log(exp(-lambda :* stimes[index,ind]:^gamm) :* Umat[,ij]) :/ (-lambda)):^(1:/gamm))
								continue
							}
							//gompertz
							if (S.modeltrans[trans,3]) {
								indices = asarray(S.transinfo,(trans,3))
								lambda = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,1]::indices[2,1]][std,] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,1]|])
								gamm = asarray(S.transinfo,(trans,dmind))[,indices[1,2]::indices[2,2]][std,] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,2]|]
								if (S.reset) newt[,ij] = log(((log(Umat[,ij]) :/ (-lambda)):*gamm):+1):/gamm
								else newt[,ij] = log(((log( exp(-lambda:/gamm:*(exp(gamm:*stimes[index,ind]):-1)) :* Umat[,ij]) :/ (-lambda)):*gamm):+1):/gamm
								continue
							}
							//llog
							if (S.modeltrans[trans,4]) {
								indices = asarray(S.transinfo,(trans,3))
								lambda = exp(-asarray(S.transinfo,(trans,dmind))[,indices[1,1]::indices[2,1]][std,] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,1]|])
								gamm = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,2]::indices[2,2]][std,] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,2]|])
								if (S.reset) newt[,ij] = (((1:/Umat[,ij]):-1):^(gamm)):/lambda
								else newt[,ij] = (((1:/(((1 :+ (lambda:*stimes[index,ind]):^(1:/gamm)):^(-1)):*Umat[,ij])):-1):^(gamm)):/lambda
								continue
							}
							//lognormal
							if (S.modeltrans[trans,5]) {
								lambda = asarray(S.transinfo,(trans,dmind))[,indices[1,1]::indices[2,1]][std,] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,1]|]
								gamm = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,2]::indices[2,2]][std,] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,2]|])
								newt[,ij] = exp(invnormal(1:-Umat[,ij]):*gamm :+ lambda)
							}
							
						}

						if (S.survsim) {	//post survival times and escape out
							st_store(.,st_local("survsim"),st_local("survsimtouse"),newt)
							exit()
						}
						
						//calculate new time and new state
						if (S.reset) newt = newt,(S.maxt:-stimes[index,ind])
						else newt = newt,J(Nsim,1,S.maxt)
						nextt = rowmin(newt)
						if (S.reset) stimes[index,ind+1] = nextt :+ stimes[index,ind]
						else stimes[index,ind+1] = nextt
						states[index,ind+1] = (nextt :== newt) * posnextstates		
						ind++
					}
				}

				if (!S.isenter & !S.prob) {
					if (S.los) {
						exit(1986)
					}
					else {
						//calculate probability of being in each state
						if (S.isstd) pt[t,] = pt[t,] :+ colsum(states[,cols(states)] :== S.statematrix) :/ S.N
						else pt[t,] = colsum(states[,cols(states)] :== S.statematrix) :/ S.N
					}
				}
				
			}
		
			if (S.isenter & !S.prob) {
				if (S.los) {
					//length of stay in each state			
					for (t=1;t<=S.obs;t++) {
						temptime = stimes
						for (s=2; s<=cols(states); s++) {
							index = selectindex(stimes[,s] :> S.predtime[t])
							temptime[index,s] = J(rows(index),1,S.predtime[t])
						}						
						temptime = temptime,J(S.N,1,S.predtime[t])
						temptime = temptime[,2..cols(temptime)] :- temptime[,1..cols(states)]
						
						temppred = J(1,S.Nstates,0)
						for (s=1; s<=S.Nstates; s++) {
							temppred[1,s] = sum((states:==s) :* temptime):/S.N
						}
						if (S.isstd) pt[t,] = pt[t,] :+ temppred
						else pt[t,] = temppred
					}
				}
				else {
					//calculate probability of being in each state
					for (t=1;t<S.obs;t++) {
						tempstate = states[,1]
						for (s=2; s<=cols(states); s++) {
							index = selectindex(stimes[,s] :<= S.predtime[t])
							if (cols(index)) tempstate[index] = states[index,s]
						}
						if (S.isstd) pt[t,] = pt[t,] :+ colsum(tempstate :== S.statematrix) :/S.N
						else pt[t,] = colsum(tempstate :== S.statematrix) :/S.N
					}
					if (S.isstd) pt[t,] = pt[t,] :+ colsum(states[,cols(states)] :== S.statematrix) :/S.N
					else pt[t,] = colsum(states[,cols(states)] :== S.statematrix) :/S.N
				}
			}

		}			
		
		if (S.prob) return(states,stimes)
		
		if (S.isstd) pt = pt:/S.K
				
		if (!S.getcis) return(pt) 
		else {
			//store them -> time points x simulations, for each current and next state combination
			if (k==1) for (tr=1;tr<=S.Nstates;tr++) asarray(A,tr,pt[,tr])								
			else for (tr=1;tr<=S.Nstates;tr++) asarray(A,tr,(asarray(A,tr),pt[,tr]))
		}
	}
	if (S.getcis) return(A)
	
}

end

