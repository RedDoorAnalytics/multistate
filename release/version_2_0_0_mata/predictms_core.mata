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

mata:

void predictms_core(					///
						`SS' S,			///
						`RS' from,		///
						`RS' at			///
					)
{

	`RS' k, std, t, ind, Nsim, trans, t2flag
	`RM' pt, states, stimes, index, posnextstates, newt, Umat, U
	`RC' b, rc

	if (S.getcis) {
		A = asarray_create("real",1)											//to hold predictions across trans and sims
		if (S.getlos) 	B = asarray_create("real",1)
		if (S.hasuser) 	C = asarray_create("real",1)
		if (S.getvisit) D = asarray_create("real",1)
	}

	// loop over simulations
	for (k=1;k<=S.M;k++) {

		S.pt[,] = J(S.obs,S.Nstates,0)											//sim specific predictions
		if (S.getlos) 	S.los[,] 	= J(S.obs,S.Nstates,0)
		if (S.getvisit) S.visit[,] 	= J(S.obs,S.Nstates,0)
		
		//std loop
		for (std=1;std<=S.K;std++) {

			for (t=1;t<=S.Nobstosim;t++) {
				
				states 	= J(S.N,1,from)											//all start in state from at time defined in enter()
				stimes 	= J(S.N,1,S.enter[t])
				ind 	= 1														//indexes the move
				done 	= J(S.N,1,0)
				Ndone	= sum(done)
				
				while (Ndone!=S.N) {

					//until they're all done, for each move need to loop over all states
					//as obs could be in any of them
					//use previous time and state as basis
					states = states,states[,ind]
					stimes = stimes,stimes[,ind]	
					
					for (j=1;j<=S.Nstates;j++) {
					
						//if there are possible next states
						if (S.Nnextstates[j]>0) {
						
							index 	= select(S.coreindex,(states[,ind]:==j) :* (done:==0))	//update index vector
							Nsim 	= rows(index)											//number of obs to simulate
							if (!Nsim) continue												//no one in current state so escape out of loop

							//simulate survival times
							newt 			= J(Nsim,S.Nnextstates[j,1],0)
							Umat 			= runiform(Nsim,S.Nnextstates[j,1])			//need to do this here to keep U draws the same from separate model approach
							posnextstates 	= asarray(S.posnextstatesj,j)				//possible next states
							
							for (ij=1; ij<=S.Nnextstates[j,1]; ij++) {
								
								trans 	= asarray(S.postrans,j)[ij]
								dm		= asarray(S.X,(trans,at))[std,]
								b 		= asarray(S.transinfo,(trans,2))[k,]'
								U 		= Umat[,ij]
								
								t2flag 	= sum(S.tscale2:==trans)  						//second timescale S.time2
								if (t2flag) t2 = S.time2[at]

								//megenreg
								if (S.modeltrans[trans,8]) {
									tmin	= smallestdouble()
									func 	= &megenreg_sim()
									if (S.reset) {
										Pf 		= predictms_megenreg_setup(b,at,Nsim,trans)
										tmax = S.maxt:-(stimes[index,ind]):+ 0.01		//plus a bit takes care of censoring in state calc below
										rc = predictms_mm_root_vec(x=J(Nsim,1,.),func,tmin,tmax,0,1000,U,Pf) 	
										newt[,ij] = x 
									}
									else {
										tmax = S.maxt :+ 0.01
										t0 = stimes[index,ind] :+ smallestdouble()		//for t0=0
										if (t2flag) {
											t0 		= t0 :+ t2						
											tmin 	= tmin :+ t2
											tmax 	= tmax :+ t2
										}
										Pf = predictms_megenreg_setup(b,at,Nsim,trans,t0)							
										rc = predictms_mm_root_vec(x=J(Nsim,1,.),func,tmin,tmax,0,1000,U,1,2,3,4,Pf,t0)	
										newt[,ij] = x 														//^fudges for t0 vector in mmroot
									}
									rmexternal(st_global("e(struct)")) 					//tidy up
								}
								//RP
								else if (S.modeltrans[trans,6]) {
									Pstpm2 	= asarray(S.transinfo,(trans,4))
									tmin	= smallestdouble()
									if (S.reset) {
										tmax = S.maxt:-(stimes[index,ind]):+ 0.01		//plus a bit takes care of censoring in state calc below
										rc = predictms_mm_root_vec(x=J(Nsim,1,.),&stpm2_sim(),tmin,tmax,0,1000,U,b,dm,at,std,Pstpm2) 	
										newt[,ij] = x 
									}
									else {
										tmax = S.maxt :+ 0.01
										t0 = stimes[index,ind] :+ smallestdouble()		//for t0=0
										if (t2flag) {
											t0 		= t0 :+ t2						
											tmin 	= tmin :+ t2
											tmax 	= tmax :+ t2
										}
										rc = predictms_mm_root_vec(x=J(Nsim,1,.),&stpm2_sim(),tmin,tmax,0,1000,U,b,dm,at,std,Pstpm2,t0)	
										newt[,ij] = x 
									}
									
								}
								//strcs
								else if (S.modeltrans[trans,7]) {
									Pstrcs 	= asarray(S.transinfo,(trans,4))
									tmin	= smallestdouble()
									if (S.reset) {
										tmax = S.maxt:-(stimes[index,ind]):+ 0.01		//plus a bit takes care of censoring in state calc below
										rc = predictms_mm_root_vec(x=J(Nsim,1,.),&strcs_sim(),tmin,tmax,0,1000,U,b,dm,at,std,Pstrcs) 	
										newt[,ij] = x 
									}
									else {
										tmax = S.maxt :+ 0.01
										t0 = stimes[index,ind] :+ smallestdouble()		//for t0=0
										if (t2flag) {
											t0 		= t0 :+ t2						
											tmin 	= tmin :+ t2
											tmax 	= tmax :+ t2
										}
										rc = predictms_mm_root_vec(x=J(Nsim,1,.),&strcs_sim(),tmin,tmax,0,1000,U,b,dm,at,std,Pstrcs,t0)	
										newt[,ij] = x 
									}
								}
								//exponential
								else if (S.modeltrans[trans,1]) {
									lambda = exp(dm * b)	
									if (S.reset) 	newt[,ij] = log(U) :/ (-lambda)
									else {
										if (t2flag) newt[,ij] = log( exp(-lambda :* (stimes[index,ind] :+ t2)) :* U) :/ (-lambda)
										else 		newt[,ij] = log( exp(-lambda :* stimes[index,ind]) :* U) :/ (-lambda)
									}
								}
								//weibull
								else if (S.modeltrans[trans,2]) {
									indices = asarray(S.transinfo,(trans,3))
									lambda 	= exp(dm[|.,indices[1,1]\.,indices[2,1]|] * b[|indices[,1]|])
									gamm 	= exp(dm[|.,indices[1,2]\.,indices[2,2]|] * b[|indices[,2]|])
									if (S.reset) 	newt[,ij] = (log(U) :/ (-lambda)):^(1:/gamm)
									else {
										if (t2flag) newt[,ij] = ((log(exp(-lambda :* (stimes[index,ind] :+ t2):^gamm) :* U) :/ (-lambda)):^(1:/gamm))
										else 		newt[,ij] = ((log(exp(-lambda :* stimes[index,ind]:^gamm) :* U) :/ (-lambda)):^(1:/gamm))
									}
								}
								//gompertz
								else if (S.modeltrans[trans,3]) {
									indices = asarray(S.transinfo,(trans,3))
									lambda 	= exp(dm[|.,indices[1,1]\.,indices[2,1]|] * b[|indices[,1]|])
									gamm 	= dm[|.,indices[1,2]\.,indices[2,2]|] * b[|indices[,2]|]
									if (S.reset) 	newt[,ij] = log(((log(U) :/ (-lambda)):*gamm):+1):/gamm
									else {
										if (t2flag) newt[,ij] = log(((log( exp(-lambda:/gamm:*(exp(gamm:*(stimes[index,ind] :+ t2)):-1)) :* U) :/ (-lambda)):*gamm):+1):/gamm
										else 		newt[,ij] = log(((log( exp(-lambda:/gamm:*(exp(gamm:*stimes[index,ind]):-1)) :* U) :/ (-lambda)):*gamm):+1):/gamm
									}
								}
								//llog
								else if (S.modeltrans[trans,4]) {
									indices = asarray(S.transinfo,(trans,3))
									lambda 	= exp(-dm[|.,indices[1,1]\.,indices[2,1]|] * b[|indices[,1]|])
									gamm 	= exp(dm[|.,indices[1,2]\.,indices[2,2]|] * b[|indices[,2]|])
									if (S.reset) 	newt[,ij] = (((1:/U):-1):^(gamm)):/lambda
									else {
										if (t2flag) newt[,ij] = (((1:/(((1 :+ (lambda:*(stimes[index,ind] :+ t2)):^(1:/gamm)):^(-1)):*U)):-1):^(gamm)):/lambda
										else 		newt[,ij] = (((1:/(((1 :+ (lambda:*stimes[index,ind]):^(1:/gamm)):^(-1)):*U)):-1):^(gamm)):/lambda
									}
								}
								//lognormal
								else {
									indices 	= asarray(S.transinfo,(trans,3))
									lambda 		= dm[|.,indices[1,1]\.,indices[2,1]|] * b[|indices[,1]|]
									gamm 		= exp(dm[|.,indices[1,2]\.,indices[2,2]|] * b[|indices[,2]|])
									newt[,ij] 	= exp(invnormal(1:-U):*gamm :+ lambda)
								}
								
								//if second timescale, need to take away again to scale back to main timescale
								//otherwise rowmin() below, will be off and will always pick the large one
								if (t2flag) newt[,ij] = newt[,ij] :- t2
							}

							if (S.survsim) {	//post survival times and escape out
								st_store(.,st_local("survsim"),st_local("survsimtouse"),newt)
								exit()
							}
							
							//calculate new time and new state
							if (S.reset) 	newt = newt,(S.maxt:-stimes[index,ind])
							else 			newt = newt,J(Nsim,1,S.maxt)
							nextt = rowmin(newt)
							
							if (S.reset) 	stimes[index,ind+1] = nextt :+ stimes[index,ind]
							else 			stimes[index,ind+1] = nextt
							states[index,ind+1] = (nextt :== newt) * posnextstates		
						
						}
					}															

					ind++	//keep here
					
					//now update done for any that are in an absorbing state
					notdoneindex = selectindex(done:==0)
					for (j=1;j<=S.Nstates;j++) {
						if (S.Nnextstates[j]==0) {
							toupdate 		= select(notdoneindex,states[notdoneindex,ind]:==j)
							Ntoupdate		= rows(toupdate)
							if (Ntoupdate) 	done[toupdate] = J(Ntoupdate,1,1)
						}
					}
					//update done for any who have reached maxt
					notdoneindex 	= selectindex(done:==0)
					toupdate 		= select(notdoneindex,stimes[notdoneindex,ind]:==S.maxt)
					Ntoupdate		= rows(toupdate)
					if (Ntoupdate) 	done[toupdate] = J(Ntoupdate,1,1)
					
					Ndone = sum(done)
					
				} //end of while loop

				if (!S.isenter) {												//fixed-horizon predictions
					if (S.getlos) {
						exit(1986)
					}
					else {
						//calculate probability of being in each state
						if (S.isstd) 	S.pt[t,] = S.pt[t,] :+ colsum(states[,cols(states)] :== S.statematrix) :/ S.N
						else 			S.pt[t,] = colsum(states[,cols(states)] :== S.statematrix) :/ S.N
					}
				}
				
			}

			if (S.isenter) {
			
				//calculate probability of being in each state
				Ntotstates = cols(states)
				for (t=1;t<S.obs;t++) {
					tempstate = states[,1]
					for (s=2; s<=Ntotstates; s++) {
						index = selectindex(stimes[,s] :<= S.predtime[t])
						if (cols(index)) tempstate[index] = states[index,s]
					}
					if (S.isstd) 	S.pt[t,] = S.pt[t,] :+ colsum(tempstate :== S.statematrix) :/S.N
					else 			S.pt[t,] = colsum(tempstate :== S.statematrix) :/ S.N
				}
				if (S.isstd) 	S.pt[t,] = S.pt[t,] :+ colsum(states[,cols(states)] :== S.statematrix) :/S.N
				else 			S.pt[t,] = colsum(states[,cols(states)] :== S.statematrix) :/S.N
				
				//length of stay in each state	
				if (S.getlos) {
					Ntotstates = cols(states)
					for (t=1;t<=S.obs;t++) {
						temptime = stimes
						for (s=2; s<=Ntotstates; s++) {
							index = selectindex(stimes[,s] :> S.predtime[t])
							temptime[index,s] = J(rows(index),1,S.predtime[t])
						}						
						temptime = temptime,J(S.N,1,S.predtime[t])
						temptime = temptime[|.,2\.,cols(temptime)|] :- temptime[|.,1\.,cols(states)|] 
						
						temppred = J(1,S.Nstates,0)
						for (s=1; s<=S.Nstates; s++) {
							temppred[1,s] = sum((states:==s) :* temptime):/S.N
						}
						if (S.isstd) 	S.los[t,] = S.los[t,] :+ temppred
						else 			S.los[t,] = temppred
					}
				}

				//prob of ever visiting each state
				if (S.getvisit) {
					Ntotstates = cols(states)
					
					for (t=1;t<=S.obs;t++) {
						temppred = J(S.N,S.Nstates,0)
						flag = stimes:<=S.predtime[t]
						for (v=1;v<=Ntotstates;v++) {
							flag2 = selectindex(flag[,v])
							if (cols(flag2) & rows(flag2)) {
								for (s=1;s<=S.Nstates;s++) {
									temppred[flag2,s] = temppred[flag2,s] :+ (states[flag2,v]:==s)
								}
							}
						}
						temppred = colsum(temppred :> 0):/S.N
						if (S.isstd) 	S.visit[t,] = S.visit[t,] :+ temppred
						else 			S.visit[t,] = temppred
					}
				}

				//user-defined mata function
				if (S.hasuser) {
					if (S.isstd) 		S.user = S.user :+ (*S.userfunc)(S)
					else  				S.user = (*S.userfunc)(S)
					if (S.Nuserflag) {
						S.Nuservars = cols(S.user)	
						S.Nuserflag = 0
					}
				}
				
			}

		} //end std			

		//standardisation
		if (S.isstd) {
			S.pt = S.pt :/ S.K
			if (S.getlos) 	S.los = S.los :/ S.K
			if (S.hasuser) 	S.user = S.user :/ S.K
		}

		//scale probs to total of 1
		S.pt = S.pt :/ rowsum(S.pt)

		//Store predictions
		if (!S.getcis) {

			asarray(S.probs,(from,at),S.pt)
			if (S.getlos) 	asarray(S.loss,(from,at),S.los)
			if (S.getvisit) asarray(S.visits,(from,at),S.visit)
			if (S.hasuser) 	asarray(S.users,(from,at),S.user)
			
		}
		else {
			
			//-> time points x simulations, for each current and next state combination
			if (k==1) {
				for (tr=1;tr<=S.Nstates;tr++) {
					asarray(A,tr,S.pt[,tr])								
				}
			}
			else {
				for (tr=1;tr<=S.Nstates;tr++) {
					asarray(A,tr,(asarray(A,tr),S.pt[,tr]))
				}
			}
			asarray(S.probs,(from,at),A)
			
			if (S.getlos) {
				if (k==1) {
					for (tr=1;tr<=S.Nstates;tr++) {
						asarray(B,tr,S.los[,tr])								
					}
				}
				else {
					for (tr=1;tr<=S.Nstates;tr++) {
						asarray(B,tr,(asarray(B,tr),S.los[,tr]))
					}
				}
				asarray(S.loss,(from,at),B)
			}
			
			if (S.getvisit) {
				if (k==1) {
					for (tr=1;tr<=S.Nstates;tr++) {
						asarray(D,tr,S.visit[,tr])								
					}
				}
				else {
					for (tr=1;tr<=S.Nstates;tr++) {
						asarray(D,tr,(asarray(D,tr),S.visit[,tr]))
					}
				}
				asarray(S.visits,(from,at),D)
			}
			
			if (S.hasuser) {
				nvars = cols(S.user)
				if (k==1) {
					for (tr=1;tr<=nvars;tr++) {
						asarray(C,tr,S.user[,tr])								
					}
				}
				else {
					for (tr=1;tr<=nvars;tr++) {
						asarray(C,tr,(asarray(C,tr),S.user[,tr]))
					}
				}
				asarray(S.users,(from,at),C)
			}
			
		}
	}

}

end

