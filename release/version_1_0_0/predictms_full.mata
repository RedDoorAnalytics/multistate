*! version 1.0.8 14may2015 MJC

/*
Notes
-> Comments marked by //&& indicate places which could be extended
-> Comments marked by //## indicates code which I don't like
-> currently coding it as a forward only transition matrix
-> should think about the rc from mm_root
-> Change to general transition matrix
*/

//now only doing 1 sim for standard predictions
//sim for each time point for dynamic predictions

/*
History
MJC 14may2015: verison 1.0.8 -
MJC 11may2015: verison 1.0.7 -
MJC 11may2015: verison 1.0.6 - fixed bug with weibull and forward approach
MJC 09may2015: version 1.0.5 - fixed bug in from() which occurred when anything but from(1) was used
							 - fixed bug in forward calculations when enter>0
							 - fixed bug that only calculated predictions to states you could go to from first state
MJC 09may2015: version 1.0.4 - clock-forward approach now the default (simulations incorporate delayed entry), reset option added to use clock-reset
							 - only reset approach allowed with streg, dist(lnormal)
MJC 06may2015: version 1.0.3 - put next states from each current state matrices etc. into arrays for speed gains
							 - synced with streg, dist(exp|gompertz|llogistic|lnormal)
							 - normal approximation for CIs synced
MJC 15apr2015: version 1.0.2 - when no ci's calculated it was using first draw from MVN, this has been fixed to be e(b)
MJC 01apr2015: version 1.0.1 - stpm2 simulation greatly improved by creating and passing struct
							 - odds and normal scales added for stpm2 models
MJC 31mar2015: version 1.0.0
*/


version 12.1

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector

mata:

void predictms_fullmodel()
{
	`SS' touse
	`RC' predtime, Nnextstates, coreindex, index
	`RS' getcis, obs, N, M, isreg, isgenreg, ispm2, Nmleqns, Nbetas, percentiles, reset, maxt
	`NM' transmat, eventsvars, simeventtimes, base, temp, X, statematrix
	
	//Core info
	reset = st_local("reset")!=""											//default clock forward, o/w reset
	touse = st_local("touse")												//touse for predictions for final predictions
	obs = strtoreal(st_local("obs"))										//# of time points to calculate probs at
	
	//time
	isenter = st_local("exit")==""
	predtime = st_data(.,st_local("timevar"),touse)
	if (isenter) {
		maxt = max(predtime)
		enter = strtoreal(st_local("enter"))
		Nobstosim = 1
	}
	else {
		maxt = strtoreal(st_local("exit"))
		enter = predtime
		Nobstosim = obs
	}

	getcis = st_local("ci")!=""												//get confidence intervals
	percentiles = st_local("normal")==""								//calculate CIs using percentiles rather than normal approximation
	
	N = strtoreal(st_local("n"))											//sample size
	if (getcis) {
		M = strtoreal(st_local("m"))										//# of draws from MVN
		A = asarray_create("real",2)										//to hold predictions across trans and sims
		st_view(bdraws,.,tokens(st_local("drawvars")),st_local("mvnind"))
	}
	else {
		M = 1
		bdraws = st_matrix(st_local("emat"))
	}
	
	transmat = st_matrix(st_local("transmatrix"))
	Nstates = cols(transmat)
	statematrix = J(N,1,1..Nstates)
	Nnextstates = rowsum(transmat:!=.)															//# of possible next states from each current state
	
	posnextstates = posnextstatesj = postrans = asarray_create("real",1)
	for (i=1;i<=Nstates;i++) {
		asarray(posnextstates,i,strtoreal(tokens(st_local("row"+strofreal(i)+"next")))')
		asarray(posnextstatesj,i,(asarray(posnextstates,i)\i))									//possible next states includng current
		asarray(postrans,i,strtoreal(tokens(st_local("row"+strofreal(i)+"trans")))')
	}
	
	from = strtoreal(tokens(st_local("from")))'
	Nstarts = rows(from)
	
	//Model details
	isreg = st_global("e(cmd2)")=="streg"
	if (isreg) {
		`RS' isweibull, isexp, isgomp, isllog, islnorm
		isweibull = st_global("e(cmd)")=="weibull"
		isexp = st_global("e(cmd)")=="ereg"
		isgomp = st_global("e(cmd)")=="gompertz"
		isllog = st_global("e(cmd)")=="llogistic"
		islnorm = st_global("e(cmd)")=="lnormal"
	}
	isgenreg = st_global("e(cmd)")=="stgenreg"
	ispm2 = st_global("e(cmd)")=="stpm2"
	
	Nmleqns = strtoreal(st_local("Nmleqns"))							//# of linear predictors from model fit
	Nbetas = strtoreal(st_local("Nparams"))								//total # of coefficients in e(b)

	//indexes to extract coefficients for each ml equation
	if (isreg & !isexp) indices = st_matrix(st_local("indices"))
			
	//Read in design matrix for all transitions built in Stata from trans#() options or get stpm2 struct
	if (ispm2) p = predictms_stpm2_fullmodel_setup()
	else X = st_matrix(st_local("dm"))	

	//at each step through the transmat, I only want to simulate survival times for those who are in the current state and so this creates an index for the appropriate rows
	//which is updated as it goes through
	coreindex = 1::N

	// loop over simulations
	for (k=1;k<=M;k++) {
	
		//get coefficients for kth sim
		bmat = bdraws[k,]'
		
		//looping each starting state listed in from()
		for (i=1;i<=Nstarts;i++) {
				
			pt = J(0,Nstates,.)
			
			for (t=1;t<=Nobstosim;t++) {
				//Do simulation for maximum requested time (if enter)
				
				//all start in state from[i] at time 0 defined in enter() or min(timevar) for exit()
				states = J(N,1,from[i])
				stimes = J(N,1,enter[t])

				ind = 1	//indexes the move
				
				for (j=from[i];j<=Nstates;j++) {
				
					//if there are possible next states
					if (Nnextstates[j]>0) {
						
						//update index vector
						index = uniqrows((states[,ind]:==j) :* coreindex)									//gets index for patients in current i th state
						Nsim = rows(index)																	//number of obs to simulate
						
						if (index[1]==0 & Nsim>1) {
							index = index[2::Nsim]															//chop off leading zero
							Nsim = Nsim - 1																	//number of obs to simulate
						}
						if (index[1]==0 & Nsim==1) continue													//no one in current state
								
						//use previous time and state as basis
						stimes = stimes,stimes[,ind]														
						states = states,states[,ind]
						
						//simulate survival times
						//transpose parameters to get nextstates as columns for newt calc.
						//survival times for each transition, and maximum fu
						if (isreg) {
						
							if (isexp) {
								lambda = exp(X[asarray(postrans,j),] * bmat)'															
								if (reset) newt = (log(runiform(Nsim,Nnextstates[j,1])) :/ (-lambda)),(maxt:-stimes[index,ind])
								else newt = (log( exp(-J(Nsim,1,lambda) :* stimes[index,ind]) :* runiform(Nsim,Nnextstates[j,1])) :/ (-lambda)),J(Nsim,1,maxt)
							}
							if (isweibull) {
								lambda = exp(X[asarray(postrans,j),indices[1,1]::indices[2,1]] * bmat[|indices[,1]|])'
								gamm = exp(X[asarray(postrans,j),indices[1,2]::indices[2,2]] * bmat[|indices[,2]|])'																			
								if (reset) newt = ((log(runiform(Nsim,Nnextstates[j,1])) :/ (-lambda)):^(1:/gamm)),(maxt:-stimes[index,ind])
								else newt = ((log(exp(-J(Nsim,1,lambda) :* stimes[index,ind]:^J(Nsim,1,gamm)) :* runiform(Nsim,Nnextstates[j,1])) :/ (-lambda)):^(1:/gamm)),J(Nsim,1,maxt)
							}
							if (isgomp) {
								lambda = exp(X[asarray(postrans,j),indices[1,1]::indices[2,1]] * bmat[|indices[,1]|])'
								gamm = (X[asarray(postrans,j),indices[1,2]::indices[2,2]] * bmat[|indices[,2]|])'
								if (reset) newt = (log(((log(runiform(Nsim,Nnextstates[j,1])) :/ (-lambda)):*gamm):+1):/gamm),(maxt:-stimes[index,ind])
								else newt = (log(((log( exp(-J(Nsim,1,lambda):/J(Nsim,1,gamm):*(exp(J(Nsim,1,gamm):*stimes[index,ind]):-1)) :* runiform(Nsim,Nnextstates[j,1])) :/ (-lambda)):*gamm):+1):/gamm),J(Nsim,1,maxt)
							}
							if (isllog) {
								lambda = exp(-X[asarray(postrans,j),indices[1,1]::indices[2,1]] * bmat[|indices[,1]|])'
								gamm = exp(X[asarray(postrans,j),indices[1,2]::indices[2,2]] * bmat[|indices[,2]|])'
								if (reset) newt = ((((1:/runiform(Nsim,Nnextstates[j,1])):-1):^(gamm)):/lambda),(maxt:-stimes[index,ind])
								else newt = ((((1:/(((1 :+ (J(Nsim,1,lambda):*stimes[index,ind]):^(1:/J(Nsim,1,gamm))):^(-1)):*runiform(Nsim,Nnextstates[j,1]))):-1):^(gamm)):/lambda),J(Nsim,1,maxt)
							}
							if (islnorm) {
								lambda = (X[asarray(postrans,j),indices[1,1]::indices[2,1]] * bmat[|indices[,1]|])'
								gamm = exp(X[asarray(postrans,j),indices[1,2]::indices[2,2]] * bmat[|indices[,2]|])'
								newt = exp(invnormal(1:-runiform(Nsim,Nnextstates[j,1])):*gamm :+ lambda),(maxt:-stimes[index,ind])
							}

							//calculate new time and new state
							nextt = rowmin(newt)
							if (reset) stimes[index,ind+1] = nextt :+ stimes[index,ind]
							else stimes[index,ind+1] = nextt
							states[index,ind+1] = (nextt :== newt) * asarray(posnextstatesj,j)			
						}

						if (ispm2) {
						
							newt = J(Nsim,Nnextstates[j,1],1)
							if (reset) {
								for (ij=1 ; ij<=Nnextstates[j,1] ; ij++) {
									for (q=1;q<=Nsim;q++){
										rc = mm_root(x=1,&stpm2_sim(),mint=1E-08,(maxt:-(stimes[index,ind])[q]),tol=0,maxit=1000,runiform(1,1),bmat,asarray(postrans,j)[ij],p) 	//lower,upper,tol,maxit,U,coefficients,trans #,pointer to struct
										newt[q,ij] = x 
									}
								}
							}
							else {
								for (ij=1 ; ij<=Nnextstates[j,1] ; ij++) {
									for (q=1;q<=Nsim;q++){
										rc = mm_root(x=1,&stpm2_sim(),mint=1E-08,maxt,tol=0,maxit=1000,runiform(1,1),bmat,asarray(postrans,j)[ij],p,(stimes[index,ind])[q]) 	//lower,upper,tol,maxit,U,coefficients,trans #,pointer to struct
										newt[q,ij] = x 
									}
								}
							}

							//calculate new time and new state
							nextt = rowmin(newt)
							ind1 = ind+1
							if (reset) stimes[index,ind1] = nextt :+ stimes[index,ind]
							else stimes[index,ind1] = nextt
							states[index,ind1] = (stimes[index,ind1] :!= maxt) :* ((nextt :== newt) * asarray(posnextstates,j)) 	//event
							states[index,ind1] = states[index,ind1] :+ (states[index,ind1] :== 0 ) :* states[index,ind]				//stay
						}
						ind++
					}
				}
			
				if (!isenter) {
					//calculate probability of being in each state
					pt = pt\(colsum(states[,cols(states)] :== statematrix) :/ N)
				}
			
			}
			//calculate probability of being in each state
			if (isenter) {
				for (t=1;t<obs;t++) {
					tempstate = states[,1]
					for (s=2; s<=cols(states); s++) {
						index = selectindex(stimes[,s] :<= predtime[t]) 
						tempstate[index] = states[index,s]
					}
					pt = pt\(colsum(tempstate :== statematrix) :/ N)
				}
				pt = pt\(colsum(states[,cols(states)] :== statematrix) :/ N)
			}

			if (!getcis) st_store(.,tokens(st_local("fromvars"+strofreal(from[i]))),touse, pt)
			else {
				//store them -> time points x simulations, for each current and next state combination
				if (k==1) for (tr=1;tr<=Nstates;tr++) asarray(A,(i,tr),pt[,tr])									//plus 1 for current
				else for (tr=1;tr<=Nstates;tr++) asarray(A,(i,tr),(asarray(A,(i,tr)),pt[,tr]))					//plus 1 for current
			}
		}

	}
	
	
	if (getcis) {
		if (percentiles) {
			`RS' med, lci, uci
			med = round(M:/2)
			lci = round(0.025:*M)
			uci = round(0.975:*M)
			for (i=1;i<=Nstarts;i++) {	
				st_view(medians,.,tokens(st_local("fromvars"+strofreal(from[i]))),touse)
				st_view(lcis,.,tokens(st_local("fromvarslci"+strofreal(from[i]))),touse)
				st_view(ucis,.,tokens(st_local("fromvarsuci"+strofreal(from[i]))),touse)
				for (tr=1;tr<=Nstates;tr++) {
					base = asarray(A,(i,tr))
					for (j=1;j<=obs;j++) {
						temp = sort(base[j,]',1)
						medians[j,tr] = temp[med,1]
						lcis[j,tr] = temp[lci,1]
						ucis[j,tr] = temp[uci,1]
					}
				}
			}
		}
		else {
			for (i=1;i<=Nstarts;i++) {	
				for (tr=1;tr<=Nstates;tr++) {
					st_view(pred=.,.,tokens(st_local("probvars_"+strofreal(from[i])+"_"+strofreal(tr))),touse)
					base = asarray(A,(i,tr))
					for (j=1;j<=obs;j++) {
						pred[j,] = quadmeanvariance(logit(base[j,])')'
					}
				}
			}
		}

	}
	
	//tidy up
	if (ispm2) rmexternal(st_local("predictms_struct"))
	
}

end









