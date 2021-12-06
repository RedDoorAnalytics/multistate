*! version 1.0.0 ?????2015 MJC

/*
Notes
*/

/*
History

*/

version 12.1

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local PS	pointer(struct predictms_sep_struct scalar) scalar

mata:

struct predictms_sep_struct
{
	`RS' M, getcis, reset, obs, N, isenter, maxt, Nobstosim, Nstates, Nstarts, los
	`RC' predtime, enter, coreindex, Nnextstates, from
	`NM' modeltrans, statematrix
	`TR' transinfo, posnextstatesj, postrans
}

`PS' predictms_sepmod_sim_setup()
{
	struct predictms_sep_struct scalar S
	pointer scalar p

	stata("tempname predictms_sep_struct")
	rmexternal(st_local("predictms_sep_struct"))
	p = crexternal(st_local("predictms_sep_struct"))

	S.getcis = st_local("ci")!=""			//get confidence intervals
	if (S.getcis) {
		S.M = strtoreal(st_local("m"))		//# of draws from MVN
	}
	else S.M = 1
	
	S.reset = st_local("reset")!=""			//default clock forward, o/w reset
	S.obs = strtoreal(st_local("obs"))		//# of time points to calculate probs at
	S.N = strtoreal(st_local("n"))			//sample size
	hasat2 = st_local("at2")!=""
	S.los = st_local("los")!=""

	//time
	S.isenter = st_local("exit")==""
	if (S.isenter) {
		S.predtime = st_data(.,st_local("timevar"),st_local("touse"))
		S.maxt = max(S.predtime)
		S.enter = strtoreal(st_local("enter"))
		S.Nobstosim = 1
	}
	else {
		S.maxt = strtoreal(st_local("exit"))
		S.enter = st_data(.,st_local("timevar"),st_local("touse"))
		S.Nobstosim = S.obs
	}

	//info for each transition
	Ntrans = strtoreal(st_local("Ntrans"))
	S.modeltrans = (tokens(st_local("cmds"))'):==J(Ntrans,1,("ereg","weibull","gompertz","llogistic","lnormal","stpm2"))
	S.transinfo = asarray_create("real",2)
	for (i=1;i<=Ntrans;i++) {
		asarray(S.transinfo,(i,1),st_matrix(st_local("dm"+strofreal(i))))													//design matrix
		if (S.getcis) asarray(S.transinfo,(i,2),st_data(.,tokens(st_local("drawvars"+strofreal(i))),st_local("mvnind")))	//parameter draws
		else asarray(S.transinfo,(i,2),st_matrix(st_local("emat"+strofreal(i))))											//parameter draws
		if (!S.modeltrans[i,6] & !S.modeltrans[i,1]) asarray(S.transinfo,(i,3),st_matrix(st_local("indices"+strofreal(i))))	//ml eqn indices
		if (S.modeltrans[i,6]) {
			stata("qui estimates restore "+st_local("modelests"+strofreal(i)))
			asarray(S.transinfo,(i,4),predictms_stpm2_sepmodel_setup(strofreal(i)))
		}
		if (hasat2) asarray(S.transinfo,(i,5),st_matrix(st_local("at2dm"+strofreal(i))))
	}
	
	//at each step through the transmat, I only want to simulate survival times for those who are in the current state and so this creates an index for the appropriate rows
	//which is updated as it goes through
	S.coreindex = 1::S.N
	transmat = st_matrix(st_local("transmatrix"))
	S.Nstates = cols(transmat)
	S.statematrix = J(S.N,1,1..S.Nstates)
	S.Nnextstates = rowsum(transmat:!=.)															//# of possible next states from each current state
	posnextstates = S.posnextstatesj = S.postrans = asarray_create("real",1)
	for (i=1;i<=S.Nstates;i++) {
		asarray(posnextstates,i,strtoreal(tokens(st_local("row"+strofreal(i)+"next")))')
		asarray(S.posnextstatesj,i,(asarray(posnextstates,i)\i))									//possible next states includng current
		asarray(S.postrans,i,strtoreal(tokens(st_local("row"+strofreal(i)+"trans")))')
	}
	
	swap((*p), S)
	return(p)
}


void predictms_sepmodel()
{
	`RS' hasat2, getcis, percentiles, difference
	`TR' res
	`SS' touse
	`PS' p
	
	touse = st_local("touse")									//touse for final predictions posting
	from = strtoreal(tokens(st_local("from")))'
	Nstarts = rows(from)
	Nstates = strtoreal(st_local("Nstates"))
	obs = strtoreal(st_local("obs"))					//<- and bove could be passed aa arguments to P_predictms_sepmod_sim_setup() instead of reading them in again within it	
	hasat2 = st_local("at2")!=""
	los = st_local("los")!=""
	if (hasat2) difference = st_local("ratio")==""
	
	getcis = st_local("ci")!=""						//get confidence intervals
	percentiles = st_local("normal")==""		//calculate CIs using percentiles rather than normal approximation
	if (getcis) {
		M = strtoreal(st_local("m"))		//# of draws from MVN
		if (percentiles) {
			med = round(M:/2)
			lci = round(0.025:*M)
			uci = round(0.975:*M)
		}	
	}
	
	p = predictms_sepmod_sim_setup()
	
	//Core
	for (i=1;i<=Nstarts;i++) {
	
		res = predictms_sepmodel_core(p,1,from[i])		//inputs: pointer to filler program, index for at desig matrix (either 1 or 5 for at() or at2())
		
		if (hasat2) {
			
			res2 = predictms_sepmodel_core(p,5,from[i])

			if (!getcis) {
				if (st_local("ratio")=="") res = res:-res2
				else res = res:/res2
			}
			else {
				if (percentiles) {
					st_view(medians,.,tokens(st_local("fromvars"+strofreal(from[i]))),touse)
					st_view(lcis,.,tokens(st_local("fromvarslci"+strofreal(from[i]))),touse)
					st_view(ucis,.,tokens(st_local("fromvarsuci"+strofreal(from[i]))),touse)
					for (tr=1;tr<=Nstates;tr++) {
						if (difference) base = asarray(res,tr) :- asarray(res2,tr)
						else base = asarray(res,tr) :/ asarray(res2,tr)
						for (j=1;j<=obs;j++) {
							temp = sort(base[j,]',1)
							medians[j,tr] = temp[med,1]
							lcis[j,tr] = temp[lci,1]
							ucis[j,tr] = temp[uci,1]
						}
					}
				}
				else {
					for (tr=1;tr<=Nstates;tr++) {
						st_view(pred=.,.,tokens(st_local("probvars_"+strofreal(from[i])+"_"+strofreal(tr))),touse)
						if (difference) {
							base = asarray(res,tr) :- asarray(res2,tr)
							if (los) {
								for (j=1;j<=obs;j++) pred[j,] = quadmeanvariance(base[j,]')'
							}
							else {
								for (j=1;j<=obs;j++) pred[j,] = quadmeanvariance(atanh(base[j,])')'
							}
						}
						else {
							base = asarray(res,tr) :/ asarray(res2,tr)
							for (j=1;j<=obs;j++) pred[j,] = quadmeanvariance(log(base[j,])')'
						}
					}
				}
			}
		}
		else {
			//post results
			if (getcis) {
				if (percentiles) {
					med = round(M:/2)
					lci = round(0.025:*M)
					uci = round(0.975:*M)
					for (i=1;i<=Nstarts;i++) {	
						st_view(medians,.,tokens(st_local("fromvars"+strofreal(from[i]))),touse)
						st_view(lcis,.,tokens(st_local("fromvarslci"+strofreal(from[i]))),touse)
						st_view(ucis,.,tokens(st_local("fromvarsuci"+strofreal(from[i]))),touse)
						for (tr=1;tr<=Nstates;tr++) {
							base = asarray(res,tr)
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
							base = asarray(res,tr)
							if (los) {
								for (j=1;j<=obs;j++) pred[j,] = quadmeanvariance(log(base[j,])')'
							}
							else {
								for (j=1;j<=obs;j++) pred[j,] = quadmeanvariance(logit(base[j,])')'
							}
						}
					}
				}
			
			}
		}
		
		if (!getcis) st_store(.,tokens(st_local("fromvars"+strofreal(from[i]))),touse, res)
		
	}
	
}

`TR' predictms_sepmodel_core(	`PS' p,			///
							`RS' dmind,		///
							`RS' from)
{

	struct predictms_sep_struct scalar S
	S = *p

	if (S.getcis) {
		transmorphic A
		A = asarray_create("real",1)		//to hold predictions across trans and sims
	}

	// loop over simulations
	for (k=1;k<=S.M;k++) {
	
		pt = J(0,S.Nstates,.)				//plus 1 for current
		for (t=1;t<=S.Nobstosim;t++) {
			
			//all start in state from[i] at time 0 defined in enter()
			states = J(S.N,1,from)
			stimes = J(S.N,1,S.enter[t])

			ind = 1	//indexes the move
			
			for (j=from;j<=S.Nstates;j++) {
			
				//if there are possible next states
				if (S.Nnextstates[j]>0) {
				
					//update index vector
					index = uniqrows((states[,ind]:==j) :* S.coreindex)									//gets index for patients in current i th state
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
					newt = J(Nsim,S.Nnextstates[j,1],0)
					for (ij=1; ij<=S.Nnextstates[j,1]; ij++) {

						trans = asarray(S.postrans,j)[ij]
						
						if (S.modeltrans[trans,6]) {
							if (S.reset) {
								for (q=1;q<=Nsim;q++){
									//plus a bit takes care of censoring in state calc below
									rc = mm_root(x=1,&stpm2_sim_sep(),mint=smallestdouble(),(S.maxt:-(stimes[index,ind])[q]+0.01),tol=0,maxit=1000,runiform(1,1),asarray(S.transinfo,(trans,2))[k,]',dmind,asarray(S.transinfo,(trans,4))) 	//lower,upper,tol,maxit,U,coefficients,pointer to struct
									newt[q,ij] = x 
								}
							}
							else {
								for (q=1;q<=Nsim;q++){
									rc = mm_root(x=1,&stpm2_sim_sep(),mint=smallestdouble(),S.maxt+0.01,tol=1e-9,maxit=1000,runiform(1,1),asarray(S.transinfo,(trans,2))[k,]',dmind,asarray(S.transinfo,(trans,4)),(stimes[index,ind])[q]) 	//lower,upper,tol,maxit,U,coefficients,pointer to struct
									newt[q,ij] = x 
								}
							}
							continue							
						}
						if (S.modeltrans[trans,1]) {
							lambda = exp(asarray(S.transinfo,(trans,dmind)) * asarray(S.transinfo,(trans,2))[k,]')
							if (S.reset) newt[,ij] = log(runiform(Nsim,1)) :/ (-lambda)
							else newt[,ij] = log( exp(-lambda :* stimes[index,ind]) :* runiform(Nsim,1)) :/ (-lambda)
							continue
						}
						if (S.modeltrans[trans,2]) {
							indices = asarray(S.transinfo,(trans,3))
							lambda = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,1]::indices[2,1]] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,1]|])
							gamm = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,2]::indices[2,2]] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,2]|])
							if (S.reset) newt[,ij] = (log(runiform(Nsim,1)) :/ (-lambda)):^(1:/gamm)
							else newt[,ij] = ((log(exp(-lambda :* stimes[index,ind]:^gamm) :* runiform(Nsim,1)) :/ (-lambda)):^(1:/gamm))
							continue
						}
						if (S.modeltrans[trans,3]) {
							indices = asarray(S.transinfo,(trans,3))
							lambda = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,1]::indices[2,1]] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,1]|])
							gamm = asarray(S.transinfo,(trans,dmind))[,indices[1,2]::indices[2,2]] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,2]|]
							if (S.reset) newt[,ij] = log(((log(runiform(Nsim,1)) :/ (-lambda)):*gamm):+1):/gamm
							else newt[,ij] = log(((log( exp(-lambda:/gamm:*(exp(gamm:*stimes[index,ind]):-1)) :* runiform(Nsim,1)) :/ (-lambda)):*gamm):+1):/gamm
							continue
						}
						if (S.modeltrans[trans,4]) {
							indices = asarray(S.transinfo,(trans,3))
							lambda = exp(-asarray(S.transinfo,(trans,dmind))[,indices[1,1]::indices[2,1]] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,1]|])
							gamm = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,2]::indices[2,2]] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,2]|])
							if (S.reset) newt[,ij] = (((1:/runiform(Nsim,1)):-1):^(gamm)):/lambda
							else newt[,ij] = (((1:/(((1 :+ (lambda:*stimes[index,ind]):^(1:/gamm)):^(-1)):*runiform(Nsim,1))):-1):^(gamm)):/lambda
							continue
						}
						if (S.modeltrans[trans,5]) {
							lambda = asarray(S.transinfo,(trans,dmind))[,indices[1,1]::indices[2,1]] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,1]|]
							gamm = exp(asarray(S.transinfo,(trans,dmind))[,indices[1,2]::indices[2,2]] * (asarray(S.transinfo,(trans,2))[k,]')[|indices[,2]|])
							newt[,ij] = exp(invnormal(1:-runiform(Nsim,1)):*gamm :+ lambda)
						}
						
					}

					//calculate new time and new state
					if (S.reset) newt = newt,(S.maxt:-stimes[index,ind])
					else newt = newt,J(Nsim,1,S.maxt)
					nextt = rowmin(newt)
					if (S.reset) stimes[index,ind+1] = nextt :+ stimes[index,ind]
					else stimes[index,ind+1] = nextt
					states[index,ind+1] = (nextt :== newt) * asarray(S.posnextstatesj,j)			
					ind++
				}
			}

			if (!S.isenter) {
				if (S.los) {
				
				}
				else {
					//calculate probability of being in each state
					pt = pt\(colsum(states[,cols(states)] :== S.statematrix) :/ S.N)
				}
			}
			
		}
	
	
		if (S.isenter) {
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
					temppred = J(1,cols(states),0)
					for (s=1; s<=cols(states); s++) {
						temppred[1,s] = sum((states:==s) :* temptime):/S.N
					}
					pt = pt\temppred
				}

			}
			else {
				//calculate probability of being in each state
				for (t=1;t<S.obs;t++) {
					tempstate = states[,1]
					for (s=2; s<=cols(states); s++) {
						index = selectindex(stimes[,s] :<= S.predtime[t])
						tempstate[index] = states[index,s]
					}
					pt = pt\(colsum(tempstate :== S.statematrix) :/ S.N)
				}
				pt = pt\(colsum(states[,cols(states)] :== S.statematrix) :/ S.N)
			}
		}

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

