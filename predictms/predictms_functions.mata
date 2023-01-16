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

`RM' predictms_gq(`RS' n)
{
	i = range(1,n,1)'
	i1 = range(1,n-1,1)'

	muzero = 2
	a = J(1,n,0)
	b = i1:/sqrt(4 :* i1:^2 :- 1)

	A= diag(a)
	for(j=1;j<=n-1;j++) {
		A[j,j+1] = b[j]
		A[j+1,j] = b[j]
	}
	symeigensystem(A,vec,nodes)
	weights = (vec[1,]:^2:*muzero)'
	weights = weights[order(nodes',1)]
	nodes = nodes'[order(nodes',1)']
	return(nodes,weights)
}

//userfunction() utilities

`RM' ms_user_prob(`SS' S, `RS' b)
{
	return(S.pt[,b])
}

`RM' ms_user_los(`SS' S, `RS' b)
{
	return(S.los[,b])
}

`RM' predictms_get_b(`SS' S, `RS' trans)
{
	if (S.hasmodels) {
		if (!S.Kind | S.novcv[trans]) 	b = asarray(S.transinfo,(trans,1))
		else 							b = asarray(S.transinfo,(trans,2))[S.Kind,]
	}
	else {
		if (!S.Kind) 					b = asarray(S.transinfo,(1))
		else 							b = asarray(S.transinfo,(2))[S.Kind,]
	}
	return(b)
}

void predictms_init_storage(`SS' S,`RS' from)
{
	if (S.getprobs) 	S.pt[,] 	= J(S.obs,S.Nstates,0)
	if (S.getlos) 	 	S.los[,] 	= J(S.obs,S.Nstates,0)
	if (S.getrmst) 	 	S.rmst[,] 	= J(S.obs,1,0)
	if (S.hasuser) 	 	S.user[,] 	= J(S.obs,S.Nuservars,0)
	if (S.getvisit)  	S.visit[,] 	= J(S.obs,S.Nstates,0)
	if (S.gethazard) 	S.hazard	= J(S.obs,S.Nnextstates[from],0)
	if (S.getsurvival) 	S.survival	= J(S.obs,S.Nnextstates[from],0)
}

//update new data draws for right censoring observations
//-> s0, starting state
//-> maxt, right censoring time

void update_rcens(`RM' newdata, `RS' s0, `RS' maxt)
{
	censindex 	= selectindex(newdata[,1]:==maxt)
	Ncens		= rows_cols(censindex)
	if (Ncens[1] & Ncens[2]) {
		newdata[censindex,2] = J(Ncens[1],1,s0)		//same state
		newdata[censindex,3] = J(Ncens[1],1,1)		//now done
	}	
}

//save simulated times and states into a Stata dataset

void predictms_sim_save(`SS' S, `RM' stimes, `RM' states)
{
	stata("preserve")
	stata("clear")
	stata("cap set obs "+strofreal(S.N))
	stata("tempvar simtouse")
	stata("qui gen byte "+st_local("simtouse")+"= _n<="+strofreal(S.N))
	Nmoves = cols(stimes)
	newstimevars = newstatevars = J(1,0,"")
	for (i=1;i<=Nmoves;i++) {
		newstimevars = newstimevars,"_time"+strofreal(i)
		newstatevars = newstatevars,"_state"+strofreal(i)		
	}
	id1 = st_addvar("int", newstatevars)
	id2 = st_addvar("double", newstimevars)
	st_store(.,id1,st_local("simtouse"),states)
	st_store(.,id2,st_local("simtouse"),stimes)
	if (S.savereplace) {
		if (!S.getcis) 	stata("qui save " + S.savestub +", replace")
		else 			stata("qui save " + S.savestub + strofreal(S.Kind)+", replace")
	}
	else {
		if (!S.getcis) 	stata("qui save " + S.savestub)
		else 			stata("qui save " + S.savestub + strofreal(S.Kind))
	}
	stata("restore")
}

void predictms_survsim_post(`RC' t)
{
	st_store(.,st_local("survsim"),st_local("survsimtouse"),t)
	exit()
}

void predictms_calc_prob(`SS' S, `RM' stimes, `RM' states)
{
	Ntotstates = cols(states)
	for (t=1;t<S.obs;t++) {
		tempstate = states[,1]
		for (s=2; s<=Ntotstates; s++) {
			index = selectindex(stimes[,s] :<= S.predtime[t])
			if (rows(index) & cols(index)) tempstate[index] = states[index,s]
		}
		if (S.standardise) 	S.pt[t,] = S.pt[t,] :+ colsum(tempstate :== S.statematrix) :/S.N
		else 			S.pt[t,] = colsum(tempstate :== S.statematrix) :/ S.N
	}
	if (S.standardise) 	S.pt[t,] = S.pt[t,] :+ colsum(states[,cols(states)] :== S.statematrix) :/S.N
	else 			S.pt[t,] = colsum(states[,cols(states)] :== S.statematrix) :/S.N
}

void predictms_calc_los(`SS' S, `RM' stimes, `RM' states)
{
	Ntotstates = cols(states)
	for (t=1;t<=S.obs;t++) {
		temptime = stimes
		for (s=2; s<=Ntotstates; s++) {
			index = selectindex(stimes[,s] :> S.predtime[t])
			if (rows(index) & cols(index)) temptime[index,s] = J(rows(index),1,S.predtime[t])
		}						
		temptime = temptime,J(S.N,1,S.predtime[t])
		temptime = temptime[|.,2\.,cols(temptime)|] :- temptime[|.,1\.,cols(states)|] 
		
		temppred = J(1,S.Nstates,0)
		for (s=1; s<=S.Nstates; s++) {
			temppred[1,s] = sum((states:==s) :* temptime):/S.N
		}
		if (S.standardise) 	S.los[t,] = S.los[t,] :+ temppred
		else 				S.los[t,] = temppred
	}
}

void predictms_calc_rmst(`SS' S, `RM' stimes, `RM' states)
{
	rmst = J(S.obs,1,0)
	for (s=1;s<=S.Nstates;s++) {
		if (S.Nnextstates[s]) rmst = rmst :+ S.los[,s]
	}
	if (S.standardise) 	S.rmst = S.rmst :+ rmst
	else 				S.rmst = rmst
}

void predictms_calc_visit(`SS' S, `RM' stimes, `RM' states)
{
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
		if (S.standardise) 	S.visit[t,] = S.visit[t,] :+ temppred
		else 				S.visit[t,] = temppred
	}
}

void predictms_calc_user(`SS' S)
{
	S.user = (*S.userfunc)(S)
	if (S.Nuserflag) {
		S.Nuservars = cols(S.user)	
		S.Nuserflag = 0
	}
}

end
