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
local PS	pointer(struct predictms_struct scalar) scalar

mata:

struct predictms_struct
{
	`RS' M, getcis, reset, obs, N, isenter, maxt, Nobstosim, Nstates, Nstarts, los, K, isstd, survsim, prob
	`RC' predtime, enter, coreindex, Nnextstates, from
	`NM' modeltrans, statematrix, bdraws, indices
	`TR' transinfo, posnextstates, posnextstatesj, postrans, Xstd
	pointer(struct predictms_stpm2 scalar) colvector pstpm2
}

`PS' predictms_setup(`RS' model)
{
	struct predictms_struct scalar S
	pointer scalar p

	stata("tempname predictms_struct")
	rmexternal(st_local("predictms_struct"))
	p = crexternal(st_local("predictms_struct"))

	S.getcis = st_local("ci")!=""					//get confidence intervals
	if (S.getcis) S.M = strtoreal(st_local("m"))	//# of draws from MVN
	else S.M = 1
	
	S.reset = st_local("reset")!=""			//default clock forward, o/w reset
	S.obs = strtoreal(st_local("obs"))		//# of time points to calculate probs at
	S.N = strtoreal(st_local("n"))			//sample size
	hasat2 = st_local("at2")!=""
	S.los = st_local("los")!=""
	S.K = strtoreal(st_local("K"))			//Number of design matrices to loop over; 1 for normal, # patients for standardised predictions
	S.isstd = S.K>1
	Ntrans = strtoreal(st_local("Ntrans"))
	S.survsim = st_local("survsim")!=""
	S.prob = st_local("prob")!=""
	
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

	S.modeltrans = (tokens(st_local("cmds"))'):==J(Ntrans,1,("ereg","weibull","gompertz","llogistic","lnormal","stpm2","strcs"))
	S.transinfo = asarray_create("real",2)

	//transition-specific info
	for (i=1;i<=Ntrans;i++) {
		
		//standardised predictions or standard design matrix
		if (!S.modeltrans[i,6] & !S.modeltrans[i,7]) {
			DM = st_matrix(st_local("dm"+strofreal(i)))
			if (S.K>1) {
				DM = J(S.K,1,DM)
				DM[,strtoreal(tokens(st_local("stdvarsindex"+strofreal(i))))] = st_data(.,tokens(st_local("stdvars"+strofreal(i))),st_local("stdtouse"))
			}
			asarray(S.transinfo,(i,1),DM)													//design matrix
		}
		//coefficient vector
		if (S.getcis) asarray(S.transinfo,(i,2),st_data(.,tokens(st_local("drawvars"+strofreal(i))),st_local("mvnind")))	//parameter draws
		else asarray(S.transinfo,(i,2),st_matrix(st_local("emat"+strofreal(i))))											//parameter draws
		//ml eqn indices
		if (!S.modeltrans[i,6] & !S.modeltrans[i,1] & !S.modeltrans[i,7]) {
			asarray(S.transinfo,(i,3),st_matrix(st_local("indices"+strofreal(i))))
		}
		//stpm2 separate model struct build
		if (S.modeltrans[i,6]) {
			if (model) stata("qui estimates restore "+st_local("modelests"+strofreal(i)))
			asarray(S.transinfo,(i,4),predictms_stpm2_setup(strofreal(i)))
		}
		//strcs separate model struct build
		if (S.modeltrans[i,7]) {
			if (model) stata("qui estimates restore "+st_local("modelests"+strofreal(i)))
			asarray(S.transinfo,(i,4),predictms_strcs_setup(strofreal(i)))
		}
		//at2() DMs for non-stpm2 models
		if (hasat2 & !S.modeltrans[i,6]) {
			DM = st_matrix(st_local("at2dm"+strofreal(i)))
			if (S.K>1) {
				DM = J(S.K,1,DM)
				DM[,strtoreal(tokens(st_local("at2stdvarsindex"+strofreal(i))))] = st_data(.,tokens(st_local("at2stdvars"+strofreal(i))),st_local("stdtouse"))
			}	
			asarray(S.transinfo,(i,5),DM)
		}
		
	}

	S.coreindex = 1::S.N
	transmat = st_matrix(st_local("transmatrix"))
	S.Nstates = cols(transmat)
	S.statematrix = J(S.N,1,1..S.Nstates)
	S.Nnextstates = rowsum(transmat:!=.)															//# of possible next states from each current state
	S.posnextstates = S.posnextstatesj = S.postrans = asarray_create("real",1)

	for (i=1;i<=S.Nstates;i++) {
		asarray(S.posnextstates,i,strtoreal(tokens(st_local("row"+strofreal(i)+"next")))')
		asarray(S.posnextstatesj,i,(asarray(S.posnextstates,i)\i))									//possible next states includng current
		asarray(S.postrans,i,strtoreal(tokens(st_local("row"+strofreal(i)+"trans")))')
	}

	swap((*p), S)
	return(p)
}


void predictms()
{
	`RS' model, hasat2, getcis, percentiles, difference, probs
	`TR' res
	`SS' touse
	`PS' p
	
	model = strtoreal(st_local("model"))
	touse = st_local("touse")									//touse for final predictions posting
	from = strtoreal(tokens(st_local("from")))'
	Nstarts = rows(from)
	Nstates = strtoreal(st_local("Nstates"))
	obs = strtoreal(st_local("obs"))					//<- and bove could be passed aa arguments to P_predictms_sepmod_sim_setup() instead of reading them in again within it	
	hasat2 = st_local("at2")!=""
	los = st_local("los")!=""
	if (hasat2) difference = st_local("ratio")==""
	probs = st_local("prob")!=""
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
	
	p = predictms_setup(model)
	
	//Core
	for (i=1;i<=Nstarts;i++) {
	
		res = predictms_core(p,1,from[i])

		if (hasat2) {

			res2 = predictms_core(p,5,from[i])

			if (!getcis) {
				if (probs) {
					nvars = ((*p).Nstates-1)*2
					res = predictms_prob(p,from[i],res,res2)
					postind = 1
					for (j=2;j<=(*p).Nstates;j++) {
						(void) st_addvar("double","prob"+strofreal(j)+"_at")
						st_store(.,"prob"+strofreal(j)+"_at",st_local("touse"),res[,postind])
						postind++
						(void) st_addvar("double","prob"+strofreal(j)+"_at2")
						st_store(.,"prob"+strofreal(j)+"_at2",st_local("touse"),res[,postind])
						postind++
					}
					exit()
				}
				else {
					if (st_local("ratio")=="") res = res:-res2
					else res = res:/res2	
				}
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

end

