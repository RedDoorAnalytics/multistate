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
local SC	string colvector
local SS	struct predictms_struct scalar
local PS	pointer(struct predictms_struct scalar) scalar
local PC	pointer(struct predictms_struct scalar) colvector
local Ps	pointer scalar

mata:

struct predictms_struct
{
	`ss' touse		//to post predictions
	`RS' method		//0 simulation, 1 aj, 2 analytic illd, ...
	`RS' nicode		//1 cr - 2 illd...
	`RM' transmat		// transition matrix
	`RS' hasmodels		//models or stacked
	`TR' probs		//to store trans. probs.
	`TR' loss		//to store length of stays
	`TR' rmsts		//store rmsts
	`TR' visits		//store visit probs
	`TR' hazards		//store hazards
	`TR' survivals		//store survivals
	`RS' getprobs		//get trans. probs.
	`RS' getlos		//get length of stays
	`RS' getrmst		//get rmst
	`RS' getvisit		//get prob visit
	`RS' gethazard		//get hazards
	`RS' getsurvival	//get survivals
	`RS' M
	`RS' getcis
	`RS' cimethod		//delta method or bootstrap
	`RS' Nbs		//number of (total) betas for DM
	`RM' VCV		//needed for DM + contrasts
	`RS' reset
	`RS' obs
	`RS' N 
	`RS' maxt
	`RS' Nobstosim
	`RS' Ntrans
	`RS' Nstates
	`RS' Nstarts
	`SC' toupdate
	`RS' K
	`RS' Kind		//contains k index, gets updated
	`RS' at			//at index, gets updated
	`RS' std 		//std index, gets updated
	`RS' standardise
	`RS' survsim
	`RS' Nats
	`RC' predtime
	`RC' enter
	`RC' coreindex
	`RC' Nnextstates
	`RM' statematrix
	`RM' indices
	`RM' tscale2
	`TR' transinfo
	`TR' posnextstates
	`TR' posnextstatesj
	`TR' postrans
	`RC' time2
	`RC' tsreset
	`RS' percentiles
	`RS' med, lci, uci
	`RS' getdiffs, getratios
	`RC' novcv		//flag transitions to ignore VCV in CI calculation
	`RS' tosave		//flag to save simulated datasets
	`ss' savestub		//save dataset stub
	`RS' savereplace	//, replace
	`RC' iscox					
	
	`RM' pt			//stores probabilities (gets updated)
	`RM' los		//stores length of stays (gets updated)
	`RM' rmst		//stores rmst
	`RM' visit		//stores visit probs
	`RM' hazard		//stores hazards
	`RM' survival		//stores survivals
	`RS' level		//CI level
	`RS' atref		//reference for diffs and ratios
	
	//stuff for userfunction
	`RS' hasuser			//logical
	`RS' Nuservars, Nuserflag	//# of vars that need to be posted
	`Ps' userfunc			//pointer to user function
	`RM' user			//store predictions
	`TR' users			//store ci predictions
	`RS' userlink			//user link function
	
	//for AJ
	`RM' varpt
	`TR' varprobs
	`RS' varaj
}

void predictms_setup(`SS' S)
{
	`RS' model
	`RM' transmat
	
	transmat 	= st_matrix(st_local("transmatrix"))
	S.Ntrans 	= strtoreal(st_local("Ntrans"))
	S.Nstates 	= cols(transmat)
	S.reset 	= st_local("reset")!=""					//default clock forward, o/w reset
	S.touse 	= st_local("touse")					//to post predictions
	S.survsim 	= st_local("survsim")!=""

	S.hasmodels	= st_local("models")!=""
	S.iscox		= (tokens(st_local("familys")):=="cox")[1]
	if (S.iscox) {
		if (st_local("latent")!="") {
			errprintf("latent not allowed with Cox models\n")
			exit(198)
		}
		if (st_local("aj")!="") {
			errprintf("aj not allowed with Cox models\n")
			exit(198)
		}
		if (st_local("ci")!="") {
			errprintf("ci not allowed with Cox models\n")
			exit(198)
		}
		if (st_local("hazard")!="" | st_local("survival")!="") {
			errprintf("hazard/survival prediction not allowed with Cox models\n")
			exit(198)
		}
                if (st_local("visit")!="") {
			errprintf("prediction type visit not supported with Cox models\n")
			exit(198)
		}
	}
	S.method 	= predictms_setup_get_method(transmat,S)
	if (S.method==1 | S.method==2) {
		if (st_local("n")!="") {
			errprintf("n() not required\n")
			exit(198)
		}
	}
	if (S.method==1) S.transmat = transmat
	
	S.obs 		= strtoreal(st_local("obs"))			//# of time points to calculate probs at
	if (S.method==0 | S.method==3) {
		if (st_local("n")=="") {
			if (st_local("ci")=="") S.N = 100000
			else			S.N = 10000
		}
		else S.N 	= strtoreal(st_local("n"))				//sample size for each simulation
	}
	
	//setup things to store results in
	S.getprobs  = st_local("probability")!=""
	if (S.getprobs) {
		S.probs 	= asarray_create("real",3)				//store results
		S.pt 		= J(S.obs,S.Nstates,.)					//gets updated
	}
	
	S.getlos 	= st_local("los")!=""
        S.getrmst 	= st_local("rmst")!=""
	if (S.getlos | S.getrmst) {
		S.loss 	= asarray_create("real",3)
		S.los 	= J(S.obs,S.Nstates,.)
	}
	if (S.getrmst) {
		S.rmsts = asarray_create("real",3)
		S.rmst 	= J(S.obs,1,.)
	}
	
	S.getvisit 	= st_local("visit")!=""
	if (S.getvisit) {
		S.visits 	= asarray_create("real",3)
		S.visit 	= J(S.obs,S.Nstates,.)
	}
	
	S.gethazard 	= st_local("hazard")!=""
	if (S.gethazard) {
		S.hazards 	= asarray_create("real",3)
	}
	
	S.getsurvival 	= st_local("survival")!=""
	if (S.getsurvival) {
		S.survivals = asarray_create("real",3)
	}
	
	S.hasuser	= st_local("userfunction")!=""
	if (S.hasuser) {
		stata("mata: pf = &"+st_local("userfunction")+"()")
		external pf
		S.userfunc = pf
		S.users = asarray_create("real",3)
		S.Nuserflag = 1
		if (st_local("userlink")=="") S.userlink = 0 
		else {
			userlink = st_local("userlink")
			if (userlink=="identity") 	S.userlink = 0 
			else if (userlink=="logit") S.userlink = 1
			else if (userlink=="log") 	S.userlink = 2
			else {
				errprintf("Unknown userlink()\n")
				exit(198)
			}
		}
	}
	
	S.toupdate 		= uniqrows(tokens(st_local("toupdate"))')
	S.standardise 	= st_local("standardise")!=""
	
	S.getcis = st_local("ci")!=""						//get confidence intervals
	predictms_novcv(S)
	if (S.getcis) {
		S.cimethod		= predictms_setup_get_ci_method(S)
		S.level 		= strtoreal(st_local("level"))
		S.percentiles 	= st_local("percentile")!=""
		if (S.cimethod) { 
			if (st_local("m")=="") 	S.M = 200
			else 					S.M = strtoreal(st_local("m"))					//# of draws from MVN
			S.varaj = st_local("varaj")!=""
			if (!S.varaj) {
				if (S.percentiles) {
					S.med = round(S.M:/2)
					sp = (100-S.level):/200
					S.lci = round(sp:*S.M)
					S.uci = round((1-sp):*S.M)
				}
			}
			else S.varprobs = asarray_create("real",2)
		}
		else {
			if (st_local("m")!="") {
				errprintf("m() only required when using bootstrapping for confidence intervals\n")
				exit(198)
			}
		}
	}
	
	S.tosave = st_local("save")!=""
	S.savestub = st_local("savestub")
	S.savereplace = st_local("savereplace")!=""
	
	//design matrix stuff
	S.Nats 			= strtoreal(st_local("Nats"))			
	S.K 			= strtoreal(st_local("K"))	//Number of design matrices to loop over; 1 for normal, # patients for standardised predictions
	S.Kind  		= 0
	
	//contrasts
	S.getdiffs	= st_local("difference")!=""
	S.getratios = st_local("ratio")!=""
	
	if (S.getdiffs | S.getratios) S.atref = strtoreal(st_local("atref"))
	
	//transition-specific resets
	if (st_local("tsreset")!="") {
		if (S.reset) {
			errprintf("reset and tsreset() can't both be used\n")
			exit(198)
		}
		tsreset = strtoreal(tokens(st_local("tsreset"))')
		S.tsreset = J(S.Ntrans,1,0)
		for (i=1;i<=S.Ntrans;i++) {
			for (j=1;j<=rows(tsreset);j++) {
				if (tsreset[j]==i) S.tsreset[i] = 1
			}
		}
		
		nrows = rows(S.tsreset)
		if (nrows>S.Ntrans) {
			errprintf("Number of tsreset() elements can't be more than number of transitions\n")
			exit(198)
		}
		//need error check for transition number checks
		
	}
	else {
		if (S.reset) 	S.tsreset = J(S.Ntrans,1,1)
		else 		S.tsreset = J(S.Ntrans,1,0)
	}
	
	//second timescale
	S.tscale2 	= J(S.Ntrans,1,0)
	S.time2		= J(S.Nats,1,0)
	if (st_local("tscale2")!="") {
		tscale2 		= strtoreal(tokens(st_local("tscale2"))')
		time2			= strtoreal(tokens(st_local("time2")))'
		if (sum(time2:==.) | sum(time2:<0)) {
			errprintf("time2() must be real > 0\n")
			exit(198)
		}
		nrows = rows(time2)
		if (nrows>1 & nrows!=S.Nats) {
			errprintf("The number of time2()s must be the same number of at#()\n")
			exit(198)
		}
		if (nrows==1 & S.Nats>1) time2 = J(S.Nats,1,time2)
		S.time2		= time2
		Nt2			= rows(tscale2)
		for (i=1;i<=Nt2;i++) S.tscale2[tscale2[i]] 	= 1
	}
	
	//time
	S.predtime 	= st_data(.,st_local("timevar"),st_local("touse"))
	S.maxt 		= max(S.predtime)
	S.enter 	= strtoreal(st_local("ltruncated"))
	S.Nobstosim = 1
	
	//transition-specific info
	
	if (S.hasmodels) {
		
		S.transinfo = asarray_create("real",2)
		//(i,1) - parameter vector
		//(i,2) - vcv or draws from MVN if cis requested
	
		for (i=1;i<=S.Ntrans;i++) {
		
			//coefficient vector
			bmat = st_matrix(st_local("emat"+strofreal(i)))
			asarray(S.transinfo,(i,1),bmat)
			
			if (S.getcis) {
				if (!S.novcv[i]) {												//must be nested else not defined
					vcvmat = st_matrix(st_local("evmat"+strofreal(i)))
					if (S.cimethod) 	asarray(S.transinfo,(i,2),predictms_mvrnorm_antithetic(bmat',vcvmat,S.M))
					else			asarray(S.transinfo,(i,2),vcvmat)
				}
			}
			
		}
		
	}
	else {
		
		S.transinfo = asarray_create("real",1)
		//(1) - parameter vector
		//(2) - vcv or draws from MVN if cis requested
		
		bmat = st_matrix(st_local("emat"))
		asarray(S.transinfo,(1),bmat)
		if (S.getcis) {
			vcvmat = st_matrix(st_local("evmat"))
			if (S.cimethod) 	asarray(S.transinfo,2,predictms_mvrnorm_antithetic(bmat',vcvmat,S.M))
			else 			asarray(S.transinfo,2,vcvmat)
		}

	}

	//extra info for core function
	if (!S.method | S.method==3) {
		S.coreindex 		= 1::S.N
		S.statematrix 		= J(S.N,1,1..S.Nstates)
	}
	S.Nnextstates 		= rowsum(transmat:!=.)			//# of possible next states from each current state
	S.posnextstates 	= asarray_create("real",1)		//possible next states
	S.posnextstatesj 	= asarray_create("real",1)		//possible next states and current
	S.postrans 		= asarray_create("real",1)

	for (i=1;i<=S.Nstates;i++) {
		asarray(S.posnextstates,i,strtoreal(tokens(st_local("row"+strofreal(i)+"next")))')
		asarray(S.posnextstatesj,i,(asarray(S.posnextstates,i)\i))									
		asarray(S.postrans,i,strtoreal(tokens(st_local("row"+strofreal(i)+"trans")))')
	}

}

void predictms_novcv(`SS' S)
{
	novcv 	= strtoreal(tokens(st_local("novcv")))'
	Nnovcv 	= rows(novcv)
	
	novcvflag = J(S.Ntrans,1,0)
	for (i=1;i<=S.Ntrans;i++) {
		for (j=1;j<=Nnovcv;j++) {
			if (novcv[j]==i) novcvflag[i] = 1
		}
	}
	S.novcv = novcvflag
}

`RM' predictms_mvrnorm_antithetic( `RC' mvec,		///
                                `RM' V,			///
                                `RS' n)	
{
	`RS' nvars
	`RM' cholV, z, res
	
	nvars = rows(mvec)
	if (nvars!=rows(V) & nvars!=cols(V)) {
		errprintf("Dimensions not valid\n")
		exit(198)
	}
	cholV = cholesky(V)
	z1 = runiform(nvars,n/2)
	z2 = z1,(1:-z1)
	z = invnormal(z2)

	//draw via x = mu + A*z
	res = J(nvars,n,.)
	for (i=1; i<=n; i++) res[.,i] = mvec + cholV*z[., i]
        res
	return(res')
}

/*
get prediction method:

0 - simulation using multinomial
1 - AJ
2 - numerical integration
3 - simulation with latent times
*/

`RS' predictms_setup_get_method(`RM' transmat, `SS' S)
{
	if (S.iscox) 				return(0)
	
	if (st_local("aj")!="") 		return(1)	//AJ estimator
	
	if (st_local("simulate")!="") {
		if (st_local("latent")!="")	return(3)	//simulation using latent times
		else				return(0)	//simulate - direct
	}
	
	if (st_local("reset")!="" | st_local("tsreset")!="") {
		if (st_local("latent")!="")     return(3)
		else 			        return(0)
	}
	
	if (st_local("singleevent")!="" | transmat==(.,1\.,.)) {
		S.nicode = 1	
		return(2)	//ni
	}
	
	if (st_local("cr")!="") {
		S.nicode = 2	//cr
		return(2)	//ni
	}
	
	Nrows = rows(transmat)
	if (sum(transmat[|2,.\Nrows,.|])==0) {
		S.nicode = 2	//cr
		return(2)	//ni
	}

	if (transmat==(.,1,2\.,.,3\.,.,.)) {
		S.nicode = 3	//illd
		return(2)	//ni
	}

	if (transmat==(.,1,2,.\.,.,.,3\.,.,.,.\.,.,.,.)) {
		S.nicode = 4	//extended illd
		return(2)	//ni
	}
	
	if (st_local("latent")!="")	return(3)
	else 				return(0)	
}

/*
get ci method:

0 - delta method
1 - bootstrap
*/

`RS' predictms_setup_get_ci_method(`SS' S)
{
	if 	(st_local("bootstrap")!="")	return(1)
	if 	(S.method==0 | S.method==3) 	return(1)
	else 	        			return(0)
}

end
