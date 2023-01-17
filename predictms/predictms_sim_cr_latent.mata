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
local Pm	pointer(struct merlin_struct scalar) scalar

mata:

`RM' predictms_sim_cr_latent(			///
                                `SS' S, 	/// predictms object
                                `RS' Nsim, 	/// # of new times to simulate
                                `RC' t0, 	/// entry time
                                `RS' s0		/// entry state
                                )
{
	`Pm' Pmerlin

	newdata		= J(Nsim,3,0)			//time,state,done
	Ntrans 		= S.Nnextstates[s0]
	newtimes 	= J(Nsim,Ntrans,.)
	
	for (s=1; s<=Ntrans; s++) {

		trans 	= asarray(S.postrans,s0)[s]
		b	= predictms_get_b(S,trans)
		logU	= log(runiform(Nsim,1))
		Pmerlin = predictms_merlin_setup(S,b,Nsim,trans,t0)

		if ((*Pmerlin).NI) {
			rc = predictms_sim_root(	
                                t=J(Nsim,1,.),		///
                                &merlin_root_ni(),	///
                                S.maxt,			///
                                0,1000,			///
                                Nsim,			///
                                t0copy=t0,		/// -> gets updated so must pass copy as it's recalled
                                logU,			///
                                1::Nsim,		///
                                Pmerlin,		///
                                S.tsreset[trans],	///
                                S.tscale2[trans],	///
                                S.time2[S.at])
		}
		else {
			if 	((*Pmerlin).familys=="exponential") {
				if (S.tsreset[trans]) {
					if (S.tscale2[trans]) 	t = merlin_simulate_exp((*Pmerlin),logU,S.maxt,J(Nobs,1,S.time2[S.at])) :- S.time2[S.at]
					else 					t = merlin_simulate_exp((*Pmerlin),logU,S.maxt)
				}
				else {
					if (S.tscale2[trans])	t = merlin_simulate_exp((*Pmerlin),logU,S.maxt,t0:+S.time2[S.at]) :- S.time2[S.at]
					else					t = merlin_simulate_exp((*Pmerlin),logU,S.maxt,t0)
				}
			}
			else if ((*Pmerlin).familys=="weibull") {
				if (S.tsreset[trans]) {
					if (S.tscale2[trans]) 	t = merlin_simulate_weibull((*Pmerlin),logU,S.maxt,J(Nobs,1,S.time2[S.at])) :- S.time2[S.at]
					else 					t = merlin_simulate_weibull((*Pmerlin),logU,S.maxt)
				}
				else {
					if (S.tscale2[trans])	t = merlin_simulate_weibull((*Pmerlin),logU,S.maxt,t0:+S.time2[S.at]) :- S.time2[S.at]
					else					t = merlin_simulate_weibull((*Pmerlin),logU,S.maxt,t0)
				}
			}
			else if ((*Pmerlin).familys=="gompertz") {
				if (S.tsreset[trans]) {
					if (S.tscale2[trans]) 	t = merlin_simulate_weibull((*Pmerlin),logU,S.maxt,J(Nobs,1,S.time2[S.at])) :- S.time2[S.at]
					else 					t = merlin_simulate_weibull((*Pmerlin),logU,S.maxt)
				}
				else {
					if (S.tscale2[trans])	t = merlin_simulate_gompertz((*Pmerlin),logU,S.maxt,t0:+S.time2[S.at]) :- S.time2[S.at]
					else					t = merlin_simulate_gompertz((*Pmerlin),logU,S.maxt,t0)
				}
			}
			else {
				rc = predictms_sim_root(	///
                                        t=J(Nsim,1,.),		///
                                        &merlin_root_ch(),	///
                                        S.maxt,			///
                                        0,1000,			///
                                        Nsim,			///
                                        t0copy=t0,		///	-> gets updated so must pass copy as it's recalled
                                        logU,			///
                                        1::Nsim,		///
                                        Pmerlin,		///
                                        S.tsreset[trans],	///
                                        S.tscale2[trans],	///
                                        S.time2[S.at])
			}
		}
		newtimes[,s] = t
		rmexternal(st_local("GML"+strofreal(trans))) 	
	}
	
	newdata[,1] = rowmin(newtimes)
	update_rcens(newdata,s0,S.maxt)	

	//new states
	eventindex 	= selectindex(newdata[,1]:<S.maxt)
	Nevents		= rows_cols(eventindex)
	if (Nevents[1] & Nevents[2]) {
		
		nextstates = asarray(S.posnextstates,s0)
		if (Ntrans==1) 	newdata[eventindex,2]	= J(Nevents[1],1,nextstates)
		else 			newdata[eventindex,2] 	= (newdata[eventindex,1]:==newtimes[eventindex,]) * nextstates

		//obs now in absorbing state are done
		for (s=1; s<=Ntrans; s++) {
			if (!S.Nnextstates[nextstates[s]]) {
				doneindex 	= select(eventindex,newdata[eventindex,2]:==nextstates[s])
				Ndone 		= rows_cols(doneindex)
				if (Ndone[1] & Ndone[2]) newdata[doneindex,3] = J(Ndone[1],1,1)
			}
		}
	}

	return(newdata)
}

function merlin_root_ni(`RC' x, `RC' t0, `RC' logU, `RC' index, `Pm' p,  ///
						`RS' reset, `RS' tscale2, `RS' time2)
{
	struct merlin_struct gml
	
	Nobs 	= rows(x)
	ch		= J(Nobs,1,0)	
	gml 	= *p
	
	//index etc. gets updated as mm_root finishes observations
	gml.N = gml.Nobs = Nobs	
	asarray(gml.xbindex,1,index)
		
	Nq 	= strtoreal(st_local("chintpoints"))
	qpw 	= predictms_gq(Nq)
	qp 		= ((x:-t0) :* J(Nobs,1,qpw[,1]') :+ x:+t0) :/ 2
	if (reset) 	qp = qp :- t0
	if (tscale2) qp = qp :+ time2
	
	chq		= J(Nobs,Nq,.)
	for (q=1; q<=Nq; q++) {
		chq[,q] = exp((*gml.Plogh)(gml,qp[,q]))
	}
	_editmissing(chq,0)		
	ch = ch :+ (x:-t0):/2 :* chq * qpw[,2]
	
	return(ch :+ logU)	
}

function merlin_root_ch(`RC' x, `RC' t0, `RC' logU, `RC' index, `Pm' p,  ///
						`RS' reset, `RS' tscale2, `RS' time2)
{
	struct merlin_struct gml
	
	Nobs 	= rows(x)
	ch		= J(Nobs,1,0)	
	gml 	= *p

	//index etc. gets updated as mm_root finishes observations
	gml.N = gml.Nobs = Nobs	
	asarray(gml.xbindex,1,index)

	mainx	= x
	if (reset) 		mainx = mainx :- t0
	if (tscale2) 	mainx = mainx :+ time2
	
	chq 	=  (*gml.Pch)(gml,mainx)
	_editmissing(chq,0)

	if (reset) {
		if (tscale2) chq0 = (*gml.Pch)(gml,J(Nobs,1,time2))
		else {
			ch = ch :+ chq
			return(ch :+ logU)
		}
	}
	else {
		if (tscale2) chq0 = (*gml.Pch)(gml,t0:+time2)
		else 				chq0 = (*gml.Pch)(gml,t0)
	}
	_editmissing(chq0,0)	
	ch = ch :+ chq :- chq0

	return(ch :+ logU)
}

`RC' merlin_simulate_exp(`gml' gml, `RC' logU, `RS' maxt, | `RC' t0)
{
	if (args()==4) 	t = t0 :- logU:/exp(merlin_util_xzb(gml))
	else 			t = -logU:/exp(merlin_util_xzb(gml))
	index = selectindex(t:>maxt)
	N = rows_cols(index)
	if (N[1] & N[2]) t[index] = J(N[1],1,maxt)
	return(t)
}

`RC' merlin_simulate_weibull(`gml' gml, `RC' logU, `RS' maxt, | `RC' t0)
{
	gam 	= merlin_util_dap(gml,1)
	lambda 	= exp(merlin_util_xzb(gml))
	if (args()==4) 	t = (t0:^gam :- logU:/lambda):^(1:/gam)
	else 			t = (-logU:/lambda):^(1:/gam)
	index = selectindex(t:>maxt)
	N = rows_cols(index)
	if (N[1] & N[2]) t[index] = J(N[1],1,maxt)
	return(t)
}

`RC' merlin_simulate_gompertz(`gml' gml, `RC' logU, `RS' maxt, | `RC' t0)
{
	gam 	= merlin_util_dap(gml,1)
	lambda 	= exp(merlin_util_xzb(gml))
	if (args()==4) 	t = log(1:-(logU :- lambda:/gam:*(exp(gam:*t0):-1)):*gam:/lambda):/gam
	else 			t = log(1:-logU:*gam:/lambda):/gam
	index = selectindex(t:>maxt)
	N = rows_cols(index)
	if (N[1] & N[2]) t[index] = J(N[1],1,maxt)
	return(t)
}

end

