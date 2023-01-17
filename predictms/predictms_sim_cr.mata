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
local PCm	pointer(struct merlin_struct scalar) colvector

mata:

`RM' predictms_sim_cr(			///
                        `SS' S, 	/// predictms object
                        `RS' Nsim, 	/// # of new times to simulate
                        `RC' t0, 	/// entry time
                        `RS' s0		/// entry state
                        )
{
	newdata	= J(Nsim,3,0)			//time,state,done
	Ntrans 	= S.Nnextstates[s0]

	//get merlin structs
	merlinP = J(Ntrans,1,NULL)
	for (s=1; s<=Ntrans; s++) {
		trans 		= asarray(S.postrans,s0)[s]
		b		= predictms_get_b(S,trans)
		merlinP[s] 	= predictms_merlin_setup(S,b,Nsim,trans,t0)
	}

	//simulate new times
	logU	= log(runiform(Nsim,1))
	rc 	= predictms_sim_root(	t=J(Nsim,1,.),			///
                                &predictms_total_ch(),	                ///
                                S.maxt,					///
                                0,1000,					///
                                Nsim,					///
                                t0,					///
                                logU,					///
                                1::Nsim,				///
                                merlinP,				///
                                Ntrans,					///
                                asarray(S.postrans,s0),	                ///
                                S.tsreset,				///
                                S.tscale2,				///
                                S.time2,				///
                                S.at)
	newdata[,1] = t
	update_rcens(newdata,s0,S.maxt)	
	
	//new states
	eventindex 	= selectindex(t:<S.maxt)
	Nevents		= rows_cols(eventindex)
	if (Nevents[1] & Nevents[2]) {
		nextstates = asarray(S.posnextstates,s0)
		newdata[eventindex,2] 	= predictms_get_newstate(t[eventindex],eventindex,merlinP,Ntrans,nextstates,rc[eventindex])

		//obs now in absorbing state are done
		for (s=1; s<=Ntrans; s++) {
			if (!S.Nnextstates[nextstates[s]]) {
				doneindex 	= select(eventindex,newdata[eventindex,2]:==nextstates[s])
				Ndone 		= rows_cols(doneindex)
				if (Ndone[1] & Ndone[2]) newdata[doneindex,3] = J(Ndone[1],1,1)
			}
		}
	}
	
	//tidy up
	for (s=1; s<=Ntrans; s++) {
		trans 		= asarray(S.postrans,s0)[s]
		rmexternal(st_local("GML"+strofreal(trans))) 	
	}
	
	return(newdata)
}

`RM' predictms_total_ch(`RC' x, `RC' t0, `RC' logU, `RC' index, `PCm' p, `RS' Ntrans, 	///
		`RC' transs, `RC' reset, `RC' tscale2, `RC' time2, `RS' at)
{
	struct merlin_struct gml
	
	Nobs 	= rows(x)
	ch	= J(Nobs,1,0)

	for (s=1; s<=Ntrans; s++) {
		
		trans 	= transs[s]
		gml 	= *p[s]
		
		//index etc. gets updated as mm_root finishes observations
		gml.N = gml.Nobs = Nobs	
		asarray(gml.xbindex,1,index)
		
		if (gml.NI) {
			Nq 	= strtoreal(st_local("chintpoints"))
			qpw 	= predictms_gq(Nq)
			qp 	= ((x:-t0) :* J(Nobs,1,qpw[,1]') :+ x:+t0) :/ 2
			if (reset[trans]) 	qp = qp :- t0
			if (tscale2[trans])     qp = qp :+ time2[at]
			
			chq = J(Nobs,Nq,.)
			for (q=1; q<=Nq; q++) {
				chq[,q] = exp((*gml.Plogh)(gml,qp[,q]))
			}
			_editmissing(chq,0)		
			ch = ch :+ (x:-t0):/2 :* chq * qpw[,2]
		}
		else {

			mainx	= x
			if (reset[trans]) 	mainx = mainx :- t0
			if (tscale2[trans]) mainx = mainx :+ time2[at]
			
			chq 	=  (*gml.Pch)(gml,mainx)
			_editmissing(chq,0)

			if (reset[trans]) {
				if (tscale2[trans]) chq0 = (*gml.Pch)(gml,J(Nobs,1,time2[at]))
				else {
					ch = ch :+ chq
					continue
				}
			}
			else {
				if (tscale2[trans]) chq0 = (*gml.Pch)(gml,t0:+time2[at])
				else 				chq0 = (*gml.Pch)(gml,t0)
			}
			_editmissing(chq0,0)
			ch = ch :+ chq :- chq0
			
		}
		
	}

	return(ch :+ logU)	
}

`RM' predictms_get_newstate(`RC' t, `RC' index, `PCm' p, `RS' Ntrans, `RM' nextstates, `RC' rc)
{
	struct merlin_struct gml
	
	Nobs 		= rows(t)
	
	if (Ntrans==1) {
		return(J(Nobs,1,nextstates))
	}
	else {
	
		probs 		= J(Nobs,Ntrans,.)
		newstate	= J(Nobs,1,.)

		//get transition-specific hazards 
		for (s=1; s<=Ntrans; s++) {
			gml = *p[s]
			gml.N = gml.Nobs = Nobs	
			asarray(gml.xbindex,1,index)
			probs[,s] = exp((*gml.Plogh)(gml,t))
		}
		
		//scale to 1
		probs = probs :/ quadrowsum(probs)
		if (min(probs)<0) {
			errprintf("A negative hazard function has been sampled\n")
			errprintf("-> try the simulate & latent options\n")
			errprintf("-> try switching to a log hazard scale model if you are using a Royston-Parmar model")
			exit(1986)
		}
		//draw transition
		for (i=1; i<=Nobs; i++) newstate[i] = rdiscrete(1,1,probs[i,])
		//return new states
		return(nextstates[newstate])	
	}
}

end
