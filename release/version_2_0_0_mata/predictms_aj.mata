*! version 1.0.0 

/*
Notes
*/

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

void predictms_aj(`SS' S, `RS' from, `RS' at)
{

	if (S.getcis) {
		A = asarray_create("real",1)					//to hold predictions across trans and sims
		if (S.getlos) {
			B = asarray_create("real",1)
		}
		if (S.hasuser) {
			C = asarray_create("real",1)
		}
	}

	ind1 = (from-1)*S.Nstates+1
	ind2 = ind1 + S.Nstates - 1
	
	// loop over simulations
	for (k=1;k<=S.M;k++) {	
	
		S.pt 	= J(S.obs,S.Nstates^2,0)				//to hold all predictions
		ch 		= J(S.obs,S.Ntrans,.)
		dch 	= ch
		
		if (S.getlos) S.los = J(S.obs,S.Nstates,0)
		
		//std loop
		for (std=1;std<=S.K;std++) {

			//get cumulative hazards for each transition
			for (trans=1;trans<=S.Ntrans;trans++) {

				dm	= asarray(S.X,(trans,at))[std,]
				b 	= asarray(S.transinfo,(trans,2))[k,]'
				
				t2flag 	= sum(S.tscale2:==trans)  						//second timescale S.time2
				if (t2flag) t2 = S.time2[at]

				//megenreg
				if (S.modeltrans[trans,8]) {
					Pmegenreg 	= predictms_megenreg_setup(b,at,S.obs,trans)
					if (t2flag) 	ch[,trans] 	= megenreg_ch(S.predtime:+t2,Pmegenreg) :- megenreg_ch(t2,Pmegenreg)
					else 			ch[,trans] 	= megenreg_ch(S.predtime,Pmegenreg)
					if (S.enter) 	ch[,trans] 	= ch[,trans] :- megenreg_ch(S.enter,Pmegenreg)
					rmexternal(st_global("e(struct)")) 					//tidy up
				}
				//RP
				else if (S.modeltrans[trans,6]) {
					Pstpm2 		= asarray(S.transinfo,(trans,4))
					if (t2flag) 	ch[,trans] 	= stpm2_ch(b,S.predtime:+t2,dm,at,std,Pstpm2) :- stpm2_ch(b,t2,dm,at,std,Pstpm2)
					else 			ch[,trans] 	= stpm2_ch(b,S.predtime,dm,at,std,Pstpm2)
					if (S.enter) 	ch[,trans] 	= ch[,trans] :- stpm2_ch(b,S.enter,dm,at,std,Pstpm2)
				}
				//strcs
				else if (S.modeltrans[trans,7]) {
					Pstrcs 		= asarray(S.transinfo,(trans,4))
					if (t2flag) 	ch[,trans] 	= strcs_ch(S.predtime:+t2,b,dm,at,std,Pstrcs) :- strcs_ch(t2,b,dm,at,std,Pstrcs)
					else 			ch[,trans] 	= strcs_ch(S.predtime,b,dm,at,std,Pstrcs)
					if (S.enter) 	ch[,trans] 	= ch[,trans] :- strcs_ch(S.enter,b,dm,at,std,Pstrcs)
				}
				//weibull
				else if (S.modeltrans[trans,2]) {	
					indices 	= asarray(S.transinfo,(trans,3))
					lambda 		= exp(dm[|.,indices[1,1]\.,indices[2,1]|] * b[|indices[,1]|])
					gamm 		= exp(dm[|.,indices[1,2]\.,indices[2,2]|] * b[|indices[,2]|])
					if (t2flag) 	ch[,trans] 	= lambda :* (S.predtime :+ t2) :^ gamm :- lambda :* t2 :^ gamm
					else 			ch[,trans] 	= lambda :* S.predtime :^ gamm
					if (S.enter) 	ch[,trans] 	= ch[,trans] :- lambda :* S.enter :^ gamm
				}
				//exp
				else if (S.modeltrans[trans,1]) {		
					lambda	 	= exp(dm * b)
					if (t2flag) 	ch[,trans] 	= lambda :* (S.predtime :+ t2) :- lambda :* t2 
					else 			ch[,trans] 	= lambda :* S.predtime 
					if (S.enter) 	ch[,trans] 	= ch[,trans] :- lambda :* S.enter
				}
				//gompertz
				else if (S.modeltrans[trans,3]) {
					indices 	= asarray(S.transinfo,(trans,3))
					lambda 		= exp(dm[|.,indices[1,1]\.,indices[2,1]|] * b[|indices[,1]|])
					gamm 		= dm[|.,indices[1,2]\.,indices[2,2]|] * b[|indices[,2]|]
					if (t2flag) 	ch[,trans] 	= lambda:/gamm :* (exp(gamm:*(S.predtime:+t2)):-1) :- lambda:/gamm :* (exp(gamm:*t2):-1)
					else 			ch[,trans] 	= lambda:/gamm :* (exp(gamm:*S.predtime):-1)
					if (S.enter) 	ch[,trans] 	= ch[,trans] :- lambda:/gamm :* (exp(gamm:*S.enter):-1)
				}
				//llog
				else if (S.modeltrans[trans,4]) {
					indices 	= asarray(S.transinfo,(trans,3))
					lambda 		= exp(-dm[|.,indices[1,1]\.,indices[2,1]|] * b[|indices[,1]|])
					gamm 		= exp(dm[|.,indices[1,2]\.,indices[2,2]|] * b[|indices[,2]|])
					if (t2flag) 	ch[,trans] 	= log(1 :+ (lambda:*(S.predtime:+t2)):^(1:/gamm)) :- log(1 :+ (lambda:*t2):^(1:/gamm))
					else 			ch[,trans] 	= log(1 :+ (lambda:*S.predtime):^(1:/gamm))
					if (S.enter) 	ch[,trans] 	= ch[,trans] :- log(1 :+ (lambda:*S.enter):^(1:/gamm))
				}
				//lognormal
				else {
					indices 	= asarray(S.transinfo,(trans,3))
					lambda 		= dm[|.,indices[1,1]\.,indices[2,1]|] * b[|indices[,1]|]
					gamm 		= exp(dm[|.,indices[1,2]\.,indices[2,2]|] * b[|indices[,2]|])
					if (t2flag) 	ch[,trans] 	= -log(1:-normal((log(S.predtime :+ t2):-lambda):/gamm)) :+ log(1:-normal((log(t2):-lambda):/gamm))
					else 			ch[,trans] 	= -log(1:-normal((log(S.predtime):-lambda):/gamm))
					if (S.enter) 	ch[,trans] 	= ch[,trans] :+ log(1:-normal((log(S.enter):-lambda):/gamm))
				}
			
			}

			//change in ch since previous timepoint
			dch[1,] = ch[1,]
			dch[|2,.\.,.|] = ch[|2,.\.,.|] :- ch[|1,.\(S.obs-1),.|]

			P = I(S.Nstates)
			
			//matrix of times x prob matrix
			probmat = J(S.obs,1,rowshape(P,1))	
			ind = 1
			for (i=1;i<=S.Nstates;i++) {
				for (j=1;j<=S.Nstates;j++) {
					if (i==j) probmat[,ind] = probmat[,ind] :- quadrowsum(dch[,asarray(S.postrans,i)])
					else {
						if (S.transmat[i,j]!=.) probmat[,ind] = dch[,S.transmat[i,j]]
					}
					ind++
				}
			}

			postP = J(S.obs,S.Nstates^2,.)
			//transition probabilities
			if (S.isenter) {
				for (i=1;i<=S.obs;i++) {
					Pi = rowshape(probmat[i,],S.Nstates)
					P = P * Pi
					postP[i,] = rowshape(P,1)
				}
			}
			else {
				/*for (i=1;i<=N;i++) {
					P = rowshape(probmat[i,],Nstates) * P
					res = rowshape(P,1) \ res 
				}*/
			}

			if (min(postP)<0) {
				errprintf("Negative probabilities, increase obs() if using a parametric model\n")
				exit(1986)
			}
			
			S.pt = S.pt :+ postP
			
			if (S.getlos) {
				hstep = S.predtime[2,1] :- S.predtime[1,1]
				los = J(S.obs,S.Nstates:^2,0)
				for (i=2;i<=S.obs;i++) {
					los[i,] = postP[1,] :+ postP[S.obs,] :+ 2 :* quadcolsum(postP[|2,.\i,.|])
				}
				los = los :* hstep :/2
				S.los = S.los :+ los[|.,ind1\.,ind2|]
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
		
		//standardisation
		if (S.isstd) {
			S.pt = S.pt :/ S.K
			if (S.getlos) 	S.los = S.los :/ S.K
			if (S.hasuser) 	S.user = S.user :/ S.K
		}

		//Store predictions	-> AJ calculates for all froms, but only storing the ones asked for
		if (!S.getcis) {

			asarray(S.probs,(from,at),S.pt[|.,ind1\.,ind2|])
			if (S.getlos) 	asarray(S.loss,(from,at),S.los)
			if (S.hasuser) 	asarray(S.users,(from,at),S.user)
			
		}
		else {

			//-> time points x simulations, for each current and next state combination
			if (k==1) {
				postind = 1
				for (tr=ind1;tr<=ind2;tr++) {
					asarray(A,postind++,S.pt[,tr])								
				}
			}
			else {
				postind = 1
				for (tr=ind1;tr<=ind2;tr++) {
					asarray(A,postind,(asarray(A,postind),S.pt[,tr]))
					postind++
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

