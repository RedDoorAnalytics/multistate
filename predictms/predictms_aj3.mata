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

void predictms_aj3(`SS' S, `RS' from, `RS' at)
{

	if (S.getcis) {
		VarCH 	= J(S.obs,S.Ntrans,.)
		S.varpt = J(S.obs,S.Nstates^2,0)
	}

		S.pt 	= J(S.obs,S.Nstates^2,0)				//to hold all predictions
		ch 		= J(S.obs,S.Ntrans,.)
		dch 	= ch
		
		//std loop
		for (std=1;std<=S.K;std++) {

			//get cumulative hazards for each transition
			for (trans=1;trans<=S.Ntrans;trans++) {

				dm	= asarray(S.X,(trans,at))[std,]
				b 	= asarray(S.transinfo,(trans,2))'
				//RP
				if (S.modeltrans[trans,6]) {
					Pstpm2 		= asarray(S.transinfo,(trans,4))
					ch[,trans] 	= stpm2_ch(b,S.predtime,dm,at,std,Pstpm2)
		
					if (S.getcis) {
						//get standard error of CH through delta method using numerical derivatives
						nbs = rows(b)
						deriv = J(S.obs,nbs,.)
						for (i=1;i<=nbs;i++) {
							hs = sqrt(epsilon(1))*b[i]
							newb = b
							newb[i] = b[i] :+ hs/2
							uh = stpm2_ch(newb,S.predtime,dm,at,std,Pstpm2)
							newb[i] = b[i] :- hs/2
							lh = stpm2_ch(newb,S.predtime,dm,at,std,Pstpm2)
							deriv[,i] = (uh :- lh):/hs
						}
						
						for (i=1;i<=S.obs;i++) VarCH[i,trans] =  deriv[i,] * asarray(S.transinfo,(trans,5)) * deriv[i,]'
					}
				}
				//strcs
				else if (S.modeltrans[trans,7]) {
					Pstrcs 		= asarray(S.transinfo,(trans,4))
					ch[,trans] 	= strcs_ch(S.predtime,b,dm,at,std,Pstrcs)
					
					if (S.getcis) {
						//get standard error of CH through delta method using numerical derivatives
						nbs = rows(b)
						deriv = J(S.obs,nbs,.)
						for (i=1;i<=nbs;i++) {
							hs = sqrt(epsilon(1))*b[i]
							newb = b
							newb[i] = b[i] :+ hs/2
							uh = strcs_ch(S.predtime,newb,dm,at,std,Pstrcs)
							newb[i] = b[i] :- hs/2
							lh = strcs_ch(S.predtime,newb,dm,at,std,Pstrcs)
							deriv[,i] = (uh :- lh):/hs
						}
						
						for (i=1;i<=S.obs;i++) VarCH[i,trans] =  deriv[i,] * asarray(S.transinfo,(trans,5)) * deriv[i,]'
					}
				}
				//weibull
				else if (S.modeltrans[trans,2]) {		

					indices 	= asarray(S.transinfo,(trans,3))
					lambda 		= exp(dm[|.,indices[1,1]\.,indices[2,1]|] * b[|indices[,1]|])
					gamm 		= exp(dm[|.,indices[1,2]\.,indices[2,2]|] * b[|indices[,2]|])
					ch[,trans] 	= lambda :* S.predtime :^ gamm
					
					if (S.getcis) {
						//get standard error of CH through delta method using numerical derivatives
						nbs = rows(b)
						deriv = J(S.obs,nbs,.)
						for (i=1;i<=nbs;i++) {
							hs = sqrt(epsilon(1))*b[i]
							newb = b
							newb[i] = b[i] :+ hs/2
							lambda 		= exp(dm[|.,indices[1,1]\.,indices[2,1]|] * newb[|indices[,1]|])
							gamm 		= exp(dm[|.,indices[1,2]\.,indices[2,2]|] * newb[|indices[,2]|])
							uh		 	= lambda :* S.predtime :^ gamm
							newb[i] = b[i] :- hs/2
							lambda 		= exp(dm[|.,indices[1,1]\.,indices[2,1]|] * newb[|indices[,1]|])
							gamm 		= exp(dm[|.,indices[1,2]\.,indices[2,2]|] * newb[|indices[,2]|])
							lh		 	= lambda :* S.predtime :^ gamm
							deriv[,i] = (uh :- lh):/hs
						}
						
						for (i=1;i<=S.obs;i++) VarCH[i,trans] =  deriv[i,] * asarray(S.transinfo,(trans,5)) * deriv[i,]'
					}
				}
				//exp
				else if (S.modeltrans[trans,1]) {		
					ch[,trans] 	= exp(dm * b) :* S.predtime 
					
					if (S.getcis) {
						//get standard error of CH through delta method using numerical derivatives
						nbs = rows(b)
						deriv = J(S.obs,nbs,.)
						for (i=1;i<=nbs;i++) {
							hs = sqrt(epsilon(1))*b[i]
							newb = b
							newb[i] = b[i] :+ hs/2
							uh		 	= exp(dm * newb) :* S.predtime 
							newb[i] = b[i] :- hs/2
							lh		 	= exp(dm * newb) :* S.predtime 
							deriv[,i] = (uh :- lh):/hs
						}
						
						for (i=1;i<=S.obs;i++) VarCH[i,trans] =  deriv[i,] * asarray(S.transinfo,(trans,5)) * deriv[i,]'
					}
				}
				
			
			}

			//change in ch since previous timepoint
			dch[1,] = ch[1,]
			dch[|2,.\.,.|] = ch[|2,.\.,.|] :- ch[|1,.\(S.obs-1),.|]

			//put variances in correct place
			if (S.getcis) {
			
				varprobmat = J(S.obs,S.Nstates^4,0)
				
				ind = 1
				for (i=1;i<=S.Nstates;i++) {
					for (j=1;j<=S.Nstates;j++) {
						for (k=1;k<=S.Nstates;k++) {
							for (p=1;p<=S.Nstates;p++) {
								if (i==k & j==p & S.transmat[i,j]!=.) {
									varprobmat[,ind] = VarCH[,S.transmat[i,j]]
								}
								ind++					
							}
						}
					}
				}
				
			}
			
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

			if (S.getcis) {
				ident 	= P
				varres 	= J(S.obs,cols(probmat),.)
				varP 	= J(S.Nstates^2,S.Nstates^2,0)
				
				vardA 				= J(S.obs,S.Nstates^4,0)
				vardA[1,] 			= varprobmat[1,]
				vardA[|2,.\.,.|] 	= varprobmat[|2,.\.,.|] :+ varprobmat[|1,.\(S.obs-1),.|]
				nrows 	= S.Nstates^2
			}
			
			
			postP = J(S.obs,S.Nstates^2,.)

			//transition probabilities
			if (S.isenter) {
				for (i=1;i<=S.obs;i++) {
					Pi = rowshape(probmat[i,],S.Nstates)
					P = P * Pi
					postP[i,] = rowshape(P,1)
					if (S.getcis) {
						varP = (Pi')#ident * varP * Pi#ident :+ ident#P * rowshape(vardA[i,],nrows) * ident#(P')
						varres[i,] = diagonal(varP)'
					}
				}
			}
			else {
				/*for (i=1;i<=N;i++) {
					P = rowshape(probmat[i,],Nstates) * P
					res = rowshape(P,1) \ res 
				}*/
			}

			S.pt = S.pt :+ postP
			S.varpt = S.varpt :+ varres
			
			if (min(postP)<0) {
				errprintf("Negative probabilities, increase obs() if using a parametric model\n")
				exit(1986)
			}
			
		}
		
		//standardisation
		if (S.isstd) {
			S.pt = S.pt :/ S.K
			if (S.getcis) S.varpt = S.varpt :/ S.K		//!!wrong
		}

		//Store predictions	-> AJ calculates for all froms, but only storing the ones asked for
		ind1 = (from-1)*S.Nstates+1
		ind2 = ind1 + S.Nstates - 1
		
		asarray(S.probs,(from,at),S.pt[|.,ind1\.,ind2|])
		if (S.getcis) asarray(S.varprobs,(from,at),S.varpt[|.,ind1\.,ind2|])
}

end

