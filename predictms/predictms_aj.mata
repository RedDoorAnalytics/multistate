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
local gml	struct merlin_struct

mata:

void predictms_aj(`SS' S, from)
{
	struct merlin_struct gml
	
	ch = dch = J(S.obs,S.Ntrans,.)

	//get cumulative hazards for each transition
	for (trans=1; trans<=S.Ntrans; trans++) {
		b		= predictms_get_b(S,trans)
		gml		= *predictms_merlin_setup(S,b,S.obs,trans)
		ch[,trans] 	= merlin_ch(S.predtime,gml)
		rmexternal(st_local("GML"+strofreal(trans)))
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

	postP = J(S.obs,S.Nstates,.)
	//transition probabilities
	for (i=1;i<=S.obs;i++) {
		Pi = rowshape(probmat[i,],S.Nstates)
		P = P * Pi
		postP[i,] = P[from,]
	}

	if (min(postP)<0) {
		errprintf("Negative probabilities, increase obs()\n")
		exit(1986)
	}

	if (S.standardise) 	S.pt = S.pt :+ postP
	else 			S.pt = postP
	
	if (S.getlos | S.getrmst) {
		hstep = S.predtime[2,1] :- S.predtime[1,1]
		los = J(S.obs,S.Nstates,0)
		for (i=2;i<=S.obs;i++) {
			los[i,] = postP[1,] :+ postP[S.obs,] :+ 2 :* quadcolsum(postP[|2,.\i,.|])
		}
		los = los :* hstep :/2
		if (S.standardise) 	S.los = S.los :+ los
		else				S.los = los
		
		if (S.getrmst) {
			rmst = J(S.obs,1,0)
			for (s=1;s<=S.Nstates;s++) {
				if (S.Nnextstates[s]) rmst = rmst :+ S.los[,s]
			}
			if (S.standardise) 	S.rmst = S.rmst :+ rmst
			else				S.rmst = rmst
		}
	}

	if (S.hasuser) predictms_calc_user(S)
	
}

`RM' merlin_ch(	`RC' x,	`gml' gml)
{
	
	//update N, Nobs and index -> used in merlin utils, and main xb()
	Nobs 	= rows(x)
	gml.N 	= gml.Nobs[gml.Nlevels,1] = Nobs
	asarray(gml.xbindex,1,asarray(gml.xbindex,1)[|1,1\Nobs,1|])

	//call CH
	if 		(gml.familys[1]=="rp") 			ch = merlin_p_rp_ch(gml,x)
	else if (gml.familys[1]=="loghazard") 	ch = merlin_p_loghazard_ch(gml,x)
	else if (gml.familys[1]=="exp") 		ch = merlin_p_exp_ch(gml,x)
	else if (gml.familys[1]=="weibull") 	ch = merlin_p_weibull_ch(gml,x)
	else if (gml.familys[1]=="gompertz") 	ch = merlin_p_gompertz_ch(gml,x)
	else if (gml.familys[1]=="ggamma") 		ch = merlin_p_ggamma_ch(gml,x)
	else if (gml.familys[1]=="user") 		ch = merlin_p_userh_ch(gml,x)

	_editmissing(ch,0)
	return(ch)
}

end

