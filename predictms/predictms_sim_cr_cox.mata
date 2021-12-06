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

`RM' predictms_sim_cr_cox(	`SS' S, 	/// predictms object
							`RS' Nsim, 	/// # of new times to simulate
							`RC' t0, 	///	entry time
							`RS' s0		/// entry state
							)
{

	newdata		= J(Nsim,3,0)				//time,state,done
	Ntrans 		= S.Nnextstates[s0]
	postrans 	= asarray(S.postrans,s0)	//possible transitions from current states
	
	noreset		= !sum(S.tsreset[postrans])
	notscale2	= !sum(S.tscale2[postrans])
	
	allreset	= sum(S.tsreset[postrans])==Ntrans
	alltscale2	= sum(S.tscale2[postrans])==Ntrans

	//get transition-specific chazards on original as estimated timescales
	merlinP = J(Ntrans,1,NULL)
	corech	= J(0,3,.)
	for (s=1; s<=Ntrans; s++) {
		trans 		= postrans[s]
		corech_s 	= predictms_merlin_cox_ch(S,trans)
		corech		= corech\corech_s,J(rows(corech_s),1,trans)
	}

	if (noreset & notscale2) {

		_sort(corech,2)
		
		//now need transition-specific ch and h at all unique failure times
		uniqt 	= sort(uniqrows(corech[,2]),1)
		Nuniqt 	= rows(uniqt)

		th = tch = J(Nuniqt,Ntrans,0)
		
		for (s=1; s<=Ntrans; s++) {
			
			trans 	= asarray(S.postrans,s0)[s]	
			chs 	= select(corech,corech[,3]:==trans)
			Nchs	= rows(chs)

			time1 	= 0
			chsub	= 0
			for (i=1;i<=Nchs;i++) {
				time2  	= chs[i,2]
				tindex	= selectindex((time1:<=uniqt) :& (uniqt:<time2))
				Ntindex = rows_cols(tindex)
				if (Ntindex[1] & Ntindex[2]) tch[tindex,s] = J(Ntindex[1],1,chsub)
				time1 	= time2
				chsub	= chs[i,1]
			}
			//final row
			tindex	= selectindex(time1:<=uniqt)
			Ntindex = rows_cols(tindex)
			if (Ntindex[1] & Ntindex[2]) tch[tindex,s] = J(Ntindex[1],1,chsub)
			
		}

		th[1,] = tch[1,]
		th[|2,.\Nuniqt,.|] = tch[|2,.\Nuniqt,.|] :- tch[|1,.\(Nuniqt-1),.|]

		for (i=1;i<=Nsim;i++) {
			timeindex		= selectindex(t0[i]:<uniqt)
			thnew 			= th[timeindex,]
			tallhnew 		= quadrowsum(thnew)
			tallchnew 		= -quadrunningsum(log(1:-tallhnew))
			newdata[i,1..2]	= predictms_cox_sample(tallchnew,rows(thnew),(uniqt[timeindex]\S.maxt))
			if (newdata[i,1]!=S.maxt) {
				if (Ntrans>1) 	newdata[i,2] = predictms_cox_sample_state(thnew[newdata[i,2],]')
				else 			newdata[i,2] = 1
			}
		}
		
	}
	
	else if (allreset) {

		_sort(corech,2)
		
		//now need transition-specific ch and h at all unique failure times
		uniqt 	= sort(uniqrows(corech[,2]),1)
		Nuniqt 	= rows(uniqt)
		
		th = tch = J(Nuniqt,Ntrans,0)
		
		for (s=1; s<=Ntrans; s++) {
			
			trans 	= asarray(S.postrans,s0)[s]	
			chs 	= select(corech,corech[,3]:==trans)
			Nchs	= rows(chs)

			time1 	= 0
			chsub	= 0
			for (i=1;i<=Nchs;i++) {
				time2  	= chs[i,2]
				tindex	= selectindex((time1:<=uniqt) :& (uniqt:<time2))
				Ntindex = rows_cols(tindex)
				if (Ntindex[1] & Ntindex[2]) tch[tindex,s] = J(Ntindex[1],1,chsub)
				time1 	= time2
				chsub	= chs[i,1]
			}
			//final row
			tindex	= selectindex(time1:<=uniqt)
			Ntindex = rows_cols(tindex)
			if (Ntindex[1] & Ntindex[2]) tch[tindex,s] = J(Ntindex[1],1,chsub)
			
		}
		
		th[1,] = tch[1,]
		th[|2,.\Nuniqt,.|] = tch[|2,.\Nuniqt,.|] :- tch[|1,.\(Nuniqt-1),.|]
		
		for (i=1;i<=Nsim;i++) {
			tallhnew 		= quadrowsum(th)
			tallchnew 		= -quadrunningsum(log(1:-tallhnew))
			newdata[i,1..2]	= predictms_cox_sample(tallchnew,rows(th),(uniqt\(S.maxt:-t0[i])))
			if (newdata[i,1]!=S.maxt) {
				if (Ntrans>1) 	newdata[i,2] = predictms_cox_sample_state(th[newdata[i,2],]')
				else 			newdata[i,2] = 1
			}
		}
		
		newdata[,1] = newdata[,1] :+ t0
		
	}
	
	else {
		
		for (i=1;i<=Nsim;i++) {
		
			icorech = corech
			
			//adjust times for reset and/or tscale2
			for (s=1; s<=Ntrans; s++) {
				
				trans 		= postrans[s]
				tindex		= selectindex(icorech[,3]:==trans)
				
				if (S.tsreset[trans]) icorech[tindex,2] = icorech[tindex,2] :+ t0[i]
				if (S.tscale2[trans]) icorech[tindex,2] = icorech[tindex,2] :+ S.time2[S.at]
				
			}
		
			_sort(icorech,2)
		
			//now need transition-specific ch and h at all unique failure times
			uniqt 	= sort(uniqrows(icorech[,2]),1)
			Nuniqt 	= rows(uniqt)
			
			th = tch = J(Nuniqt,Ntrans,0)
			
			for (s=1; s<=Ntrans; s++) {
				
				trans 	= asarray(S.postrans,s0)[s]	
				chs 	= select(icorech,icorech[,3]:==trans)
				Nchs	= rows(chs)

				time1 	= 0
				chsub	= 0
				for (c=1;c<=Nchs;c++) {
					time2  	= chs[c,2]
					tindex	= selectindex((time1:<=uniqt) :& (uniqt:<time2))
					Ntindex = rows_cols(tindex)
					if (Ntindex[1] & Ntindex[2]) tch[tindex,s] = J(Ntindex[1],1,chsub)
					time1 	= time2
					chsub	= chs[c,1]
				}
				//final row
				tindex	= selectindex(time1:<=uniqt)
				Ntindex = rows_cols(tindex)
				if (Ntindex[1] & Ntindex[2]) tch[tindex,s] = J(Ntindex[1],1,chsub)
				
			}
			
			th[1,] = tch[1,]
			th[|2,.\Nuniqt,.|] = tch[|2,.\Nuniqt,.|] :- tch[|1,.\(Nuniqt-1),.|]
		
		
			timeindex		= selectindex(t0[i]:<uniqt)
			thnew 			= th[timeindex,]
			tallhnew 		= quadrowsum(thnew)
			tallchnew 		= -quadrunningsum(log(1:-tallhnew))
			newdata[i,1..2]	= predictms_cox_sample(tallchnew,rows(thnew),(uniqt[timeindex]\S.maxt))
			if (newdata[i,1]!=S.maxt) {
				if (Ntrans>1) 	newdata[i,2] = predictms_cox_sample_state(thnew[newdata[i,2],]')
				else 			newdata[i,2] = 1
			}
		}		
	
	}

	update_rcens(newdata,s0,S.maxt)
	update_events(newdata,S,s0)
	return(newdata)
}

/*
get CH for a Cox model
-> need to build the at rather than replace vars, since baseline hazard must be based on original dataset
*/

`RM' predictms_merlin_cox_ch(`SS' S, `RS' trans)
{
	strtrans 	= strofreal(trans)
	at 			= S.at
	std 		= S.std

	stata("qui preserve")

	if (S.hasmodels) stata("cap estimates restore "+st_local("modelests"+strtrans))

	//get baseline hazard before any dataset changes 
	stata("tempvar baseh0"+strtrans)
	stata("qui predict "+st_local("baseh0"+strtrans)+", basehazard")

	if (st_local("survsim")=="") {

		//start by replacing allvars with zeroes (if not in stdvars)
		allvars 	= tokens(st_global("e(allvars)"))
		y			= tokens(st_global("e(response1)"))
		Nallvars	= cols(allvars)
		zerovars 	= J(1,0,"")
		
		if (S.standardise) {
			stdvars = tokens(st_local("at"+strofreal(at)+"stdvars"+strtrans))
			for (i=1; i<=Nallvars; i++) {
				flag = 0
				for (j=1;j<=cols(stdvars);j++) {
					if (stdvars[j]==allvars[i]) flag = 1
				}
				if (!flag) zerovars = zerovars,allvars[i]
			}
		}
		else zerovars = allvars

		for (i=1;i<=cols(zerovars);i++) {
			if (!sum(zerovars[i]:==y)) stata("qui replace "+zerovars[i]+" = 0 if e(sample)")
		}

		//need to replace any variables with their at#()
		ats 		= tokens(st_local("at"+strofreal(at)))'
		Natstodo 	= rows(ats)/2
		i = 1
		j = 2
		for (a=1;a<=Natstodo;a++) {
			stata("qui replace "+ats[i]+" = "+ats[j]+" if e(sample)")
			i = i+2
			j = j+2
		}

		//standardisation
		if (S.standardise) {
			st_view(poststdvars=.,.,stdvars,st_local("stdtouse"))
			for (i=1;i<=cols(stdvars);i++) {
				stata("qui replace "+stdvars[i]+" = "+strofreal(poststdvars[std,i])+" if e(sample)")
			}
		}

		//if stacked model, update *_trans# variables
		if (!S.hasmodels) {

			for (i=1;i<=S.Ntrans;i++) {
				ti = "_trans"+strofreal(i)
				
				if (i==trans) 	stata("qui replace "+ti+" = 1 if e(sample)")
				else 			stata("qui replace "+ti+" = 0 if e(sample)")
				
				//rebuild any *_transi variables
				if (st_local("toupdate")!="") {
					for (k=1;k<=rows(S.toupdate);k++) {
						stata("cap replace "+S.toupdate[k]+ti+"= "+S.toupdate[k]+"*"+ti+" if e(sample)")
					}
				}
			}
			
		}

	}	

	stata("tempvar ch"+strtrans+" chuse"+strtrans)
	stata("qui predict "+st_local("ch"+strtrans)+", chazard passbaseh("+st_local("baseh0"+strtrans)+")")
	stata("qui gen byte "+st_local("chuse"+strtrans)+"= e(sample) & "+st_global("e(failure1)"))
	vars 	= st_local("ch"+strtrans),y[1]
	ch 		= st_data(.,vars,st_local("chuse"+strtrans))

	stata("qui restore")
	return(ch)
}

`RR' predictms_cox_sample(`RC' ch, `RS' N, `RC' times)
{
	cdf 				= 1:-exp(-ch)
	probs 				= J(N+1,1,.)
	probs[1] 			= cdf[1]
	probs[|2,1\N,1|] 	= cdf[|2,1\N,1|] :- cdf[|1,1\(N-1),1|]
	probs[N+1] 			= exp(-ch[N])
	row 				= rdiscrete(1,1,probs)
	time 				= times[row]
	if (time>times[N+1]) time = times[N+1]
	return(time,row)
	
}

`RS' predictms_cox_sample_state(`RC' hazards)
{
	hazards = hazards :/ sum(hazards)
	return(rdiscrete(1,1,hazards))
}

//update new data draws for events i.e. new states
//-> s0, starting state

void update_events(`RM' newdata, `SS' S, `RS' s0)
{
	eventindex 	= selectindex(newdata[,1]:<S.maxt)
	Nevents		= rows_cols(eventindex)
	if (Nevents[1] & Nevents[2]) {
	
		nextstates 				= asarray(S.posnextstates,s0)
		newdata[eventindex,2] 	= nextstates[newdata[eventindex,2],]

		//obs now in absorbing state are done
		for (s=1; s<=S.Nnextstates[s0]; s++) {
			if (!S.Nnextstates[nextstates[s]]) {
				doneindex 	= select(eventindex,newdata[eventindex,2]:==nextstates[s])
				Ndone 		= rows_cols(doneindex)
				if (Ndone[1] & Ndone[2]) newdata[doneindex,3] = J(Ndone[1],1,1)
			}
		}
	}
}

end
