*! version 1.0.0 ?????2015 MJC

version 12.1

local ss 	string scalar
local RS	real scalar
local NM	numeric matrix
local RM 	real matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local PS	pointer(struct predictms_struct scalar) scalar
local SS	struct predictms_struct scalar
local PC	pointer colvector
local gml 	struct merlin_struct scalar
local cgml 	struct merlin_struct colvector
local Pcm	pointer(struct merlin_struct scalar) colvector

mata:

// illness-death
void predictms_analytic_illd(`SS' S, `RS' from, `Pcm' Pmerlin, `RS' Nobs)
{
	t = S.predtime
	
	if (S.getprobs) {
		
		if (from==1) 	pred = predictms_analytic_illd_prob1(Pmerlin,t,S.enter,S.chips)
		else 		pred = predictms_analytic_illd_prob2(Pmerlin,t,S.enter)
		
		if (S.standardise) 	S.pt = S.pt :+ pred
		else 			S.pt = pred
	}
	
	if (S.getlos | S.getrmst) {
		
		Nqp 	= S.chips
		gq 	= predictms_gq(Nqp)
		qp	= (t:-S.enter) :/ 2 :* J(Nobs,1,gq[,1]') :+ (t:+S.enter) :/2
		if (from==1) { 
			pred	= J(Nobs,3,0)
			for (q=1; q<=Nqp; q++) {
				pred = pred :+ predictms_analytic_illd_prob1(Pmerlin,qp[,q],S.enter,S.chips) :* gq[q,2]
			}
		}
		else {
			pred	= J(Nobs,3,0)
			for (q=1; q<=Nqp; q++) {
				pred = pred :+ predictms_analytic_illd_prob2(Pmerlin,qp[,q],S.enter) :* gq[q,2]
			}
		}
		pred = pred :* (t:-S.enter) :/ 2
		
		if (S.standardise) 	S.los = S.los :+ pred
		else 				S.los = pred
		
	}
		
	if (S.getrmst) {
		if (S.standardise) 	S.rmst = S.rmst :+ pred[,1] :+ pred[,2]
		else 				S.rmst = pred[,1] :+ pred[,2]
	}

	if (S.hasuser) predictms_calc_user(S)	
		
}

`RM' predictms_analytic_illd_prob1(`Pcm' Pmerlin, `RC' t, `RS' enter, `RS' chips)
{
	`gml' gml1, gml2, gml3
	gml1 = *Pmerlin[1]
	gml2 = *Pmerlin[2]
	gml3 = *Pmerlin[3]

	Nobs 	= rows(t)
	result 	= J(Nobs,3,.)

	//1
	ch = (*gml1.Pch[1])(gml1,t) :+ (*gml2.Pch[1])(gml2,t)
 	_editmissing(ch,0)
	result[,1] = exp(-ch)

	//2
	Ngq 	= chips
	gq 	= predictms_gq(Ngq)
	qw	= (t:-enter) :/ 2
	qp	= qw :* J(Nobs,1,gq[,1]') :+ (t:+enter):/2
	pred 	= J(Nobs,1,0)
	
	for (q=1; q<=Ngq; q++) {					
		res = (*gml1.Plogh[1])(gml1,qp[,q]) :- (*gml1.Pch[1])(gml1,qp[,q])
		res = res :- (*gml2.Pch[1])(gml2,qp[,q])
		res = res :- (*gml3.Pch[1])(gml3,t,qp[,q]) :+ (*gml3.Pch[1])(gml3,qp[,q],qp[,q])
		pred = pred :+ exp(res :+ log(gq[q,2]))
	}
	pred = pred :* qw
	
	if (missing(pred)==Nobs) merlin_error("Error in numerical integration; try the -simulate- option")
	
 	_editmissing(pred,0)
	result[,2] = pred
	
	//delayed entry
	if (enter) {
		result[,(1,2)] = result[,(1,2)] :/ exp(-(*gml1.Pch[1])(gml1,J(Nobs,1,enter)))
		result[,(1,2)] = result[,(1,2)] :/ exp(-(*gml2.Pch[1])(gml2,J(Nobs,1,enter)))
	}
	
	//3
	result[,3] = 1 :- result[,1] :- result[,2]
	
	return(result)
}

`RM' predictms_analytic_illd_prob2(`Pcm' Pmerlin, `RC' t, `RS' enter)
{
	`gml' gml
	gml = *Pmerlin[3]
	
	Nobs 		= rows(t)
	pred 		= J(Nobs,3,0)
	ch		= (*gml.Pch[1])(gml,t)
	_editmissing(ch,0)
	pred[,2] 	= exp(-ch)
	
	if (enter) {
		pred[,2] = pred[,2] :/ exp(-(*gml.Pch[1])(gml,J(Nobs,1,enter)))
	}
	pred[,3]	= 1:-pred[,2]
	
	return(pred)
}


// illness-death, standardised
void predictms_analytic_illd_std(`SS' S, `RS' from, `Pcm' Pmerlin, `RS' Nobs)
{
	`RC' t
	`RS' Nt
	`gml' mlob
	`cgml' gml
	gml = J(3,1,mlob)
	
	t = S.predtime
	Nt = rows(t)
	gml[1] = *Pmerlin[1]
	gml[2] = *Pmerlin[2]
	gml[3] = *Pmerlin[3]
	
	if (S.getprobs) {
		if (from==1) {
			for (i=1;i<=Nt;i++) {
				S.pt[i,] = predictms_analytic_illd_p1_std(gml,t[i],S.enter,S.chips)
			}
		}
		else {
			for (i=1;i<=Nt;i++) {
				S.pt[i,] = predictms_analytic_illd_p2_std(gml,t[i],S.enter)
			}
		}
	}
	
	if (S.getlos | S.getrmst) {
		
		Nqp 	= S.chips
		gq 	= predictms_gq(Nqp)
		qp	= (t:-S.enter) :/ 2 :* J(Nobs,1,gq[,1]') :+ (t:+S.enter) :/2
		if (from==1) { 
			pred	= J(Nobs,3,0)
			for (q=1; q<=Nqp; q++) {
				pred = pred :+ predictms_analytic_illd_prob1(Pmerlin,qp[,q],S.enter,S.chips) :* gq[q,2]
			}
		}
		else {
			pred	= J(Nobs,3,0)
			for (q=1; q<=Nqp; q++) {
				pred = pred :+ predictms_analytic_illd_prob2(Pmerlin,qp[,q],S.enter) :* gq[q,2]
			}
		}
		pred = pred :* (t:-S.enter) :/ 2
		
		if (S.standardise) 	S.los = S.los :+ pred
		else 				S.los = pred
		
	}
		
	if (S.getrmst) {
		if (S.standardise) 	S.rmst = S.rmst :+ pred[,1] :+ pred[,2]
		else 				S.rmst = pred[,1] :+ pred[,2]
	}

	if (S.hasuser) predictms_calc_user(S)	
		
}

`RM' predictms_analytic_illd_p1_std(`cgml' gml, `RS' t, `RS' enter, `RS' chips)
{
	Nobs 	= rows(t)
	result 	= J(Nobs,3,.)

	//1
	ch = (*gml[1].Pch[1])(gml[1],t) :+ (*gml[2].Pch[1])(gml[2],t)
 	_editmissing(ch,0)
	result[,1] = exp(-ch)

	//2
	Ngq 	= chips
	gq 	= predictms_gq(Ngq)
	qw	= (t:-enter) :/ 2
	qp	= qw :* J(Nobs,1,gq[,1]') :+ (t:+enter):/2
	pred 	= J(Nobs,1,0)
	
	for (q=1; q<=Ngq; q++) {					
		res = (*gml[1].Plogh[1])(gml[1],qp[,q]) :- (*gml[1].Pch[1])(gml[1],qp[,q])
		res = res :- (*gml[2].Pch[1])(gml[2],qp[,q])
		res = res :- (*gml[3].Pch[1])(gml[3],t,qp[,q]) :+ (*gml[3].Pch[1])(gml,qp[,q],qp[,q])
		pred = pred :+ exp(res :+ log(gq[q,2]))
	}
	pred = pred :* qw
	
	if (missing(pred)==Nobs) merlin_error("Error in numerical integration; try the -simulate- option")
	
 	_editmissing(pred,0)
	result[,2] = pred
	
	//delayed entry
	if (enter) {
		result[,(1,2)] = result[,(1,2)] :/ exp(-(*gml[1].Pch[1])(gml[1],J(Nobs,1,enter)))
		result[,(1,2)] = result[,(1,2)] :/ exp(-(*gml[2].Pch[1])(gml[2],J(Nobs,1,enter)))
	}
	
	//3
	result[,3] = 1 :- result[,1] :- result[,2]
	return(result)
}

`RM' predictms_analytic_illd_p2_std(`cgml' gml, `RC' t, `RS' enter)
{	
	Nobs 		= rows(t)
	pred 		= J(Nobs,3,0)
	ch		= (*(gml[3].Pch[1]))(gml,t)
	_editmissing(ch,0)
	pred[,2] 	= exp(-ch)
	
	if (enter) {
		pred[,2] = pred[,2] :/ exp(-(*gml[3].Pch[1])(gml,J(Nobs,1,enter)))
	}
	pred[,3]	= 1:-pred[,2]	
	return(pred)
}
end
