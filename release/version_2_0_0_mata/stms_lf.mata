*! version 1.0.0 ?????2014 MJC


//!! main touse for reading stuff in
// could add cov indicators in setup

version 12.1

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix
local TR 	transmorphic
mata:

struct stms_main {
	`RS' Ntrans, Nobs, hasbhaz
	real colvector t, d, t0, bhazard
	real rowvector delentry
	`NM' tousevars, modeltrans
	`TR' touseindex, t0index
}

void stms_setup()
{
	struct stms_main scalar S
	pointer scalar p

	//get tempname and create struct
	stata("tempname stms_struct")
	rmexternal(st_local("stms_struct"))
	p = crexternal(st_local("stms_struct"))

	S.t = st_data(.,"_t")
	S.d = st_data(.,"_d")
	S.delentry = strtoreal(tokens(st_local("delentry")))
	S.Nobs = rows(S.d)
	S.Ntrans = strtoreal(st_local("Ntrans"))
	S.hasbhaz = st_local("bhazard")!=""
	if (S.hasbhaz) S.bhazard = st_data(.,st_local("bhazard"))
	S.modeltrans = (tokens(st_local("models"))'):==J(S.Ntrans,1,("exp","weibull","gompertz","llogistic","lnormal","gamma","stpm2"))
	
	S.touseindex = asarray_create("real",1)
	for (i=1;i<=S.Ntrans;i++) {
		asarray(S.touseindex,i,selectindex(st_data(.,st_local("touse"+strofreal(i)))))
	}
		
	if (sum(S.delentry)) {
		S.t0 = st_data(.,"_t0")
		S.t0index = asarray_create("real",1)
		for (i=1;i<=S.Ntrans;i++) {
			if (S.delentry[1,i]) {
				asarray(S.t0index,i,selectindex(st_data(.,st_local("t0touse"+strofreal(i)))))
			}
		}
	}

	/*
	S.orthog = st_local("rmatopt")!=""
	if (S.orthog) S.rmatrix = st_matrix(st_local("rmat"))
	S.knots = strtoreal(tokens(st_local("ln_bknots")))
	S.Nsplines = strtoreal(st_local("Nsplines"))	
	*/
	
	/*
	*/
	
	//Done 	
	swap((*p), S)
}

void stms_lf(	transmorphic scalar M, 		///
				real scalar todo, 			///
				real rowvector b, 			///
				real colvector lnfj, 		///
				real matrix S, 				///
				real matrix H)
{
	pointer(struct stms_main scalar) scalar pS
	struct stms_main scalar A
	real colvector index, t0index
	
	pS = &moptimize_util_userinfo(M,1)
	A = *pS
	
	lnfj = J(A.Nobs,1,0)
	ind = 1

	for (i=1; i<=A.Ntrans; i++) {

		index = asarray(A.touseindex,i)
		if (A.delentry[1,i]) t0index = asarray(A.t0index,i)

		if (A.modeltrans[i,1]) {
			ln_lambda = moptimize_util_xb(M,b,ind++)
			if (rows(ln_lambda)==1) ln_lambda = J(A.Nobs,1,ln_lambda)

			if (A.hasbhaz)  lnfj[index] = A.d[index] :* log(exp(ln_lambda[index]) :+ A.bhazard[index]) :- exp(ln_lambda[index]) :* A.t[index]
			else lnfj[index] = A.d[index] :* (ln_lambda[index] :+ log(A.t[index])) :- exp(ln_lambda[index]) :* A.t[index]
			if (A.delentry[1,i]) lnfj[t0index] = lnfj[t0index] :+ exp(ln_lambda[t0index]) :* A.t0[t0index]
			continue
		}
		else if (A.modeltrans[i,2]) {
			ln_lambda = moptimize_util_xb(M,b,ind++)
			if (rows(ln_lambda)==1) ln_lambda = J(A.Nobs,1,ln_lambda)
			gam = exp(moptimize_util_xb(M,b,ind++))
			if (rows(gam)==1) gam = J(A.Nobs,1,gam)
			
			if (A.hasbhaz) lnfj[index] = A.d[index] :* log(exp(ln_lambda[index]) :* gam[index] :* A.t[index]:^(gam[index]:-1) :+ A.bhazard[index]) :- exp(ln_lambda[index]) :* A.t[index] :^ gam[index]
			else lnfj[index] = A.d[index] :* (ln_lambda[index] :+ log(gam[index]) :+ (gam[index]:-1):*log(A.t[index]) :+ log(A.t[index])) :- exp(ln_lambda[index]) :* A.t[index] :^ gam[index]
			if (A.delentry[1,i]) lnfj[t0index] = lnfj[t0index] :+ exp(ln_lambda[t0index]) :* A.t0[t0index] :^ gam[t0index]
			continue
		}
		else if (A.modeltrans[i,3]) {
			ln_lambda = moptimize_util_xb(M,b,ind++)
			if (rows(ln_lambda)==1) ln_lambda = J(A.Nobs,1,ln_lambda)
			gam = moptimize_util_xb(M,b,ind++)
			if (rows(gam)==1) gam = J(A.Nobs,1,gam)
			
			if (A.hasbhaz) lnfj[index] = A.d[index] :* log(exp(ln_lambda[index]) :* exp(gam[index] :* A.t[index]) :+ A.bhazard[index]) :- exp(ln_lambda[index]) :* (1:/gam[index]) :* (exp(gam[index]:*A.t[index]):-1)
			else lnfj[index] = A.d[index] :* (ln_lambda[index] :+ log(gam[index]) :+ gam[index]:*A.t[index] :+ log(A.t[index])) :- exp(ln_lambda[index]) :* (1:/gam[index]) :* (exp(gam[index]:*A.t[index]):-1)
			if (A.delentry[1,i]) lnfj[t0index] = lnfj[t0index] :+ exp(ln_lambda[t0index]) :* (1:/gam[t0index]) :* (exp(gam[t0index]:*A.t0[t0index]):-1)
			continue
		}
		else if (A.modeltrans[i,4]) {
			lambda = exp(-moptimize_util_xb(M,b,ind++))
			if (rows(lambda)==1) lambda = J(A.Nobs,1,lambda)
			gam = exp(moptimize_util_xb(M,b,ind++))
			if (rows(gam)==1) gam = J(A.Nobs,1,gam)
			
			if (A.hasbhaz) lnfj[index] = A.d[index] :* log( A.bhazard[index] :* (1 :+ (lambda[index]:*A.t[index]):^(1:/gam[index])):^(-1)     :+  (lambda[index]:^(1:/gam[index]):* A.t[index]:^(1:/gam[index]:-1)):/(gam[index]:*(1:+(lambda[index]:*A.t[index]):^(1:/gam[index])):^2)) :- (1:-A.d[index]) :* log(1 :+ (lambda[index]:*A.t[index]):^(1:/gam[index]))
			else lnfj[index] = A.d[index] :* log(A.t[index]) :+ A.d[index] :* log((lambda[index]:^(1:/gam[index]):* A.t[index]:^(1:/gam[index]:-1)):/(gam[index]:*(1:+(lambda[index]:*A.t[index]):^(1:/gam[index])):^2)) :- (1:-A.d[index]) :* log(1 :+ (lambda[index]:*A.t[index]):^(1:/gam[index]))
			if (A.delentry[1,i]) lnfj[t0index] = lnfj[t0index] :+ log(1 :+ (lambda[t0index]:*A.t0[t0index]):^(1:/gam[t0index]))
			continue
		}
		else if (A.modeltrans[i,5]) {
			lambda = moptimize_util_xb(M,b,ind++)
			if (rows(lambda)==1) lambda = J(A.Nobs,1,lambda)
			gam = exp(moptimize_util_xb(M,b,ind++))
			if (rows(gam)==1) gam = J(A.Nobs,1,gam)
			
			if (A.hasbhaz) lnfj[index] = A.d[index] :* log( A.bhazard[index] :* (1 :- normal((log(A.t[index]):-lambda[index]):/gam[index]))  :+  1:/(A.t[index]:*gam[index]:*sqrt(2:*pi())) :* exp( -1 :/ (2:* gam[index]:^2) :* (A.t[index]:-lambda[index]):^2)) :+ (1:-A.d[index]) :* log(1 :- normal((log(A.t[index]):-lambda[index]):/gam[index]))
			else lnfj[index] = A.d[index] :* log(A.t[index]) :+ A.d[index] :* log(1:/(A.t[index]:*gam[index]:*sqrt(2:*pi())) :* exp( -1 :/ (2:* gam[index]:^2) :* (A.t[index]:-lambda[index]):^2)) :+ (1:-A.d[index]) :* log(1 :- normal((log(A.t[index]):-lambda[index]):/gam[index]))
			if (A.delentry[1,i]) lnfj[t0index] = lnfj[t0index] :- log(1 :- normal((log(A.t0[t0index]):-lambda[t0index]):/gam[t0index]))
			continue
		}
		else if (A.modeltrans[i,6]) {
			mu = moptimize_util_xb(M,b,ind++)
			if (rows(mu)==1) mu = J(A.Nobs,1,mu)
			sigma = exp(moptimize_util_xb(M,b,ind++))
			if (rows(sigma)==1) sigma = J(A.Nobs,1,sigma)
			kapp = moptimize_util_xb(M,b,ind++)
			if (rows(kapp)==1) kapp = J(A.Nobs,1,kapp)

			gam = abs(kapp[index]):^(-2)
			z = sign(kapp[index]) :* ( log(A.t[index]) :-mu[index]  ):/sigma[index]
			u = gam :* exp(abs(kapp[index]):*z)
			pdf1 = (gam:^gam) :/ (sigma[index] :* A.t[index] :* sqrt(gam):*gamma(gam)) :* exp(z :* sqrt(gam) :- u)
			pdf2 = 1:/ (sigma[index] :* A.t[index] :* sqrt(2:* pi())) :* exp(-(z:^2):/2)
			surv1 = 1:-gammap(gam,u)
			surv2 = 1:- normal(z)
			surv3 = gammap(gam,u)

			if (A.hasbhaz) lnfj[index] = A.d[index] :* log( (kapp[index]:!=0) :* pdf1 :+ (kapp[index]:==0) :* pdf2 :+ A.bhazard[index] :* ((kapp[index]:>0) :* surv1 :+ (kapp[index]:==0) :* surv2 :+ (kapp[index]:<0) :* surv3)) :+ (1:-A.d[index]) :* log((kapp[index]:>0) :* surv1 :+ (kapp[index]:==0) :* surv2 :+ (kapp[index]:<0) :* surv3)
			else lnfj[index] = A.d[index] :* log(A.t[index]) :+ A.d[index] :* log((kapp[index]:!=0) :* pdf1 :+ (kapp[index]:==0) :* pdf2) :+ (1:-A.d[index]) :* log((kapp[index]:>0) :* surv1 :+ (kapp[index]:==0) :* surv2 :+ (kapp[index]:<0) :* surv3)
			if (A.delentry[1,i]) {
				gam0 = abs(kapp[t0index]):^(-2)
				z0 = sign(kapp[t0index]) :* ( log(A.t0[t0index]) :-mu[t0index]  ):/sigma[t0index]
				u0 = gam0 :* exp(abs(kapp[t0index]):*z0)
				surv10 = 1:- gammap(gam0,u0)
				surv20 = 1:- normal(z0)
				surv30 = gammap(gam0,u0)
				lnfj[t0index] = lnfj[t0index] :- log((kapp[t0index]:>0) :* surv10 :+ (kapp[t0index]:==0) :* surv20 :+ (kapp[t0index]:<0) :* surv30)
			}
			continue
		}
		else if (A.modeltrans[i,7]) {
			xb = moptimize_util_xb(M,b,ind++)
			dxb = moptimize_util_xb(M,b,ind++)
			if (A.delentry[1,i]) s0xb = moptimize_util_xb(M,b,ind++)

			if (A.hasbhaz)  lnfj[index] = A.d[index] :* (log(dxb[index]) :+ xb[index] :+ A.bhazard[index]) :- exp(xb[index])
			else lnfj[index] = A.d[index] :* (log(dxb[index]) :+ xb[index]) :- exp(xb[index])
			if (A.delentry[1,i]) lnfj[t0index] = lnfj[t0index] :+ exp(s0xb[t0index])
		}

	}

	if (todo==0) return

}
end

exit
