
clear
set seed 249857
set obs 10000
gen trt = runiform()>0.5
gen age = rnormal(50,5)
gen agec = age - 50

//simulate from time since entry and (centred) attained age timescales
survsim stime died , maxt(5)  cov(trt -0.5) ///
	hazard(0.1:*1.2:*{t}:^0.2 :* exp(0.1 :* (agec :+ {t})))

//true model
merlin (stime trt rcs(stime, df(1) offset(agec)) 	///
		, family(weibull, failure(died)) timevar(stime)) , 
est store m1

cap range time 0 10 100

//bench predictions from merlin
predict _prob_m1 , survival at(trt 1 agec 10) timevar(time)
predict _prob_m2 , cif at(trt 1 agec 10) timevar(time)
predict _los_m1 , rmst at(trt 1 agec 10) timevar(time)
predict _los_m2 , timelost at(trt 1 agec 10) timevar(time)

//predictms - numerical integration
predictms , singleevent models(m1) 	///
		at1(trt 1 agec 10)  	///
		timevar(time) 		///
		prob los 			
rename _prob_at1* prob_at1*
rename _los_at1* los_at1*

assert (_prob_m1 - prob_at1_1_1)<1e-8 in 1/100
assert (_prob_m2 - prob_at1_1_2)<1e-8 in 1/100
assert (_los_m1 - los_at1_1_1)<1e-8 in 2/100	//missings in 1 at t=0
assert (_los_m2 - los_at1_1_2)<1e-8 in 2/100

//predictms using simulation
predictms , singleevent models(m1) 	///
		at1(trt 1 agec 10)  	///
		timevar(time) 		///
		prob los simulate 	///
		n(1000000)	

assert (_prob_m1 - _prob_at1_1_1)<1e-3 in 1/100
assert (_prob_m2 - _prob_at1_1_2)<1e-3 in 1/100
assert (_los_m1 - _los_at1_1_1)<1e-2 in 2/100
assert (_los_m2 - _los_at1_1_2)<1e-2 in 2/100

//==========================================================================//
//standardised predictions 

drop if _n>1000

predict sstd1, surv standardise timevar(time)

predictms , singleevent models(m1) 	///
		at1()  			///
		timevar(time) 		///
		prob standardise
		
assert (sstd1 - _prob_at1_1_1)<1e-8 in 1/100
		
predictms , singleevent models(m1) 	///
		at1()  			///
		timevar(time) 		///
		surv standardise ci
		
assert (sstd1 - _survival_at1_1_2)<1e-8 in 1/100
