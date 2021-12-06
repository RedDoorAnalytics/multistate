
clear
set seed 249857
set obs 10000
gen trt = runiform()>0.5
gen age = rnormal(50,5)
gen agec = age - 50

//simulate from time since entry and (centred) attained age timescales
survsim stime died , maxt(5)  cov(trt -0.5) ///
	hazard(0.1:*1.2:*{t}:^0.2 :* exp(0.1 :* (agec :+ {t})))

//truth
merlin (stime trt rcs(stime, df(1) offset(agec)) 	///
					, family(weibull, failure(died)) timevar(stime)) , 
est store m1

cap range time 0 10 100

predict _prob_m1 , survival at(trt 1 agec 10) timevar(time)
predict _prob_m2 , cif at(trt 1 agec 10) timevar(time)
predict _los_m1 , rmst at(trt 1 agec 10) timevar(time)
predict _los_m2 , timelost at(trt 1 agec 10) timevar(time)

//galahad using simulation
galahad , 	survival models(m1) ///
			at1(trt 1 agec 10)  ///
			timevar(time) 		///
			los 				///
			n(100000)
			
