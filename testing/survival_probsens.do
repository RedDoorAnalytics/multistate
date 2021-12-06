//example showing how to incorporate a dummy covariate 
//can vary the value of the covariate and its coefficient

clear
set seed 249857
set obs 1000
gen trt = runiform()>0.5

survsim stime died , maxt(5) cov(trt -0.5) ///
	dist(weib) lambda(0.1) gamma(1.5) 

//create dummy variable
gen dummy = 1	
	
//Use @0 to constrain the parameter of dummy to be 0 -> for later manipulation
merlin (stime trt dummy@0 , family(weib, failure(died))) , 
est store m1

//timevar to make predictions at
cap range time 0 5 100

//predict survival using original model
galahad , 	survival models(m1) 		///
			at1(trt 1 dummy 1)  		///
			timevar(time) 				///
			los 						///
			n(100000)
//rename it to compare to a second call
rename _prob_at1_1_1 prob1

//adapt coefficient of dummy so it has a beneficial effect on survival
mat b = e(b)
mat b2 = b
mat b2[1,2] = -0.5	//change coefficient of dummy
mat list b2
erepost b = b2		//ssc install erepost
est store m1		//restore the model object to pass to galahad

//galahad using simulation
galahad , 	survival models(m1) 		///
			at1(trt 1 dummy 1)  		///
			timevar(time) 				///
			los 						///
			n(100000)
			
//show the effect
scatter _prob_at1_1_1 prob1 time
