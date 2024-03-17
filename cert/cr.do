local drive /Users/Michael/My Drive/software/multistate

use "`drive'/data/multistate_example",clear
set seed 98775
keep if runiform()<0.2
// qui {
	msset, id(pid) states(rfi osi) times(rf os)  cr
	mat tmat = (.,1,2\.,.,3\.,.,.) 
	stset _stop, enter(_start) failure(_status==1) scale(12) 
	tab size, gen(sz)

	stmerlin hormon age if _trans1==1, dist(rp) df(3) 
	est store m1

	stmerlin hormon age if _trans2==1, dist(rp) df(3)
	est store m2
// }

// merlin predict check

merlin 	(_t hormon age if _trans1==1, family(rp, df(3) failure(_d)))	///
	(_t hormon age if _trans2==1, family(rp, df(3) failure(_d)))	///
	, 

cap drop tvar
cap range tvar 0 5 200
predict cif1, cif at(age 55) zeros outcome(1) timevar(tvar)
predict cif2, cif at(age 55) zeros outcome(2) timevar(tvar)	

predictms , cr models(m1 m2) 		///
	prob				///
	at1(age 55) 			///
	timevar(tvar) 

assert abs(cif1-_prob_at1_1_2)<1e-07 in 1/200
assert abs(cif2-_prob_at1_1_3)<1e-07 in 1/200


//==========================================================================//
// standardise


predictms , cr models(m1 m2) 		///
	prob				///
	at1(age 55) 			///
	timevar(tvar) 			///
	 standardise 
rename _prob* prob*

predictms , cr models(m1 m2) 		///
	prob				///
	at1(age 55) 			///
	timevar(tvar) 			///
	 standardise aj

assert abs(_prob_at1_1_1 - prob_at1_1_1) <1e-03 in 1/100

keep if runiform()<0.5
cap drop tvar2
cap range tvar2 0 5 10
cap drop prob*
predictms , cr models(m1 m2) 		///
	prob				///
	at1(age 55) 			///
	timevar(tvar2) 			///
	 standardise 
rename _prob* prob*

predictms , cr models(m1 m2) 		///
	prob				///
	at1(age 55) 			///
	timevar(tvar2) 			///
	 standardise simulate
	 
assert abs(_prob_at1_1_1 - prob_at1_1_1) <1e-03 in 1/10

