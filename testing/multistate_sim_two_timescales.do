clear

local N = 30000
local tmax = 10
set obs `N'

gen sex = runiform()>0.5
gen double age = rnormal(50,7)

//time to recurrence
//-> main timescale
survsim rectime, dist(weibull) lambda(0.15) gamma(0.8) cov(age 0.01 sex -0.5)

//time to death
//-> age as the timescale
local ld  0.1*exp(-0.5 * sex)
local gd = 1.2
gen double dtime2 = ((log(exp(-`ld' * age^`gd') * runiform()) / (-`ld'))^(1/`gd'))
gen dtime = dtime - age

gen admin = `tmax'

egen double rmtime = rowmin(rectime dtime admin)
egen double dmtime = rowmin(rectime dtime admin)

gen ri = rectime==rmtime
gen di = dtime==dmtime

//death after relapse
local ld 0.1*exp(-0.5*sex + 0.01*age)
local gd = 1
replace dmtime = ((log(exp(-`ld' * rmtime^`gd') * runiform()) / (-`ld'))^(1/`gd')) if ri==1

replace dmtime = `tmax' if dmtime>=`tmax' & ri==1
replace di = 0 if dmtime==`tmax'

replace di = 1 if ri==1 & dmtime<`tmax'

gen id = _n

msset , id(id) times(rmtime dmtime) states(ri di) cov(age sex)
mat tmat = r(transmatrix)

stset _stop, enter(_start) f(_status)

streg age sex if _trans1==1, dist(weib) nohr
est store m1

streg age sex if _trans3==1, dist(weib) nohr
est store m3

gen stop2 = _stop + age
gen start2 = _start + age
stset stop2, enter(start2) f(_status)

streg sex if _trans2==1, dist(weib) nohr
est store m2

range tvar 0 10 100
predictms , transmat(tmat) models(m1 m2 m3) at(age 50) tscale2(2) time2(50) timevar(tvar) graph




