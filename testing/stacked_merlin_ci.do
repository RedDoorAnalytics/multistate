
//local drive Z:/
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

do ./build/buildmlib.do
mata mata clear

local drive /Users/Michael/Documents/reddooranalytics/products/multistate
cd "`drive'"
adopath ++ "."
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./stms"
adopath ++ "./msboxes"
adopath ++ "./msaj"
adopath ++ "./graphms"
clear all

tr:do ./build/buildmlib.do
mata mata clear

clear all


use "./data/readmission_MSM.dta", clear

cap tab _trans, gen(_trans)

mat tmat=  (.,1,.,.,.,2   \ ///
.,.,3,.,.,4   \ ///
.,.,.,5,.,6   \ ///
.,.,.,.,7,8   \ ///
.,.,.,.,.,9   \ ///
.,.,.,.,.,.)

mat li tmat
cap range tt 2 5 4 

stset _stop , enter(_start) failure(_status=1) scale(365.25)

constraint 1 _b[_cmp_1_10_1:_cons]=_b[_cmp_1_11_1:_cons]	
constraint 2 _b[_cmp_1_11_1:_cons]=_b[_cmp_1_12_1:_cons]
constraint 3 _b[_cmp_1_12_1:_cons]=_b[_cmp_1_13_1:_cons]

constraint 4 _b[_cmp_1_10_2:_cons]=_b[_cmp_1_11_2:_cons]	
constraint 5 _b[_cmp_1_11_2:_cons]=_b[_cmp_1_12_2:_cons]
constraint 6 _b[_cmp_1_12_2:_cons]=_b[_cmp_1_13_2:_cons]

constraint 7 _b[_cmp_1_10_3:_cons]=_b[_cmp_1_11_3:_cons]	
constraint 8 _b[_cmp_1_11_3:_cons]=_b[_cmp_1_12_3:_cons]
constraint 9 _b[_cmp_1_12_3:_cons]=_b[_cmp_1_13_3:_cons]

merlin (_t  _trans1 _trans3 _trans5 _trans7         /// 1-4 
_trans2 _trans4 _trans6 _trans8 _trans9 /// 5-9  
_trans1#rcs(_t,  df(3) event log orthog) /// 10 (_1 _2 _3)
_trans3#rcs(_t,  df(3) event log orthog) /// 11 (_1 _2 _3)
_trans5#rcs(_t,  df(3) event log orthog) /// 12 (_1 _2 _3)
_trans7#rcs(_t,  df(3) event log orthog) /// 13 (_1 _2 _3)
, nocons constraint(1 2 3 4 5 6 7 8 9  ) /// 
family(loghazard,  failure(_d) ltruncated(_t0)))

set seed 155
predictms, transmatrix(tmat) timevar(tt) probability  ///
        latent n(5000) ci m(20) 

