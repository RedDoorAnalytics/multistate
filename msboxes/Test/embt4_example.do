// using msboxes with ebmt4

local drive c:\
//local drive z:\
//local drive /Users/Michael/Documents

cd "`drive'/multistate/multistate"
adopath ++ "./msset"
adopath ++ "./predictms"

clear all

use "data/ebmt4",clear

mat tmat = (.,1,2,.,3,4\.,.,.,5,6,7\.,.,.,8,9,10\.,.,.,.,11,12\.,.,.,.,.,.\.,.,.,.,.,.)
mat list tmat

msset, id(id) states(recs aes recaes rels srvs) times(rec ae recae rel srv) transmatrix(tmat)

mat list r(transmatrix)
mat list r(freqmatrix)
mat tmat = r(transmatrix)
