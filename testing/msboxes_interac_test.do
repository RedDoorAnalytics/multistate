
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

do ./build/buildmlib.do
mata mata clear

local drive /Users/Michael/Documents
cd "`drive'/multistate/multistate"
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

clear
use http://fmwww.bc.edu/repec/bocode/m/multistate_example
msset, id(pid) states(rfi osi) times(rf os)
matrix tmat = r(transmatrix)
	 
//Specify timevar
range tt 0 500 101

// Run the msboxes command
msboxes, transmatrix(tmat) id(pid) xvalues(0.2 0.7 0.45) yvalues(0.7 0.7 0.2) ///
 statenames("Surgery" "Relapse" "Dead") ///
 transnames("h1" "h2" "h3") freqat(tt) scale(365.25)  ///
 interactive //jsonpath("$N\4.Json_files")
