local drive /Users/Michael/Documents
//cd "`drive'/multistate/multistate"
// cd "Z:/My Documents/Multistate repository/multistate"
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./stms"
adopath ++ "./msboxes"
adopath ++ "./msaj"
clear all

do ./build/buildmlib.do
mata mata clear


* Set global tolerance
global tol 1e-10


di as text "msaj 1. Options"
do ./cert/msaj_1_options.do

di as text "msaj 2. Errors"
do ./cert/msaj_2_errors.do

di as text "msaj 3. General"
do ./cert/msaj_3_general.do

di as text "msaj 4. HAI"
do ./cert/msaj_4_hai.do

di as text "msaj 5. EBMT"
do ./cert/msaj_5_ebmt.do

di as text "msaj 6. EBMT CR"
do ./cert/msaj_6_ebmt_cr.do

di as text "msaj 7. SIMBI"
do ./cert/msaj_7_simbi.do


di in yellow "msaj checks done, no errors"
macro drop tol
