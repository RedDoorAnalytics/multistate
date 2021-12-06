
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

do ./build/buildmlib.do
mata mata clear

clear all

di as text "Running Markov tests"

do ./cert/markov.do

di as text "Running semi-Markov tests"

do ./cert/semi-markov.do

di as text "Running msaj tests"

do ./cert/msaj_allcerts.do


di in yellow "Done, no errors"
