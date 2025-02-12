//============================================================================//
// cert. script for multistate

local drive /Users/michael/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"
clear all

do ./build/buildmlib.do
mata mata clear

local drive /Users/michael/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software
cd "`drive'/multistate"
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

di as text "Running survival tests"

do ./cert/survival.do

di as text "Running CR tests"

do ./cert/cr.do

di as text "Running Markov tests"

do ./cert/markov.do

di as text "Running semi-Markov tests"

do ./cert/semi-markov.do

di as text "Running msaj tests"

do ./cert/msaj_allcerts.do


di in yellow "Done, no errors"
