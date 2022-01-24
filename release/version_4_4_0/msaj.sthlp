{smcl}
{* *! version 1.0.0 02/03/2020}{...}
{vieweralsosee "[multistate] multistate" "help multistate"}{...}
{vieweralsosee "[multistate] msboxes" "help msboxes"}{...}
{vieweralsosee "[multistate] msaj" "help msaj"}{...}
{vieweralsosee "[multistate] predictms" "help predictms"}{...}
{vieweralsosee "[multistate] graphms" "help graphms"}{...}
{vieweralsosee "[merlin] stmerlin" "help stmerlin"}{...}
{vieweralsosee "[merlin] merlin" "help merlin"}{...}
{viewerjumpto "Syntax" "msaj##syntax"}{...}
{viewerjumpto "Description" "msaj##description"}{...}
{viewerjumpto "Options" "msaj##options"}{...}
{viewerjumpto "Examples" "msaj##examples"}{...}

{title:Title}

{p2colset 5 13 16 2}{...}
{p2col :{hi:msaj} {hline 2}}Aalen-Johansen estimates of transition probabilities{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}{cmd:msaj} {ifin} {cmd:,} [{it:options}]

{marker options}{...}
{synoptset 29 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt by(varname)}}calculate by {it:varname}{p_end}
{synopt :{opt transm:atrix(matname)}}name of transition matrix{p_end}
{synopt :{opt cr}}shortcut for competing risks{p_end}
{synopt :{opt from(#)}}starting state from which predictions are made{p_end}
{synopt :{opt lt:runcated(#)}}time at which predictions are made from{p_end}
{synopt :{opt exit(#)}}time at which predictions are made until{p_end}
{synopt :{opt ci}}calculate confidence intervals (for probability estimates only){p_end}
{synopt :{opt se}}calculate standard errors (for probability estimates only){p_end}
{synopt :{opt los}}calculate length of stay estimates{p_end}
{synoptline}
{p 4 6 2}


{title:Description}

{pstd}
{cmd:msaj} calculates the Aalen-Johansen estimates of the transition probabilities. Probability estimates are output only for 
observed event times, i.e. where _d == 1 (within the boundaries of {opt ltruncated(#)} and {opt exit(#)}). Length of stay estimates 
are output for all observed times (within the boundaries of {opt ltruncated(#)} and {opt exit(#)}). Before using {cmd:msaj}, you 
should use {cmd:msset} and then {cmd:stset}. See the example below.
{p_end}

{pstd}
Probability estimates are stored in new variables P_AJ_*, where star represents the numbered states in the multi-state model. {p_end}

{pstd}
Standard errors are stored in new variables P_AJ_*_se, when the {opt se} option is specified. {p_end}

{pstd}
Confidence intervals are stored in new variables P_AJ_*_lci and P_AJ_*_uci, when the {opt ci} option is specified. {p_end}

{pstd}
Length of stay estimates are stored in new variables LOS_AJ_*, when the {opt los} option is specified. {p_end}

{phang}
{cmd:msaj} is part of the {helpb multistate} package.
{p_end}


{title:Options}

{phang} 
{opt by(varname)} will calculate estimates separately by {it:varname}. Variable must be numeric. Estimates are not calculated 
for {it:varname} equal to missing.
{p_end}

{phang}
{opt transmatrix(matname)} specifies the transition matrix used in the multi-state model that was fitted. The matrix must 
have missing entries on the diagonal. Transitions must be numbered as an increasing sequence of integers from 1,...,K. 
This transition matrix should be the same as that used/produced by {cmd: msset}. {opt transmatrix(matname)} or {opt cr} 
must be specified.
{p_end}

{phang}
{opt cr} states that it is a competing risks analysis, i.e. all transitions are in the first row of the transition matrix. 
This means that you do not need to specify the {opt transmatrix(matname)} option.
{p_end}

{phang}
{opt from(#)} specifies the state you wish to make predicitons from (starting state). Default is state 1. The state must not 
be an absorbing state. 
{p_end}

{phang}
{opt ltruncated(#)} specifies the time at which predictions are made from. It is time s when estimating the probability P(s,t). 
The default is time 0. This must be strictly before the exit time.
{p_end}

{phang}
{opt exit(#)} specifies the time at which predictions are made until. Default is the last observed event time, i.e. where 
{cmd:_d == 1}. When the {opt by(varname)} option is also specified, the default exit time is the last observed event time 
from any group, excluding when {it:varname} is missing. This default is used for all groups and so the estimates may be 
extroplated past the last event time for some groups. Exit time must be strictly after the ltruncated time. Exit time must be 
on or before the last observed time in the data.
{p_end}

{phang}
{opt ci} calculates confidence intervals for the transition probabilities using Greenwood standard errors. See Andersen 
et al (1993) and Allignol et al (2011) for details on how the standard errors were calculated. Estimates of confidence 
intervals are truncated to ensure they lie in the [0,1] boundary.
{p_end}

{phang} 
{opt se} calculates Greenwood standard errors for the transition probabilities. See Andersen et al (1993) and Allignol 
et al (2011) for details on how the standard errors were calculated. 
{p_end}

{phang}
{opt los} calculates length of stay estimates. These are calculated at all observed times (within the boundaries of 
{opt ltruncated(#)} and {opt exit(#)}). They are calculated by integrating the area under the step probability function.
{p_end}


{title:Examples}

{pstd}Load example dataset:{p_end}
{phang}{stata "use http://fmwww.bc.edu/repec/bocode/m/multistate_example":. use http://fmwww.bc.edu/repec/bocode/m/multistate_example}{p_end}

{pstd}{cmd:msset} the data:{p_end}
{phang}{stata "msset, id(pid) states(rfi osi) times(rf os)":. msset, id(pid) states(rfi osi) times(rf os)}{p_end}

{pstd}Store the transition matrix:{p_end}
{phang}{stata "mat tmat = r(transmatrix)":. mat tmat = r(transmatrix)}{p_end}

{pstd}{cmd:stset} the data using the variables created by {cmd:msset}{p_end}
{phang}{stata "stset _stop, enter(_start) failure(_status=1)":. stset _stop, enter(_start) failure(_status=1)}{p_end}

{pstd}Calculate transition probabilities using {cmd:msaj}{p_end}
{phang}{stata "msaj, transmat(tmat) ci":. msaj, transmat(tmat) ci}{p_end}

{pstd}Probability in State 1 (alive) with confidence intervals {p_end}
{phang}{stata "line P_AJ_1* _t, sort connect(stairstep stairstep stairstep)":. line P_AJ_1* _t, sort connect(stairstep stairstep stairstep)}{p_end}

{pstd}Probability in State 2 (recurrence) with confidence intervals {p_end}
{phang}{stata "line P_AJ_2* _t, sort connect(stairstep stairstep stairstep)":. line P_AJ_2* _t, sort connect(stairstep stairstep stairstep)}{p_end}

{pstd}Probability in State 3 (dead) with confidence intervals {p_end}
{phang}{stata "line P_AJ_3* _t, sort connect(stairstep stairstep stairstep)":. line P_AJ_3* _t, sort connect(stairstep stairstep stairstep)}{p_end}

{pstd}Calculate transition probabilities from state 2 at time 50 until time 150 using {cmd:msaj}{p_end}
{phang}{stata "cap drop P_AJ*":. cap drop P_AJ*}{p_end}
{phang}{stata "msaj, transmat(tmat) from(2) ltruncated(50) exit(150) ci":. msaj, transmat(tmat) from(2) ltruncated(50) exit(150) ci}{p_end}

{pstd}Probability in State 2 (recurrence) from state 2 at time 50 until time 150 {p_end}
{phang}{stata "line P_AJ_2* _t, sort connect(stairstep stairstep stairstep)":. line P_AJ_2* _t, sort connect(stairstep stairstep stairstep)}{p_end}

{pstd}Calculate length of stay estimates using {cmd:msaj} {p_end}
{phang}{stata "cap drop P_AJ*":. cap drop P_AJ*}{p_end}
{phang}{stata "msaj, transmat(tmat) los":. msaj, transmat(tmat) los}{p_end}

{pstd}Length of stay in State 1 (alive){p_end}
{phang}{stata "line LOS_AJ_1 _t, sort":. line LOS_AJ_1 _t, sort}{p_end}

{pstd}Length of stay in State 2 (recurrance){p_end}
{phang}{stata "line LOS_AJ_2 _t, sort":. line LOS_AJ_2 _t, sort}{p_end}


{title:Authors}

{pstd}
Micki Hill, University of Leicester, UK.
({browse "mailto:mh594@leicester.ac.uk":mh594@leicester.ac.uk})

{pstd}
Paul Lambert, University of Leicester, UK and Karolinska Institutet, Stockholm, Sweden.
({browse "mailto:paul.lambert@leicester.ac.uk":paul.lambert@leicester.ac.uk})

{pstd}
Michael Crowther, Red Door Analytics and Karolinska Institutet, Stockholm, Sweden.
({browse "michael@reddooranalytics.se":michael@reddooranalytics.se})


{title:References}

{phang}
Allignol A, Schumacher M, Beyersmann J. Empirical transition matrix of multi-state models: The etm package. {it:Journal of Statistical Software} 2011;38(4):1-15.

{phang}
Andersen PK, Borgan O, Gill RD, Keiding N. {it:Statistical models based on counting processes}. Springer Series in Statistics 1993. 

{phang}
Putter H, Fiocco M, Geskus RB. Tutorial in biostatistics: competing risks and multi-state models. {it:Statistics in Medicine} 2007;26:2389-2430.{p_end}
