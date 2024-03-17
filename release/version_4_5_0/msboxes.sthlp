{smcl}
{* *! version 0.2.0}{...}
{vieweralsosee "[multistate] multistate" "help multistate"}{...}
{vieweralsosee "[multistate] msset" "help msset"}{...}
{vieweralsosee "[multistate] msaj" "help msaj"}{...}
{vieweralsosee "[multistate] predictms" "help predictms"}{...}
{vieweralsosee "[multistate] graphms" "help graphms"}{...}
{vieweralsosee "[merlin] stmerlin" "help stmerlin"}{...}
{vieweralsosee "[merlin] merlin" "help merlin"}{...}
{viewerjumpto "Syntax" "msboxes##syntax"}{...}
{viewerjumpto "Description" "msboxes##description"}{...}
{viewerjumpto "Options" "msboxes##options"}{...}
{viewerjumpto "Examples" "msboxes##examples"}{...}
{title:Title}

{p2colset 5 16 19 2}{...}
{p2col :{hi:msboxes} {hline 2}}Simple plot to summarise states and transitions in a multi-state model{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}{cmd:msboxes} [{cmd:,} {it:options}]

{synoptset 29 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt boxh:eight(#)}}height of boxes{p_end}
{synopt :{opt boxw:idth(#)}}width of boxes{p_end}
{synopt :{opt id(varname)}}name of subject ID variable{p_end}
{synopt :{opt grid}}add a grid to the plot{p_end}
{synopt :{opt staten:ames(string)}}List of names of states{p_end}
{synopt :{opt transn:ames(string)}}List of names of transitions{p_end}
{synopt :{opt transm:at}}name of transition matrix{p_end}
{synopt :{opt yran:ge(numlist)}}range of y-axis{p_end}
{synopt :{opt xran:ge(numlist)}}range of x-axis{p_end}
{synopt :{opt ysize(#)}}y size of plot{p_end}
{synopt :{opt xsize(#)}}x size of plot{p_end}
{synopt :{opt yval:ues(numlist)}}y values of the centre of each box{p_end}
{synopt :{opt xval:ues(numlist)}}x values of the centre of each box{p_end}
{synopt :{opt freqat(varname)}}time points at which to calculate frequencies in each state; see details{p_end}
{synopt :{opt scale(#)}}scales the time variable specified in {cmd:freqat()}{p_end}
{synopt :{opt interactive}}create a JSON file for interactive plotting{p_end}
{synopt :{opt jsonpath(string)}}path to save the JSON file; default is current directory{p_end}
{synoptline}
{p 4 6 2}


{marker description}{...}
{title:Description}

{pstd}
{cmd:msboxes} is a simple descriptive tool to summarise data to be used in a multistate model. It will plot boxes for each 
state and summarise the number at risk at the start and end of follow-up. Transitions are denoted by arrows between the 
boxes and show the number of subjects that transition between the different states. Before using {cmd:msboxes} you should 
use {helpb msset} and then {helpb stset}. See the examples below.{p_end}

{pstd}
By default the boxes are plotted on a (0,1) (0,1) grid and the user must give sensible values for the centre of each box. 
There are simple rules on how to join the boxes with arrows and where to place the text.
{p_end}

{pstd}
{cmd:msboxes} can produce JSON files for use with interactive multi-state model visualisations.
{p_end}

{phang}
{cmd:msboxes} is part of the {helpb multistate} package.
{p_end}


{marker options}{...}
{title:Options}

{phang}
{opt boxheight(#)} Height of the boxes to denote states. Default is 0.3{it:varname}. 

{phang}
{opt boxwidth(#)} Width of the boxes to denote states. Default is 0.2{it:varname}. 

{phang}
{opt grid} Add a grid to the plot. Useful when trying to get the plot to look nice by showing where it my be a good 
idea to move the boxes. 

{phang}
{opt id(varname)} Name of id variable.

{phang}
{opt statenames(string)} The names of each state. These should be given as a list of names, for example
{cmd:statenames("Healthy" "Diseased" "Dead before disease" "Dead with disease")}.

{phang}
{opt transnames(string)} The name of each transition. These should be given as a list of names, for example 
{cmd:transnames("Transition 1" " Transition 2" " Transition 3")}.

{phang}
{opt transmatrix(matrix)}  specifies the transition matrix used in the multi-state model that was fitted. This must be an 
upper triangular matrix (with diagonal and lower triangle elements coded missing). Transitions must be numbered as an 
increasing sequence of integers from 1,...,K. This transition matrix should be the same as that used/produced by {helpb msset}.

{phang}
{opt yrange(numlist)} gives the range of the y-axis. By default this is (0 1).

{phang}
{opt xrange(numlist)} gives the range of the x-axis. By default this is (0 1).

{phang}
{opt ysize(#)} gives the ysize of the plot. Default taken from the default scheme.(0 1).

{phang}
{opt xsize(#)} gives the xsize of the plot. Default taken from the default scheme.(0 1).

{phang}
{opt yvalues(numlist)} gives the y locaion of the centre of each box for each state. It should be a {it:numlist} of length {it:K}, where {it:K} is the number of states. 

{phang}
{opt xvalues(numlist)} gives the x locaion of the centre of each box for each state. It should be a {it:numlist} of length {it:K}, where {it:K} is the number of states. 

{phang}
{opt freqat(varname)} calculates the number of individuals being in each state at times specified in the provided {cmd:varname}. 
It also calculates the number of individuals that have made each transition by time {it:t}. {cmd:freqat()} takes as argument a 
variable that supplies the time points of interest. The calculated measures are returned in {cmd:r(frequencies)}. In the case that 
options {cmd:interactive} and {cmd:jsonpath()} are used, this matrix is automatically saved in the produced JSON file.

{phang}
{opt scale(#)} rescale time values provided at {cmd:freqat()}. For example, if {cmd:freqat()} provides time points on a daily timescale 
but {cmd:_start}, {cmd:_stop} are on a year timescale, we can to provide {cmd:scale(365.25)}. Default value is 1.

{phang}
{opt interactive} create a JSON file named {cmd:"msboxes"}. This jsonfile can be used as input to the MSMplus web tool for 
the interactive exploration of the multi-state graph.

{phang}
{opt jsonpath(string)} If option {cmd:interactive} is defined, the user has to specify the position in the folder that the 
json file will be saved by providing the pathway. For example, jsonpath("C:\documents\multistate").


{title:Remarks}

{pstd}
This is a fairly basic implementation and aims to give a quick summary of the data that feeds into the multi-state model. 
It is not aimed to give publication quality figures. 
You may need some trial and error of the location and size of the boxes.
{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
This dataset contains information on 2982 patients with breast cancer. Baseline is defined as time of surgery, and patients 
can experience relapse, relapse then death, or death with no relapse. Time of relapse is stored in {cmd:rf}, with event indicator 
{cmd:rfi}, and time of death is stored in {cmd:os}, with event indicator {cmd:osi}.
{p_end}

{pstd}
{bf:1. Illness death model:}
{p_end}

{cmd:    . use http://fmwww.bc.edu/repec/bocode/m/multistate_example}
{cmd:    . msset, id(pid) states(rfi osi) times(rf os)}
{cmd:    . matrix tmat = r(transmatrix)}
{cmd:    . msboxes, transmat(tmat) id(pid) xvalues(0.2 0.7 0.45) yvalues(0.7 0.7 0.2) ///}
{cmd:    >          statenames("Surgery" "Relapse" "Dead")}
{phang}{it:({stata msboxes_examples 1:click to run})}{p_end}

{pstd}
{bf:2. Extended illness death model. Separate deaths before/after recurrence.}
{p_end}

{cmd:    . use http://fmwww.bc.edu/repec/bocode/m/multistate_example}
{cmd:    . matrix tmat = (.,1,2,. \ .,.,.,3 \ .,.,.,. \ .,.,.,.)}
{cmd:    . matrix list tmat}
{cmd:    . msset, id(pid) states(rfi osi osi) times(rf os os) transmatrix(tmat)}
{cmd:    . msboxes, transmatrix(tmat) id(pid)                  ///}
{cmd:    >          xvalues(0.2 0.7 0.2 0.7)                   ///}
{cmd:    >          yvalues(0.7 0.7 0.2 0.2)                   ///}
{cmd:    >          statenames(Surgery Relapse Dead Dead)      ///}
{cmd:    >          boxheight(0.2) yrange(0.09 0.81) ysize(3)}
{phang}{it:({stata msboxes_examples 2:click to run})}{p_end}


{title:Authors}	

{pstd}
Paul Lambert, University of Leicester, UK.
({browse "mailto:paul.lambert@leicester.ac.uk":paul.lambert@leicester.ac.uk})

{pstd}
Nikolaos Skourlis, Karolinska Institutet.


{title:Acknowledgement}

{phang}
This is based on the R command boxes written by Bendix Carstensen.
{p_end}
