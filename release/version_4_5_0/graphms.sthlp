{smcl}
{vieweralsosee "[multistate] multistate" "help multistate"}{...}
{vieweralsosee "[multistate] msset" "help msset"}{...}
{vieweralsosee "[multistate] msboxes" "help msboxes"}{...}
{vieweralsosee "[multistate] msaj" "help msaj"}{...}
{vieweralsosee "[multistate] predictms" "help predictms"}{...}
{viewerjumpto "Syntax" "graphms##syntax"}{...}
{viewerjumpto "Description" "graphms##description"}{...}
{viewerjumpto "Options" "graphms##options"}{...}
{viewerjumpto "Examples" "graphms##examples"}{...}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{hi:graphms} {hline 2}}create a stacked transition probability graph{p_end}
{p2colreset}{...}

{phang}
{cmd:graphms} is part of the {helpb multistate} package, for use after calculating transition probabilities with {helpb predictms}
{p_end}

{marker syntax}{...}
{title:Syntax}

{phang2}
{cmd: graphms} {cmd:,} [{it:options}]

{synoptset 29 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opth from(numlist)}}starting state(s) for predictions{p_end}
{synopt:{opth at(numlist)}}calculate predictions for the defined {cmd:at#()} patterns; see details{p_end}
{synopt:{opth time:var(varname)}}time variable used for the x-axis{p_end}
{synopt:{opth Nstates(#)}}defines the number of potential states{p_end}
{synoptline}

{phang}
All options are optional...if they are not specified, {cmd:graphms} will obtain the appropriate information from the 
{cmd:return list} of {cmd:predictms} (if it is still in memory).
{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:graphms} creates a stacked transition probability plot, following the estimation of a multi-state model 
(using {helpb stmerlin} or {helpb merlin}) and prediction of such quantities using {helpb predictms}. {cmd:graphms} uses 
the {cmd:_prob_at} variables creates by {cmd:predictms} to create the plot.
{p_end}

{phang}
{cmd:graphms} is part of the {helpb multistate} package. 
{p_end}
	

{marker options}{...}
{title:Options}

{phang}
{opt from(numlist)} define the starting state for all predictions. Multiple starting states can be defined, 
which will create multiple graphs.

{phang}
{opt at(numlist)} defines the {cmd:at#()} patterns to create plots for. If {cmd:at(1)} is specified, then plots will 
be creates for the {cmd:_prob_at1_*} predictions. Multiple {cmd:at()}s can be specified.

{phang}
{opt timevar(varname)} variable to use as the x-axis, i.e. the time points at which the {cmd:_prob*} predictions 
were calculated at.

{phang}
{opt Nstates(#)} defines the number of potential states.


{marker examples}{...}
{title:Example 1:}

{pstd}
This dataset contains information on 2982 patients with breast cancer. Baseline is defined as time of surgery, and patients can experience 
relapse, relapse then death, or death with no relapse. Time of relapse is stored in {cmd:rf}, with event indicator {cmd:rfi}, and time of death 
is stored in {cmd:os}, with event indicator {cmd:osi}.
{p_end}

{pstd}Load example dataset:{p_end}
{phang}{stata "use http://fmwww.bc.edu/repec/bocode/m/multistate_example":. use http://fmwww.bc.edu/repec/bocode/m/multistate_example}{p_end}

{pstd}{helpb msset} the data (from the {cmd:multistate} package):{p_end}
{phang}{stata "msset, id(pid) states(rfi osi) times(rf os)":. msset, id(pid) states(rfi osi) times(rf os)}{p_end}

{pstd}Store the transition matrix:{p_end}
{phang}{stata "mat tmat = r(transmatrix)":. mat tmat = r(transmatrix)}{p_end}

{pstd}stset the data using the variables created by {cmd:msset}{p_end}
{phang}{stata "stset _stop, enter(_start) failure(_status=1)":. stset _stop, enter(_start) failure(_status=1)}{p_end}

{pstd}We fit separate Weibull models, so a fully stratified model, also allowing transition specific age effects:{p_end}

{phang}{stata "stmerlin age if _trans1==1, distribution(weibull)":. stmerlin age if _trans1==1, distribution(weibull)}{p_end}
{phang}{stata "estimates store m1":. estimate store m1}{p_end}

{phang}{stata "stmerlin age if _trans2==1, distribution(weibull)":. stmerlin age if _trans2==1, distribution(weibull)}{p_end}
{phang}{stata "estimates store m2":. estimate store m2}{p_end}

{phang}{stata "stmerlin age if _trans3==1, distribution(weibull)":. stmerlin age if _trans3==1, distribution(weibull)}{p_end}
{phang}{stata "estimates store m3":. estimate store m3}{p_end}

{pstd}Calculate transition probabilities for a patient with age 50:{p_end}
{phang}{stata "predictms, transmatrix(tmat) models(m1 m2 m3) at1(age 50)":. predictms, transmatrix(tmat) models(m1 m2 m3) at1(age 50)}{p_end}

{pstd}Create a stacked plot of transition probabilities:{p_end}
{phang}{stata "graphms":. graphms}{p_end}


{title:Author}

{pstd}Michael J. Crowther{p_end}
{pstd}Red Door Analytics{p_end}
{pstd}Stockholm, Sweden{p_end}
{pstd}E-mail: {browse "mailto:michael@reddooranalytics.se":michael@reddooranalytics.se}{p_end}

{phang}
Please report any errors you may find.{p_end}


{title:References}

{phang}
Crowther MJ. Extended multivariate generalised linear and non-linear mixed effects models. 
{browse "https://arxiv.org/abs/1710.02223":https://arxiv.org/abs/1710.02223}
{p_end}

{phang}
Crowther MJ. merlin - a unified framework for data analysis and methods development in Stata. 
{browse "https://arxiv.org/abs/1806.01615":https://arxiv.org/abs/1806.01615}
{p_end}

{phang}
Crowther MJ, Lambert PC. Parametric multi-state survival models: flexible modelling allowing transition-specific distributions with 
application to estimating clinically useful measures of effect differences. {it: Statistics in Medicine} 2017;36(29):4719-4742.
{p_end}

{phang}
de Wreede LC, Fiocco M, Putter H. mstate: An R Package for the Analysis of Competing Risks and Multi-State Models. 
{it:Journal of Statistical Software} 2011;38:1-30.
{p_end}

{phang}
Putter H, Fiocco M, Geskus RB. Tutorial in biostatistics: competing risks and multi-state models. 
{it:Statistics in Medicine} 2007;26:2389-2430.
{p_end}

{phang}
Weibull CE, Lambert PC, Eloranta S, Andersson TML, Dickman PW, Crowther MJ. A multi-state model incorporating 
estimation of excess hazards and multiple time scales. {it: Submitted.}
{p_end}


