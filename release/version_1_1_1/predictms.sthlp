{smcl}
{* *! version 1.0.0 ?????2012}{...}
{vieweralsosee "streg" "help streg"}{...}
{vieweralsosee "stpm2" "help stpm2"}{...}
{vieweralsosee "msset" "help msset"}{...}
{viewerjumpto "Syntax" "predictms##syntax"}{...}
{viewerjumpto "Description" "predictms##description"}{...}
{viewerjumpto "Options" "predictms##options"}{...}
{viewerjumpto "Examples" "predictms##examples"}{...}
{title:Title}

{p2colset 5 18 18 2}{...}
{p2col :{hi:predictms} {hline 2}}predictions from a multi-state survival model{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{phang2}
{cmd: predictms} {cmd:,} {opt transmatrix(varname)} [{it:options}]


{synoptset 29 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opth gen(stub)}}stub for generated new variables; default is {it:pred}{p_end}
{synopt:{opth transm:atrix(matname)}}transition matrix{p_end}
{synopt:{opth models(namelist)}}list of estimates stored for # transition{p_end}
{synopt:{opt reset}}use clock-reset approach{p_end}
{synopt:{opth from(numlist)}}starting state(s) for predictions{p_end}
{synopt:{opth obs(#)}}number of time points to calculate predictions at between {cmd:mint()} and {cmd:maxt()}{p_end}
{synopt:{opth mint(#)}}minimum time at which to calculate predictions{p_end}
{synopt:{opth maxt(#)}}maximum time at which to calculate predictions{p_end}
{synopt:{opth time:var(varname)}}time points at which to calculate predictions{p_end}
{synopt:{opth enter(#)}}time that observations enter model, default 0, for forward predictions{p_end}
{synopt:{opth exit(#)}}time that observations exit the model, for fixed horizon predictions{p_end}
{synopt:{opth n(#)}}sample size of simulated dataset{p_end}
{synopt:{opth m(#)}}number of simulation repetitions for calculating confidence intervals{p_end}
{synopt:{opt ci}}calculate confidence intervals of predictions{p_end}
{synopt:{opt normal}}calculate confidence intervals using a normal approximation{p_end}
{synopt:{opth seed(#)}}set the simulation seed{p_end}
{synopt:{opth l:evel(#)}}calculate confidence intervals at specific level, default is 95{p_end}
{synopt:{opt los}}calculate length of stay{p_end}
{synopt:{opt graph}}produce stacked transition probability graphs{p_end}
{synopt:{opth graphopts(options)}}pass options to twoway{p_end}
{synopt:{opth at(options)}}calculate predictions at covariate patterns{p_end}
{synopt:{opth at2(options)}}calculate contrast predictions at covariate patterns{p_end}
{synopt:{opt ratio}}calculate ratios of predictions{p_end}
{synopt:{opt surv:ival}}standard single outcome survival analysis{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:predictms} calculates a variety of predictions from a multi-state survival model, including 
transition probabilities and length of stay (restricted mean survival times, partitioned by cause). Predictions can be made 
at user-specified covariate patterns. Differences and ratios of predictions across two covariate patterns can also be calculated. 
Confidence intervals for all quantities can be obtained. Simulation is used to calculate all quantities.
{p_end}

{pstd}
Survival model fitting commands supported include {cmd:stpm2}, {cmd:strcs}, and all those in {cmd:streg}, except {cmd:distribution(gamma)}.
{p_end}

{pstd}
The user must provide the transition matrix used in the fitted nodel, through the {bf:transmatrix()} option.
Default predictions assume all subjects start in state {bf:_from = 1}, at time {bf:_start = 0}.
{p_end}

{pstd}
{cmd:predictms} creates the following variables:
{p_end}

	{bf:_time}            times at which predictions are calculated
	{bf:pred_#1_#2}       variables which contain calculated prediction (from state #1 to state #2)
	{bf:pred_#1_#2_lci}   variables which contain calculated lower confidence interval of prediction (from state #1 to state #2)
	{bf:pred_#1_#2_uci}   variables which contain calculated upper confidence interval of prediction (from state #1 to state #2)

{pstd}
{cmd:predictms} drops any variables called {cmd:pred_*} first.

{pstd}
The stub {it:pred} of the generated variables can be overidden by use of the {cmd:gen}() option.

{pstd}
Factor variables must not be used in fitted models. You must create your own dummy variables by specifying them in the {bf:covariates()} option of {bf:msset}, or by making sure they end in {bf:_trans#}.
	
{marker options}{...}
{title:Options}

{phang}
{opt gen(stub)} stub for the new variables created to hold predictions. Default is {it:pred}.

{phang}
{opt transmatrix(matname)} specifies the transition matrix used in the multi-state model 
that was fitted. This must be an upper triangular matrix (with diagonal and lower 
triangle elements coded missing). Transition must be numbered as an increasing sequence 
of integers from 1,...,K. If obtaining predictions after a model fitted to the stacked data created 
by {cmd:msset}, this this transition matrix must be the same as that used/produced by {cmd:msset}.

{phang}
{opt models(namelist)} specifies the names of the {cmd:estimates store} objects containing the estimates of the model fitted for transition 1, 2, 3, ..., for example,

{phang2}
{cmd: predictms} {cmd:,} {opt transmatrix(tmat)} {opt models(m1 m2 m3)}

{phang2}
When {bf:models()} is not specified, it is assumed a model has been fitted to the stacked data created by {bf:msset}.

{phang}
{opt reset} use the clock-reset approach, i.e. a semi-Markov model, in the simulation framework when calculating 
predictions, default is clock-forward using delayed entry

{phang}
{opt from(numlist)} define the starting state for all observations, default is state 1. 
Multiple starting states can be defined, which will calculate all possible predictions for each starting state.

{phang}
{opt obs(#)} the number of time points at which to calculate predictions, equally spaced between {cmd:mint()} and {cmd:maxt()}, with a default of 100

{phang}
{opt mint(#)} minimum time at which to calculate predictions. If {cmd:timevar()} is not specified, then a default time variable will 
be created, using {cmd:mint()}, {cmd:maxt()} and {cmd:obs()}.

{phang}
{opt maxt(#)} maximum time at which to calculate predictions. If {cmd:timevar()} is not specified, then a default time variable will 
be created, using {cmd:mint()}, {cmd:maxt()} and {cmd:obs()}.

{phang}
{opt timevar(varname)} variable which contains time points at which to calculate predictions, which overrides the default.

{phang}
{opt enter(#)} defines the time that observations enter the model, default 0, for forward predictions.

{phang}
{opt exit(#)} defines the time that observations exit the model, for fixed horizon predictions.

{phang}
{opt n(#)} samples size of each simulated dataset, default is 100,000 unless {cmd:ci} is specified, then it is 10,000. 
Accuracy increases the higher the sample.

{phang}
{opt m(#)} number of simulation repetitions for calculating confidence intervals.

{phang}
{opt ci} calculate confidence intervals of predictions. Default method uses normal approximation of simulated repetitions to calculate the standard error.

{phang}
{opt normal} calculate confidence intervals based on a normal approximation, instead of the default perecentile based.

{phang}
{opt seed(#)} sets the simulation seed.

{phang}
{opt level(#)} confidence interval level, default {cmd:level(95)}.

{phang}
{opt los} calculate length of stay in each state, i.e. restricted mean survival time partitioned by states.

{phang}
{opt graph} produce stacked transition probability plots.

{phang}
{opt graphopts(twoway_options)} pass options to the {cmd:twoway} command when using the {cmd:graph} option.

{phang}
{opt at(vn # ...)} calculate predictions at specified covariate patterns, e.g. {cmd:at(female 1 age 45)}. Only for use with the {cmd:model#()} syntax.

{phang}
{opt at2(vn # ...)} calculate a contrast prediction at specified covariate patterns, e.g. {cmd:at(female 1 age 55)}. For use with {cmd:at()}. Default prediction is 
differences in transition probabilities. See {cmd:ratio} and {cmd:los}.

{phang}
{opt ratio} calculate the ratio of the predictions, specified with {cmd:at()} and {cmd:at2()}, instead of the default difference.

{phang}
{opt survival} indicates that you have fitted a standard single event survival analysis model. This corresponds to a transition matrix of
 (.,1\.,.). {bf:transmatrix()} does not need to be specified. All predictions from {bf:predictms} are valid.


{marker examples}{...}
{title:Example 1:}

{pstd}
This dataset contains information on 2982 patients with breast cancer. Baseline is defined as time of surgery, and patients can experience 
relapse, relapse then death, or death with no relapse. Time of relapse is stored in {cmd:rf}, with event indicator {cmd:rfi}, and time of death 
is stored in {cmd:os}, with event indicator {cmd:osi}.
{p_end}

{pstd}Load example dataset:{p_end}
{phang}{stata "use http://fmwww.bc.edu/repec/bocode/m/multistate_example":. use http://fmwww.bc.edu/repec/bocode/m/multistate_example}{p_end}

{pstd}{cmd:msset} the data:{p_end}
{phang}{stata "msset, id(pid) states(rfi osi) times(rf os)":. msset, id(pid) states(rfi osi) times(rf os)}{p_end}

{pstd}Store the transition matrix:{p_end}
{phang}{stata "mat tmat = r(transmatrix)":. mat tmat = r(transmatrix)}{p_end}

{pstd}stset the data using the variables created by {cmd:msset}{p_end}
{phang}{stata "stset _stop, enter(_start) failure(_status=1)":. stset _stop, enter(_start) failure(_status=1)}{p_end}

{pstd}Fit a Weibull model, allowing a separate baseline (stratified), but the same effect of age across transitions, assuming transition 1 as the reference:{p_end}
{phang}{stata "streg age _trans2 _trans3, dist(weibull) ancillary(_trans2 _trans3)":. streg age _trans2 _trans3, dist(weibull) ancillary(_trans2 _trans3)}{p_end}

{pstd}Calculate transition probabilities for a patient with age 50:{p_end}
{phang}{stata "predictms, transmatrix(tmat) at(age 50)":. predictms, transmatrix(tmat) at(age 50)}{p_end}

{pstd}Alternatively, we fit separate Weibull models, also allowing transition specific age effects:{p_end}

{phang}{stata "streg age if _trans1==1, dist(weibull)":. streg age if _trans1==1, dist(weibull)}{p_end}
{phang}{stata "estimates store m1":. estimate store m1}{p_end}

{phang}{stata "streg age if _trans2==1, dist(weibull)":. streg age if _trans2==1, dist(weibull)}{p_end}
{phang}{stata "estimates store m2":. estimate store m2}{p_end}

{phang}{stata "streg age if _trans3==1, dist(weibull)":. streg age if _trans3==1, dist(weibull)}{p_end}
{phang}{stata "estimates store m3":. estimate store m3}{p_end}

{pstd}Calculate transition probabilities for a patient with age 50:{p_end}
{phang}{stata "predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50)":. predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50)}{p_end}

{pstd}Calculate length of stay for a patient with age 50:{p_end}
{phang}{stata "predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50) los":. predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50) los}{p_end}

{pstd}Calculate the difference in transition probabilities for a patient with age 50 compared to a patient aged 60:{p_end}
{phang}{stata "predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50) at2(age 60)":. predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50) at2(age 60)}{p_end}

{pstd}Calculate the ratio of transition probabilities for a patient with age 50 compared to a patient aged 60:{p_end}
{phang}{stata "predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50) at2(age 60) ratio":. predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50) at2(age 60) ratio}{p_end}

{pstd}Calculate the ratio of length of stay for a patient with age 50 compared to a patient aged 60, with confidence intervals:{p_end}
{phang}{stata "predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50) at2(age 60) los ratio ci":. predictms, transmatrix(tmat) models(m1 m2 m3) at(age 50) at2(age 60) los ratio ci}{p_end}


{title:Author}

{pstd}Michael J. Crowther{p_end}
{pstd}Department of Health Sciences{p_end}
{pstd}University of Leicester{p_end}
{pstd}E-mail: {browse "mailto:michael.crowther@le.ac.uk":michael.crowther@le.ac.uk}{p_end}

{phang}
Please report any errors you may find.{p_end}


{title:References}

{phang}
Crowther MJ and Lambert PC. Parametric multi-state survival models: flexible modelling allowing transition-specific distributions with application to estimating clinically useful measures of effect differences (Submitted).
{p_end}

{phang}
de Wreede LC, Fiocco M and Putter H. mstate: An R Package for the Analysis of Competing Risks and Multi-State Models. {it:Journal of Statistical Software} 2011;38:1-30.
{p_end}

{phang}
Putter H, Fiocco M and Geskus RB. Tutorial in biostatistics: competing risks and multi-state models. {it:Statistics in Medicine} 2007;26:2389-2430.
{p_end}

