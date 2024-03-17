{smcl}
{* *! version 1.0.0}{...}
{vieweralsosee "[multistate] multistate" "help multistate"}{...}
{vieweralsosee "[multistate] msset" "help msset"}{...}
{vieweralsosee "[multistate] msboxes" "help msboxes"}{...}
{vieweralsosee "[multistate] msaj" "help msaj"}{...}
{vieweralsosee "[multistate] graphms" "help graphms"}{...}
{vieweralsosee "[merlin] stmerlin" "help stmerlin"}{...}
{vieweralsosee "[merlin] merlin" "help merlin"}{...}
{viewerjumpto "Syntax" "predictms##syntax"}{...}
{viewerjumpto "Description" "predictms##description"}{...}
{viewerjumpto "Options" "predictms##options"}{...}
{viewerjumpto "Examples" "predictms##examples"}{...}
{title:Title}

{p2colset 5 18 18 2}{...}
{p2col :{hi:predictms} {hline 2}}predictions from a multi-state model{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{phang2}
{cmd: predictms} {cmd:,} {help predictms##transopt:{it:transmatrix}} 
{help predictms##statistic:{it:statistic}} 
[{help predictms##model_spec:{it:model_spec}} 
{help predictms##pred_spec:{it:prediction_spec}} 
{help predictms##options:{it:options}}]


{marker transopt}{...}
{synoptset 29 tabbed}{...}
{synopthdr:transition matrix}
{synoptline}
{synopt:{opth transm:atrix(matname)}}transition matrix{p_end}
{synopt:{opt singleevent}}shorthand to define a {cmd:transmatrix()} for a single event survival model{p_end}
{synopt:{opt cr}}shorthand to define a {cmd:transmatrix()} for a competing risks model{p_end}
{synoptline}
{pstd}One of the above options must be specified

{marker statistic}{...}
{synoptset 29 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{synopt:{opt prob:ability}}calculate transition probabilities{p_end}
{synopt:{opt los}}calculate length of stay in each state{p_end}
{synopt:{opt rmst}}calculate restricted mean survival time{p_end}
{synopt:{opt visit}}probability of ever visiting each state{p_end}
{synopt:{opt haz:ard}}transition-specific hazard function{p_end}
{synopt:{opt surv:ival}}transition-specific survival function{p_end}
{synopt:{cmd: {ul:userf}unction(}{it:func_name}{cmd:)}}user-defined Mata function for bespoke predictions; see details{p_end}
{synoptline}
{pstd}At least one of the above options must be specified

{marker model_spec}{...}
{synoptset 29 tabbed}{...}
{synopthdr:model specification}
{synoptline}
{syntab:{it:}}
{synopt:{opth models(namelist)}}list of estimates stored for all transition hazards{p_end}
{synopt:{opt reset}}use clock-reset approach{p_end}
{synopt:{opt tscale2(numlist)}}transition models on a second timescale{p_end}
{synopt:{opt time2(numlist)}}time to add to main timescale{p_end}
{synoptline}

{marker pred_spec}{...}
{synoptset 29 tabbed}{...}
{synopthdr:prediction specification}
{synoptline}
{synopt:{opth from(numlist)}}starting state(s) for predictions{p_end}
{synopt:{opth lt:runcated(#)}}starting time, i.e. time at which the starting state(s) {cmd:from()} are entered{p_end}
{synopt:{opth time:var(varname)}}time points at which to calculate predictions{p_end}
{synopt:{opth mint(#)}}minimum time at which to calculate predictions{p_end}
{synopt:{opth maxt(#)}}maximum time at which to calculate predictions{p_end}
{synopt:{opth obs(#)}}number of time points to calculate predictions at between {cmd:mint()} and {cmd:maxt()}{p_end}
{synopt:{opt at#(at_list)}}calculate predictions (and contrasts) at covariate patterns{p_end}
{synopt:{opt stand:ardise}}calculates standardised/population-averaged predictions; see details{p_end}
{synopt:{opt standif(condition)}}restricts the observations that are standardised over; see details{p_end}
{synopt:{opt diff:erence}}calculate differences between predictions{p_end}
{synopt:{opt ratio}}calculate ratios of predictions{p_end}
{synopt:{opth atref:erence(#)}}specifies the reference prediction for {cmd:difference} and {cmd:ratio} contrasts{p_end}
{synopt:{opt userl:ink(string)}}link function used in calculation of confidence intervals for {cmd:userfunction()}; default {cmd:identity}{p_end}
{synopt:{opt out:sample}}for out of sample predictions; see below.{p_end}

{synopt:{it:confidence intervals}}{p_end}
{synopt:{opt ci}}calculate confidence intervals of predictions{p_end}
{synopt:{opth l:evel(#)}}calculate confidence intervals at specific level, default is 95{p_end}
{synopt:{opt boot:strap}}calculate confidence intervals using the parametric bootstrap{p_end}
{synopt:{opth m(#)}}number of bootstrap repetitions for calculating confidence intervals{p_end}
{synopt:{opt perc:entile}}calculate confidence intervals using percentiles of the bootstrap repetitions{p_end}
{synopt:{opt novcv(numlist)}}transition models that should be assumed are not estimated with uncertainty; see details{p_end}
{synoptline}

{marker model_spec}{...}
{synoptset 29 tabbed}{...}
{synopthdr:options}
{synoptline}
{synopt:{opt aj}}use the Aalen-Johansen estimator for transition probabilities; see details{p_end}
{synopt:{opt sim:ulate}}calculate predictions using large-sample simulation; see details{p_end}
{synopt:{opt latent}}use latent times as the method of simulation; see details{p_end}
{synopt:{opth n(#)}}sample size of simulated dataset{p_end}
{synopt:{opth chintpoints(#)}}number of Gauss-Legendre quadrature points; see details{p_end}
{synopt:{opt save(name, [replace])}}save each simulated dataset used to calculate predictions and confidence intervals; see details{p_end}
{synopt:{opth seed(#)}}set the simulation seed{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:predictms} calculates a variety of predictions from a multi-state survival model, including 
transition probabilities, length of stay (restricted mean survival times in each state), the probability of ever visiting 
each state, and more. Predictions are made at user-specified covariate patterns. Differences and ratios of predictions across 
covariate patterns can also be calculated. Standardised (study population-averaged) predictions can be obtained. Confidence 
intervals for all quantities are available. Numerical integration, simulation or the Aalen-Johansen estimator (with 
parametric estimates of the cumulative hazard) are used to calculate all quantities. User-defined predictions can also be 
calculated by providing a user-written {cmd:Mata} function, to provide complete flexibility.
{p_end}

{pstd}
The transition models must be estimated using the {helpb stmerlin} or {helpb merlin} commands, with the current range of 
supported models including:
{p_end}

{phang2}- exponential proportional hazards model{p_end}
{phang2}- Weibull proportional hazards model{p_end}
{phang2}- Gompertz proportional hazards model{p_end}
{phang2}- generalised gamma accelerated failure time model{p_end}
{phang2}- log normal accelerated failure time model{p_end}
{phang2}- log logistic accelerated failure time model{p_end}
{phang2}- piecewise-exponential proportional hazards model{p_end}
{phang2}- Royston-Parmar flexible parametric model{p_end}
{phang2}- general log hazard scale model{p_end}
{phang2}- general log cumulative hazard scale model{p_end}
{phang2}- general additive hazard model{p_end}
{phang2}- user-defined/custom hazard model{p_end}

{pstd}
The user (usually) must provide the transition matrix used in the fitted model, through the {bf:transmatrix()} option. 
The transition matrix may be cyclic or acyclic. Default predictions assume all subjects start in state {bf:from(1)}, at time 
{cmd:ltruncated(0)}.
{p_end}

{pstd}
{cmd:predictms} creates the following variables:
{p_end}

{phang2}
	{bf:_time} times at which all predictions are calculated
{p_end}

{pstd}	
If {cmd:probability} is requested, then {cmd:predictms} creates the following variables:
{p_end}	
{phang2}
{bf:_prob_at{it:i}_{it:a}_{it:b}}       transition probability for {cmd:at{it:i}()} (from state {it:a} to state {it:b})
{p_end}

{phang2}
	{bf:_prob_at{it:i}_{it:a}_{it:b}_lci}   lower confidence interval of transition probability for {cmd:at{it:i}()} (from state {it:a} to state {it:b})
{p_end}
	
{phang2}
	{bf:_prob_at#_{it:a}_{it:b}_uci}   upper confidence interval of transition probability for {cmd:at{it:i}()} (from state {it:a} to state {it:b})
{p_end}

{pstd}	
If {cmd:los} is requested, then {cmd:predictms} also creates the following variables:
{p_end}

{phang2}
{bf:_los_at{it:i}_{it:a}_{it:b}}       length of stay in state {it:b} for {cmd:at{it:i}()} (given they started from state {it:a})

{phang2}
	{bf:_los_at{it:i}_{it:a}_{it:b}_lci}   lower confidence interval of the length of stay in state {it:b} for {cmd:at{it:i}()} (given they started from state {it:a})
	
{phang2}
	{bf:_los_at{it:i}_{it:a}_{it:b}_uci}   upper confidence interval of the length of stay in state {it:b} for {cmd:at{it:i}()} (given they started from state {it:a})
	
{pstd}	
If {cmd:rmst} is requested, then {cmd:predictms} also creates the following variables:
{p_end}

{phang2}
{bf:_rmst_at{it:i}_{it:a}}       total time spent in non-absorbing states for {cmd:at{it:i}()} (given they started from state {it:a})

{phang2}
	{bf:_rmst_at{it:i}_{it:a}_lci}   lower confidence interval of the total time spent in non-absorbing states for {cmd:at{it:i}()} (given they started from state {it:a})
	
{phang2}
	{bf:_rmst_at{it:i}_{it:a}_uci}   upper confidence interval of the total time spent in non-absorbing states for {cmd:at{it:i}()} (given they started from state {it:a})
	
{pstd}	
If {cmd:visit} is requested, then {cmd:predictms} also creates the following variables:
{p_end}

{phang2}
{bf:_visit_at{it:i}_{it:a}_{it:b}}       probability of ever visiting state {it:b} for {cmd:at{it:i}()} (given they started from state {it:a})

{phang2}
	{bf:_visit_at{it:i}_{it:a}_{it:b}_lci}   lower confidence interval of the probability of ever visiting state {it:b} for {cmd:at{it:i}()} (given they started from state {it:a})
	
{phang2}
	{bf:_visit_at{it:i}_{it:a}_{it:b}_uci}   upper confidence interval of the probability of ever visiting state {it:b} for {cmd:at{it:i}()} (given they started from state {it:a})

{pstd}	
If {cmd:hazard} is requested, then {cmd:predictms} also creates the following variables:
{p_end}

{phang2}
{bf:_hazard_at{it:i}_{it:a}_{it:b}}       predicted hazard function for the transition going from state {it:a} to {it:b} for {cmd:at{it:i}()}

{phang2}
	{bf:_hazard_at{it:i}_{it:a}_{it:b}_lci}   lower confidence interval of the predicted hazard function for the transition going from state {it:a} to {it:b} for {cmd:at{it:i}()}
	
{phang2}
	{bf:_hazard_at{it:i}_{it:a}_{it:b}_uci}   upper confidence interval of the predicted hazard function for the transition going from state {it:a} to {it:b} for {cmd:at{it:i}()}

{pstd}	
If a {cmd:userfunction()} is provided, then {cmd:predictms} also creates the following variables:
{p_end}

{phang2}
	{bf:_user_at{it:i}_{it:a}_{it:b}}       returned values for column {it:b} of the user-function for {cmd:at{it:i}()} (given they started from state {it:a})

{phang2}
	{bf:_user_at{it:i}_{it:a}_{it:b}_lci}   lower confidence interval of the returned values for column {it:b} of the user-function for {cmd:at{it:i}()} (given they started from state {it:a})
	
{phang2}
	{bf:_user_at{it:i}_{it:a}_{it:b}_uci}   upper confidence interval of the returned values for column {it:b} of the user-function for {cmd:at{it:i}()} (given they started from state {it:a})
	
{pstd}	
If a {cmd:difference} is requested, then {cmd:predictms} also creates some or all of the following variables:
{p_end}

{phang2}
{bf:_diff_prob_at{it:i}_{it:a}_{it:b}} difference in transition probabilities between {cmd:at{it:i}()} and that specified in the reference {cmd:at{it:j}()}

{phang2}
{bf:_diff_los_at{it:i}_{it:a}_{it:b}} difference in length of stay between {cmd:at{it:i}()} and that specified in the reference {cmd:at{it:j}()}

{phang2}
{bf:_diff_user_at{it:i}_{it:a}_{it:b}} difference in the user-function between {cmd:at{it:i}()} and that specified in the reference {cmd:at{it:j}()}

{phang2}
{bf:*_lci} and {bf:*_uci} are returned as appropriate.

{pstd}	
If a {cmd:ratio} is requested, then {cmd:predictms} also creates some or all of the following variables:
{p_end}

{phang2}
{bf:_ratio_prob_at{it:i}_{it:a}_{it:b}} ratio in transition probabilities between {cmd:at{it:i}()} and that specified in the reference {cmd:at{it:j}()}

{phang2}
{bf:_ratio_los_at{it:i}_{it:a}_{it:b}} ratio in length of stay between {cmd:at{it:i}()} and that specified in the reference {cmd:at{it:j}()}

{phang2}
{bf:_ratio_user_at{it:i}_{it:a}_{it:b}}	ratio in the user-function between {cmd:at{it:i}()} and that specified in the reference {cmd:at{it:j}()}

{phang2}
{bf:*_lci} and {bf:*_uci} are returned as appropriate.
	
{pstd}
All appropriate combinations of {it:a} and {it:b} are returned. {cmd:predictms} first drops any variables it may have calculated 
from a previous call (e.g. {cmd:_prob_*}).

{phang}
{cmd:predictms} is part of the {helpb multistate} package. Further 
details here: {bf:{browse "https://www.mjcrowther.co.uk/software/multistate":mjcrowther.co.uk/software/multistate}}
{p_end}
	

{marker options}{...}
{title:Options}

{dlgtab:Transition matrix}

{phang}
{opt transmatrix(matname)} specifies the transition matrix which governs the multi-state model 
that was fitted. Transitions must be numbered as an increasing sequence 
of integers from 1,...,K, from left to right, top to bottom of the matrix. Reversible transitions are allowed. 

{phang}
{opt singleevent} indicates that you have fitted a standard single event survival analysis model. This corresponds to a 
transition matrix of (.,1\.,.). {bf:transmatrix()} does not need to be specified.

{phang}
{opt cr} indicates that you have fitted a competing risks model, and is a useful way of avoiding having to specify a 
{bf:transmatrix()}. For use with {cmd:models()}, as the number of competing risks corresponds to the number of model 
objects. For example, for a model with two competing events:

{phang2}
{cmd:. matrix define tmat = (.,1,2\.,.,.\.,.,.)}

{dlgtab:Statistic}

{phang}
{opt probability} calculate transition probabilities. If {cmd:aj} has been used then I recommend at least {cmd:obs(500)} 
is used in your {cmd:timevar()}.

{phang}
{opt los} calculate (restricted) length of stay in each state, i.e. restricted mean survival time for each transient state. 
This is the integral of the transition probabilities across follow-up time. If {cmd:aj} has been used then length of stay 
is calculated using numerical integration of the transition probabilities, and I recommend at least {cmd:obs(500)} 
is used in your {cmd:timevar()}.

{phang}
{opt rmst} calculate (restricted) mean survival time, i.e. the total time spent in non-absorbing/transient states. 
If {cmd:aj} has been used then restricted mean survival time is calculated using numerical integration of the 
transition probabilities, and I recommend at least {cmd:obs(500)} is used in your {cmd:timevar()}.

{phang}
{opt visit} calculate the probability of ever visiting each state within the time interval defined from {cmd:mint()} and 
the prediction time. Can currently only be calculated using the large-sample simulation method.

{phang}
{opt hazard} calculate the predicted hazard function for each transition from the {cmd:from()} state(s).

{phang}
{opt survival} calculate the predicted survival function for each transition from the {cmd:from()} state(s).

{phang}
{opt userfunction(func_name)} defines a Mata function which returns a {it: real matrix}, to calculate user-defined quantities. 
The function must be of the form:

{p 8 12 2}
{cmd:{it: real matrix} func_name(S)}
{p_end}
{p 8 12 2}
{bf:{c -(}}
{p_end}
{p 12 12 2}
{cmd:...}
{p_end}
{p 12 12 2}
{cmd:{it:some code}}
{p_end}
{p 12 12 2}
{cmd:...}
{p_end}
{p 12 12 2}
{cmd:pred = {it:some more code}}
{p_end}
{p 12 12 2}
{cmd:return(pred)}
{p_end}
{p 8 12 2}
{bf:{c )-}}
{p_end}

{phang2}
where S is a transmorphic object which should not be changed. It is passed to utility functions, which give access to the 
transition probabilites and/or length of stays, for example:

{p 8 12 2}
{cmd:{it: real matrix} func_name(S)}
{p_end}
{p 8 12 2}
{bf:{c -(}}
{p_end}
{p 12 12 2}
{cmd:p1 = ms_user_prob(S,1)}
{p_end}
{p 12 12 2}
{cmd:p2 = ms_user_prob(S,2)}
{p_end}
{p 12 12 2}
{cmd:p3 = ms_user_prob(S,3)}
{p_end}
{p 12 12 2}
{cmd:pred = p1,p2,p3}
{p_end}
{p 12 12 2}
{cmd:return(pred)}
{p_end}
{p 8 12 2}
{bf:{c )-}}
{p_end}

{phang2}
which would give us our transition probabilities for each state. Length of stays are accessed using {cmd:ms_user_los()} in the 
same way, but the {cmd:los} option must also be specified.

{dlgtab:Model specification}

{phang}
{opt models(namelist)} specifies the names of the {helpb estimates store} objects containing the estimates of the model 
fitted for transition 1, 2, 3, ..., for example,

{phang2}
{cmd:. matrix tmat = (.,1,2\.,.,3\.,.,.)}
{p_end}
{phang2}
{cmd:. predictms} {cmd:, transmatrix(tmat)} {opt models(m1 m2 m3)}
{p_end}

{phang2}If {cmd:models()} is not specified, then {cmd:predictms} assumes you have fitted a single, stacked, multi-state model to 
your dataset.

{phang}
{opt reset} tells {cmd:predictms} that the clock-reset approach, i.e. a Markov renewal multi-state model has been fitted. The 
default assumes clock-forward. When {cmd:reset} is specified, then the time at which predictions are made, stored in 
{cmd:_time}, represent the main time scale, time since entry to the {cmd:from()} state(s), which may be left truncated if 
{cmd:ltruncated()} is greater than 0.

{phang}
{opt tscale2(numlist)} specifies any transition models (using the same index as specified in {cmd:transmatrix()}), which 
are modelled on a secondary timescale (different to the main one), enabling transition-specific timescales (Weibull et al. 
{it:Under review.}). This can be used for example when modelling age as the timescale for some transitions, and time since 
diagnosis for others. 

{phang}
{opt time2(numlist)} specifies the value to be added to the main timescale at entry, for each of the transition models 
specified in {cmd:tscale2()}. If age is the second timescale, then this would be age at baseline. Each element of 
{cmd:time2()} corresponds to each {cmd:at#()}. If only one {cmd:time2()} is specified then it will be assumed for all 
{cmd:at#()}s.

{dlgtab:Prediction specification}

{phang}
{opt from(numlist)} define the starting state for all predictions (state {it:a} in the variable descriptions above), default 
is state 1. Multiple starting states can be defined, which will calculate all possible predictions for each starting state.

{phang}
{opt ltruncated(#)} defines the starting time, i.e. the time at which the starting state(s) {cmd:from()} are entered, with a default 
of {cmd:ltruncated(0)}

{phang}
{opt timevar(varname)} variable which contains time points at which to calculate predictions, which overrides the default 
created by using {cmd:mint()}, {cmd:maxt()} and {cmd:obs()}. Predictions are made at these timepoints, conditional on 
starting in state(s) {cmd:from()} at time {cmd:ltruncated()}.

{phang}
{opt mint(#)} minimum time at which to calculate predictions, with a default of {cmd:mint(0)}. If {cmd:timevar()} is not 
specified, then a default time variable will be created, using {cmd:mint()}, {cmd:maxt()} and {cmd:obs()}.

{phang}
{opt maxt(#)} maximum time at which to calculate predictions. If {cmd:timevar()} is not specified, then a default time 
variable will be created, using {cmd:mint()}, {cmd:maxt()} and {cmd:obs()}.

{phang}
{opt obs(#)} the number of time points at which to calculate predictions, equally spaced between {cmd:mint()} and {cmd:maxt()}, 
with a default of 100 (500 when using {cmd:aj}).

{phang}
{opt at#(vn # ...)} calculates predictions at specified covariate patterns, e.g. {cmd:at1(female 1 age 55)}. 
Specifying multiple {cmd:at#()}s means only one call of {cmd:predictms} has to be made to calculate many predictions. 
Note that any covariates not specified in {cmd:at#()} are set to {cmd:0}.

{phang}
{opt standardise} calculates standardised/population-averaged predictions. For each {it:statistic} requested, the prediction 
will be averaged over all covariate patterns found either in your dataset, or in the observations identified using the
{cmd:standif()} condition. For each observation, any covariates not specified in {cmd:at#()} will take their observed values - 
this is done for all observations, and the prediction is then averaged. Predictions can then be interpreted as marginal 
predictions with respect to the observed covariate distributions that have been averaged over. A standardised prediction 
can be very computationally intenstive to obtain.

{phang}
{opt standif(condition)} restricts the observations that are standardised over, for example {cmd:standif(_n<100)}, or 
{cmd:standif(trt==1)}

{phang}
{opt difference} calculate the difference between predictions, across different {cmd:at#()} specifications. 
See {cmd:atreference()}.

{phang}
{opt ratio} calculate the ratio of the predictions, across different {cmd:at#()} specifications.

{phang}
{opt atreference(#)} specifies the reference {cmd:at#()} for calculating prediction contrasts. Default is {cmd:atref(1)}, meaning 
that each prediction is contrasted to that calculated using {cmd:at1()}.

{phang}
{opt userlink(link_name)} specifies the link function used when calculating confidence intervals with the normal approximation 
applied to the {cmd:userfunction()}. Options include {cmd:log}, {cmd:logit}, and the default {cmd:identity}. The link function 
is applied, and the mean and standard error are then calculated on the transformed scale to calculate confidence intervals, 
and then transformed back. This ensures all predictions are within the desired range (e.g. 0 and 1 for probabilities).

{phang}
{opt outsample} specifies predictions are being made out of sample, which suppresses checks that variables specified in 
{cmd:at#()}s are in your current dataset. This is of use when transition {cmd:models()} come from different datasets. As 
this suppresses some error checks, users should be extra careful when specifying their {cmd:at#()}s.

{phang}
{opt ci} requests confidence intervals of all predictions. The delta method or parametric bootstrap are used, depending on the model 
specification and the requested statistics. The most efficient method will be chosen by default. The number of bootstrap draws from 
the estimated coefficient vectors, and associated variance-covariance matrix, is specified in the {cmd:m()} option. When bootstrapping, 
the default method uses a normal approximation of simulated repetitions to calculate the standard error, see also {cmd:percentile}. 

{phang}
{opt level(#)} confidence interval level, default {cmd:level(95)}.

{phang}
{opt bootstrap} forces {cmd:predictms} to use the parametric bootstrap to calculate confidence intervals, rather than the 
delta method, which is the default for survival, competing risks & illness-death models.

{phang}
{opt m(#)} number of parametric bootstrap repetitions for calculating confidence intervals; default is {cmd:m(200)}. I 
recommend using an increasing number to check a sufficiently large {cmd:m()} has been used.

{phang}
{opt percentile} calculate confidence intervals based on centiles of the predictions across the {cmd:m(#)} sets, instead of 
the default normal based calculation.

{phang}
{opt novcv(numlist)} specifies transitions (using the same index as specified in {cmd:transmatrix()}) that when calculating 
confidence intervals on predictions (using {cmd:ci}), are assumed to be estimated free of uncertainty, i.e. draws will not be 
made from the multivariate normal centred on the parameter estimates and with VCV the estimated VCV from that transition, but 
will use the estimated parameter vector every time. 

{dlgtab:Options}

{phang}
{opt aj} specifies a hybrid Aalen-Johansen estimator be used to calculate predictions, instead of either numerical integration or 
large-sample simulation. This method uses the non-parametric formula for the AJ estimator, but uses parametric estimates of the 
transition-specific cumulatve hazard functions from the fitted {cmd:models()}, or stacked model, as the ingredients, rather 
than the non-parametric Nelson-Aalen estimates. This approach assumes constant hazards between timepoints, and hence more 
timepoints should be used when calculating predictions, i.e. at least {cmd:obs(100)}. This option is only valid with 
Markov models so can't currently be used with a {cmd:reset} model. The {cmd:aj} method is substantially faster than the 
large-sample simulation approach.

{phang}
{opt simulate} forces {cmd:predictms} to use large-sample simulation to calculate predictions, rather than numerical integration 
which is the default for survival, competing risks, and illness-death settings. Any other transition matrix structure will use 
simulation by default.

{phang}
{opt latent} specifies that latent times are used as the method of simulation, rather than the default method of Beyersmann et al. (2009). 
The latent times approach conducts a series of competing risks simulations, simulating an event time from each cause-specific hazard 
function, and taking the minimum as the observed event, compared to Beyersmann's method which simulated one event time from the 
total hazard, and a multinomial draw to determine the event that occured, based on the cause-specific hazard contributions. Both 
methods will give essentially identical results, but the Beyersmann method is generally more computationally efficient.

{phang}
{opt n(#)} defines the number of individual simulated trajectories through the multi-state model which are used to calculate predictions 
under the large-sample simulation approach. The default is {cmd:n(100,000)} unless {cmd:ci} is specified, then it is {cmd:n(10,000)}. 
Accuracy increases with a higher sample, reducing Monte Carlo error, but increasing computation time.

{phang}
{opt chintpoints(#)} defines the number of Gauss-Legendre quadrature points used to calculate any analytically intractible integrals. 
Numerical integration is required in the default prediction method for survival, competing risks and illness-death settings, and for 
simulating transition times from spline and complex time-dependent models. The defaukt is {bf:chintpoints(30)}. It is good 
practice to increase this number to ensure the predictions are reliable. Accuracy increases with a higher number, but increasing 
computation time.

{phang}
{opt save(name, [replace])} saves the simulated trajectories dataset(s) generated when using large-sample simulation to calculate predictions and 
confidence intervals. The dataset used to calculate point estimates will be saved in {it:name}.dta, and each simulated dataset generated in the 
{it:m}th bootstrap sample will be saved in {it:name{bf:i}}.dta. If the datasets already exist, you will need to add the {cmd:replace} 
option to overwrite.{p_end}

{phang}
{opt seed(#)} sets the simulation seed. When using large-sample simulation to obtain predictions, it is important to set the seed 
to ensure reproducibility.


{marker examples}{...}
{title:Examples}

{phang}
{bf:Example 1: A separate transition-specific illness-death model}
{p_end}

{pstd}
This dataset contains information on 2982 patients with breast cancer. Baseline is defined as time of surgery, and patients 
can experience relapse, relapse then death, or death with no relapse. Time of relapse is stored in {cmd:rf}, with event 
indicator {cmd:rfi}, and time of death is stored in {cmd:os}, with event indicator {cmd:osi}.
{p_end}

{pstd}Load example dataset:{p_end}
{cmd:    . use http://fmwww.bc.edu/repec/bocode/m/multistate_example}

{pstd}{helpb msset} the data (from the {cmd:multistate} package):{p_end}
{cmd:    . msset, id(pid) states(rfi osi) times(rf os)}

{pstd}Store the transition matrix:{p_end}
{cmd:    . mat tmat = r(transmatrix)}

{pstd}stset the data using the variables created by {cmd:msset}{p_end}
{cmd:    . stset _stop, enter(_start) failure(_status=1)}

{pstd}We fit separate Weibull models, so a fully stratified model, also allowing transition specific age effects:{p_end}

{cmd:    . stmerlin age if _trans1==1, distribution(weibull)}
{cmd:    . estimate store m1}

{cmd:    . stmerlin age if _trans2==1, distribution(weibull)}
{cmd:    . estimate store m2}

{cmd:    . stmerlin age if _trans3==1, distribution(weibull)}
{cmd:    . estimates store m3}

{pstd}Calculate transition probabilities for a patient with age 50:{p_end}
{cmd:    . predictms, transmatrix(tmat) models(m1 m2 m3) probability at1(age 50)}

{pstd}Create a stacked transition probability plot:{p_end}
{cmd:    . graphms}

{pstd}Calculate transition probabilities and length of stay for a patient with age 50:{p_end}
{cmd:    . predictms, transmatrix(tmat) models(m1 m2 m3) at1(age 50) probability los}

{pstd}Calculate the difference in transition probabilities for a patient with age 60 compared to a patient aged 50:{p_end}
{cmd:    . predictms, transmatrix(tmat) models(m1 m2 m3) probability at1(age 50) at2(age 60) difference}

{pstd}Calculate the ratio of transition probabilities for a patient with age 60 compared to a patient aged 50:{p_end}
{cmd:    . predictms, transmatrix(tmat) models(m1 m2 m3) probability at1(age 50) at2(age 60) ratio}

{pstd}Calculate differences and ratios of length of stay and transition probabilities for a patient with age 60 compared to 
a patient aged 50, with confidence intervals:{p_end}
{cmd:    . predictms, transmatrix(tmat) models(m1 m2 m3) probability at1(age 50) at2(age 60) los diff ratio ci}

{phang}
{bf:Example 2: A stacked, proportional, illness-death model}
{p_end}

{pstd}
Here we replicate the above example, but this time using a single {cmd:merlin} model, and using interactions between covariates 
and transition-specific indicators to allow for transition-specific effects.
{p_end}

{pstd}Load example dataset:{p_end}
{cmd:    . use http://fmwww.bc.edu/repec/bocode/m/multistate_example}

{pstd}{helpb msset} the data (from the {cmd:multistate} package):{p_end}
{cmd:    . msset, id(pid) states(rfi osi) times(rf os)}

{pstd}Store the transition matrix:{p_end}
{cmd:    . mat tmat = r(transmatrix)}

{pstd}stset the data using the variables created by {cmd:msset}{p_end}
{cmd:    . stset _stop, enter(_start) failure(_status=1)}

{pstd}We fit proportional Weibull models, with a common age effect:{p_end}
{cmd:    . merlin (_t age _trans2 _trans3 , family(weibull, failure(_d) ltruncated(_t0)))}

{pstd}We fit proportional Weibull models, with a transition-specific age effects:{p_end}
{cmd:    . merlin (_t age#_trans1 age#_trans2 age#_trans3 _trans2 _trans3 , family(weibull, failure(_d) ltruncated(_t0)))}

{pstd}All predictions are as above, but we no longer need the {cmd:models()} option. Calculate transition probabilities 
for a patient with age 50:{p_end}
{cmd:    . predictms, transmatrix(tmat) probability at1(age 50)}

{pstd}Create a stacked transition probability plot:{p_end}
{cmd:    . graphms}

{pstd}Calculate transition probabilities and length of stay for a patient with age 50:{p_end}
{cmd:    . predictms, transmatrix(tmat) at1(age 50) probability los}

{pstd}Calculate the difference in transition probabilities for a patient with age 60 compared to a patient aged 50:{p_end}
{cmd:    . predictms, transmatrix(tmat) probability at1(age 50) at2(age 60) difference}

{pstd}Calculate the ratio of transition probabilities for a patient with age 60 compared to a patient aged 50:{p_end}
{cmd:    . predictms, transmatrix(tmat) probability at1(age 50) at2(age 60) ratio}

{pstd}Calculate differences and ratios of length of stay and transition probabilities for a patient with age 60 compared 
to a patient aged 50, with confidence intervals:{p_end}
{cmd:    . predictms, transmatrix(tmat) probability at1(age 50) at2(age 60) los diff ratio ci}


{title:Authors}

{phang}
Michael J. Crowther (1,*), Paul C. Lambert (2,3)
{p_end}

{phang}
(1) Red Door Analytics, Stockholm, Sweden
{p_end}
{phang}
(2) Department of Medical Epidemiology and Biostatistics, Karolinska Institutet, Sweden
{p_end}
{phang}
(3) Biostatistics Research Group, Department of Health Sciences, University of Leicester, UK
{p_end}
{phang}
(*) michael@reddooranalytics.se
{p_end}

{phang}
Please report any errors you may find.{p_end}


{title:References}

{phang}
Beyersmann J, Latouche A, Buchholz A, Schumacher M. Simulating competing risks data in survival analysis. 
{it: Statistics in Medicine} 2009;28(6):956-971.
{p_end}

{phang}
Crowther MJ, Lambert PC. Parametric multi-state survival models: flexible modelling allowing transition-specific distributions with 
application to estimating clinically useful measures of effect differences. {it: Statistics in Medicine} 2017;36(29):4719-4742.
{p_end}

{phang}
Crowther MJ. Extended multivariate generalised linear and non-linear mixed effects models. 
{browse "https://arxiv.org/abs/1710.02223":https://arxiv.org/abs/1710.02223}
{p_end}

{phang}
Crowther MJ. merlin - a unified framework for data analysis and methods development in Stata. 
{browse "https://arxiv.org/abs/1806.01615":https://arxiv.org/abs/1806.01615}
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
Weibull CE, Lambert PC, Eloranta S, Andersson TM-L, Dickman PW, Crowther MJ. A multi-state model incorporating 
estimation of excess hazards and multiple time scales. {it: Statistics in Medicine} 2021; (In Press).
{p_end}


