{smcl}
{* *! version 1.0.0}{...}
{vieweralsosee "[multistate] msset" "help msset"}{...}
{vieweralsosee "[multistate] msboxes" "help msboxes"}{...}
{vieweralsosee "[multistate] msaj" "help msaj"}{...}
{vieweralsosee "[multistate] predictms" "help predictms"}{...}
{vieweralsosee "[multistate] graphms" "help graphms"}{...}
{vieweralsosee "[merlin] merlin" "help merlin"}{...}
{viewerjumpto "Description" "multistate##description"}{...}
{title:Title}

{p2colset 5 19 19 2}{...}
{p2col:{helpb multistate} {hline 2}}Multi-state survival analysis{p_end}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:multistate} provides a set of commands, described below, for multi-state survival analysis. This includes data 
preparation tools, obtaining predictions from general continuous time multi-state survival models, both Markov and 
semi-Markov, and plotting utilities. Transition hazard models must be estimated using the {helpb stmerlin} or 
{helpb merlin} commands.
{p_end}

{pstd}
There are a number of commands in the {cmd:multistate} package, including:

{phang2}
{helpb msset} is a data preparation tool which converts a dataset from wide (one observation per subject, multiple time 
and status variables) to long (one observation for each transition of which a subject is at risk).

{phang2}
{helpb msboxes} creates a descriptive plot of the multi-state process through the transition matrix and 
numbers at risk.

{phang2}
{helpb msaj} calculates the non-parametric Aalen-Johansen estimates of transition probabilities, and the length of stay in each 
state.

{phang2}
{helpb predictms} calculates a variety of predictions from a Markov or semi-Markov multi-state survival model, including transition 
probabilities, length of stay (restricted mean time in each state), the probability of ever visiting each state and transition specific 
hazard and survival functions. Predictions are made at user-specified covariate patterns. Differences and ratios of predictions across 
covariate patterns can also be calculated. Standardised (study population-averaged) predictions can be obtained. Confidence intervals 
for all quantities are available. User-defined predictions can also be calculated by providing a user-written Mata function, to provide 
complete flexibility. {helpb predictms} can be used with a general transition matrix (cyclic or acyclic), and allows the use of 
transition-specific timescales.

{phang2}
{helpb graphms} creates stacked transition probability plots, following a {helpb predictms} call. 


{title:Website}

{pstd}
Visit the {bf:{browse "https://www.mjcrowther.co.uk/software/multistate":multistate homepage}} for more information on tutorials, 
short courses, and the package version history.


{title:Authors}

{phang}
Michael J. Crowther (1,2,*), Micki Hill (3), Paul C. Lambert (2,3)
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
Weibull CE, Lambert PC, Eloranta S, Andersson TM-L, Dickman PW, Crowther MJ. A multi-state model incorporating 
estimation of excess hazards and multiple time scales. {it: Statistics in Medicine} 2021; (In Press).
{pstd}

