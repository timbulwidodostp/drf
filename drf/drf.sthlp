{smcl}
{* 30March2013}{...}
{cmd:help drf}{right: ({browse "http://www.stata-journal.com/article.html?article=st0352":SJ14-3: st0352})}
{hline}

{title:Title}

{p2colset 5 12 14 5}{...}
{p2col :{hi:drf} {hline 2}}Semiparametric estimation of dose-response functions{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 11 2}
{cmd:drf}
{varlist} 
{ifin}
{weight}{cmd:,}
{opt outcome(varname)}
{opt treatment(varname)}
{cmd:cutpoints(}{it:varname}{cmd:)}
{cmd:index(}{it:string}{cmd:)}
{cmd:nq_gps(}{it:#}{cmd:)}
{cmd:method(}{it:type}{cmd:)}
[{cmd:gps} {cmd:family(}{it:familyname}{cmd:)}
{cmd:link(}{it:linkname}{cmd:)}
{cmd:vce(}{it:vcetype}{cmd:)}
{cmd:nolog(}{it:#}{cmd:)}
{cmd:search}
{cmd:common(}{it:#}{cmd:)}
{cmdab:numov:erlap(}{it:#}{cmd:)}
{cmd:test_varlist(}{it:varlist}{cmd:)}
{cmd:test(}{it:type}{cmd:)}
{cmd:flag(}{it:#}{cmd:)}
{cmd:tpoints(}{it:vector}{cmd:)}
{cmd:npoints(}{it:#}{cmd:)}
{cmdab:nper:centiles(}{it:#}{cmd:)}
{cmd:det}
{cmd:delta(}{it:#}{cmd:)}
{cmdab:band:width(}{it:#}{cmd:)}
{cmd:nknots(}{it:#}{cmd:)}
{cmd:knots(}{it:#}{cmd:)}
{cmdab:stand:ardized}
{cmd:degree1(}{it:#}{cmd:)}
{cmd:degree2(}{it:#}{cmd:)}
{cmd:nknots1(}{it:#}{cmd:)}
{cmd:nknots2(}{it:#}{cmd:)}
{cmd:knots1(}{it:#}{cmd:)}
{cmd:knots2(}{it:#}{cmd:)}
{cmdab:add:itive}
{cmdab:estop:ts(}{it:string}{cmd:)}]

{pstd}
{cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s are allowed;
see {help weight}.


{title:Description}

{pstd}
{cmd:drf} i) estimates the generalized propensity score (GPS) under various
alternative parametric assumptions; ii) imposes the common support condition
as defined in Flores et al. (2012) and assesses the balance of covariates
after adjusting for the estimated GPS; and iii) estimates the dose-response
function using the estimated GPS by using either the nonparametric
inverse-weighting (IW) estimator developed in Flores et al. (2012) or a new
set of semiparametric estimators based on penalized spline techniques.


{title:Options}

    {title:Required}

{phang}
{cmd:outcome(}{it:varname}{cmd:)} specifies that {it:varname} is the outcome
variable.

{phang}
{cmd:treatment(}{it:varname}{cmd:)} specifies that {it:varname} is the
treatment variable.

{phang}
{cmd:cutpoints(}{it:varname}{cmd:)} divides the range or set of the possible
treatment values, Y, into intervals within which the balancing properties of
the GPS are checked using a "blocking on the GPS" approach.  {it:varname} is a
variable indicating to which interval each observation belongs.  This option
is required unless {cmd:flag()} is set to {cmd:0} (see below).

{phang}
{cmd:index(}{it:string}{cmd:)} specifies the representative point of the
treatment variable at which the GPS must be evaluated within each treatment
interval specified in {cmd:cutpoints()}.  {it:string} identifies either the
mean ({it:string} = {cmd:mean}) or a percentile ({it:string} = {cmd:p1}, ...,
{cmd:p100}) of the treatment.  This is used when checking the balancing
properties of the GPS using a "blocking on the GPS" approach.  This option is
required unless {cmd:flag()} is set to {cmd:0} (see below).

{phang}
{cmd:nq_gps(}{it:#}{cmd:)} specifies that for each treatment interval defined
in {cmd:cutpoints()}, the values of the GPS evaluated at the representative
point {cmd:index()} have to be divided into {it:#} (1 <= {it:#} <= 100)
intervals, defined by the quantiles of the GPS evaluated at the representative
point {cmd:index()}.  This is used when checking the balancing properties of
the GPS using a "blocking on the GPS" approach.  This option is required
unless {cmd:flag()} is set to {cmd:0} (see below).

{phang}
{cmd:method(}{it:type}{cmd:)} specifies the {it:type} of approach to be used
to estimate the dose-response function.  The approaches are bivariate
penalized splines ({it:type} = {cmd:mtpspline}), bivariate penalized radial
splines ({it:type} = {cmd:radialpspline}), or IW kernel ({it:type} =
{cmd:iwkernel}).

    {title:Global options}

{phang}
{cmd:gps} stores the estimated generalized propensity score in the
{cmd:gpscore} variable that is added to the dataset.  This option must not be
specified when running the {cmd:bootstrap}.

{phang}
{cmd:family(}{it:familyname}{cmd:)} specifies the distribution to be used to
estimate the GPS.  The available distributional families are Gaussian (normal)
({cmd:family(}{cmd:gaussian}{cmd:)}), inverse Gaussian
({cmd:family(}{cmd:igaussian}{cmd:)}), Gamma
({cmd:family(}{cmd:gamma}{cmd:)}), and Beta ({cmd:family(}{cmd:beta}{cmd:)}).
The default is {cmd:family(gaussian)}.  The Gaussian, inverse Gaussian, and
Gamma distributional families are fit using {cmd:glm}, and the Beta
distribution is fit using {cmd:betafit}.

{pstd}
The following four options are for the {cmd:glm} command, so they can be
specified only when the Gaussian, inverse Gaussian, or Gamma distribution is
assumed for the treatment variable.

{phang}
{cmd:link(}{it:linkname}{cmd:)} specifies the link function for the Gaussian,
inverse Gaussian, and Gamma distributional families.  The available links are
{cmd:link(identity)}, {cmd:link(log)}, and {cmd:link(pow)}, and the default is
the canonical link for the {cmd:family()} specified (see help for {cmd:glm}
for further details).

{phang}
{cmd:vce(}{it:vcetype}{cmd:)} specifies the type of standard error reported
for the GPS estimation when the Gaussian, inverse Gaussian, or Gamma
distribution is assumed for the treatment variable.  {it:vcetype} may be
{cmd:oim}, {cmdab:r:obust}, {cmdab:cl:uster} {it:clustvar}, {cmd:eim},
{cmd:opg}, {cmdab:boot:strap}, {cmdab:jack:knife}, {cmd:hac} {it:kernel},
{cmd:jackknife1} (see help for {cmd:glm} for further details).

{phang}
{cmd:nolog(}{it:#}{cmd:)} is a flag ({it:#} = 0,1) that suppresses the
iterations of the algorithm toward eventual convergence when running the
{cmd:glm} command.  The default is {cmd:nolog(0)}.

{phang}
{cmd:search} searches for good starting values for the parameters of the
generalized linear model used to estimate the generalized propensity score
(see {manhelp glm R} for further details).
    
    {title:Overlap options}

{phang}
{opt common(#)} is a flag ({it:#} = 0,1) that restricts the inference to the
subsample satisfying the common support condition when it is implemented
({it:#} = 1).  The default is {cmd:common(1)}.

{phang}
{opt numoverlap(#)} specifies that the common support condition is imposed by
dividing the sample into {it:#} groups according to {it:#} quantiles of the
treatment distribution.  By default, the sample is divided into 5 groups,
cutting at the 20th, 40th, 60th, and 80th percentiles of the distribution if
{cmd:common(1)}.

    {title:Balancing property assessment options}

{phang}
{cmd:test_varlist(}{it:varlist}{cmd:)} specifies that the balancing property
must be assessed for each variable in {it:varlist}.  The default
{cmd:test_varlist()} consists of all the variables used to estimate the GPS.

{phang}
{cmd:test(}{it:type}{cmd:)} allows users to specify whether the balancing
property is to be assessed using a "blocking on the GPS" approach employing
either standard two-sided t tests ({cmd:test(t_test)}) or Bayes factors
({cmd:test(Bayes_factor)}) or using a model-comparison approach with a
likelihood-ratio test {cmd:test(}{cmd:L_like}{cmd:)}.

{pmore}
The "blocking on the GPS" approach, based on standard two-sided t tests
provides the values of the test statistics before and after adjusting for the
GPS for each pretreatment variable included in
{cmd:test_varlist()} and for each prefixed treatment
interval specified in {cmd:cutpoints()}.  Specifically, let p be the number of
control variables in {cmd:test_varlist()}, and let H be the number of
treatment intervals specified in {cmd:cutpoints()}.  Then the program
calculates and shows p x H values of the test statistic before and after
adjusting for the GPS, where the adjustment is done by dividing the values of
the GPS evaluated at the representative point {opt index()} into the
number of intervals specified in {opt nq_gps()}.

{pmore}
The model-comparison approach uses a likelihood-ratio test to compare an
unrestricted model for T_i, including all the covariates plus the GPS (up to a
cubic term), with a restricted model that sets the coefficients of all
covariates to 0.  By default, both the "blocking on the GPS" approach and the
model-comparison approach are applied.

{phang}
{cmd:flag(}{it:#}{cmd:)} allows the user to specify that {cmd:drf} estimates
the GPS without performing the balancing test.  The default is {cmd:flag(1)},
which means that the balancing property must be assessed.

    {title:Dose-response function options}

{phang}
{cmd:tpoints(}{it:vector}{cmd:)} indicates that the the dose-response function
is evaluated at each level of the treatment in {it:vector}.  By default, the
{cmd:drf} program creates a vector with jth element equal to the jth observed
treatment value.  This option cannot be used with {cmd:npoints()} or
{cmd:npercentiles()} (see below).

{phang}
{cmd:npoints(}{it:#}{cmd:)} indicates that the dose-response function is
evaluated at each level of the treatment belonging to a set of evenly spaced
values t_0, t_1, ..., t_{it:#} that cover the range of the observed treatment.
This option cannot be used with {cmd:tpoints()} (see above) or
{cmdab:npercentiles()} (see below).

{phang}
{cmdab:npercentiles(}{it:#}{cmd:)} indicates that the dose-response function
is evaluated at each level of the treatment corresponding to the percentiles
t_q0, t_q1, ..., t_q{it:#} of the treatment's empirical distribution.
This option cannot be used with {cmd:tpoints()} or {cmd:npoints()} (see
above).

{phang}
{cmd:det} displays more detailed output on the dose-response function
estimation.  When {cmd:det} is not specified, the program displays only the
chosen dose-response function estimator: {cmd:method(radialpspline)},
{cmd:method(mtpspline)}, or {cmd:method(mtpspline)}.

{phang}
{cmd:delta(}{it:#}{cmd:)} specifies that {cmd:drf} also estimate the
treatment-effect function u(t + {it:#}) - u(t).  The default is
{cmd:delta(0)}, which means that {cmd:drf} estimates only the dose-response
function, u(t).

    {title:Options for the IW kernel estimator ({cmd:iwkernel})}

{phang}
{opt bandwith(#)} specifies the bandwidth to be used.  By default, the global
bandwidth is chosen using the automatic procedure described in Fan and Gijbels
(1996).  This procedure estimates the unknown terms in the optimal global
bandwidth by using a global polynomial of order p + 3, where p is the order of
the local polynomial fitted.

    {title:Options for the radial penalized spline estimator ({cmd:radialpspline})}

{phang}
{opt nknots(#)} specifies the number of knots to be selected in the
two-dimensional space of the treatment variable and the GPS.  The default
{cmd:max(20, min(}n{cmd:/4, 150)))}, where n is the number of unique (T_i,
R_i) (Ruppert, Wand, and Carroll 2003).  When this option is specified, the
subroutines {cmd:radialpspline} and {cmd:spacefill} (Bia and Van Kerm 2014)
are called.  This option cannot be used with the {cmd:knots()} option (see
below).

{phang}
{cmd:knots(}{it:numlist}{cmd:)} specifies the list of knots for the treatment
and the GPS variable.  This option cannot be used with the {cmd:nknots()}
option (see above).

{phang}
{cmd:standardized} implies that the {cmd:spacefill} algorithm standardizes the
treatment variable and the GPS variables before selecting the knots.  The
knots are chosen using the standardized variables.

    {title:Options for the tensor-product penalized spline estimator ({cmd:mtpspline})}

{phang}
{cmd:degree1(}{it:#}{cmd:)} specifies the power of the treatment variable
included in the penalized spline model.  The default is {cmd:degree1(1)}.

{phang}
{cmd:degree2(}{it:#}{cmd:)} specifies the power of the GPS included in the
penalized spline model.  The default is {cmd:degree2(1)}.

{phang}
{cmd:nknots1(}{it:#}{cmd:)} specifies the number ({it:#}) of knots for the
treatment variable.  The location of the K_kth knot is defined as
{(k+1)/({it:#}+2)}th sample quantile of the unique T_i for k = 1, ..., {it:#}.
The default is {cmd:nknots1(max(5, min(}n{cmd:/4, 35)))}, where n is the
number of unique T_i (Ruppert, Wand, and Carroll 2003).  This option cannot be
used with the {cmd:knots1()} option (see below).

{phang}
{cmd:nknots2(}{it:#}{cmd:)} specifies the number ({it:#}) of knots for the
GPS.  The location of the K_kth knot is defined as (k+1)/({it:#}+2)th sample
quantile of the unique R_i, for k = 1, ..., {it:#}.  The default is
{cmd:nknots2(max(5, min(}n{cmd:/4, 35)))}, where n is the number of unique R_i
(Ruppert, Wand, and Carroll 2003).  This option cannot be used with the
{cmd:knots2()} option (see below).

{phang}
{cmd:knots1(}{it:#}{cmd:)} specifies the list of knots for the treatment
variable.  This option cannot be used with the {cmd:nknots1()} option (see
above).

{phang}
{cmd:knots2(}{it:#}{cmd:)} specifies the list of knots for the GPS.  This
option cannot be used with the {cmd:nknots2()} option (see above).

{phang}
{cmd:additive} allows users to implement penalized splines using the additive
model without including the product terms.

    {title:Mutual options for the tensor-product and radial penalized spline estimators}

{pstd}
Mutual options for the tensor-product and radial penalized spline estimators
involve either the {cmd:mtpspline} subroutine or the {cmd:radialpspline}
subroutine, depending on which estimator is used.

{phang}
{opt estopts(string)} specifies all the possible options allowed when running
the {cmd:xtmixed} models to fit penalized spline models (see {helpb xtmixed}
for further details).


{title:Remarks} 

{pstd}
Please remember to use the {cmdab:update query} command before running this
program to make sure you have an up-to-date version of Stata installed.
Otherwise, this program might not run properly.

{pstd}
The treatment has to be continuous.

{pstd}
Make sure that the variables in {it:varlist} do not contain missing values.


{title:Examples}

{phang}
{cmd:. drf agew ownhs owncoll male tixbot workthen yearm1 yearm2 yearm3 yearm4, outcome(year6) treatment(prize)} {cmd:gps flag(1) tpoints(tp) cutpoints(cut)} {cmd:index(p50) nq_gps(5)} {cmd:numoverlap(3)} {cmd:method(}{cmd:radialpspline}{cmd:)} 
{cmd:family(}{cmd:gaussian}{cmd:)} {cmd:link(log) nknots(10) det delta(1)}

{phang}
{cmd:. drf agew ownhs owncoll male tixbot workthen yearm1 yearm2 yearm3 yearm4, outcome(year6) treatment(prize) gps flag(0) tpoints(tp) numoverlap(3) method(iwkernel) family(gaussian) link(log) nknots(10) det}
 

{title:References}

{phang}
Bia, M., and P. Van Kerm. 2014.
{browse "http://www.stata-journal.com/article.html?article=st0353":Space-filling location selection}.
{it:Stata Journal} 14: 605-622.

{phang}
Fan, J., and I. Gijbels. 1996.
{it:Local Polynomial Modelling and Its Applications}.
New York: Chapman and Hall/CRC.

{phang}
Flores, C. A., A. Flores-Lagunes, A. Gonzalez, and T. C. Neumann. 2012. 
Estimating the effects of length of exposure to instruction in a training
program: The case of job corps. {it:Review of Economics and Statistics} 
94: 153-171.

{phang}
Ruppert, D., M. P. Wand, and R. J. Carroll. 2003. 
{it:Semiparametric Regression}. Cambridge: Cambridge University Press.


{title:Authors}

{pstd}Michela Bia{p_end}
{pstd}CEPS/INSTEAD{p_end}
{pstd}Esch-Sur-Alzette, Luxembourg{p_end}
{pstd}{browse "mailto:michela.bia@ceps.lu":michela.bia@ceps.lu}{p_end}

{pstd}Carlos A. Flores{p_end}
{pstd}Department of Economics{p_end}
{pstd}California Polytechnic State University{p_end}
{pstd}San Luis Obispo, CA{p_end}
{pstd}{browse "cflore32@calpoly.edu":cflore32@calpoly.edu}{p_end}

{pstd}Alfonso Flores-Lagunes{p_end}
{pstd}Department of Economics{p_end}
{pstd}State University of New York, Binghamton{p_end}
{pstd}Binghamton, NY{p_end}
{pstd}{browse "aflores@binghamton.edu":aflores@binghamton.edu}{p_end}

{pstd}Alessandra Mattei{p_end}
{pstd}Department of Statistics, Informatics, Applications, "Giuseppe Parenti"{p_end}
{pstd}University of Florence{p_end}
{pstd}Florence, Italy{p_end}
{pstd}{browse "mailto:mattei@disia.unifi.it":mattei@disia.unifi.it}{p_end}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 14, number 3: {browse "http://www.stata-journal.com/article.html?article=st0352":st0352}{p_end}
