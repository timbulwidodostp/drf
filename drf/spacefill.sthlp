{smcl}
{* 02Apr2014}{...}
{cmd:help spacefill}{right: ({browse "http://www.stata-journal.com/article.html?article=st0353":SJ14-3: st0353})}
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{hi:spacefill} {hline 2}}Perform space-filling location selection{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:spacefill}
{varlist} 
{ifin}
{weight}
[{cmd:,}
{cmdab:nd:esign}{cmd:(}{it:#}{cmd:)} {cmdab:design0}{cmd:(}{it:varlist}{cmd:)}
{cmdab:fix:ed}{cmd:(}{it:varname}{cmd:)}
{cmdab:exclude}{cmd:(}{it:varname}{cmd:)}
{cmd:p}{cmd:(}{it:#}{cmd:)}
{cmd:q}{cmd:(}{it:#}{cmd:)}
{cmdab:nn:frac}{cmd:(}{it:#}{cmd:)}
{cmdab:nnp:oints}{cmd:(}{it:#}{cmd:)}
{cmdab:nr:uns}{cmd:(}{it:#}{cmd:)}
{cmdab:stand:ardize}
{cmd:standardize2}
{cmd:standardize3}
{cmdab:spher:icize}
{cmd:ranks}
{cmdab:gen:erate}{cmd:(}{it:newvar}{cmd:)}
{cmdab:genmark:er}{cmd:(}{it:newvar}{cmd:)}
{cmdab:nover:bose}]

{p 4 8 2}
{cmd:aweight}s, {cmd:fweight}s, and {cmd:iweight}s are
allowed; see {help weight}.

{pstd}{it:varlist} and the {cmd:if} or {cmd:in} qualifier identify the data
from which the optimal subset is selected.


{title:Description}

{pstd}{cmd:spacefill} performs space-filling location selection: it
selects n points (a "design") from a candidate set of size N that best
"cover" the data according to a geometric coverage criterion.  This is
useful in site selection problems, but also in various nonparametric
estimation procedures, for example, when selecting (multivariate) knot
locations in spline regression analysis or choosing grid points.

{pstd}{cmd:spacefill} is an implementation of Royle and Nychka's (1998)
point-swapping algorithm.

{pstd}{cmd:spacefill} options allow forced inclusion or exclusion of
particular observations, user-specified initial design, and automatic
standardization of location coordinates.


{title:Options}

{phang}
{opt ndesign(#)} specifies n, the size of the design.  The default is
{cmd:ndesign(4)}.

{phang}
{opt design0(varlist)} identifies a set of initial designs identified by
observations with nonzero {it:varlist}.  If multiple variables are passed, one
optimization is performed on each initial design, and the selected design is
the one with best coverage.

{phang}
{opt fixed(varname)} identifies observations that are included in all designs
when {it:varname} is nonzero.

{phang}
{opt exclude(varname)} identifies observations excluded from all designs when
{it:varname} is nonzero.

{phang}
{opt p(#)} specifies a scalar value for the distance parameter for calculating
the distance of each location to the design; for example, p=-1 gives harmonic
mean distance, and p=-infinity gives the minimum distance.  The default is
{cmd:p(-5)}, as recommended in Royle and Nychka (1998).

{phang}
{opt q(#)} specifies a scalar value for the parameter q.  The default is
{cmd:q(1)} (the arithmetic mean).

{phang}
{opt nnfrac(#)} specifies the fraction of data to consider as nearest
neighbors in the point-swapping iterations.  Limiting checks to the nearest
neighbors improves speed but does not guarantee convergence to the best
design; therefore, setting {cmd:nruns(}{it:#}{cmd:)} is recommended.  The
default is {cmd:nnfrac(0.50)}.

{phang}
{opt nnpoints(#)} specifies the number of nearest neighbors
considered in the point-swapping iterations.  Limiting checks to nearest
neighbors improves speed. {cmd:nnfrac(}{it:#}{cmd:)} and
{cmd:nnpoints(}{it:#}{cmd:)} are mutually exclusive.

{phang}
{opt nruns(#)} sets the number of independent runs performed
on alternative random initial designs.  The selected design is the one with
best coverage across the runs.  The default is {cmd:nruns(5)}.

{phang}
{cmd:standardize} standardizes all variables in {it:varlist} to zero mean and
unit standard deviation (SD) before calculating distances between
observations.

{phang}
{cmd:standardize2} standardizes all variables in {it:varlist} to
zero mean and SD before calculating distances between observations, with
an estimator of the SD as 0.7413 times the interquartile range.

{phang}
{cmd:standardize3} standardizes all variables in {it:varlist} to zero median
and SD before calculating distances between observations, with an estimator of
the SD as 0.7413 times the interquartile range.

{phang}
{cmd:sphericize} transforms all variables in {it:varlist} into zero mean, SD,
and zero covariance using a Cholesky decomposition of the variance-covariance
matrix before calculating distances between observations.

{phang}
{cmd:ranks} transforms all variables in {it:varlist} into their (fractional)
ranks and uses distances between these observation ranks in each dimension to
evaluate distances between observations.

{phang}
{opt generate(newvar)} specifies the names for new variables containing the
locations of the best design points.  If one variable is specified, it is used
as a {it:stubname}; otherwise, the number of new variable names must match the
number of variables in {it:varlist}.

{phang}
{opt genmarker(newvar)} specifies the name of a new binary variable equal to
one for observations selected in the best design and zero otherwise.

{phang}
{cmd:noverbose} suppresses output display.


{title:Dependencies}

{pstd}
Ben Jann's {helpb moremata} (available from the Statistical Software
Components archive) is required for options {cmd:standardize2},
{cmd:standardize3}, and {cmd:ranks}.


{title:Examples}

{phang}{cmd:. clear}{p_end}
{phang}{cmd:. set obs 16}{p_end}
{phang}{cmd:. range lon -95 -80 16}{p_end}
{phang}{cmd:. range lat 36 46 11}{p_end}
{phang}{cmd:. fillin lon lat}{p_end}
{phang}{cmd:. gen byte sample = 0}{p_end}
{phang}{cmd:. save gridlatlon.dta , replace}{p_end}
{phang}{cmd:. clear}{p_end}
{phang}{cmd:. insheet using ozone2.txt}{p_end}
{phang}{cmd:. keep lat lon}{p_end}
{phang}{cmd:. gen byte sample = 1}{p_end}
{phang}{cmd:. append using gridlatlon}{p_end}
{phang}{cmd:. spacefill lon lat, ndesign(10) nnpoints(40)}{p_end}
{phang}{cmd:. spacefill lon lat, ndesign(20) nruns(10) nnfrac(0.3) gen(Best_Des) genmark(set1)}{p_end}
{phang}{cmd:. sysuse auto, clear}{p_end}
{phang}{cmd:. spacefill weight trunk, ndesign(6) standardize nruns(10) nnfrac(1) genmark(selected) p(-10) q(5)}{p_end}
{phang}{cmd:. twoway (scatter weight trunk) (scatter weight trunk if selected==1, msize(vlarge) mlabel(make))}{p_end}


{title:Reference}

{phang}Royle, J. A., and D. Nychka. 1998.  An algorithm for the construction
of spatial coverage designs with implementation in SPLUS.
{it:Computers and Geosciences} 24: 479-488.{p_end}


{title:Acknowledgments}

{pstd}This research is part of the project "Estimation of direct and
indirect causal effects using semiparametric and nonparametric methods"
project, which is supported by the Luxembourg "Fonds National de la
Recherche", cofunded under the Marie Curie Actions of the European
Commission (FP7-COFUND).  Van Kerm acknowledges funding for the project
"Information and Wage Inequality", which is supported by the Luxembourg
Fonds National de la Recherche (contract C10/LM/785657).{p_end}


{title:Authors}

{phang}Michela Bia{p_end}
{phang}CEPS/INSTEAD{p_end}
{phang}Esch-sur-Alzette, Luxembourg{p_end}
{phang}{browse "mailto:michela.bia@ceps.lu":michela.bia@ceps.lu}{p_end}

{phang}Philippe Van Kerm{p_end}
{phang}CEPS/INSTEAD, Luxembourg{p_end}
{phang}Esch-sur-Alzette, Luxembourg{p_end}
{phang}{browse "mailto:philippe.vankerm@ceps.lu":philippe.vankerm@ceps.lu}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 14, number 3: {browse "http://www.stata-journal.com/article.html?article=st0353":st0353}
{p_end}
