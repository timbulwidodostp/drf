* ! v1.0.3, 2014-05-7 Bia, Mattei, Flores, Flores-Lagunes: second production release! (SJ submission) 
	*(–- seed option excluded, after some trials putting a seed from outside seem to suffice for making ///
	* the spacefill replicable); U and gps vars changed into temporary variables.
* ! v1.0.2, 2014-02-25 Bia: third production release! (SJ submission) (–- nolog option inserted) and a bug in glm models fixed.
* ! v1.0.1, 2014-02-25 Bia, Mattei, Flores, Flores-Lagunes: second production release! (SJ submission) (–- seed option inserted)
* ! v1.0.0, 2013-03-30 Bia, Mattei, Flores, Flores-Lagunes: first production release! (SJ submission)
	* March 2013
	* Program: drf
	* INPUT 
	* covariates

	* Required option:
	* outcome(varname) = Name of the outcome 
	* t(varname) 	     = Name of the treatment variable
	* gps  = Name that contains the estimated gpscore when gps is "activated"
	* predict(string)  = Name of the variable which contains the hat treatment values
	* drf(string) = Name of the variable which contains the estimated values of the dose-response function

	* predict(string)
	* theta(string)
	* gps         
	* dose(namelist)
	
	program define drf, eclass
	version 11
	
	#delimit ;
	syntax varlist [if] [in] [fweight iweight pweight],
	outcome(varname)
	treatment(varname)        
	method(string)
	[  
	gps
	family(string)
	link(string)
	vce(string) 
	nolog(int 0)
	search
	common(int 1)
	NUMOVerlap(int 5)
	tpoints(namelist)
	npoints(numlist)
	NPERcentiles(numlist)
	cutpoints(varname numeric)
	index(string)
	nq_gps(numlist)
	test_varlist(varlist)
	test(string) 
	flag(int 1)
	det                        
	delta(real 0)
	/// /*OPTIONS FOR IW-Kernel*/
	BANDwidth(numlist)
	/// /*OPTIONS FOR MTPSPLINE*/
	degree1(numlist)              ///
	degree2(numlist)               ///
    nknots1(numlist)              ///
	nknots2(numlist)              ///
    knots1(numlist sort)         ///
	knots2(numlist sort)         ///
	ADDitive      	          ///	
	/// /*OPTIONS FOR RADIALPSPLINE*/
	nknots(numlist)   ///
	knots(numlist sort)   ///
	STANDardized		///
	/// /*OPTIONS FOR BOTH MTPSPLINE AND RADIALPSPLINE*/
	ESTOPts(str asis)            ///
	*filename(string) /*PER IL LOG FILE DOVE SALVARE I RISULTATI DEL BILANCIAMENTO?*/
	]
	;
	#delimit cr
	
	*counting the number of observations
	qui count
	local obs = r(N) 
	
	
	/*If weights are specified, create local variable*/

	if "`weight'" != ""{
	tempvar wei 
	qui gen double `wei' `exp'
	local w [`weight' = `wei']
	}

	tokenize `varlist'
	marksample touse
 
 
	local k: word count `varlist'
	tempvar gpscore
	
	*confirm new variable `gpscore'

	if("`method'" != "mtpspline"){
		if("`degree1'" != ""  | "`degree2'" != "" | "`nknots1'" != "" | "`nknots2'" != "" | ///
			`"`knots1'"'!=""  | `"`knots2'"'!=""  | `"`additive'"' != `""'){
			di as error "An option for mtpspline has been selected, but the estimation method is not mtpspline"
			exit
		}
	}
	
	if("`method'" != "radialpspline"){
		if("`nknots'" != "" | `"`knots'"'!= `""'  | `"`standardized'"' != `""'){
			di as error "An option for radialspline has been selected, but the estimation method is not radialspline"
			exit
		}
	}	
	
	


	di in ye _newline(1)   "******************************************************"
	di in ye	     	   "Algorithm to estimate the generalized propensity score "
	di in ye	     	   "****************************************************** "


	di _newline(3) "Estimation of the propensity score "


	
	if "`family'" == "poisson" |  "`family'" =="nbinomial"  | "`family'" =="binomial" {
	di as error "Poisson, binomial and Negative binomial distributions are not supported"
	exit 
	}
	
	tempvar etreat
	
if("`family'" == "igaussian" | "`family'" == "gamma"){
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" == "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'") 
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" == ""  {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  vce("`vce'") 
	}
	else if "`vce'"=="" & ("`nolog'" =="1") & "`search'" == "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  nolog 
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" == "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  vce("`vce'") nolog 
	}
	else if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" != "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'") search
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" != ""  {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  vce("`vce'") search
	}
	else if "`vce'"=="" & ("`nolog'" =="1") & "`search'" != "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  nolog  search
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" != "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  vce("`vce'") nolog  search
	}
sca scale_parameter = e(phi) 
}


if("`family'" == "" | "`family'" == "gaussian"){
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" == "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'") 
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" == ""  {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  vce("`vce'") 
	}
	else if "`vce'"=="" & ("`nolog'" =="1") & "`search'" == "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  nolog 
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" == "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  vce("`vce'") nolog 
	}
	
	else if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" != "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'") search
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" != ""  {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  vce("`vce'") search
	}
	else if "`vce'"=="" & ("`nolog'" =="1") & "`search'" != "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  nolog search
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" != "" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse', ///
	family("`family'") link("`link'")  vce("`vce'") nolog  search
	}
sca scale_parameter = e(phi) 
}

	
if("`family'" == "beta"){
	if("`link'" !=""){
	di as err "The link function is not allowed using beta density"
	exit
	}
	if("`vce'" !="" & "`vce'"!="robust"){
	di as err "Beta only supports vce(robust) option (Huber/White/sandwich estimator of variance) "
	exit
	}
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") {
	betafit `treatment' [`weight'`exp'] if `touse', mu(`varlist')  
	}
	if "`vce'"=="robust" & ("`nolog'" =="0"|"`nolog'" ==" "){
	betafit `treatment' [`weight'`exp'] if `touse', mu(`varlist') robust 
	}
	if "`vce'"=="" & ("`nolog'" =="1") {
	betafit `treatment' [`weight'`exp'] if `touse', mu(`varlist') nolog  
	}
	if "`vce'"=="robust" & ("`nolog'" =="1"){
	betafit `treatment' [`weight'`exp'] if `touse', mu(`varlist') robust nolog
	}
sca scale_parameter = r(est) 
}
qui predict double `etreat'  if e(sample) & `touse'		


	tempname theta1
	if "`family'" == "" | "`family'"== "gaussian" {
		qui gen `theta1' = `etreat'
		sca theta2 = scale_parameter
		qui gen double `gpscore' = normalden(`treatment', `theta1', sqrt(theta2))  if `touse'
	}

	 
	if "`family'"== "igaussian" {
		qui gen `theta1' = `etreat'
		*see GLM models and extensions for the IG pronbability density distribution (Stata press-Hardin and Hilbe 2007) pg.105-109
		** computing sigma2 (pg. 109)
		sca theta2 = sqrt(scale_parameter)
		** computing the IG probability density distribution (pg. 105)
		igaussianden `gpscore' `treatment' `theta1' theta2 if `touse'
	}


	if "`family'"== "gamma" {
		sca theta2 = 1/scale_parameter
		qui gen `theta1' = `etreat'/theta2
		qui gen double `gpscore' = gammaden(theta2, `theta1', 0, `treatment') if `touse'
		*di theta2
		*gen b =`theta1'
	}


	if "`family'"== "beta" {
		tempname theta2
		qui gen `theta1' = (`etreat')*scale_parameter
		qui gen `theta2' = (1-`etreat')*scale_parameter
		qui gen double `gpscore' = betaden(`theta1',`theta2',`treatment') if `touse'
	}
	

tempvar overlap
qui gen `overlap'=1 if `touse'


if("`common'" == " " | "`common'" == "1") {
	
	
	forvalues i = 1/`numoverlap'{
	
	local perc = (100/`numoverlap')*`i'
	
	
		if `i' < `numoverlap'{
		tempvar Q`i'
		quietly egen `Q`i''=pctile(`treatment') if `touse', p(`perc')
		}
			
		tempvar G`i' TT`i' gps_G`i' 
		
		qui gen `G`i'' = 0
	
		if `i' ==1{
		quietly replace `G`i'' =1 if `treatment' < `Q`i'' & `touse'
			}
	
		if `i' > 1 & `i' < `numoverlap'{ 
		local befperc = `i'-1
		quietly replace `G`i''=1 if  `treatment' >= `Q`befperc''   &  `treatment' < `Q`i'' & `touse'
		}
	
		if `i'==`numoverlap'{
		local befperc =`i'-1
		quietly replace `G`i''=1 if `treatment'>=`Q`befperc'' & `touse'
		}
		
		
	*** Get median value of `treatment' for each of the `numoverlap' groups, 
	*** and calculate the GPS, f(t|xi), at each of these values for all units 
	*** (i.e., `numoverlap' gps for each individual).
 
	 qui egen `TT`i''=pctile(`treatment') if `G`i''==1 & `touse', p(50)
	 qui sum `TT`i''
	 qui replace `TT`i'' = r(mean) if `TT`i''==. & `touse'
	
	if "`family'" == "" | "`family'"== "gaussian" {
	qui gen double `gps_G`i'' = normalden(`TT`i'', `theta1', sqrt(theta2))  if `touse'
	}
	
	if "`family'"== "igaussian" {
	igaussianden `gps_G`i'' `TT`i'' `theta1' theta2 if `touse' 
	}
	
	if "`family'"== "gamma" {
	qui gen double `gps_G`i''= gammaden(theta2, `theta1',0,`TT`i'') if `touse' 
	}

	if "`family'"== "beta" {
	qui gen double `gps_G`i'' = betaden(`theta1',`theta2',`TT`i'') if `touse'
	}
	
	/*
	if "`family'"== "binomial" {
	qui gen double `gpscore' = binomialp(1,1,`theta1') if `touse'
	}
	*/
	

	* Generate all minima and maxima used in overlap condition.
	tempvar  mx_G`i'_0 mx_G`i'_1 mn_G`i'_0 mn_G`i'_1
	
	
	quietly sum `gps_G`i'' if `G`i''==0
	quietly gen `mx_G`i'_0'= r(max)
	quietly gen `mn_G`i'_0'= r(min)
		
		
	quietly sum `gps_G`i'' if `G`i''==1
	quietly gen `mx_G`i'_1'= r(max)
	quietly gen `mn_G`i'_1'= r(min)


	* We get the maximum of the minimums and the minimum of the maximum for each of the
	* numoverlap comparisons. These are the values that bind, the others can be ignored.   
	tempvar max_min_`i' min_max_`i'
	qui gen     `max_min_`i''= 0   // Gives the maximum of the minima. Only one of the two lines below holds.
	qui replace `max_min_`i''= `mn_G`i'_0' if  `mn_G`i'_0' >= `mn_G`i'_1'   // Set equal to avoid ties. If equal does not matter.
	qui replace `max_min_`i''= `mn_G`i'_1' if  `mn_G`i'_1' > `mn_G`i'_0'
		
	qui gen `min_max_`i''=0   // Gives the minimum of the maxima. Only one of the two lines below holds.
	qui	replace `min_max_`i''= `mx_G`i'_0' if  `mx_G`i'_0' <= `mx_G`i'_1'   // Set equal to avoid ties. If equal does not matter.
	qui	replace `min_max_`i''= `mx_G`i'_1' if  `mx_G`i'_1' < `mx_G`i'_0'
	
	*** Impose overlap condition as the intersection of the groups.  
	qui replace `overlap' = 0 if (`gps_G`i'' >=  `min_max_`i'') | (`gps_G`i'' <=  `max_min_`i'')
	} /*End loop over i= 1/`numoverlap'*/
	
	
	di as text "Note: The common support condition is imposed" 
	
	qui count if `touse' & `overlap'
	local drop_obs = `obs' - r(N)

	di in ye _newline(1) "***************************************************************** "
	di                   "`drop_obs' observations are dropped after imposing common support "
	di                   "***************************************************************** "
	
	} /*End loop over common==1*/


label var `gpscore' "drf_gpscore"
sum `gpscore' if `touse' & `overlap'==1, det 
* gpscore gpscore2 gpscore3 (that we will use below) cannot be tempvars
* if we want to display them not as local vars 
* when running the glm models
qui gen double drf_gpscore = `gpscore'
forvalues p =2/3{
qui gen double drf_gpscore`p' = `gpscore'^`p'
}



di in ye _newline(1) "******************************************** "
di                   "End of the algorithm to estimate the gpscore "
di                   "******************************************** "


if("`flag'" == " " | "`flag'" == "1"){

	
if ("`test'"=="" | "`test'"=="t_test" | "`test'"=="Bayes_factor")  { /* Begin of test t_test BF_test */	



di in ye _newline(1)   "*****************************************************************"
di    			"Beginning of the assessment of the balancing property of the GPS"
di             		"*****************************************************************"

*******************************************************************
*We split the treatment range in sub - intervals
*The bounds of each treatment sub-interval are given by "cutpoints"
*******************************************************************

tempvar broken_t
qui xtile `broken_t' = `treatment' if `touse' & `overlap'==1, cutpoints(`cutpoints')
qui tab `broken_t' if `touse' & `overlap'==1, gen(`broken_t') 
local nblock_t = r(r)


di in ye _newline(1) "*******************************************************************************************"
di                   "The set of the potential treatment values (range of T) is divided into `nblock_t' intervals"
di 		     "*******************************************************************************************"


if("`index'" == "mean"){
	local i = 1
	while(`i' <= `nblock_t'){
		tempvar mean_t_`i' gpscore_`i'
		qui sum `treatment' if `broken_t'`i' ==1 & `touse' & `overlap'
		qui gen `mean_t_`i'' = r(mean)


		if "`family'" == "" | "`family'"== "gaussian" {
			qui gen double `gpscore_`i'' = normalden(`mean_t_`i'', `theta1', sqrt(theta2))  if `touse' & `overlap'
		}
		
		if "`family'" == "igaussian" {
			igaussianden `gpscore_`i'' `mean_t_`i'' `theta1' theta2 if `touse' & `overlap'
		}

		if "`family'" == "gamma" {
			qui gen double `gpscore_`i'' = gammaden(theta2, `theta1', 0, `mean_t_`i'') if `touse' & `overlap'
		}
			
		if "`family'"== "beta" {
		qui gen double `gpscore_`i'' = betaden(`theta1',`theta2',`mean_t_`i'') if `touse' & `overlap'
		}
	local i = `i' + 1
	}
}
		
		
		
	foreach x of numlist 1/100{
	 if("`index'" == "p`x'"){
	  local i = 1
	   while(`i' <= `nblock_t'){
		tempvar p`x'_t_`i' gpscore_`i' 
		qui egen `p`x'_t_`i'' = pctile(`treatment')   if  `broken_t'`i' ==1 & `touse' & `overlap', p(`x')
		qui sum `p`x'_t_`i'' if `touse' & `overlap'
		qui replace `p`x'_t_`i''  = r(mean)
		
		
		if "`family'" == "" | "`family'"== "gaussian" {
			qui gen double `gpscore_`i'' = normalden(`p`x'_t_`i'', `theta1', sqrt(theta2))  if `touse' & `overlap'
		}
		
		if "`family'" == "igaussian" {
			igaussianden `gpscore_`i'' `p`x'_t_`i'' `theta1' theta2 if `touse' & `overlap'
		}

		
		if "`family'" == "gamma" {
			*di theta2
			qui gen double `gpscore_`i'' = gammaden(theta2, `theta1', 0, `p`x'_t_`i'') if `touse' & `overlap'
		}	
		
		
		if "`family'" == "beta" {
		qui gen double `gpscore_`i'' = betaden(`theta1',`theta2',`p`x'_t_`i'') if `touse' & `overlap'
		}
		
	local i = `i' + 1
	}
  }
}
*********************************************************************
*We split the gpscore in nq_gps sub-classes in each treatment class 
*********************************************************************



	di in ye _newline(1)   "************************************************************************ "
	di    			"The values of the gpscore evaluated at the representative point of each "
	di			"treatment interval are divided into `nq_gps' intervals			 "
	di                     "*************************************************************************"

local i = 1
while(`i' <= `nblock_t'){
	tempvar broken_gps_`i'  
	qui xtile  `broken_gps_`i''   = `gpscore_`i''  if  `broken_t'`i' ==1 & `touse' & `overlap', n(`nq_gps')

	local j = 1
	while(`j' <= `nq_gps'){
		tempvar min_`i'`j'  
		tempvar max_`i'`j'
		qui sum `gpscore_`i''  if `broken_gps_`i'' == `j' & `touse' & `overlap'
		qui gen `max_`i'`j'' = r(max)
		local j = `j' + 1
	}

	qui replace `broken_gps_`i'' = 1  if  `gpscore_`i'' <= `max_`i'1'  & `broken_gps_`i'' == . & `touse' & `overlap'
	local j = 2
	while(`j' <= `nq_gps'){
		local k = `j'-1
		qui replace `broken_gps_`i'' = `j'  if `gpscore_`i'' > `max_`i'`k'' & `gpscore_`i'' <= `max_`i'`j''  & `broken_gps_`i'' == . & `touse' & `overlap'
		local j = `j' + 1
	}

	local i = `i' + 1
}

 


	di in ye _newline(1) "***********************************************************"
	di                   "Summary statistics of the distribution of the GPS evaluated" 
	di                   "at the representative point of each treatment interval"      
	di                   "***********************************************************"


	local i = 1
	while(`i' <= `nblock_t'){
	qui gen gps_`i' = `gpscore_`i'' 
	sum gps_`i'  if `overlap'
	drop gps_`i'
	local i = `i' + 1
	}


if !("`test'"=="" | "`test'"=="t_test" | "`test'"=="Bayes_factor" | "`test'"=="L_like") {
	di as error "Balancing Test `test' not recognized"
	exit 198
}


local test_varlist `test_varlist'
if ("`test_varlist'"=="") {
	local test_varlist `varlist'
	}

	
	di in ye _newline(2) "************************************************************************************************"
	di                   " Test that the mean of the pre-treatment variables is not different between units who belong to "
	di					 " a particular treatment interval and units who belong to all other treatment intervals       "
	di                   "************************************************************************************************"

local i = 1

while(`i' <= `nblock_t'){ /*Begin "while" on treatment classes*/
	foreach x of varlist `test_varlist'{
		tempvar ddiff_`x'_`i' vvar_diff_`x'_`i'  NNt_`x'_`i' NNc_`x'_`i'
		qui gen `ddiff_`x'_`i'' = 0
		qui gen `vvar_diff_`x'_`i'' = 0 
		qui gen `NNt_`x'_`i'' = 0
		qui gen `NNc_`x'_`i'' = 0 
}
 
 

foreach x of varlist `test_varlist'{
		qui ttest `x', by(`broken_t'`i') 
		qui replace `ddiff_`x'_`i'' = `ddiff_`x'_`i'' + ((r(N_1) + r(N_2))/_N)*(r(mu_1) - r(mu_2))
		qui replace `vvar_diff_`x'_`i'' = `vvar_diff_`x'_`i''  +  (((r(N_1) + r(N_2))/_N)^2)*(r(se)^2)
		qui replace `NNt_`x'_`i'' = `NNt_`x'_`i''  + ((r(N_1) + r(N_2))/_N)*r(N_2)
		qui replace `NNc_`x'_`i'' = `NNc_`x'_`i''  + ((r(N_1) + r(N_2))/_N)*r(N_1) 

		tempvar sse_diff_`x'_`i'  tt_diff_`x'_`i'  BFF_`x'_`i'
		qui gen `sse_diff_`x'_`i'' = sqrt(`vvar_diff_`x'_`i'')
		qui gen `tt_diff_`x'_`i'' = `ddiff_`x'_`i''/`sse_diff_`x'_`i''
		qui gen `BFF_`x'_`i'' = (3/2)*sqrt((`NNt_`x'_`i''*`NNc_`x'_`i'')/(`NNt_`x'_`i'' + `NNc_`x'_`i''))*((1+ ((`tt_diff_`x'_`i'')^2)/((`NNt_`x'_`i'' + `NNc_`x'_`i''-2)))^(-0.5*(`NNt_`x'_`i'' + `NNc_`x'_`i'')))
 

  }
 local i =`i' +1
 }


if ("`test'"=="" | "`test'"=="t_test") {/*Begin "if" t - test*/
tempvar tt_max
qui gen `tt_max'= 0

local i = 1
while(`i' <= `nblock_t'){ /*Begin "while" on treatment classes*/
	qui sum `treatment' if `broken_t'`i'==1 & `touse' & `overlap'
	local ttmin = r(min)
	local ttmax = r(max)
	
		di ""
		di as text    "Treatment Interval No `i' - [`ttmin', `ttmax']"   
		di ""
		di as text    "               Mean "          "       Standard   " 
		di as text    "               Difference"     "  Deviation   "   "t-value"   
		di ""

	tempvar tt_max`i'
	qui gen `tt_max`i''= 0

	foreach x of varlist `test_varlist'{
			di as text %12s abbrev("`x'            ",12) "  " as result %7.0g `ddiff_`x'_`i'' "      " as result %7.0g `sse_diff_`x'_`i''      "    " as result %7.0g `tt_diff_`x'_`i''  
			di ""    
	

	qui replace `tt_max`i''= abs(`tt_diff_`x'_`i'') if  abs(`tt_diff_`x'_`i'')>`tt_max`i''
  }
 qui replace `tt_max'= `tt_max`i'' if  `tt_max' < `tt_max`i''

 local i =`i' +1
 }/*End "while" on treatment classes*/
} /*End t_test*/

if ("`test'"=="Bayes_factor") { /*Begin "if" Bayes-Factor*/

tempvar BFF_min
qui gen `BFF_min'= .


local i = 1
while(`i' <= `nblock_t'){ /*Begin "while" on treatment classes*/
qui sum `treatment' if `broken_t'`i'==1 & `touse' & `overlap'
local ttmin = r(min)
local ttmax = r(max)
di ""
 di as text    "Treatment Group No `i' - [`ttmin', `ttmax']"   
 di ""
 di as text    "               Mean "          "       Standard   " 
 di as text    "               Difference"     "  Deviation   "   "Bayes-Factor"   
 di ""

tempvar BFF_min`i'
qui gen `BFF_min`i''= .

foreach x of varlist `test_varlist'{

	di as text %12s abbrev("`x'            ",12) "  " as result %7.0g `ddiff_`x'_`i'' "      " as result %7.0g `sse_diff_`x'_`i''      "    " as result %7.0g `BFF_`x'_`i''  
	di ""    


qui replace `BFF_min`i''= `BFF_`x'_`i'' if  `BFF_`x'_`i'' < `BFF_min`i''
}

qui replace `BFF_min'= `BFF_min`i'' if  `BFF_min' > `BFF_min`i''

local i =`i' +1
} /*End "while" on treatment classes*/

}/*End "if" BF*/

	di in ye _newline(2) "************************************************************************************"
	di                   "Test that the conditional mean of the pre-treatment variables given the generalized "
	di			   		 "propensity score is not different between units who belong to a particular treatment"
	di			   		 "interval and units who belong to all other treatment intervals					  "
	di 				     " (Mean differences adjusted by the GPS)											  "
	di                   "************************************************************************************"


local i = 1
while(`i' <= `nblock_t'){ /*Begin "while" on treatment classes*/
	foreach x of varlist `test_varlist'{
		tempvar diff_`x'_`i' var_diff_`x'_`i'  Nt_`x'_`i' Nc_`x'_`i'
		qui gen `diff_`x'_`i'' = 0
		qui gen `var_diff_`x'_`i'' = 0 
		qui gen `Nt_`x'_`i'' = 0
		qui gen `Nc_`x'_`i'' = 0 
}

local j = 1
while(`j' <= `nq_gps'){ /*Begin "while" on gpscore sub-classes in each treatment class*/
	foreach x of varlist `test_varlist'{
		qui ttest `x' if  `broken_gps_`i'' ==`j' , by(`broken_t'`i') 
		qui replace `diff_`x'_`i'' = `diff_`x'_`i'' + ((r(N_1) + r(N_2))/_N)*(r(mu_1) - r(mu_2))
		qui replace `var_diff_`x'_`i'' = `var_diff_`x'_`i''  +  (((r(N_1) + r(N_2))/_N)^2)*(r(se)^2)
		qui replace  `Nt_`x'_`i'' = `Nt_`x'_`i''  + ((r(N_1) + r(N_2))/_N)*r(N_2)
		qui replace  `Nc_`x'_`i'' = `Nc_`x'_`i''  + ((r(N_1) + r(N_2))/_N)*r(N_1) 
}
local j = `j' +1
}/*End "while" on gpscore sub-classes in each treatment class*/
foreach x of varlist `test_varlist'{
	tempvar se_diff_`x'_`i'  t_diff_`x'_`i'  BF_`x'_`i'
	qui gen `se_diff_`x'_`i'' = sqrt(`var_diff_`x'_`i'')
	qui gen `t_diff_`x'_`i'' = `diff_`x'_`i''/`se_diff_`x'_`i''
	qui gen `BF_`x'_`i'' = (3/2)*sqrt((`Nt_`x'_`i''*`Nc_`x'_`i'')/(`Nt_`x'_`i'' + `Nc_`x'_`i''))*((1+ ((`t_diff_`x'_`i'')^2)/((`Nt_`x'_`i'' + `Nc_`x'_`i''-2)))^(-0.5*(`Nt_`x'_`i'' + `Nc_`x'_`i'')))
}
local i =`i' +1
}/*End "while" on treatment classes*/


if ("`test'"=="" | "`test'"=="t_test") {/*Begin "if" t - test*/
tempvar t_max
qui gen `t_max'= 0


local i = 1
while(`i' <= `nblock_t'){ /*Begin "while" on treatment classes*/
	qui sum `treatment' if `broken_t'`i'==1 & `touse' & `overlap'
	local tmin = r(min)
	local tmax = r(max)
	
		di ""
		di as text    "Treatment Interval No `i' - [`tmin', `tmax']"   
		di ""
	    di as text    "               Mean "          "       Standard   " 
		di as text    "               Difference"     "  Deviation   "   "t-value"   
		di ""

	tempvar t_max`i'
	qui gen `t_max`i''= 0

	foreach x of varlist `test_varlist'{
		
		di as text %12s abbrev("`x'            ",12) "  " as result %7.0g `diff_`x'_`i'' "      " as result %7.0g `se_diff_`x'_`i''      "    " as result %7.0g `t_diff_`x'_`i''  
		di ""    


	qui replace `t_max`i''= abs(`t_diff_`x'_`i'') if  abs(`t_diff_`x'_`i'')>`t_max`i''
  }
qui replace `t_max'= `t_max`i'' if  `t_max' < `t_max`i''

local i =`i' +1
}/*End "while" on treatment classes*/


} /*End "if" t_test*/


if ("`test'"=="Bayes_factor") { /*Begin "if" Bayes-Factor*/

tempvar BF_min
qui gen `BF_min'= .


local i = 1
while(`i' <= `nblock_t'){ /*Begin "while" on treatment classes*/
qui sum `treatment' if `broken_t'`i'==1 & `touse' & `overlap'
local tmin = r(min)
local tmax = r(max)

 di ""
 di as text    "Treatment Gruop No `i' - [`tmin', `tmax']"   
 di ""
 di as text    "               Mean "          "       Standard   " 
 di as text    "               Difference"     "  Deviation   "   "Bayes-Factor"   
 di ""

tempvar BF_min`i'
qui gen `BF_min`i''= .

foreach x of varlist `test_varlist'{

 di as text %12s abbrev("`x'            ",12) "  " as result %7.0g `diff_`x'_`i'' "      " as result %7.0g `se_diff_`x'_`i''      "    " as result %7.0g `BF_`x'_`i''  
 di ""    

qui replace `BF_min`i''= `BF_`x'_`i'' if  `BF_`x'_`i'' < `BF_min`i''
}

qui replace `BF_min'= `BF_min`i'' if  `BF_min' > `BF_min`i''

local i =`i' +1
} /*End "while" on treatment classes*/


}/*End "if" BF*/


} /* End test t_test BF_test */



************************************************************************
* Begin of the L_likelihood Test for Unrestricted and Restricted Model *
************************************************************************

if ("`test'"=="" | "`test'"=="L_like") {

	di in ye _newline(2) "**********************************************************"
	di					 "Log-Likelihood test for Unrestricted and Restricted Model "
	di 					 "**********************************************************"



	di in ye _newline(2) "****************************************************"
	di                   "              Unrestricted Model                    "
	di					 " link(E[T]) = GPSCORE + GPSCORE^2  + GPSCORE^3  + X "
	di                   "****************************************************"


if("`family'" == "igaussian" | "`family'" == "gamma"){
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" =="" {
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" ==""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") 
	}
	else if "`vce'"=="" & ("`nolog'" =="1") & "`search'" ==""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog 
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" =="" {
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'") nolog 
	}
	else if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" !="" {
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") search
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" !=""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'")  search
	}
	else if "`vce'"=="" & ("`nolog'" =="1") & "`search'" !=""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog search
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" !="" {
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'") nolog  search
	}
}


if("`family'" == "" | "`family'" == "gaussian"){
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" =="" {
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") 
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" ==""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'")
	}
	else if "`vce'"=="" & "`nolog'" =="1" & "`search'" ==""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  nolog
	} 
	else if "`vce'"!="" & "`nolog'" =="1" & "`search'" ==""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'") nolog
	} 	
	else if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" !="" {
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") search
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" !=""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'") search
	}
	else if "`vce'"=="" & "`nolog'" =="1" & "`search'" !=""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  nolog search
	} 
	else if "`vce'"!="" & "`nolog'" =="1" & "`search'" !=""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'") nolog search
	} 
} 



	
if("`family'" == "beta"){
	if("`link'" !=""){
	di as err "The link function is not allowed using beta density"
	exit
	}
	if("`vce'" !="" & "`vce'"!="robust"){
	di as err "Beta only supports vce(robust) option (Huber/White/sandwich estimator of variance) "
	exit
	}	
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" "){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(drf_gpscore drf_gpscore2 drf_gpscore3 `varlist')  
	}
	if "`vce'"=="robust" & ("`nolog'" =="0"|"`nolog'" ==" "){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(drf_gpscore drf_gpscore2 drf_gpscore3 `varlist') robust 
	}
	if "`vce'"=="" & ("`nolog'" =="1"){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(drf_gpscore drf_gpscore2 drf_gpscore3 `varlist')  nolog
	}
	if "`vce'"=="robust" & ("`nolog'" =="0"|"`nolog'" ==" "){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(drf_gpscore drf_gpscore2 drf_gpscore3 `varlist') robust nolog
	}
}

	
mat LL_U1 = e(ll)
estimates store U1


di in ye _newline(2) "********************************************************"
di                   " Restricted Model: Pretreatment variables are excluded  "
di					 "     link(E[T]) = GPSCORE + GPSCORE^2  + GPSCORE^3      "
di                   "********************************************************"



if("`family'" == "igaussian" | "`family'" == "gamma"){
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" ==""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")
	}
	if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" ==""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'")  
	}
	if "`vce'"=="" & ("`nolog'" =="1")& "`search'" ==""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog  
	}
	else if "`vce'"!="" & ("`nolog'" =="1")& "`search'" ==""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") nolog 
	}
	else if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" !=""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") search
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" !=""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'")  search
	}
	if "`vce'"=="" & ("`nolog'" =="1")& "`search'" !=""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog  search
	}
	else if "`vce'"!="" & ("`nolog'" =="1")& "`search'" !=""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") nolog search
	}
}



if("`family'" == "" | "`family'" == "gaussian"){
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" == ""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") 
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ")& "`search'" == ""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") 
	}
	else if "`vce'"=="" & ("`nolog'" =="1")& "`search'" == ""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog 
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" == ""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") nolog
	}
	else if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" != ""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") search
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ")& "`search'" != ""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'")  search
	}
	else if "`vce'"=="" & ("`nolog'" =="1")& "`search'" != ""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog search
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" != ""{
	glm `treatment' drf_gpscore drf_gpscore2 drf_gpscore3 [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") nolog search
	} 
} 



	
if("`family'" == "beta"){
	if("`link'" !=""){
	di as err "The link function is not allowed using beta density"
	exit
	}
	if("`vce'" !="" & "`vce'"!="robust"){
	di as err "Beta only supports vce(robust) option (Huber/White/sandwich estimator of variance) "
	exit
	}
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" "){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(drf_gpscore drf_gpscore2 drf_gpscore3)  
	}
	if "`vce'"=="robust" & ("`nolog'" =="0"|"`nolog'" ==" "){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(drf_gpscore drf_gpscore2 drf_gpscore3) robust 
	}
	if "`vce'"=="" & ("`nolog'" =="1"){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(drf_gpscore drf_gpscore2 drf_gpscore3)  nolog
	}
	if "`vce'"=="robust" & ("`nolog'" =="1"){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(drf_gpscore drf_gpscore2 drf_gpscore3) robust nolog
	}
}


	
mat LL_R1 = e(ll)
estimates store R1

qui lrtest U1 R1

mat Tstat1 = r(chi2)
mat Pval1 = r(p)
mat NofR1 = r(df)

di in ye _newline(2) "**********************************************************"
di                   " Restricted Model: GPS terms are excluded (link(E[T]) = X)"
di                   "**********************************************************"



if("`family'" == "igaussian" | "`family'" == "gamma"){
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" =="" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") 
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ")& "`search'" =="" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'") 
	}
	else if "`vce'"=="" & ("`nolog'" =="1")& "`search'" ==""{
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog 
	}
	else if "`vce'"!="" & ("`nolog'" =="1")& "`search'" ==""{
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") nolog 
	}
	else if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" ") & "`search'" !="" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") search
	}
	else if "`vce'"!="" & ("`nolog'" =="0"|"`nolog'" ==" ")& "`search'" !="" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'") search
	}
	else if "`vce'"=="" & ("`nolog'" =="1") & "`search'" !=""{
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog search
	}
	else if "`vce'"!="" & ("`nolog'" =="1")& "`search'" !=""{
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") nolog  search
	}
}



if("`family'" == "" | "`family'" == "gaussian") {
	if "`vce'"=="" & ("`nolog'" ==" "|"`nolog'" =="0") & "`search'" =="" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") 
	}
	else if "`vce'"!="" & ("`nolog'" ==" "|"`nolog'" =="0") & "`search'" ==""{
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'") 
	}
	else if "`vce'"=="" & ("`nolog'" =="1") & "`search'" ==""{
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog 
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" =="" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") nolog 
	} 
	else if "`vce'"=="" & ("`nolog'" ==" "|"`nolog'" =="0") & "`search'" !="" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") search
	}
	else if "`vce'"!="" & ("`nolog'" ==" "|"`nolog'" =="0") & "`search'" !=""{
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") vce("`vce'") search
	}
	else if "`vce'"=="" & ("`nolog'" =="1") & "`search'" !=""{
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'") nolog search
	}
	else if "`vce'"!="" & ("`nolog'" =="1") & "`search'" !="" {
	glm `treatment' `varlist' [`weight'`exp'] if `touse' & `overlap', ///
	family("`family'") link("`link'")  vce("`vce'") nolog  search
	} 
} 



	
if("`family'" == "beta"){
	if("`link'" !=""){
	di as err "The link function is not allowed using beta density"
	exit
	}
	if("`vce'" !="" & "`vce'"!="robust"){
	di as err "Beta only supports vce(robust) option (Huber/White/sandwich estimator of variance) "
	exit
	}
	if "`vce'"=="" & ("`nolog'" =="0"|"`nolog'" ==" "){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(`varlist')  
	}
	if "`vce'"=="robust" & ("`nolog'" =="0"|"`nolog'" ==" "){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(`varlist') robust 
	}
	if "`vce'"=="" & ("`nolog'" =="1"){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(`varlist')  nolog
	}
	if "`vce'"=="robust" & ("`nolog'" =="1"){
	betafit `treatment' [`weight'`exp'] if `touse' & `overlap', mu(`varlist') robust nolog
	}
}


	
mat LL_R2 = e(ll)
estimates store R2

qui lrtest U1 R2

mat Tstat2 = r(chi2)
mat Pval2 = r(p)
mat NofR2 = r(df)

*sca N_LLike = e(N)

/*
mat Mat_LLike = LL_U1\LL_R1\Tstat1\Pval1\NofR1\LL_U1\LL_R2\Tstat2\Pval2\NofR2\N_LLike
mat coln Mat_LLike = "Lrtest"
mat rown Mat_LLike = "Unrestricted" "Restricted X" "T-Stat. X" "p-value X " "Restrictions X" ///
"Unrestricted" "Restricted GPS" "T-Stat.GPS" "p-value GPS" "Restrictions GPS"  "N"
*/

mat empty_mat=.
mat LR_TEST = LL_U1, empty_mat, empty_mat,empty_mat\LL_R1, Tstat1, Pval1, NofR1\LL_R2,Tstat2,Pval2,NofR2
mat coln LR_TEST = "Lrtest" "T-Statistics" "p-value" "Restrictions" 
mat rown LR_TEST = "Unrestricted" "Covariates X" "GPS terms" 

di in ye _newline(2) "********************************************************************"
di                   " 			    Likelihood-ratio tests:	                 	     	  "
di                   " Comparison between the unrestricted model and the restricted models"
di                   "********************************************************************"
mat li LR_TEST
di
di as text "Number of observations = " e(N)

mat drop LR_TEST LL_U1 LL_R1 Tstat1 Pval1 NofR1 LL_U1 LL_R2 Tstat2 Pval2 NofR2 
drop _est*

} /*End of the L_likelihood Test*/



di in ye _newline(2) "***********************************************************"
di                   " End of the assesment of the balancing property of the GPS "
di                   "***********************************************************"

} /*End "if" flag*/

/*END BALANCING TEST*/






di in ye _newline(2) "****************"
di                   " DRF estimation "
di                   "****************"


* BEGIN OF THE DOSE-RESPONSE COMPUTATION


if (("`npoints'"!="" & "`tpoints'"!="") | ("`npoints'"!="" & "`npercentiles'"!="") | ("`tpoints'"!="" & "`npercentiles'"!="")) {	
		di as error "You cannot specify the options tpoints and npoints and npercentiles simultaneously"
		exit 198
	}

tempvar treatment_values  tag id_treat 
tempname tvector mat_aux 

qui gen `treatment_values'  = .
if ("`npoints'"=="" & "`tpoints'"=="" &  "`npercentiles'"=="") {
qui duplicates report  `treatment'
local ntreat_values = r(unique_value) 
sort `treatment'
qui replace `treatment_values' = `treatment' if `treatment'[_n-1]!=`treatment'[_n]

sort  `treatment_values' 
mkmat `treatment_values', matrix(`mat_aux')
matrix def `tvector' = `mat_aux'[1..`ntreat_values',1]
}

if ("`npoints'"!="" & "`tpoints'"=="" & "`npercentiles'"=="") {	
tempvar max_t min_t 
qui sum `treatment'
qui gen `max_t'= r(max)
qui gen `min_t'= r(min)

mkmat `min_t' if _n==1, matrix(`tvector')

local i =1
while `i'<=`npoints'{
qui replace `treatment_values' = `min_t' + `i'*((`max_t'-`min_t')/`npoints') ///
				if `i'< `npoints' | (`i'== `npoints' & `treatment_values' <=`max_t')  
qui replace `treatment_values' = `max_t'  /// 
                if  (`i'== `npoints' & `treatment_values' > `max_t')
mkmat `treatment_values' if _n==1, matrix(`mat_aux')
matrix  `tvector'= `tvector' \ `mat_aux'
local i = 	`i' + 1
}
}

if ("`npoints'"=="" & "`tpoints'"=="" & "`npercentiles'"!="") {	
tempvar Tper Tper_raw

pctile  `Tper_raw' = `treatment', nquantiles(`npercentiles')

qui duplicates report  `Tper_raw'
local nper = r(unique_value) - 1
if `nper' < `npercentiles'-1{
di as text "There are duplicates in the percentiles: They are discarded"
}

sort `Tper_raw'
qui gen `Tper' = `Tper_raw' if `Tper_raw'[_n-1]!=`Tper_raw'[_n]


sort  `Tper' 
mkmat `Tper', matrix(`mat_aux')
matrix def `tvector' = `mat_aux'[1..`nper',1]
}

if ("`npoints'"=="" & "`tpoints'"!="" & "`npercentiles'"=="") {	
matrix def   `tvector' = `tpoints'
}
	
*matrix list `tvector'

global J = rowsof(`tvector')
local j = 1
while `j' <= $J {
tempvar Tpoint_`j'  RTpoint_`j'
qui gen `Tpoint_`j''  = el(`tvector',`j',1) if `touse' & `overlap'==1

if "`family'" == "" | "`family'"== "gaussian" {
qui gen double `RTpoint_`j'' = normalden(`Tpoint_`j'', `etreat', sqrt(theta2))  if `touse' & `overlap'==1
}

if "`family'"== "gamma" {
qui gen double `RTpoint_`j'' = gammaden(theta2, `theta1', 0, `Tpoint_`j'') if `touse' & `overlap'==1
}

if "`family'"== "igaussian" {
igaussianden `RTpoint_`j'' `Tpoint_`j'' `theta1' theta2 if `touse' & `overlap'==1
}

if "`family'"== "beta" {
qui gen double `RTpoint_`j'' =  betaden(`theta1',`theta2',`Tpoint_`j'') if `touse' & `overlap'==1
}


if (`delta' > 0) & ("`method'" != "`iwkernel'") {
tempvar Tpointplus_`j'  RTpointplus_`j'
qui gen `Tpointplus_`j''  = el(`tvector',`j',1) + `delta'

if "`family'" == "" | "`family'"== "gaussian" {
qui gen double `RTpointplus_`j'' = normalden(`Tpointplus_`j'', `etreat', sqrt(theta2))  if `touse' & `overlap'==1
}

if "`family'"== "gamma" {
qui gen double `RTpointplus_`j'' = gammaden(theta2, `theta1', 0, `Tpointplus_`j'') if `touse' & `overlap'==1
}

if "`family'"== "igaussian" {
qui igaussianden `RTpointplus_`j'' `Tpointplus_`j'' `theta1' theta2 if `touse' & `overlap'==1
}

if "`family'"== "beta" {
qui gen double `RTpointplus_`j'' =  betaden(`theta1',`theta2',`Tpointplus_`j'') if `touse' & `overlap'==1
}


} /*End if over derivative*/

local j = `j'+1
}



if "`method'" != "iwkernel" &  "`method'" !="mtpspline" & "`method'" != "radialpspline" {
	di as error "Wrong specification of the estimation method"
	exit 
	}
    
	
if `"`det'"' == `""'{ 
	 local q quietly	 
	}
	
	
if("`method'" == "iwkernel"){
* 6. GET ESTIMATE OF THE DRF USING SEMIPARAMETRIC WEIGHTING APPROACH.
*** a). Get bandwidth. We use Fan and Gijbels rule of thumb. 
tempvar treatment2 treatment3 treatment4 treatment5
forvalues pow =2/5{
qui gen double drf_`treatment'`pow' = `treatment'^`pow'
}

di in ye _newline(2) "IW Kernel estimator"

tempvar d2m2
`q' reg `outcome' `treatment' drf_`treatment'2 drf_`treatment'3 drf_`treatment'4 if `touse' & `overlap'==1

sca sigma2=e(rss)/e(df_r)

* Get estimate of sum of squares of second derivatives    
qui gen `d2m2'= ( (2*_b[drf_`treatment'2]) + (6*_b[drf_`treatment'3]*`treatment') + ///
            (12*_b[drf_`treatment'4]*(drf_`treatment'2)))^2 if `touse' & `overlap'==1
qui sum `d2m2'
scalar sumd2m2=r(sum) 
* Get range of Treatment
qui sum `treatment' if `touse' & `overlap'==1
scalar ranget=r(max)-r(min)
* Get bandwidth. we write it as scalar.  

if ("`bandwidth'" == "") {
	scalar bw= ( (sigma2*ranget) / (2*sqrt(_pi)*sumd2m2) )^(0.2) 
	disp("Bandwidth used in weighting")
	sca list bw
}
else {
	scalar bw = `bandwidth'
	disp("Bandwidth used in weighting")
	sca li bw
}


scalar drop sigma2 sumd2m2 ranget
drop drf_`treatment'2 drf_`treatment'3 drf_`treatment'4 drf_`treatment'5

*** b). Get DRF using inverse by the GPS based on local linear regression.
** Weighting estimator is equivalent to a local linear regression using a "new kernel" which 
** equals the original one now divided by the gps at t: f(t|xi). 

tempname DRF
local JJ = rowsof(`tvector')
local JJ_JJ = 2*`JJ'
matrix `DRF'=J(1,`JJ_JJ',0) // Matrix where to save DRF for weighting

forvalues j = 1/`JJ' {
tempvar weights Ttemp
quietly gen `weights'= exp( -0.5 * (((`Tpoint_`j''-`treatment')/bw))^2 ) / `RTpoint_`j''  //Generate weights. Normal kernel.
quietly gen `Ttemp'= `treatment'-`Tpoint_`j'' if `touse' & `overlap'==1
quietly reg `outcome' `Ttemp' [pweight=`weights'] if `touse' & `overlap'==1 // Run local linear regression
matrix `DRF'[1,`j']=_b[_cons]
}

if(`delta' > 0){
* 9. GET ESTIMATE OF the dDRF using the same approach implemented for the psplines and the radialpsplines. 

tempname DRFplus dDRFw
local JJ = rowsof(`tvector')
matrix `DRFplus'=J(1,`JJ',0)
matrix `dDRFw'=J(1,`JJ',0)   // Matrix where to save the difference between DRF and DRFplus using weighting
forvalues j = 1/`JJ' {

tempvar dweights dTtemp
quietly gen `dweights'= exp( -0.5 * (((`Tpointplus_`j''-`treatment')/bw))^2 ) / `RTpointplus_`j''  //Generate weights. Normal kernel.
quietly gen `dTtemp'= `treatment'-`Tpointplus_`j'' if `touse' & `overlap'==1
quietly reg `outcome' `dTtemp' [pweight=`dweights'] if `touse' & `overlap'==1 // Run local linear regression
matrix `DRFplus'[1,`j']=_b[_cons]


matrix `dDRFw'[1,`j']=`DRFplus'[1,`j']-`DRF'[1,`j']
local jj=`JJ' + `j'
matrix `DRF'[1,`jj'] = el(`dDRFw',1,`j')
}
}/*END IF OVER 	`"`derivative'"' != `""'*/
matrix def DRF_MAT = `DRF'

} /*END IF OVER "`method'" == "iwkernel"*/
	

if("`method'" == "mtpspline" | "`method'" == "radialpspline"){
tempname bb KK U
tempname DRF DRFps DRFps_tplus dDRFps
local JJ = rowsof(`tvector')
matrix `DRFps'=J(1,`JJ',0)  
matrix `DRFps_tplus'=J(1,`JJ',0)  
matrix `dDRFps'=J(1,`JJ',0)  
local JJ_JJ = 2*`JJ'
matrix `DRF'=J(1,`JJ_JJ',0)  

if("`method'" == "mtpspline"){

if("`degree1'" == ""){ 
local degree1 = 1
}

if("`degree2'" == ""){ 
local degree2 = 1
}

if("`nknots1'"==""){
local nknots1 =-1
}

if("`nknots2'"==""){
local nknots2 =-1
}

di in ye _newline(2) "Multivariate penalized spline estimator"
`q' mtpspline `outcome' `treatment' drf_gpscore if `touse' & `overlap'==1, ///
        degree1(`degree1') degree2(`degree2')               ///
        nknots1(`nknots1') nknots2(`nknots2')              ///
		knots1(`knots1')   knots2(`knots2')         	///
		estopts(`estopts') `det'  `additive'  	       

	sca npar = e(k_f)
	matrix def `bb' = e(b)
	matrix def `bb' = `bb'[1,1..npar]
	matrix def `U' = e(Uvars)
	mat def `KK' = e(Knots)
	qui count if `touse' & `overlap'==1
	local nobs = r(N) 
		
	qui count
	local obsoverlap = r(N)-`nobs' + 1
	matrix def `U' = `U'[`obsoverlap'...,1...]
	local nspl = colsof(`U')
	qui svmat `U', names(`U'_)
	*sum U*
	
	local a = colsof(`KK')
	local b = rowsof(`KK')
	local nspline = 2*(`a'+`b') + `a'*`b'
	
	if (`nspline' != `nspl' ){
	di in yellow "Warning message: Some error terms u cannot be estimated" 
	}

	
	* Uvars times splines (see pg. 240 Ruppert et al. (2003))

 
    tempname cost
	qui count if `touse' & `overlap'==1
	local nobs = r(N) 
	mat def `cost' =J(`nobs', 1,1)
		
	*tempvar yhat 	
	forvalues j = 1/`JJ' {
	mksplinevars `Tpoint_`j'' `RTpoint_`j'' if `touse' & `overlap'==1,  knots(`KK') degree1(`degree1') degree2(`degree2') `additive' `standardized'
	tempname Powersj  Powersj_cost Spline
	mat def `Powersj' = e(Powers)
    mat def `Spline' = e(Spline)
	local nspline = colsof(`Spline')
	
	mat def `Spline' = `Spline'[1...,1..`nspl']
	
	
	/*After getting spline and U variables do not introduce the overlap==1  condition. It's already included! */
	
	qui svmat `Spline', names(Z)

	tempvar ZU
	
	forv spl = 1/`nspl'{
	tempvar Zu`spl'
	qui gen `Zu`spl'' = Z`spl'*`U'_`spl'
	drop Z`spl' 
	if `spl'==1{
	qui gen `ZU' = `Zu`spl''
	}
	else{
	qui replace `ZU' = `ZU' + `Zu`spl''
	}
	}

	mat def `Powersj_cost' = `Powersj',`cost'
	tempname xbeta
	mat def `xbeta' = `Powersj_cost'*`bb'' /*`bb'' -> transponse of `bb'  */
	svmat `xbeta', names(xb)
	tempvar yhat
	qui gen `yhat' = xb1 + `ZU'
	drop xb1
	qui sum `yhat'
	matrix `DRFps'[1,`j'] = r(mean) 
	matrix `DRF'[1,`j'] =  el(`DRFps',1,`j')

	
	if(`delta' > 0){
	
	mksplinevars `Tpointplus_`j'' `RTpointplus_`j'' if `touse' & `overlap'==1,  knots(`KK') degree1(`degree1') degree2(`degree2') `additive'  `standardized'
	tempname PowersPlusj  PowersPlusj_cost SplinePlus
	mat def `PowersPlusj' = e(Powers)
    mat def `SplinePlus' = e(Spline)
	mat def `SplinePlus' = `SplinePlus'[1...,1..`nspl']
	
	/*After getting spline and U variables do not introduce the overlap==1  condition. It's already included! */
	
	qui svmat `SplinePlus', names(Zplus)
	
	tempvar ZUplus
	
	forv spl = 1/`nspl'{
	tempvar ZuPlus`spl'
	qui gen `ZuPlus`spl'' = Zplus`spl'*`U'_`spl'
	drop Zplus`spl' 
	if `spl'==1{
	qui gen `ZUplus' = `ZuPlus`spl''
	}
	else{
	qui replace `ZUplus' = `ZUplus' + `ZuPlus`spl''
	}
	}
	
	mat def `PowersPlusj_cost' = `PowersPlusj',`cost'
	tempname xbetaplus 
	mat def `xbetaplus' = `PowersPlusj_cost'*`bb'' /*l'apice in fondo significa trasposto --> non togliere!*/
	svmat `xbetaplus'
	
	
	tempvar yhatplus 
	qui gen `yhatplus' = `xbetaplus' + `ZUplus'
	
	qui sum `yhatplus'
	matrix `DRFps_tplus'[1,`j'] = r(mean) 
	matrix `dDRFps'[1,`j'] = (`DRFps_tplus'[1,`j'] -`DRFps'[1,`j'])
	local jj=`JJ' + `j'
	matrix `DRF'[1,`jj'] = el(`dDRFps',1,`j')
	
	} /*END if(`delta' > 0)*/

	} /*END LOOP FOR j = 1...number of tpoints*/

drop `U'_*
matrix def DRF_MAT = `DRF'	

}/*END IF OVER "`method'" == "mtpspline" */
	
if("`method'" == "radialpspline"){
if("`nknots'"==""){
 local nknots =-1
}
di  in ye _newline(2) "Radial penalized spline estimator"
`q' radialpspline `outcome' `treatment' drf_gpscore if `touse' & `overlap'==1, ///
        nknots(`nknots')   knots(`knots')  `standardized'  estopts(`estopts') `det' 
		
	sca npar = e(k_f)
	matrix def `bb' = e(b)
	matrix def `bb' = `bb'[1,1..npar]
	matrix def `KK' = e(RadKnots)	
	matrix def `U' = e(Uvars)
	
	* if some error terms U cannot be estimated (e.g., because of a lack of observations) 
	* the number of knots saved in `KK' matrix will be always > than the number 
	* of U estimated by using the xtmixed procedure. 
	* In this regard we introduce a warning messagge to inform the final user.
	local nRadK = rowsof(`KK')
	local nU = colsof(`U')
	
	if(`nU' < `nRadK'){
		di in yellow "Warning message: Some error terms u cannot be estimated. The number of knots is `nRadK'. Reduce the number of knots specified."
	}
	
	qui count if `touse' & `overlap'==1
	local nobs = r(N) 

	qui count
	local obsoverlap = r(N)-`nobs' + 1
	matrix def `U' = `U'[`obsoverlap'...,1...]
	
	qui svmat `U', names(`U'_)

	tempname cost
	qui count if `touse' & `overlap'==1
	local nobs = r(N) 
	mat def `cost' =J(`nobs', 1,1)

	
	forvalues j = 1/`JJ' {
	mkradialpspline `Tpoint_`j'' `RTpoint_`j'' if `touse' & `overlap'==1,  knots(`KK') 
	tempname Z0 Xmodel
	mat def `Z0' = e(Z0)
	local nspl = colsof(`Z0')
	
	
	/*After getting spline and U variables do not introduce the overlap==1  condition. It's already included! */
	
	qui svmat `Z0', names(Z)
	tempvar ZU
	
	forv spl = 1/`nspl'{
	tempvar Zu`spl'
	qui gen `Zu`spl'' = Z`spl'*`U'_`spl'
	qui drop Z`spl' 
		if `spl'==1{
		qui gen `ZU' = `Zu`spl''
		}
		else{
		qui replace `ZU' = `ZU' + `Zu`spl''
		}
	}
	
	tempname X
	mkmat `Tpoint_`j'' `RTpoint_`j''  if `touse' & `overlap'==1, matrix(`X')
	mat def `Xmodel' = `X',`cost'
	
	tempname xbeta
	mat def `xbeta' = `Xmodel'*`bb'' /*l'apice in fondo significa trasposto --> non togliere!*/
	
	svmat `xbeta', names(xb)
	tempvar yhat
	qui gen `yhat' = xb1 + `ZU'
	drop xb1
	qui sum `yhat'
	matrix `DRFps'[1,`j'] = r(mean)
	matrix `DRF'[1,`j'] =  el(`DRFps',1,`j')
	if(`delta' > 0){
	mkradialpspline `Tpointplus_`j'' `RTpointplus_`j'' if `touse' & `overlap'==1,  knots(`KK') 
	tempname Z0plus Xmodelplus
	mat def `Z0plus' = e(Z0)
	
	local nspl = colsof(`Z0plus')
	
	
	/*After getting spline and U variables do not introduce the overlap==1  condition. It's already included! */
	
	qui svmat `Z0plus', names(Zplus)
	tempvar ZUplus
	
	forv spl = 1/`nspl'{
	tempvar Zuplus`spl'
	qui gen `Zuplus`spl'' = Zplus`spl'*`U'_`spl'
	drop Zplus`spl' 
	if `spl'==1{
	qui gen `ZUplus' = `Zuplus`spl''
	}
	else{
	qui replace `ZUplus' = `ZUplus' + `Zuplus`spl''
	}
	}
	
	tempname Xplus
	mkmat `Tpointplus_`j'' `RTpointplus_`j''  if `touse' & `overlap'==1, matrix(`Xplus')
	mat def `Xmodelplus' = `Xplus',`cost'
	
	tempname xbetaplus
	mat def `xbetaplus' = `Xmodelplus'*`bb'' /* `bb'' is the transpose of `bb' */
	
	svmat `xbetaplus', names(xbplus)
	tempvar yhatplus
	qui gen `yhatplus' = xbplus1 + `ZUplus'
	drop xbplus1
	qui sum `yhatplus'
	matrix `DRFps_tplus'[1,`j'] = r(mean) 
	matrix `dDRFps'[1,`j'] = (`DRFps_tplus'[1,`j'] -`DRFps'[1,`j']) 
	local jj=`JJ' + `j'
	matrix `DRF'[1,`jj'] = el(`dDRFps',1,`j')
	} /*END if(`delta > 0')*/	
} /*END LOOP FOR j = 1...number of tpoints*/
drop `U'_*
matrix def DRF_MAT = `DRF'


}/*END IF OVER "`method'" == "radialpspline"*/	

}/*END IF OVER "`method'" == "mtpspline" | "`method'" == "radialpspline"*/	


*label var `gpscore' "Estimated generalized propensity score"
drop drf_gpscore drf_gpscore2 drf_gpscore3

if(`"`gps'"' != `""') {
	qui gen gpscore = `gpscore'
}

ereturn post DRF_MAT, esample(`touse') properties("b")
ereturn local properties "b"
*ereturn local cmd "drf"
ereturn local cmd "bootstrap"

end

******************************************************************************************
******************************************************************************************
******************************************************************************************

*Utilities
program define igaussianden, rclass
version 11.0

args igd x mu sigma2

if `sigma2'<=0{
di as err "The scale parameter must be positive"
exit
}

quietly gen `igd' = (1/sqrt(2*_pi*(`x'^3)*`sigma2'))* ///		
			exp(-((`x'- `mu')^2)/(2*(`mu'^2)*`sigma2'*`x')) 
end
		
		
* Additive model: 
* Y -> outcome
* T -> Treatment
* GPS -> gpscore
* Y = beta0 + f(T) + f(GPS) + e_t

* Knots option: 
* 1) knots list
* 2) number of knots -> K_k =  {(k+1)/(K+2)}th sample quantile of the unique x_i's for k = 1, ...K  where K = number of knots
* 3) No number, no list of knots -> K = max(5, min(1/4*number of unique x_i's, 35 )) -> See Ruppert(2002) and then the knots are located using rule in 2).

* Note that we define knots and nknots and Degree for both the treatment and GPS variables: 
* degree1, nknots1 and knots1 refer to the treatment variable
* degree2, nknots2 and knots2 refer to the GPS 

*No use of option starting with "no"





