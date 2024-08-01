* v1.1.2, 2014-04-09, PVK (SJ resubmission)
* v1.1.1, 2014-04-02, PVK 
*         - more and relabelled standardization options (dependency on moremata, include wvar in standadirzation)
*         - bug fix on 1/distance to power p calculation with distnace==0
* v1.1.0, 2014-03-01, Bia IQR option added
* v1.0.0, 2013-07-30, first production release! (SJ submission)
* v0.2.0, 2013-07-23, noverbose
* v0.1.1, 2013-04-04, PVK complete redesign
*         - fixed and exclude option
*         - weighted calculation
* v0.1.0, 2013-04-03, PVK complete redesign
*         - speeding up formula
*         - improved standardize
*         - improved saving
* v0.0.3, 2011-11-20, Michela Bia, changing the generate option into a dummy variale and checking the program
* v0.0.3, 2011-11-15, Michela Bia, fixing the 'standardized variables' option and checking the program
* v0.0.2, 2011-09-26, Michela Bia, adding the 'standardized variables' option and checking the program
* v0.0.1, 2008-12-09, Philippe Van Kerm, Stata implementation of the space-filling algorithm described
* in J. Andrew Royle & Douglas Nychka, 1998, 'An Algorithm for the Construction of Spatial Coverage Designs with Implementation
* in Splus', Computers and Geosciences, vol. 24, Number 5, pp. 479-488.

          
program def spacefill , rclass sortpreserve
    version 9.2
    syntax varlist(min=1 numeric) [if] [in] [aw iw fw] ,   ///
        [                         /// 
        NDesign(real -1)           /// number of design points (overridden if design0 is specified)
        DESIGN0(varlist numeric)  /// set of initial designs (identified by varname>0)	      
        FIXed(varname numeric)    /// identify locations included in all designs (identified by varname>0)	      
        EXCLUDE(varname numeric)  /// identify locations excluded from all designs (identified by varname>0)	      
        p(real -5)  q(real 1)     /// parameters 
        NRuns(real 5)             /// number of runs (overridden if design0 is specified)
        NNPoints(real -1)         /// number of nearest neighbours 
        NNfrac(real -1)           /// fraction of data in nn 
	      STANDardize		            /// standardizes variables to mean zero and unit SD
	      STANDARDIZE2	            /// standardizes variables to mean zero and unit SD (as estimated by 0.7413*IQR)
	      STANDARDIZE3	            /// standardizes variables to median zero and unit SD (as estimated by 0.7413*IQR)
	      SPHERicize		            /// standardizes variable matrix to mean zero, unit diagonal covariance matrix
	      RANKS     		            /// transform variables to fractional ranks
        GENerate(string)          /// if one word: then create new variables using string as stub
			                            ///   otherwise create new variables with selected locations 
        GENMARKer(string)         /// generate new dummy variable marking selected locations 
        noVERbose                 /// deletes any output 
        ]  
   
    // --- parse options: check validity
    if (`nruns'<=0) {
      di as error "nruns() invalid: number of runs must be positive"
      exit 198
    }
    if (`nnfrac'>1) { 
      di as error "nnfrac() invalid: fraction of nearest neighbours must in [0,1]"
      exit 198
    }
    if ((`ndesign'>0)+("`design0'"!="")==2) {
      di as error "ndesign() and design0() are mutually exclusive"
      exit 198
    }  
    if (`ndesign'<=0) {
      loc ndesign 4 
    }
    
    if (("`standardize'"!="")+("`standardize2'"!="")+("`standardize3'"!="")+("`sphericize'"!="")+("`ranks'"!="")>1) {
      di as error "Input standardization options are mutually exclusive"
      exit 198
    }  
    if (("`standardize2'"!="")+("`standardize3'"!="")+("`ranks'"!="")>0) {
      capt findfile lmoremata.mlib
      if _rc {
        di as error "Ben Jann's -moremata- is required for options standardize2, standardize3 and ranks; type {stata ssc install moremata}"
        exit 499
      }
    }  
    
    
    // --- parse options : set defaults
    if ((`nnfrac'<0)+(`nnpoints'<0)==2) loc nnfrac 0.50
    
    if ("`generate'"!="") {
      gettoken genword genopts : generate , parse(",")
      if (`:word count `genword'' != `:word count `varlist'') {
        if (`:word count `genword''==1) {
          loc generate 
          forv i=1/`: word count `varlist'' {
            loc generate "`generate' `genword'`i'" 
          }
          loc generate "`generate' `genopts'" 
        }
        else {
          di as error "Number of names in generate() does not match number of variables in varlist"
          exit 198
        }    
      }  
      _parse_genvar generate `generate'
    }  
	 
    if ("`genmarker'"!="") {
      _parse_genvar genmarker `genmarker'
    }  
    
    
    // --- set sample marker
    marksample touse , zeroweight

    // --- handle weights
    tempvar wvar 
    if ("`weight'`exp'"=="") {
      qui gen byte `wvar' = 1 if `touse'
    }
    else {
      qui gen byte `wvar' `exp'  if `touse'
    }
    
    // --- : check designs
    if ("`design0'"!="") {
      loc ndesign
      loc nruns 0
      foreach var of varlist `design0' {
        loc ++nruns
        qui count if `var'>0 & !mi(`var') & `touse'
        if ("`ndesign0'"!="") {
          if (r(N)!=`ndesign') {
            di as error "Number of design points varies across initial designs"
            exit 198
          }
        }
        loc ndesign = r(N)
      }
      if (`ndesign'==0) {
        di as error "No design points identified in initial design(s)"
        exit 198
      }
    }
      
    // --- check size and set nearest neighbour counts
    qui count if `touse'
    loc N = r(N)
    if (`N'<=`ndesign') {
      di as error "Number of design points must be smaller than observations (i.e. candidate points)"
      exit 198
    }
    if (`nnfrac'>0) {
      loc nnpoints = ceil(`nnfrac'*(`N'-`ndesign'))
    }
    loc nnpoints = min((`N'-`ndesign'),`nnpoints')
       
    
    // --- check fixed 
    tempvar tmpfixed
    if ("`fixed'"!="") {
      qui count if (`fixed'>0 & !mi(`fixed'))  & `touse'
      loc nfixed = r(N)
      di as text "`nfixed' fixed design points (`fixed'>0)"
      if ( r(N)>=`ndesign') {
        di as error "Number of fixed design points exceeds total design size (ndesign())"
        exit 198
      }
      if ("`design0'"!="") {
        foreach var of varlist `design0' {
          qui count if (`fixed'>0 & !mi(`fixed')) & (`var'==0 | mi(`var'))  & `touse'
          if (r(N)>0) {
            di as error "Fixed locations in `fixed' are missing from initial design `var'"
            exit 198
          }
        }
      }
      qui gen byte `tmpfixed' = `fixed' if `touse'
      qui replace `tmpfixed' = 0 if mi(`tmpfixed') & `touse'  // set missings to zero
    }
    else {
      qui gen byte `tmpfixed' = 0 if `touse'
      loc nfixed 0
    }

    
    // --- check exclude 
    tempvar tmpexcluded
    if ("`exclude'"!="") {
      qui count if (`exclude'>0 & !mi(`exclude'))  & `touse'
      loc nexcluded = r(N)
      di as text "`nexcluded' points excluded from designs (`exclude'>0)"
      qui count if `touse'
      if ( `ndesign' >= (r(N)-`nexcluded')) {
        di as error "Too many excluded points for the requested design size"
        exit 198
      }
      if ("`design0'"!="") {
        foreach var of varlist `design0' {
          qui count if (`exclude'>0 & !mi(`exclude')) & (`var'>0 & !mi(`var'))  & `touse'
          if (r(N)>0) {
            di as error "Excluded locations found in initial design set `var'"
            exit 198
          }
        }
      }
      if ("`fixed'"!="") {
        qui count if (`exclude'>0 & !mi(`exclude')) & (`fixed'>0 & !mi(`var'))  & `touse'
        if (r(N)>0) {
          di as error "Excluded locations also found in fixed design set"
            exit 198
        }
      }
      qui gen byte `tmpexcluded' = `exclude' if `touse'
      qui replace `tmpexcluded' = 0 if mi(`tmpexcluded') & `touse'  // set missings to zero
    }
    else {
      qui gen byte `tmpexcluded' = 0 if `touse'
      loc nexcluded 0
    }
        
    
    // --- initialize locations matrix (with transformations)
    tempvar n
    qui gen `n' = _n if `touse'
    loc transform = ("`standardize'"!="") + 2*("`standardize2'"!="") + 3*("`standardize3'"!="") + 4*("`sphericize'"!="") + 5*("`ranks'"!="") 
    tempname L
    mata: `L' = spacefill_locations_init("`varlist'", "`n'","`touse'", "`wvar'", `transform')    // read the locations from data and apply transformation
  
    
    // --- runs 
    tempname Cpq
    
    loc maxCpq = 10e13
    
	if ("`verbose'"=="noverbose")  loc qui qui
	
    forv i=1/`nruns' {
      `qui' di as text "Run " as res "`i'"  _col(8) _c
      tempvar design`i' 
      qui gen byte `design`i'' = 0   if `touse'    // in the future: set to 1 for fixed points
      if ("`design0'"!="") {
        `qui' mata: spacefill_run( `L' , `ndesign' , `nnpoints', `p' , `q' , "`design`i''",  "`tmpfixed'",  "`tmpexcluded'", "`wvar'", "`touse'", "`:word `i' of `design0''" )  
      }  
      else  {
        `qui' mata: spacefill_run( `L' , `ndesign' , `nnpoints', `p' , `q' , "`design`i''",  "`tmpfixed'",  "`tmpexcluded'", "`wvar'", "`touse'" )  
      } 
      * here `design`i'' contains 0 and 1 with ones in the selected design 
      *  `Cpq'    contains the distance statistic with final design
      tempname Cpq`i' 
      sca `Cpq`i'' = `Cpq'
  	  if (scalar(`Cpq`i'') < `maxCpq' ) {
        loc maxi `i'
        loc maxCpq = `Cpq`i''
      } 
  
      `qui' di _col(45) as text "(Cpq = " as res %14.2f `Cpq`i'' as text ")"
    } 
      
    tempname Dmat
    qui mkmat `varlist' if `design`maxi''==1 , matrix(`Dmat')

    // --- create new variables
    if ("`genmarker'"!="") {
      qui gen byte `genmarker' = `design`maxi''
    }
           
    if ("`generate'"!="") {
      qui svmat `Dmat' 
		  loc i 0
      foreach var of local generate {
        loc ++i
        rename `Dmat'`i' `var'  		  
		  }
		}  
		
    // --- return values    
    return local varlist "`varlist'"
    if ("`generate'"!="") {
      return local generate "`generate'"
    }  
    if ("`genmarker'"!="") {
      return local genmarker "`genmarker'"
    }  
    if ("`fixed'"!="") {
      return local fixed "`fixed'"
    }  
    if ("`exclude'"!="") {
      return local exclude "`exclude'"
    }  
    return scalar N = `N'
    return scalar ndesign = `ndesign'
    return scalar nfixed = `nfixed'
    return scalar nexcluded = `nexcluded'
    return scalar Cpq = `maxCpq'
    return scalar nn = `nnpoints'
    return scalar p = `p'
    return scalar q = `q'
    
    return matrix Best_Design = `Dmat'
    
        
end    

version 9.2 
mata: 
  real matrix spacefill_locations_init(
          string scalar varnames, string scalar n, string scalar touse,
          string scalar wvar,           
          real scalar transform
          )
  {
    real matrix L, w, mw, m, v
    
    L = st_data( .  , tokens(varnames) ,  touse )
    w = st_data( .  , tokens(wvar) ,  touse )
    
    if (transform==1) {
      mv = meanvariance(L,w)
      L = ( L :- (mv[1,.]) ) :/ sqrt(diagonal(mv[|2,1 \ .,.|])') 
    }
    if (transform==2) {
      mv = mean(L ,w)
      L = ( L :- (mv[1,.]) ) :/ (mm_iqrange(L)*(0.7413))
    }
    if (transform==3) {
      L = ( L :- mm_median(L,w) ) :/ (mm_iqrange(L,w)*(0.7413))
    }
    if (transform==4) {
      m = mean(L,w)
      v = variance(L,w)
      L =  ( L :- m ) * cholesky(invsym(v))  
      //Alternative: L = (d:-m) * (matpowersym(v,-0.5))
    }
    if (transform==5) {
      L = mm_ecdf(L,w)
    }
    
    return(L) 
  }    

  real matrix spacefill_distance(real matrix A,b) 
  {    
   real matrix euclid
   real scalar i, j
   euclid = J(rows(A),rows(b),.)
   for (i=1; i<=rows(A); i++) {
     for (j=1; j<=rows(b); j++) {
       euclid[i,j] = sqrt(rowsum(abs(A[i,]:-b[j,]):^2))
     } 
   }
   return(euclid)
  }
  
  
  void spacefill_run(
          real matrix L , real scalar n , real scalar m, 
          real scalar p, real scalar q, 
          string scalar designvar, 
          string scalar fixedvar, 
          string scalar excludevar, 
          string scalar wvar,           
          string scalar touse , 
          | string scalar inidesignvar
          )
  {
    real matrix W, indx, fixed, excluded, oinit, design, candidates, dist, w, r, xi, yj, rnew, wnew, tmp
    real scalar nfixed, nexcluded, N, Cpq, anyswapped, toswap, maxdCpq, Cpqnew
    
    // initiqlize results vector
    st_view(dvec = ., . , tokens(designvar) , touse )

    // read in weight vectr
    W = st_data( .  , tokens(wvar) ,  touse )
        
    // handle fixed design points
    st_view(fixcol = ., . , tokens(fixedvar) , touse )
    indx = (1 :: rows(fixcol))
    fixed = select(indx , fixcol:>0) 
    nfixed = rows(fixed) 

    // handle excluded points
    st_view(excol = ., . , tokens(excludevar) , touse )
    indx = (1 :: rows(excol))
    excluded = select(indx , excol:>0) 
    nexcluded = rows(excluded) 
        
    N = rows(L) - n
    
    if (args()<11) {
      // randomize starting index sets into design and candidates
      oinit = revorder( order((fixcol,excol,runiform((N+n),1)),(1,2,3)) ) 
      design = oinit[((nfixed+nexcluded+1) :: (n+nexcluded))]
      candidates = oinit[((n+nexcluded+1) :: (N+n))]
    }
    else {
      // take starting index sets from inidesign ==0 and >0
      st_view( init = . , . , tokens(inidesignvar) , touse )
      oinit = revorder( order((init,fixcol,excol),(1,2,3)) ) 
      design = oinit[((nfixed+1) :: (n-nfixed))]
      //design = oinit[(1 :: (n-nfixed))]
      candidates = oinit[((n+nexcluded+1) :: (N+n))]
    }
    
    dist = spacefill_distance( L[(candidates \ excluded),.] , L[(fixed \ design),.] ) // Euclidian distance between each candidate point to each design point (Nxn)
    dist = editmissing(dist:^p , 10e20)           // Euclidian distance raised to power p 
    w = W[(candidates \ excluded)] 
    r = rowsum(dist)         // sum over all design points for each candidate point
    Cpq = (colsum(w:*(r:^(q/p))))^(1/q)  // aggregation over all candidate points  

    anyswapped = 1
    while (anyswapped==1) {
      anyswapped = 0
      printf(".")
      displayflush()
      for (i=1; i<=(n-nfixed); i++) {   // foreach column
        toswap = 0
        maxdCpq = 0
        xi = L[design[i],.]  // ith row in design index vector and pick it from L
        orderveci = order(dist[(1::(N-nexcluded)),.],i) // permutation vector of dis matrix in order if ith col (the last obs in dist are excluded)
        for (j=1; j<=min((m,(N-nexcluded))); j++) {
          yj = L[candidates[orderveci[j]],.]  // jth line from ordering vector, so take the resulting indexc from candiate set
        
          // create test distance:
          rnew = r :- editmissing(spacefill_distance( L[(candidates \ excluded),.] , xi ):^p, 10e20)  :+ editmissing(spacefill_distance( L[(candidates \ excluded),.] , yj ):^p , 10e20)   
          rnew[orderveci[j]] =  rowsum(   spacefill_distance( xi , L[(fixed \ design),.] ):^p  )  :+  editmissing(spacefill_distance( xi , yj ):^p , 10e20)    
          wnew = W[(candidates \ excluded)] 
          wnew[orderveci[j]] = W[design[i]]
          Cpqnew = (colsum(wnew:*(rnew:^(q/p))))^(1/q)   
          // check reduction:
          if ( (Cpq - Cpqnew)>maxdCpq ) {
            maxdCpq = Cpq - Cpqnew
            toswap = orderveci[j]
          }  
        }
	      if (toswap>0) {
          // printf("Swap %f and %f", i, toswap)
          anyswapped = 1
          
          // do swap:
          tmp = candidates[toswap]
          candidates[toswap] = design[i]
          design[i] = tmp
          // reset distance and row vector (needed because sdwapped out point can be swapped in again)
          dist = spacefill_distance( L[(candidates \ excluded),.] , L[(fixed \ design),.] ) 
          dist = editmissing(dist:^p , 10e20)    
          w = W[(candidates \ excluded)] 
          r = rowsum(dist)  
          Cpq = (colsum(w:*(r:^(q/p))))^(1/q)
        }

       // printf("%8.5g," , Cpq)   
      }
    }
        
    // finally fill with ones the design variable and save smallest distance
    st_numscalar(st_local("Cpq"), Cpq)     
    dvec[(fixed \ design),1] =  J(n,1,1) 
    
  }    // end spacefill_run
        

      
    
end
    




*! v1.0.0, 2008-12-03, Philippe Van Kerm, Parse options of the form -genstuff(varname , replace)-
program define _parse_genvar
    version 9.2
    gettoken locname 0 : 0
    gettoken genvar genopts : 0 , parse(",")
    gettoken  comma genopts : genopts  , parse(",")
    if (trim("`genopts'")=="")  confirm new variable `genvar'
    else {
      if (trim("`genopts'")=="replace") cap drop `genvar'
      else {
        di as error "`0'   invalid"
        exit 198
      }
    }
    c_local `locname' "`genvar'"
end

exit



  
