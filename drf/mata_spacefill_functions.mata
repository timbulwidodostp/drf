*! v 1.0.0 Authors: M. Bia and P. Van Kerm, Jan, 2014
version 9.2 
mata: 
  real matrix spacefill_locations_init(
          string scalar varnames, string scalar n, string scalar touse,
          real scalar transform
          )
  {
    L = st_data( .  , tokens(varnames) ,  touse )
    if (transform==1) {
      mv = meanvariance(L)
      L = ( L :- (mv[1,.]) ) :/ sqrt(diagonal(mv[|2,1 \ .,.|])') 
    }
    if (transform==2) {
      // to do
    }
    return(L) 
  }    

  real matrix spacefill_distance(real matrix A,b) 
  {    
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
    dist = editmissing(dist:^p , 0)           // Euclidian distance raised to power p 
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
          rnew = r :- editmissing(spacefill_distance( L[(candidates \ excluded),.] , xi ):^p, 0)  :+ editmissing(spacefill_distance( L[(candidates \ excluded),.] , yj ):^p , 0)   
          rnew[orderveci[j]] =  rowsum(   spacefill_distance( xi , L[(fixed \ design),.] ):^p  )  :+  editmissing(spacefill_distance( xi , yj ):^p , 0)    
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
          dist = editmissing(dist:^p , 0)    
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
* syntax: _parse_genvar genstuff `genstuff'
* reads the content of a string of the form "genstuff(blah)" or "genstuff(blah, replace)"
* if -replace- found, then drops any existing var called -blah-
* otherwise checks if -blah- does not exist 
* then returns -genstuff- with the -, replace- dropped
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
