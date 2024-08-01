*! v 1.0.0 Authors: M. Bia and A. Mattei. Jan, 2014
version 9.2
mata: 

  void mm_radial(
          string scalar varnames1, string scalar varnames2, string scalar selectobs1,
		  string scalar selectobs2
           )

    {

    X  = st_data( .  , tokens(varnames1) ,  selectobs1 )
    DK = st_data( .  , tokens(varnames2) ,  selectobs2 )
	   

    Z = _mm_dist(X,DK) 
	
	Omega = _mm_dist(DK,DK) 
	
	
    st_matrix(st_local("Zmat"), Z)    
    st_matrix(st_local("Omegamat"), Omega)    
        
}

 real matrix _mm_dist(real matrix A,b) 
  {    
   Dist = J(rows(A),rows(b),.)

   for (i=1; i<=rows(A); i++) {
     for (j=1; j<=rows(b); j++) {
       Dist[i,j] = rowsum(abs(A[i,]:-b[j,]):^2)*ln(sqrt(rowsum(abs(A[i,]:-b[j,]):^2)))
     } 
   }
	   return(Dist)
  }



  end

  

exit
	
