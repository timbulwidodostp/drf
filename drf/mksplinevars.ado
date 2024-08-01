* ! v1.0.0, 2013-03-30 Bia, Mattei: first production release! (SJ submission)
program mksplinevars, eclass
version 11.0

syntax varlist(min=2 max=2 numeric)  [if] [in], knots(namelist) [ ///
		Degree1(int 1) Degree2(int 1) additive ///
		]
		
		
    if `degree1'<0 {
        di as err "degree()<0 not allowed for the treatment variable"
        exit 198
    }
	 if `degree2'<0 {
        di as err "degree()<0 not allowed for the GPS"
        exit 198
    }

	
		
	marksample touse
    gettoken var1 var2 : varlist
  
  		
	// degree1(): generate T^2, T^3, ...
    forv i = 1/`degree1' {
        if `i'==1 {
            local powers1 "`var1'"
            continue
        }
        tempvar power1`i'
        qui gen `power1`i'' = `var1'^`i' if `touse'
        local powers1 "`powers1' `power1`i''"
    }
	local pow1: word count `powers1'
	
	// degree2(): generate GPS^2, GPS^3, ...
    forv i = 1/`degree2' {
        if `i'==1 {
            local powers2 "`var2'"
            continue
        }
        tempvar power2`i'
        qui gen `power2`i'' = `var2'^`i' if `touse'
        local powers2 "`powers2' `power2`i''"
    }

	local pow2: word count `powers2'
	
	local c = colsof(`knots')
	sca nc1 = 0
	sca nc2 = 0
	forv j = 1/`c'{
	if el(`knots', 1,`j') <. {
	sca nc1 = nc1 + 1 
	}
	if el(`knots', 2,`j') <.{
	sca nc2 = nc2 + 1 
	}
	}
		

	local K1 = nc1
		// - generate splines for treatment variable
		if `K1'>0 {
			local vartype1 = cond(`degree1'==0, "byte", "")
			forv i = 1/`K1' {
				tempvar splineT`i' 
				local splineT "`splineT' `splineT`i''"
				local k1 = el(`knots',1,`i')
				qui gen `vartype1' `splineT`i'' = cond(`var1' - `k1' > 0, (`var1' - `k1')^`degree1', 0) if `touse'
				if `degree1'==0 {
					qui replace `splineT`i'' = 1 if `var1'==`k1' & `touse' // include lower boundary
				} /* end of ->  if degree1 == 0 */
				
					
					if  `"`additive'"' == `""' {
					tempname SplineT_GPS
					local p2=0
					foreach x of varlist `powers2'{
					local p2 = `p2' +1 
					tempvar splineT_GPS`i'`p2'
						local splineT_GPS "`splineT_GPS' `splineT_GPS`i'`p2''"
						qui gen `splineT_GPS`i'`p2'' = `splineT`i''*`x' if `touse'
					}
						
					}
			} /* end of loop over K1 -> numero nodi di T */
							
		} /* end of -> if `K1' > 0*/

		local K2 = nc2
	    tempname SplineGPS
			// - generate splines for GPS 
			if `K2'>0 {	
				local vartype2 = cond(`degree2'==0, "byte", "")
				forv i = 1/`K2' {
					tempvar splineGPS`i' 
					local splineGPS "`splineGPS' `splineGPS`i''"
					local k2 = el(`knots',2,`i')
					qui gen `vartype2' `splineGPS`i'' = cond(`var2' - `k2' > 0, (`var2' - `k2')^`degree2', 0) if `touse'
					if `degree2'==0 {
						qui replace `splineGPS`i'' = 1 if `var2'==`k2' & `touse' // include lower boundary
					} /* end of ->  if degree2 == 0 */
					
					
					if  `"`additive'"' == `""' {		
					    tempname SplineGPS_T
						local p1=0
						foreach x of varlist `powers1'{
						local p1 = `p1' +1 
						tempvar splineGPS_T`i'`p1'
						local splineGPS_T "`splineGPS_T' `splineGPS_T`i'`p1''"
						qui gen `splineGPS_T`i'`p1'' = `splineGPS`i''*`x' if `touse'
							}
						
					}
				} /* end of loop over K2 -> numero nodi di GPS */
				
				*mat li `knotmat2' 
				
			} /* end of -> if `K2' > 0*/		

				if  `"`additive'"' == `""'{	
				local p1=0
				local p2=0
					foreach treat of varlist `splineT'{
					local p1 = `p1' + 1
					foreach gps  of varlist `splineGPS'{
					local p2 = `p2' + 1
					tempvar splineGPS_splineT`p1'`p2'
					local splineGPS_splineT "`splineGPS_splineT' `splineGPS_splineT`p1'`p2''"
					qui gen `splineGPS_splineT`p1'`p2'' = `treat'*`gps' if `touse'
					}
					}
					
				local p1=0
				local p2=0				
					foreach treat of varlist `powers1'{
					local p1 = `p1' + 1
					foreach gps   of varlist `powers2'{
					local p2 = `p2' + 1
					tempvar powers1_powers2`p1'`p2'
					local powers1_powers2 "`powers1_powers2' `powers1_powers2`p1'`p2''"
					qui gen `powers1_powers2`p1'`p2'' = `treat'*`gps' if `touse'
					}
					}
			
			}
	
	/*CHECKS --> CHECKS DONE: IT WORKS
	di "POW 1"
	sum `powers1'
	di "POW 2"
	sum `powers2'
	 if  `"`additive'"' == `""'{	
	di "POW 1 & POW 2"
	sum `powers1_powers2' 
	}
	*/


	/*
 tempname Spline
 if  `"`additive'"' == `""'{	
 mkmat `splineT' `splineT_GPS' `splineGPS' `splineGPS_T'  `splineGPS_splineT' if `touse', matrix(`Spline')
 }
 else{
 mkmat `splineT' `splineGPS' if `touse', matrix(`Spline')
 }

 tempname Powers
 if  `"`additive'"' == `""'{	
 mkmat `powers1' `powers2' `powers1_powers2' if `touse', matrix(`Powers')
 }
 else{
 mkmat  `powers1' `powers2' if `touse', matrix(`Powers')
 }

			
ereturn matrix Spline = `Spline'	
ereturn matrix Powers = `Powers'	
*/


/*
 if  `"`additive'"' == `""'{	
	global Spline "`splineT' `splineT_GPS' `splineGPS' `splineGPS_T'  `splineGPS_splineT'"
	}
else{
	global Spline "`splineT'  `splineGPS'"
}

	foreach nome of global Spline{
	mkmat `nome' if `touse', matrix(`nome')
	ereturn matrix `nome' = `nome'
	}
*/


 tempname Spline
 if  `"`additive'"' == `""'{	
 mkmat `splineT' `splineT_GPS' `splineGPS' `splineGPS_T'  `splineGPS_splineT' if `touse', matrix(`Spline')
 }
 else{
 mkmat `splineT' `splineGPS' if `touse', matrix(`Spline')
 }


 tempname Powers
 if  `"`additive'"' == `""'{	
 mkmat `powers1' `powers2' `powers1_powers2' if `touse', matrix(`Powers')
 }
 else{
 mkmat  `powers1' `powers2' if `touse', matrix(`Powers')
 }

			
	
ereturn matrix Spline = `Spline'		
ereturn matrix Powers = `Powers'	
end
