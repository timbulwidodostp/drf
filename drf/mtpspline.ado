* ! v1.0.0, 2013-03-30 Bia, Mattei: first production release! (SJ submission)
program define mtpspline, eclass
* mtpspline = Multidimensional Truncated Power Pspline 

version 11

syntax varlist(min=3 max=3 numeric) [if] [in], [ ///
        Degree1(int 1)               ///
		Degree2(int 1)               ///
        nknots1(int -1)              ///
		nknots2(int -1)              ///
        Knots1(numlist sort)         ///
		Knots2(numlist sort)         ///
		ESTOPts(str asis)            ///
		det                        ///
		additive      	             ///
       ]

		 
	 if `"`det'"' == `""'{ 
	 local q quietly	 
	 }
	 
	 if `nknots1'>=0 & `"`knots1'"'!="" {
        di as err "nknots() and knots() not both allowed for the treatment variable"
        exit 198
    }
	 if `nknots2'>=0 & `"`knots2'"'!="" {
        di as err "nknots() and knots() not both allowed for the GPS"
        exit 198
    }
    if `degree1'<0 {
        di as err "degree()<0 not allowed for the treatment variable"
        exit 198
    }
	 if `degree2'<0 {
        di as err "degree()<0 not allowed for the GPS"
        exit 198
    }


		
	marksample touse
    gettoken depv controls : varlist
	gettoken var1 var2 : controls
    sort `touse' `var1' 
    tempvar unique1 
    qui gen `unique1' = `var1' if `var1'[_n]!=`var1'[_n-1] & `touse'
	sort `touse' `var2'
	tempvar unique2
    qui gen `unique2' = `var2' if `var2'!=`var2'[_n-1] & `touse'

		
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

	// - determine nknots for the treatment
        if `"`knots1'"'=="" {
            if `nknots1'<0 {
                qui count if `unique1'<. &  `touse'
                local U1 = r(N)
                local K1 = max(5, min(int(`U1'/4), 35))
				/*
                if `K1'<2 {
                    di as err "not enough distinct values" //=> linear only in the treatment 
                    exit 198
					}
				*/
				} /*end  if `nknots1' < 0 */
		 	else local K1 `nknots1'
			*di `K1'
			if `K1'>0 {
                forv i = 1/`K1' {
                    local pct1 "`pct1' `=(`i'+1)/(`K1'+2)*100'"
                    // Ruppert et al. (2003; eq. 5.8 pg 126) give rule (i+1)/(k+2); 
					}
                _pctile `unique1' if `touse', p(`pct1')
				forv i = 1/`K1' {
					scalar r1`i' = r(r`i')
					*di r1`i'
				}
				}
				
			} /*end of if knots1 == "" */
			else { 
				qui sum `var1' if `touse'
				local i 0
				foreach k1 of local knots1 {
					if (`k1'<r(min)) | (`k1'>r(max)) {
						di as err "(`k1' not in treatment range. knot will be dropped)"
						continue
					}
				local ++i
				scalar r1`i' = `k1'
				*di r1`i'
			} /* end foreach k1 */
		    
			if `i'==0 {
				di as err "no valid knots"
				exit 198
			}
			local K1 `i'
		   } /* end else -> knots1 != "" */	
	
		
			
			// - determine nknots for the GPS
			if `"`knots2'"'=="" {
            if `nknots2'<0 {
                qui count if `unique2'<. & `touse'
                local U2 = r(N)
                local K2 = max(5, min(int(`U2'/4), 35))
                /*
				if `K2'<2 {
                    di as err "not enough distinct values" //=> linear only in the GPS
                    exit 198
					}
				*/
				} /*end  if `nknots2' < 0 */
				
				else local K2 `nknots2'
				di `K2'
				if `K2'>0 {
                forv i = 1/`K2' {
                    local pct2 "`pct2' `=(`i'+1)/(`K2'+2)*100'"
                    // Ruppert et al. (2003; eq. 5.8 pg 126) give rule (i+1)/(k+2); 
					}
				  _pctile `unique2' if `touse', p(`pct2')
				  forv i = 1/`K2' {
						scalar r2`i' = r(r`i')
						*di r2`i'
						}
				  }
				
			} /*end of if knots2 == "" */
			else { 
				qui sum `var2' if `touse'
				local i 0
				foreach k2 of local knots2 {
					if (`k2'<r(min)) | (`k2'>r(max)) {
						di as err "(`k2' not in GPS range. knot will be dropped)"
						continue
					}
					local ++i
					scalar r2`i' = `k2'
					*di r2`i'
				} /* end foreach k2 */
		    
			    if `i'==0 {
					di as err "no valid knots"
					exit 198
				}
				local K2 `i'
	
			} /* end else -> knots2 != "" */	
			
			

			// - generate splines for treatment variable
			if `K1'>0 {
				tempname knotmat1
				mat `knotmat1' = J(1,`K1', 0)
				local vartype1 = cond(`degree1'==0, "byte", "")
				forv i = 1/`K1' {
					tempvar splineT`i' 
					local splineT "`splineT' `splineT`i''"
					local k1 = r1`i'
					mat `knotmat1'[1,`i'] = `k1'
					local knottick1 "`knottick1' `k1'"
					qui gen `vartype1' `splineT`i'' = cond(`var1' - `k1' > 0, (`var1' - `k1')^`degree1', 0) if `touse'
					if `degree1'==0 {
						qui replace `splineT`i'' = 1 if `var1'==`k1' & `touse' // include lower boundary
					} /* end of ->  if degree1 == 0 */
					
					
					if `"`additive'"' == `""'{
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
			
						
			// - generate splines for GPS 
			if `K2'>0 {
				tempname knotmat2
				mat `knotmat2' = J(1,`K2', 0)
				local vartype2 = cond(`degree2'==0, "byte", "")
				forv i = 1/`K2' {
					tempvar splineGPS`i' 
					local splineGPS "`splineGPS' `splineGPS`i''"
					local k2 = r2`i'
					mat `knotmat2'[1,`i'] = `k2'
					local knottick2 "`knottick2' `k2'"
					qui gen `vartype2' `splineGPS`i'' = cond(`var2' - `k2' > 0, (`var2' - `k2')^`degree2', 0) if `touse'
					if `degree2'==0 {
						qui replace `splineGPS`i'' = 1 if `var2'==`k2' & `touse' // include lower boundary
					} /* end of ->  if degree2 == 0 */
					
					
				if `"`additive'"' == `""'{
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
			
				tempvar uvars
			 if `"`additive'"' != `""'{	
			`q'  xtmixed `depv'  `powers1' `powers2' if `touse' ///
                 || _all:`splineT', cov(identity) noconstant || _all: `splineGPS', cov(identity) noconstant `estopts'
				 qui predict `uvars'* if `touse', reffects 
            }		
			else{
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

			
			    `q'  xtmixed `depv' `powers1' `powers2' `powers1_powers2'  if `touse' ///
                 || _all:`splineT', cov(identity) noconstant || _all: `splineGPS', cov(identity) noconstant ///
				 || _all:`splineGPS_T', cov(identity) noconstant || _all: `splineT_GPS', cov(identity) noconstant ///
				 || _all:`splineGPS_splineT', cov(identity)  noconstant `estopts'
				
				 if e(k_res) == 0 {
				 di as err "Standard-error calculation failed using tensor products, switch to additive model "
				 }

			     /*  `q'  xtmixed `depv' `powers1' `powers2' `powers1_powers2'  if `touse' ///
                 || _all:`splineT' `splineGPS' `splineGPS_T' `splineT_GPS' `splineGPS_splineT', cov(identity)  noconstant `estopts' 
				 */
				  qui predict `uvars'* if `touse', reffects
				 }
		
	
		
		tempname knotmat
		local c1 = colsof(`knotmat1')
		local c2 = colsof(`knotmat2')
		if(`c1'>= `c2'){
		mat `knotmat' = J(2,`c1',.)
		}
		else{
		mat `knotmat' = J(2,`c2',.)
		}
		
		forv j= 1/`c1'{
		mat `knotmat'[1, `j'] = el(`knotmat1',1, `j')
		}
		
		forv j= 1/`c2'{
		mat `knotmat'[2, `j'] = el(`knotmat2',1,`j')
		}
		
		tempname Uvars
		mkmat `uvars'*, matrix(`Uvars')
		
		ereturn matrix Knots = `knotmat'
		ereturn matrix Uvars = `Uvars'	
	
			
		
end
