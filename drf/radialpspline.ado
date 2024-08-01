* ! v1.0.1, 2014-05-7 Bia: seed option excluded. From outside makes the spacefill replicable anyway.
* ! v1.0.1, 2014-04-25 Bia: seed option inserted.
* ! v1.0.0, 2013-03-30 Bia, Mattei: first production release! (SJ submission)
program define radialpspline, eclass
* radialpspline = Radial Pspline 

version 9.2

syntax varlist(min=3 numeric) [if] [in], [ ///
        NKnots(int -1)               ///
		Knots(namelist)              ///
		ESTOPts(str asis)            ///
		STANDardized		       ///
		det                        ///
		]   
		
	 marksample touse 
	 if `"`det'"' == `""'{ 
	 local q quietly	 
	 }
	 
	 if `nknots'>=0  & `"`knots'"'!="" {
        di as err "nknots() and knots() not both allowed for the treatment variable"
        exit 198
    }
	      
     if  `"`knots'"'!="" {
		matrix def   DKnots = `knots'
		local K   = rowsof(DKnots)
		local nvar = colsof(DKnots)
		if `nvar' != 2 {
        di as err "dim(knots) must be equal to K X 2 (K=Numbers of Knots, 2= T + GPS)"
        exit 198
		}
    }
	
	gettoken depv controls : varlist
	gettoken var1 var2 : controls
	sort `touse' `var1' `var2'
	
    tempvar unique
	
    qui gen     `unique' =   1  if  (`var1'!=`var1'[_n-1] | `var2'!=`var2'[_n-1])  & `touse' 
	qui replace `unique' =   0  if `unique'==. & `touse'  

	
		
	// Determine nknots 
	if `"`knots'"' == ""{
			if `nknots'<0 {
                qui count if `unique'<.
                local U = r(N)
                local K = max(20, min(int(`U'/4), 150))

				} /*end  if `nknots' < 0 */
				
				else local K `nknots'				
				
				if `K'>0 {
				
				if(`"`standardized'"'!=`""'){
				local standardized "`standardized'"
				local standardize "`standardize'"
				}
				else{
				local standardized " " 
				local standardize " " 
				}
				spacefill `var1' `var2' if `touse', ndesign(`K') `standardize' 
				matrix def DKnots = r(Best_Design)
				
				}	/*End loop in `K'>0 */
				
				
	} /* End if over knots*/
	
	*mat li DKnots
	svmat DKnots
	local nvar_knots= colsof(DKnots)
	foreach v of numlist 1/ `nvar_knots'{
	tempvar var_`v'_knots
	local varknots "`varknots' `var_`v'_knots'" 
	qui gen `var_`v'_knots' = DKnots`v'
	drop DKnots`v'
	}
	
	tempvar selectobs1 selectobs2
	qui gen byte `selectobs1' = 0
	qui replace  `selectobs1' = 1 if `touse'
	
	qui gen byte `selectobs2' = 0
	qui replace  `selectobs2' = 1 if _n<= `K'
	
	
    tempname Zmat  Omegamat Zk Omega  
    mata: mm_radial("`controls'", "`varknots'", "`selectobs1'", "`selectobs2'")
	mat `Zk' = `Zmat'
	mat `Omega' = `Omegamat'
	
	foreach j of numlist 1/`K'{
	mat `Omega'[`j',`j']=0
	}
	
	*mat list `Omega'
	
	tempname U w V sqrt_w
	matrix svd `U' `w' `V' = `Omega'
	mat `sqrt_w' = J(1,`K',.)
	foreach j of numlist 1/`K'{
	mat `sqrt_w'[1,`j'] = sqrt(el(`w',1,`j'))
	}
	tempname sqrt_Omega
	matrix `sqrt_Omega' =  `U'*diag(`sqrt_w')*`V''
	
	tempname Z
	matrix `Z' = `Zk'*inv(`sqrt_Omega')
	
	
	svmat `Z'
	local nvar_Z = colsof(`Z')
	foreach v of numlist 1/ `nvar_Z'{
	tempvar Zcontrols`v'
	local Zcontrols "`Zcontrols' `Zcontrols`v''" 
	qui gen `Zcontrols`v'' = `Z'`v'
	drop `Z'`v'
	}
	
	tempvar uvars
	`q'  xtmixed `depv' `controls' if `touse' ///
                 || _all:`Zcontrols', cov(identity) noconstant `estopts'
	
	
	qui predict `uvars'* if `touse', reffects
	tempname Uvars 
	qui mkmat `uvars'*, matrix(`Uvars')
	
	ereturn matrix Uvars = `Uvars'	
	ereturn matrix Sqrt_Omega = `sqrt_Omega'
	ereturn matrix RadKnots = DKnots 

end	


