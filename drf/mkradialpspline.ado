* ! v1.0.0, 2013-03-30 Bia, Mattei: first production release! (SJ submission)
program define mkradialpspline, eclass
* radialpspline = Radial Pspline 

version 11

syntax varlist(min=2 max=2 numeric) [if] [in],  Knots(namelist)         
	 marksample touse
	 
	 if `"`det'"' == `""'{ 
	 local q quietly	 
	 }
	 
	     
		matrix def   DKnots = `knots'
		local K   = rowsof(DKnots)
		local nvar = colsof(DKnots)
		if `nvar' != 2 {
        di as err "dim(knots) must be equal to K X 2 (K=Numbers of Knots, 2= T + GPS)"
        exit 198
		}
    
	
	gettoken var1 var2 : varlist
	sort `touse' `var1' `var2'
	
	local controls "`controls' `var1' `var2'"
   
	
	svmat DKnots
	local nvar_knots= colsof(DKnots)
	foreach v of numlist 1/ `nvar_knots'{
	tempvar var_`v'_knots
	local varknots "`varknots' `var_`v'_knots'" 
	qui gen `var_`v'_knots' = DKnots`v'
	drop DKnots`v'
	}
	
	tempvar selectobs1 selectobs2
	`q' gen byte `selectobs1' = 0
	`q' replace  `selectobs1' = 1 if `touse'
	
	`q' gen byte `selectobs2' = 0
	`q' replace  `selectobs2' = 1 if _n<= `K'
	
	
    tempname Zmat  Omegamat Z Omega  
    mata: mm_radial("`controls'", "`varknots'", "`selectobs1'", "`selectobs2'")
	mat `Z' = `Zmat'
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
	
	tempname ZZ0
	matrix `ZZ0' = `Z'*inv(`sqrt_Omega')
	
	
	ereturn matrix Z0 =  `ZZ0'
	

end

