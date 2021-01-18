*!version 1.0 Helmut Farbmacher (January 2021)
* 
***********************************************

prog heckman_lasso, eclass
version 16

syntax varlist, seldep(varlist) [twostep notpen(varlist) tuning(string) robust cluster(passthru) verbose lassofirst simulation]		

gettoken lhs rhs : varlist
loc selind: list notpen | rhs

// checks and defaults
if `:word count `seldep''!=1 {
	dis "{red}Only one variable is allowed in the seldep-option. You specified: {res}`seldep' {err}"
	exit 103
}
if "`tuning'"=="" local tuning "cv"
else if !inlist("`tuning'", "ftest","cv","cvse")  {
    dis "{res}`tuning' {err}is not an implemented tuning"
    exit 198
}

// display more results
if "`verbose'"=="" {
	local qui qui
}

// use Lasso estimates to derive weights for adaptive Lasso
if "`lassofirst'"=="" {
	local unpenalized unpenalized
}

// correcting for sample selection
tempvar xb lambda
qui probit `seldep' `selind'
qui predict `xb', xb
qui gen `lambda'=normalden(`xb')/normal(`xb')

// partialling out
foreach var of varlist `selind' {
	tempvar `var'_raw
	qui gen ``var'_raw'=`var'
	tempvar `var'_res
	qui reg `var' `lambda' if `seldep'
	qui predict ``var'_res', resid
	qui replace `var'=``var'_res'
}

dis ""
dis "{txt}--------------------------------------------"
dis "{txt}Heckman with unknown exclusion restrictions "
dis "{txt}--------------------------------------------"


`qui' lasso linear `lhs' `selind' if `seldep', selection(adaptive, `unpenalized')
`qui' est store lassomodel
`qui' lassoknots
mat B=r(table)
local rowsofB=rowsof(B)

foreach var of varlist `selind' {
	qui replace `var'=``var'_raw'
}

if "`simulation'"!="" {		//Calculate all three Heckman Post-Lasso regressions at once
	*Ftest
	`qui' heckman `lhs' `selind', select(`seldep'=`selind') `twostep' `robust' `cluster'
	`qui' est store testing
	`qui' testparm `selind', equation(`lhs')
	local pvalue=r(p)
	qui sum `lhs' if `seldep'
	local rN=r(N)
	local r=1
	while (`pvalue'<(0.10/ln(`rN')) & `r'<=`rowsofB') {
		local steps = B[`r',1]
		`qui' est restore lassomodel	
		`qui' lassoselect id=`steps'
		local nonzero=e(k_nonzero_sel)
		local allvars=e(k_allvars)
		if `nonzero'!=`allvars' {
			local active_ftest=e(allvars_sel)
			local testset: list selind-active_ftest
			`qui' est restore testing
			*`qui' test `testset'
			`qui' testparm `testset', equation(`lhs')
			local pvalue=r(p)
		}
		else {
			local active_ftest `selind'
		}
		local r=`r'+1	
	}
	*CV
	`qui' est restore lassomodel
	local lam=e(lambda_cv)
	lassoselect lambda=`lam'
	local nonzero=e(k_nonzero_sel)
	if `nonzero'!=0 {
		local active_cv=e(allvars_sel)
	}
	else {
		local active_cv
	}
	*CVSE
	`qui' est restore lassomodel
	local lam=e(lambda_serule)
	lassoselect lambda=`lam'
	local nonzero=e(k_nonzero_sel)
	if `nonzero'!=0 {
		local active_cvse=e(allvars_sel)
	}
	else {
		local active_cvse
	}	
}
else {
	if "`tuning'"=="ftest" {
		
		dis ""
		dis "{txt}-----------------------------------------------------------"
		dis "{txt}Post-Lasso-OLS Heckman model (F-test as stopping criterion)"
		dis "{txt}-----------------------------------------------------------"

		`qui' heckman `lhs' `selind', select(`seldep'=`selind') `twostep' `robust' `cluster'
		`qui' est store testing
		`qui' testparm `selind', equation(`lhs')
		local pvalue=r(p)
		qui sum `lhs' if `seldep'
		local rN=r(N)
		local r=1
		while (`pvalue'<(0.10/ln(`rN')) & `r'<=`rowsofB') {
			local steps = B[`r',1]
			`qui' est restore lassomodel	
			`qui' lassoselect id=`steps'
			local nonzero=e(k_nonzero_sel)
			local allvars=e(k_allvars)
			if `nonzero'!=`allvars' {
				local active=e(allvars_sel)
				local testset: list selind-active
				`qui' est restore testing
				*`qui' test `testset'
				`qui' testparm `testset', equation(`lhs')
				local pvalue=r(p)
			}
			else {
				local active `selind'
			}
			local r=`r'+1	
		}
	}

	if "`tuning'"=="cv" {
		
		dis ""
		dis "{txt}-----------------------------------------------------------------------"
		dis "{txt}Post-Lasso-OLS Heckman model (CV without SE-rule as stopping criterion)"
		dis "{txt}-----------------------------------------------------------------------"

		local lam=e(lambda_cv)
		lassoselect lambda=`lam'
		local nonzero=e(k_nonzero_sel)
		if `nonzero'!=0 {
			local active=e(allvars_sel)
		}
		else {
			local active
		}
	}

	if "`tuning'"=="cvse" {
		
		dis ""
		dis "{txt}--------------------------------------------------------------------"
		dis "{txt}Post-Lasso-OLS Heckman model (CV with SE-rule as stopping criterion)"
		dis "{txt}--------------------------------------------------------------------"

		local lam=e(lambda_serule)
		lassoselect lambda=`lam'
		local nonzero=e(k_nonzero_sel)
		if `nonzero'!=0 {
			local active=e(allvars_sel)
		}
		else {
			local active
		}	
	}
}

if "`notpen'"!="" {
	if "`simulation'"!="" {	
		local active_ftest: list notpen | active_ftest
		local active_cv: list notpen | active_cv
		local active_cvse: list notpen | active_cvse
	}
	else {
		local active: list notpen | active
	}
}

qui local coll `s(collinear)'
	qui	_rmcoll `selind', ///
		`constan' `coll' 
	local selind `r(varlist)'	

if "`simulation'"!="" {
	foreach p in ftest cv cvse {
		heckman `lhs' `active_`p'', select(`seldep'=`selind') `twostep' `robust' `cluster'
		local exclusions_`p': list selind-active_`p'
		local b_x_`p'=[y]_b[x]
		local se_x_`p'=[y]_se[x]
		local b_lambda_`p'=[/mills]_b[lambda]
		local se_lambda_`p'=[/mills]_se[lambda]
	}
	eret scalar nuexclusions_ftest=`:word count `exclusions_ftest''
	eret local exclusions_ftest = "`exclusions_ftest'"
	eret scalar nuexclusions_cv=`:word count `exclusions_cv''
	eret local exclusions_cv = "`exclusions_cv'"
	eret scalar nuexclusions_cvse=`:word count `exclusions_cvse''
	eret local exclusions_cvse = "`exclusions_cvse'"
	foreach p in ftest cv cvse {
		eret local b_x_`p'=`b_x_`p''
		eret local se_x_`p'=`se_x_`p''
		eret local b_lambda_`p'=`b_lambda_`p''
		eret local se_lambda_`p'=`se_lambda_`p''
	}
}
else {
	heckman `lhs' `active', select(`seldep'=`selind') `twostep' `robust' `cluster'
	local exclusions: list selind-active
	dis "{txt}Found exclusion restrictions:" "{res} `exclusions'" 
	eret scalar nuexclusions=`:word count `exclusions''
	eret local cmd2 "heckman_lasso"
	eret local exclusions = "`exclusions'"
}

end

