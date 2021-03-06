# heckman_lasso
 Adaptive post-Lasso Heckman regression

## Description
 `heckman_lasso` estimates a post-Lasso Heckman regression where the exclusion restrictions are determined in a data-driven way.
    Lasso estimations are performed using the built-in Stata command `lasso`, which requires Stata 16. Post-Lasso Heckman estimation is performed using the built-in Stata command `heckman`. For general information about adaptive Lasso see Zou (2006).
    
## Installing `heckman_lasso`
 Way 1: Use the Stata module `github` to install `heckman_lasso`. Type the following code in Stata:
 
 1. ```{js}
    net install github, from("https://haghish.github.io/github/")
    ```
 2. ```{js}
    github install farbmacher/heckman_lasso
    ```

 Way 2: Copy the `heckman_lasso.ado` and `heckman_lasso.sthlp` files into your personal ado folder or the current working directory.

## Example

 Perform a post-Lasso Heckman regression with `y` the outcome (`ds` indicating observations for which `y` is observed) 
 and `x1, x2,...` a set of exogenous variables containing both control variables and potential exclusion restrictions

 1. ```{js}
     heckman_lasso y x1 x2, seldep(ds) twostep
     ```
 
 Let `z1, z2,...` be a set of additional exogenous control variables, which should always be included in the model, then
 
 2. ```{js}
     heckman_lasso y x1 x2 z1 z2, seldep(ds) twostep notpen(z1 z2)
     ```
        
## References
 * Farbmacher, H. (2021): Selection Models with Data-Driven Exclusion Restrictions in Managerial Economics, Discussion Paper.
 * Zou H. (2006): The Adaptive Lasso and Its Oracle Properties, *Journal of the American Statistical
             Association* 101, 1418-1429.
