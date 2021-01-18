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
 Let `y` be the outcome (with `ds` indicating observations where `y` is observed) and `x1, x2,...` be a set of exogenous variables containing both control variables and potential exclusion restrictions, then the post-Lasso Heckman regression would be
        
 1. ```{js}
    heckman_lasso y x1 x2, seldep(ds)
    ```
        
## References
 * Farbmacher, H. (2021): Sample Selection Models with Unknown Exclusion Restrictions, Discussion Paper.
 * Zou H. (2006): The Adaptive Lasso and Its Oracle Properties, *Journal of the American Statistical
             Association* 101, 1418-1429.
