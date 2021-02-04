{smcl}
{cmd:help heckman_lasso}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col :{hi:heckman_lasso} {hline 2} Adaptive post-Lasso Heckman regression}{p_end}
{p2colreset}{...}

{title:Syntax}

	{opt heckman_lasso} {depvar} {indepvars}{cmd:,} {cmd:seldep(}{var}{cmd:)} [{cmd:notpen(}{vars}{cmd:)} {it:options}]
	
{pstd}

{synoptset 20}{...}
{synopthdr}
{synoptline}
{synopt :{opt seldep(varname)}}specifies the selection indicator; seldep=1 if y observed, otherwise seldep=0 [required]{p_end}
{synopt :{opt notpen(varlist)}}specifies a list of variables that should, in any case, be included in the final post-Lasso Heckman regression {p_end}
{synopt :{opt lassofirst}}use Lasso estimates to derive weights for adaptive Lasso {p_end}
{synopt :{opt twostep}}performs two-step Heckman regression instead of Maximum Likelihood {p_end}
{synopt :{opt robust}}calculates robust SEs; works only with Maximum Likelihood {p_end}
{synopt :{opt cluster(clustvar)}}calculates cluster-robust SEs; works only with Maximum Likelihood {p_end}
{synopt :{opt verbose}}display additional results{p_end}
{synoptline}

{title:Description}

{pstd} {cmd:heckman_lasso} estimates a post-Lasso Heckman regression where the exclusion restrictions are determined in a data-driven way.
Lasso estimations are performed using the built-in Stata command {manhelp lasso R:lasso linear}, which requires Stata 16. Post-Lasso Heckman estimation 
is performed using the built-in Stata command {manhelp heckman R:heckman}. For general information about adaptive Lasso see Zou (2006).

{title:Example}

{pstd}Let {opt y} be the outcome (with {opt ds} indicating observations for which {opt y} is observed) and {opt x1, x2,...} be a set of exogenous variables 
containing both control variables and potential exclusion restrictions, then the post-Lasso Heckman regression would be{p_end}

{phang2}{cmd:. heckman_lasso y x1 x2, seldep(ds)}{p_end}

{title:Author}

{pstd}Helmut Farbmacher{p_end}
{pstd}Munich Center for the Economics of Aging (MEA){p_end}
{pstd}Max Planck Society, Germany{p_end}
{pstd}farbmacher@mea.mpisoc.mpg.de{p_end}

{title:Reference}

{psee}Farbmacher, H. (2021): {it:Sample Selection Models with Unknown Exclusion Restrictions}, Discussion Paper.{p_end}
{psee}Zou, H. (2006): {it:The Adaptive Lasso and Its Oracle Properties}, Journal of the American Statistical Association 101, 1418-1429.{p_end}

