# Issues, Future Directions, and Ideas

1. It might be useful to have a function to check whether a vine (including the array) is valid (like in `CopulaModel`, but allowing for truncated vines).

2. It might be useful to somehow add functionality so that `U-V` asymmetric copulas can be used:
	* In addition to the 4 possible reflections, there's also the `U-V` flip for those -- meaning that there's 8 copula families in total. Are the latter 4 really needed?
	* Involves being careful about what variable comes first when specifying a vine edge, and getting the copula right when the vine array is moved around (such as re-leafing).

3. Package should probably be renamed "rvine". Or, perhaps notation should switch to "gvine", or "regvine".

4. The vignette needs to be improved:
	* It contains too much extraneous information that would only be useful for developers.
	* (This might follow from the previous downfall) it reads more like a pile of information, as opposed to a meaningful "story" that communicates what a user would want to use the package for. 
	* Get the vignette working so that, upon typing `vignette('copsupp')` in R, the vignette appears.

5. `fit.rvine()` outputs a different vine array depending on the order of the variables input into the `vars` argument. It shouldn't.

6. `qcondrvine()` would be handy to have in some situations.

7. Problems with `RVineCopSelect()` in `VineCopula` (may have to write my own version):
	* Forces you to use a pre-defined set of bivariate copula families.
	* Doesn't let you choose the copula model matrix. It would even be nice to be able to choose _parts_ of the copula model matrix too.
	* When you request a "bvtcop" to be fit (via familyset = 2), it can instead fit a "bvncop" (returning family = 1) even though you didn't want one to be fit.
	* (Also with `BiCopSelect()`) Puts negative parameters on 90- or 270-degree rotated copulas, but those models actually have positive parameters. 
	* (Also with `BiCopSelect()`) Can't seem to fit the BB6 copula properly. See this example: 

    udat <- CopulaModel::rbb6(10000, c(3,3))
    lapply(VineCopula::BiCopSelect(1-udat[, 1], udat[, 2], familyset = 28), identity)

8. A fitting procedure for vines that's similar to that of the marginal should be included. So, parts of the vine might be specified, parts of the marginals might be specified, and maybe even some parameters ought to be shared amongst some edges.
	* Would it make sense to add components to an `"rvine"` object, such as data and corresponding measures of fit like AIC/BIC?

9. Perhaps a function `edge_possibilities(i, j, G, rvine = TRUE)` or something would be useful to see what variables can possibly go in `G[i,j]`, based on the parts of `G` already specified. This way I can "fill-in" a vine array if need be, but more importantly so that, when building a new layer, I can choose variables so that the result is still a regular vine (so that I can finish coding the `fitseq_rvine()` function).

10. In `pcondrvine()`, I should change 'vbls' to 'cond' (or some other name), so that it's the conditioning set -- it's just more natural.

11. Multivariate integration should be sped up (in, for example, `pcondrvine()`). Use gauss quadrature or something (keep in mind the bounds of integration are probably always contained in [0,1]).

12. `fitlayer_cnqr()` should allow for partial specification of parameters.

13. If you want to partially specify parts of a vine, like in `layeropts()` for example, it might be useful to have a function that returns the desired matrix (or whatever) with NULLs and NAs in their appropriate places. Something like `specifycop([2,3] = c("frk", "bvncop"))`.

14. `fitrvine_basic()`, and probably many others, throw an error when the data is a data frame (must be matrix). Fix this to allow for data frames.

15. It would be nice if `summary.rvine()` combined the copmat and cparmat by making the entries like "frk(0.8)" and "bvtcop(0.7, 4)" (this way you can see all parameters in case there are more than 1)

16. Speed-up multivariate integration by using a more appropriate method than nested "integrate".
	* cubature package.
	* Or, check the suite of packages for this purpose on [this page](https://cran.r-project.org/web/views/NumericalMathematics.html) of CRAN task view.

17. Instead of making flipped copulas by appending "u", "v", and "r", perhaps make functionals like `pcopu()` and `qcondcopu()` that take a function (like `pfrk()` and `qcondgum()`) and returns the same type of function, but flipped accordingly. This way, when defining their own, people don't have to make these flipped versions. Also, it makes programming easier.