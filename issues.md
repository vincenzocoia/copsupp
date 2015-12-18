# Issues, Future Directions, and Ideas

1. It might be useful to have a function to check whether a vine (including the array) is valid (like in `CopulaModel`, but allowing for truncated vines).

2. It might be useful to somehow add functionality so that `U-V` asymmetric copulas can be used:
	* In addition to the 4 possible reflections, there's also the `U-V` flip for those -- meaning that there's 8 copula families in total. Are the latter 4 really needed?
	* Involves being careful about what variable comes first when specifying a vine edge, and getting the copula right when the vine array is moved around (such as re-leafing).

3. Package should probably be renamed "rvine".

4. The vignette needs to be improved:
	* It contains too much extraneous information that would only be useful for developers.
	* (This might follow from the previous downfall) it reads more like a pile of information, as opposed to a meaningful "story" that communicates what a user would want to use the package for. 
	* Get the vignette working so that, upon typing `vignette('copsupp')` in R, the vignette appears.

5. `fit.rvine()` outputs a different vine array depending on the order of the variables input into the `vars` argument. It shouldn't.

6. `qcondrvine()` would be handy to have in some situations.

7. Problems with `RVineCopSelect()` in `VineCopula` (may have to write my own version):
	* Forces you to use a pre-defined set of bivariate copula families.
	* Doesn't let you choose the copula model matrix. It would even be nice to be able to choose _parts_ of the copula model matrix too.
	* (Also with `BiCopSelect()`) Puts negative parameters on 90- or 270-degree rotated copulas, but those models actually have positive parameters. 
	* (Also with `BiCopSelect()`) Can't seem to fit the BB6 copula properly. See this example: 

    udat <- CopulaModel::rbb6(10000, c(3,3))
    lapply(VineCopula::BiCopSelect(1-udat[, 1], udat[, 2], familyset = 28), identity)

