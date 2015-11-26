# Issues, Future Directions, and Ideas

1. The package needs a consistent naming convention:
	* `relabel.varray()` and `truncvarray()` are inconsistent -- the latter is missing a period (but that's because if there's a period, `roxygen` (I believe) wants to make it an S3 object.
	* `rvinesubset()` is not of the form `verb.object()` like the others tend to be.
	* When is it appropriate to use `"varray"` or `"rvine"`? 
		* And should `"vine"` be used instead of `"rvine"`, for simplicity? 
	* May want to give a name to the "convenient vine arrays" which put the labels on top. Would need to rename `Atocon()` and `contoA()` too to reflect that.

2. It might be useful to have a function to check whether a vine array is valid (like in `CopulaModel`, but allowing for truncated vines).

3. It might be useful to somehow add functionality so that `U-V` asymmetric copulas can be used:
	* In addition to the 4 possible reflections, there's also the `U-V` flip for those -- meaning that there's 8 copula families in total. Are the latter 4 really needed?
	* Involves being careful about what variable comes first when specifying a vine edge, and getting the copula right when the vine array is moved around (such as re-leafing).

4. As more functionality is developed, it might be useful to work with "vine objects", which would contain the vine array, copula model matrix, and copula parameter matrix instead of working with each separately. They should be specified in that order, with the possibility that "upstream" features (like the copula parameter matrix, or both copula and parameter matrices) are unspecified and therefore unknown.
	* `reform.copmat()` would not be needed anymore.

5. It might be useful to change the format of the copula families so that **number of parameters** and **parameter space** can be extracted (so that the user doesn't have to look it up all the time). Here are some ideas:
	* Could add that info as 'attributes' for each function. 
		* Problem: too many functions to do this too!
	* Instead of having a bunch of functions for each family under a naming convention, make a single object that is a list of such functions. This list could also contain this additional info. 
		* For example, `frk$pcop` and `frk$qcond` for the cdf and conditional quantile functions, and `frk$npar` and `frk$parspace` for number of parameters and parameter space.
		* Problem: That's a lot that would need changing.
	* Have objects come with the `copsupp` package that contains this information. Specifically:
		* `cparspacecop()` can be a function that takes the parameter vector and return `TRUE` if it's in the parameter space copula "cop", and `FALSE` otherwise.
		* `ncparcop` or `nparcop` can be integer objects carrying the number of parameters family "cop" accepts. Or maybe put all that info in one list called `ncpar` to be used like `ncpar$cop`. Or just have it as a function that takes no arguments, like `ncparcop()`. Or have one function that accepts the copula name as an argument.

6. `fvinesim()` does not work when an independence copula is in the mix.
	* Solution: trick `CopulaModel` by substituting the independence copula with a copula family with a parameter such that it equals the independence copula.