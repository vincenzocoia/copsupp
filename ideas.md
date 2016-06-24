# Ideas

1. It may be useful to consider making an object of type "copula" (actually a copula _family_?). I'd have to be clever in what it looks like.
	* Query/action functions:
		* Check whether a copula exists.
		* Get parameter spaces.
		* Extract name of copula (still use "frk" over "frank", it's shorter).
		* Computation functions like `dcop()`, `pcop()`, `qcondcop()`, etc, with arguments `u,v,copula`.
		* Rotate/flip?
	* Construction function `copula()` (allow users to make their own):
		* Specify formulas for density, cdf, conditional, etc.
		* Specify parameter space.
	* Benefits:
		* Reduces a slew of functions (pcop, rcop, qcondcop, pcondcop12, etc) to one object.
		* Can handle copula flips and rotations more easily. 
	* Could consider multivariate copulas. Then `"rvine"` object can be a special case, so it would have class `c("rvine", "copula")`. 

2. When writing about the package, differentiate between the variable labels (which are arbitrary), and the order (which are integers 1,...,d). Otherwise, always labelling with 1,...,d adds confusion: hard to know when to look at 1,...,d as the order or just arbitrary labels. (So, sometimes the labels correspond to the variable names, or "column number" of the data if they're integers, or they may refer to the order that they're introduced into the vine, in which case they'd be 1:ncol(A). I find the former much more intuitive)