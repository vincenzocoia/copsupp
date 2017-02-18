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
	* NOTE: Add an additional "parameter", being the "reflection" and/or "permutation". This forms a group called the dihedral group, which consists of 8 elements. Which element the copula family takes can be considered a parameter, albeit with a finite base.
		* Make a function that does such group theory computations. Allow the user to input reflections and permutations, and it will compute the proper element.
		* When describing a copula family, like "igcop", its "parameter space" should be expanded to indicate how many unique elements in the dihedral group there are. My first feeling is that a copula family can either have 1, 4, or 8 unique elements forming the group. 
	* Possibly make an object that describes a copula _family_, where one can describe the parameter space (again, parameter here means in the traditional sense, but also the elements of the dihedral group that are unique). 
		* Or, perhaps just make this an entry in the 'bicop' object? The reason I don't like that is because the family is just that: a family. A set of copulas. The family should therefore describe the parameter space, which a simple character string can't do. I think it's a little unwieldy to carry the parameter space along with the copula family name in the copula object, because it gives the impression that the copula family _name_ and the _parameter space_ are not intricately linked, when really they are. 
	* I don't need to decide between, for example, `dcop(u, v, cop="frk")` and `dfrk(u, v)` -- just include both!

2. When writing about the package, differentiate between the variable labels (which are arbitrary), and the order (which are integers 1,...,d). Otherwise, always labelling with 1,...,d adds confusion: hard to know when to look at 1,...,d as the order or just arbitrary labels. (So, sometimes the labels correspond to the variable names, or "column number" of the data if they're integers, or they may refer to the order that they're introduced into the vine, in which case they'd be 1:ncol(A). I find the former much more intuitive)

3. How about making a function to extract the copula from an edge? For example, get the copula linking (1,2)|3.

4. Additional functionality:
	* I'm currently working on my PhD analysis, and it would be useful if I could extract the AIC from a fitted vine.