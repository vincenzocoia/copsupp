# Issues

1. The package needs a consistent naming convention:
	* `relabel.varray()` and `truncvarray()` are inconsistent -- the latter is missing a period (but that's because if there's a period, `roxygen` (I believe) wants to make it an S3 object.
	* `rvinesubset()` is not of the form `verb.object()` like the others tend to be.
	* When is it appropriate to use `"varray"` or `"rvine"`? 
		* And should `"vine"` be used instead of `"rvine"`, for simplicity? 

2. A vignette is needed to demonstrate how one might use this package.