# countprop 1.1.1

* Added argument `lr.penalty` to functions `mleLR()` and `mlePath()`
allow penalization of the covariance matrix
on both the additive log-ratio and centered log-ratio scales.

* Added argument `lr` to both `logitNormalVariation()` and `naiveVariation()`
to allow calculation of metrics of proportionality on either
the ALR or CLR scale.

* Added function `convertSigma()` to allow easy conversion
between covariance matrices on ALR and CLR scales.

* Function `pluginVariation()` is now deprecated as its results will
now be equivalent to `naiveVariation()`.

# countprop 1.0.1

* Initial CRAN submission.

