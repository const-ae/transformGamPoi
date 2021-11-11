# transformGamPoi 1.1.x

* Implement faster scaling with size factors for acosh and log-based transformations


# transformGamPoi 0.1.x

* Add clipping functionality to `residual_transform()`
* Add check against `residual_type` argument in `transformGamPoi()` 
* Fix bug in `acosh_transform()` related to sparse input and `on_disk = FALSE`
* Change default of `overdispersion_shrinkage` to `TRUE` if `overdispersion = TRUE` 
for `acosh_transform()` and `shifted_log_transform()`

# transformGamPoi 0.1.0

* Initial release of transformGamPoi on GitHub https://github.com/const-ae/transformGamPoi
