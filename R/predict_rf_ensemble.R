#' Make a prediction from a ranger::ranger object emulating an enemble instance
#' 
#' @param model_object a tidymodels workflow object. For `ranger` models, must have been fit with `keep.inbag = TRUE`
#' @param new_data data to predict, must have a column named `"parameter_seed"` in addition to the model variables
predict_tidy_ensemble <- function(model_object, new_data) {
 # if(inherits(model_object, "bundle")) model_object <- unbundle(model_object)
  stopifnot("parameter_seed" %in% colnames(new_data))
  suppressWarnings(preds <- ranger:::predict.ranger(model_object, data = new_data, type = "se"))
  #preds_se <- predict(model_object, new_data = new_data, type = "conf_int")
  if(any(is.nan(preds$se))) stop("Provide more new_data to estimate standard errors")
  preds_out <- numeric(nrow(new_data))
  for (i in seq_along(preds_out)) {
    set.seed(new_data$parameter_seed[i])
    preds_out[i] <- rnorm(1, mean = preds$predictions[i], 
                       sd = preds$se[i])
  }
  
  preds_out
}

