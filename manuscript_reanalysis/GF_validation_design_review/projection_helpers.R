# Helpers for the six-predictor exploratory projection. No model fitting here.
gf_predictors <- c('prec_08','prec_01','bio_15','bio_3','tmin_10','bio_19')

# Apply ONLY when extracting raw calendar-month rasters, never to aligned tables.
calendar_month <- function(aligned_month, latitude) {
  stopifnot(length(aligned_month)==1L, aligned_month %in% 1:12,
            all(is.finite(latitude)), all(abs(latitude)<=90))
  ifelse(latitude < 0, (aligned_month + 5L) %% 12L + 1L, aligned_month)
}

check_environment <- function(x) {
  stopifnot(all(gf_predictors %in% names(x)))
  x <- as.data.frame(x[,gf_predictors,drop=FALSE])
  if (!all(vapply(x,is.numeric,logical(1))) || any(!is.finite(as.matrix(x))))
    stop('All six environmental predictors must be numeric and finite.')
  x
}

# Outside any marginal range is a transparent extrapolation warning, not a
# complete multivariate novelty test. Nearest distance adds continuous context;
# no arbitrary distance cutoff is used to classify habitat suitability.
climate_novelty <- function(training, query) {
  training <- unique(check_environment(training))
  query <- check_environment(query)
  lo <- vapply(training,min,numeric(1)); hi <- vapply(training,max,numeric(1))
  outside <- sweep(as.matrix(query),2,lo,'<') | sweep(as.matrix(query),2,hi,'>')
  center <- colMeans(training); scale <- vapply(training,sd,numeric(1))
  if(any(!is.finite(scale) | scale<=0)) stop('Constant training predictor.')
  ztrain <- sweep(sweep(as.matrix(training),2,center,'-'),2,scale,'/')
  zquery <- sweep(sweep(as.matrix(query),2,center,'-'),2,scale,'/')
  nearest <- vapply(seq_len(nrow(query)),function(i)
    sqrt(min(rowSums(sweep(ztrain,2,zquery[i,],'-')^2))),numeric(1))
  data.frame(outside_training_range=rowSums(outside)>0,
             n_predictors_outside=rowSums(outside),
             predictors_outside=apply(outside,1,function(z) paste(gf_predictors[z],collapse=';')),
             nearest_standardized_climate_distance=nearest)
}

transform_offset <- function(model, current, future) {
  current <- check_environment(current); future <- check_environment(future)
  stopifnot(nrow(current)==nrow(future))
  a <- as.matrix(predict(model,current)); b <- as.matrix(predict(model,future))
  if(!setequal(colnames(a),gf_predictors) || !setequal(colnames(b),gf_predictors))
    stop('Unexpected GF transformation dimensions.')
  a <- a[,gf_predictors,drop=FALSE]; b <- b[,gf_predictors,drop=FALSE]
  if(any(!is.finite(a)) || any(!is.finite(b))) stop('Nonfinite GF transformation.')
  sqrt(rowSums((a-b)^2))
}
