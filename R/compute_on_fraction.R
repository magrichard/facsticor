#
# Authors: Magali Richard, CNRS
# mag.richardt@ens-lyon.fr
#
#---------------------------------------------
#' Perfom annova on the data
#'
#' This function works on FL1H, FSC and SSC data acquired on 96-well plates
#' @param data The data normalized
#' @param threshold The ON_OFF threshold
#' @return An data table with on fractions
#' @export

add_on_fraction = function (data, threshold) {
  data$anno$on_threshold = threshold

    idx = 1:length(data$gdat)
  data$anno$on_cells = NA
  data$anno$off_cells = NA
  for (i in idx) {
    data$anno[i, ]$on_cells = length(which(data$ngdat[[i]]$FL1.H > data$anno[i, ]$on_threshold))
    data$anno[i, ]$off_cells = length(which(data$ngdat[[i]]$FL1.H <= data$anno[i, ]$on_threshold))
  }
  data$anno$on_frac = data$anno$on_cells/data$anno$nb_gcells
return(data)
}

# add_on_fraction = function (data) {
#   data$anno$ngoff_threshold = NA
#   data$anno$goff_threshold = NA
#   bgs = sort(unique(data$anno[data$anno$time == 0, ]$bg))
#   NA_data = setdiff(unique(data$anno$bg), bgs)
#   print(paste("the following data are NA because there is no 0 time point", NA_data, sep=" "))
#   for (bg in bgs ){
#     idx = which(data$anno$bg == bg & data$anno$time == 0)
#     ngoff_threshold = max(sapply(idx, function(i){quantile(data$ngdat[[i]]$FL1.H, 0.95)}))
#     data$anno[data$anno$bg == bg, ]$ngoff_threshold = ngoff_threshold
#     goff_threshold = max(sapply(idx, function(i){quantile(data$gdat[[i]]$FL1.H, 0.95)}))
#     data$anno[data$anno$bg == bg, ]$goff_threshold = goff_threshold
#   }
#   idx = 1:length(data$gdat)
#   data$anno$ngon_cells = NA
#   data$anno$ngoff_cells = NA
#   data$anno$gon_cells = NA
#   data$anno$goff_cells = NA
#   for (i in idx) {
#     data$anno[i, ]$ngon_cells = length(which(data$ngdat[[i]]$FL1.H > data$anno[i, ]$ngoff_threshold))
#     data$anno[i, ]$ngoff_cells = length(which(data$ngdat[[i]]$FL1.H <= data$anno[i, ]$ngoff_threshold))
#     data$anno[i, ]$gon_cells = length(which(data$gdat[[i]]$FL1.H > data$anno[i, ]$goff_threshold))
#     data$anno[i, ]$goff_cells = length(which(data$gdat[[i]]$FL1.H <= data$anno[i, ]$goff_threshold))
#   }
#   data$anno$ngon_frac = data$anno$ngon_cells/data$anno$nb_gcells
#   data$anno$gon_frac = data$anno$gon_cells/data$anno$nb_gcells
#   return(data)
# }

#---------------------------------------------
#' Perfom annova on the data
#'
#' This function works on FL1H, FSC and SSC data acquired on 96-well plates
#' @param data The data normalized
#' @param on_cells Cell used to modelized on cells
#' @param off_cells Cell used to modelized off cells
#' @return An data table with the number of cell on and the number of cell off
#' @export

get_probas = function (data,  mu_on, var_on, mu_off, var_off) {
  compute_proba = lapply(1:length(data$ngdat), function(idx) {
    print(paste(idx, "/", length(data$ngdat), "computing..."))
    data$anno$key = paste(data$anno$strain, data$anno$time, data$anno$gal, sep="_")
    cells = data$ngdat[[idx]]$FL1.H
    probas = lapply(1:length(cells), function(i) {
      on_proba = pnorm(cells[i],  mu_on,  var_on)
      off_proba = 1- pnorm(cells[i],  mu_off,  var_off)
      return(cbind(on_proba,off_proba))
    })
    probas = do.call(rbind, probas)
    probas = data.frame(lapply(data.frame(probas, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
    proba.threshold = 1e-05
    nb_on = sum(probas$on_proba > probas$off_proba & probas$on_proba > proba.threshold)
    nb_off = sum(probas$off_proba > probas$on_proba & probas$off_proba > proba.threshold)
    loglike = sum(log(rowSums(probas)))
    normloglike = loglike/length(cells)
    fraction_on = nb_on/length(cells)
    return(list(nb_on= nb_on, nb_off= nb_off, nloglike_on_off = normloglike, fraction_on = fraction_on))
  })
  compute_proba = do.call(rbind, compute_proba)
  compute_proba = data.frame(lapply(data.frame(compute_proba, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
  data$anno = cbind(data$anno, compute_proba)
  return(data)
}

#---------------------------------------------
#' Perfom annova on the data
#'
#' This function works on FL1H, FSC and SSC data acquired on 96-well plates
#' @param data The data normalized
#' @param eps eps
#' @param  mu_on The mean of ON cells
#' @param  sd_on The standart deviation of ON cells
#' @param  mu_off The standart deviation of OFF cells
#' @param  sd_on The variance of OFF cells
#' @return An data table with the number of cell on and the number of cell off
#' @export
#'

# get_threshold = function(c, eps, step, iter, mu_on, var_on, mu_off, var_off) {
#   threshold =lapply(1:iter, function(i) {
#     if ((1-pnorm(c, mu_off, var_off)) <  (pnorm(c, mu_on, var_on))){
#     c = c - step
#     } else if ((1-pnorm(c, mu_off, var_off)) <  (pnorm(c, mu_on, var_on))){
#     c = c + step
#     } else if ((1-pnorm(c, mu_off, var_off)) >  (pnorm(c, mu_on, var_on)-eps) & (1-pnorm(c, mu_off, var_off)) <  (pnorm(c, mu_on, var_on)+eps)){
#     print(paste(c, "is the winner", sep =' '))
#     }
#     return(c)
#   })
#   return(threshold[[iter]])
# }

get_threshold = function(t, eps,  mu_on, sd_on, mu_off, sd_off) {
  diff = function(t) {
    proba_on = pnorm(t, mu_on, sd_on)
    proba_off = 1-pnorm(t, mu_off, sd_off)
    return(proba_on-proba_off)
  }
  # if (diff(t) > 0 & diff(t) > eps) {
	if (diff(t) > 0 & abs(diff(t)) > eps) {
    t = mean(c(t, mu_off))
    get_threshold(t, eps,  mu_on, sd_on, mu_off, sd_off)
  # } else if (diff(t) < 0 & diff(t) > eps) {
	} else if (diff(t) < 0 & abs(diff(t)) > eps) {
    t = mean(c(t, mu_on))
      get_threshold(t, eps,  mu_on, sd_on, mu_off, sd_off)
    }  else if ( abs(diff(t)) < eps) {
      print(diff(t))
  return (t)
    }
}


#---------------------------------------------
#'Fit ON gaussian on the data
#'
#' This function works on FL1H, FSC and SSC data acquired on 96-well plates
#' @param on_cells The ID of well in which cells are supposed to be on
#' @return An data table with the number of cell on and the number of cell off
#' @export
#'

fit_ON_gaussian = function(on_cells){
  mixmdl = mixtools::normalmixEM(on_cells)
  plot(mixmdl,which=2)
  return(cbind(mu_on=mixmdl$mu[2], sigma_on=mixmdl$sigma[2]))
}


#---------------------------------------------
#'Compute BIC and loglike value of gaussian models
#'
#' This function uses the Model-Based Clustering function Mclust to extract the optimal BIC value and the corresponding log-likelihood for a 1-gaussian model and a 2-gaussian model.
#' @param data The data table
#' @return An data table with BIC and loglik for each model
#' @importFrom  mclust Mclust
#' @importFrom  mclust mclustBIC
#' @export
#'

compute_model_values = function (data) {
  data$anno$gau_1_BIC = NA
  data$anno$gau_2_BIC = NA
  data$anno$gau_1_loglik = NA
  data$anno$gau_2_loglik = NA
  idx = 1:length(data$gdat)
  for (i in idx) {
	  x = data$ngdat[[i]]$FL1.H
    data$anno[i, ]$gau_1_BIC = mclust::Mclust(x, G=1)$bic
    data$anno[i, ]$gau_2_BIC =  mclust::Mclust(x, G=2)$bic
    data$anno[i, ]$gau_1_loglik =  mclust::Mclust(x, G=1)$loglik
    data$anno[i, ]$gau_2_loglik =  mclust::Mclust(x, G=2)$loglik
	
  }
return(data)
}
