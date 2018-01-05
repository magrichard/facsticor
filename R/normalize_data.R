#
# Authors: Florent Chuffart, INSERM and Magali Richard, CNRS
# florent.chuffart@univ-grenoble-alpes.fr
# mag.richardt@ens-lyon.fr
#
#---------------------------------------------
#' Perfom annova on the data
#'
#' This function works on FL1H, FSC and SSC data acquired on 96-well plates
#' @param data A particular element of the list data$gdat (coresponds to wells present in all plates for normalization)
#' @param model The model to be used
#' @return An annova model
#' @examples à finaliser
#' load("data/data.rda")
#' idx_to_be_keep = data$anno$rep%in%2:5 & data$anno$strain != "GY6" & TRUE
#' model_data = filter_data(data, !idx_to_be_keep)
#' model_p.s.t = mean_gfl1~plate+strain+time
# 'm_annov = perform_anova(data=model_data, model=model_p.s.t)
#' @export

perform_anova = function(data, model) {
  mean_gfl1 = sapply(1:length(data$gdat), function(i) {
    mean(data$gdat[[i]]$FL1.H)
  })
  median_gfl1 = sapply(1:length(data$gdat), function(i) {
    median(data$gdat[[i]]$FL1.H)
  })
  foo = data$anno[,c("plate", "well","time", "date", "strain", "rep", "bg","gal")]
  foo$mean_gfl1 = mean_gfl1
  foo$median_gfl1 = median_gfl1
  foo$plate = as.factor(foo$plate)
  foo$well = as.factor(foo$well)
  foo$time = as.factor(foo$time)
  foo$date = as.factor(foo$date)
  foo$strain = as.factor(foo$strain)
  foo$gal = as.factor(foo$gal)
  foo$bg = as.factor(foo$bg)
  foo$rep = as.factor(foo$rep)
  m = aov(model, data=foo)
  print(shapiro.test(residuals(m)))
  print(summary(m))
  t = model.tables(m)
  print(t)
  print(foo)
  return(m)
}

#---------------------------------------------
#' Normalize the data according to the annova model
#'
#' This function works on FL1H, FSC and SSC data acquired on 96-well plates
#' @param data The data to be corrected
#' @param offsets The offsets to be used
#' @return A normalized data list
#' @examples à finaliser
#' load("data/data.rda")
#' final_data = filter_data(data, data$anno$rep != 1)
#' offsets_plate = model.tables(m_annov)$tables$plate
#' final_data = adjust_effect(final_data, offsets_plate, "plate")
#' @export

adjust_effect = function(data, offsets, offset_key) {
  data$ngdat = data$gdat
  for (i in 1:length(data$gdat)) {
    offset_key_val = data$anno[[offset_key]][i]
    offset = offsets[[offset_key_val]]
    data$ngdat[[i]]$FL1.H = data$gdat[[i]]$FL1.H - offset
  }
  foo = lapply(data$ngdat, function(l){
    list(
      median_ngfl1 = median(l$FL1.H),
      median_ngfsc = median(l$FSC.H),
      median_ngssc = median(l$SSC.H),
      sd_ngfl1 = sd(l$FL1.H),
      sd_ngfsc = sd(l$FSC.H),
      sd_ngssc = sd(l$SSC.H),
      nb_ngcells = length(l$FL1.H)
    )
  })
  foo = do.call(rbind, foo)
  foo = data.frame(lapply(data.frame(foo, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
  data$anno = cbind(data$anno, foo)
  return(data)
}
