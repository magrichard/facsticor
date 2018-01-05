#
# Authors: Florent Chuffart, INSERM and Magali Richard, CNRS
# florent.chuffart@univ-grenoble-alpes.fr
# mag.richard@ens-lyon.fr
#
#---------------------------------------------
#' Export final data
#'
#' This function works on FL1H, FSC and SSC data acquired on 96-well plates
#' @param data The data to be exported
#' @return A data tab
#' @export

export_data = function (data ){
  meanFL1lin = NULL
  varFL1lin  = NULL
  meanFL1log = NULL
  varFL1log  = NULL
  on_frac  = NULL
  ncells  = NULL
  meanFSC = NULL
  meanSSC = NULL
  on_thresholdlin = NULL
  on_thresholdlog = NULL
  for (i in 1:dim(data$anno)[1]){
      print(i)
    cur = data$ngdat
    linfl1 = 10^(4 * cur[[i]]$FL1.H / 1024) * 1
    on_thresholdlin[i] = 10^(4 * data$anno$on_threshold[i] / 1024) * 1
    on_thresholdlog[i] = data$anno$on_threshold[i]
    on_frac[i]   = data$anno$on_frac[i]
    meanFL1lin[i]   = mean(linfl1)
    meanFL1log[i] = mean(cur[[i]]$FL1.H)
    varFL1lin[i]    = var(linfl1)
    varFL1log[i]    = var(cur[[i]]$FL1.H)
    ncells[i]    = dim(cur[[i]])[1]
    meanFSC[i]   = mean(cur[[i]]$FSC.H)
    meanSSC[i]   = mean(cur[[i]]$SSC.H)
  }
  export <- cbind(data$anno[,c("date","time", "gal", "strain", "bg")], meanFL1lin, varFL1lin, on_frac, on_thresholdlin, meanFL1log, varFL1log, on_frac, on_thresholdlog, ncells, meanFSC, meanSSC)
  date=unique(data$anno$date)
  write.table(export, file = paste("outputs/",date,"_gal3_priors.tab", sep = ""), quote = FALSE, row.names=FALSE)
  return(export)
}



# export_data = function (data){
#   meanFL1 = NULL
#   varFL1  = NULL
#   ncells  = NULL
#   meanFSC = NULL
#   meanSSC = NULL
#   nb_on = NULL
#   nb_off = NULL
#   fraction_on  = NULL
#   for (i in 1:dim(data$anno)[1]){
#     cur = data$ngdat
#     linfl1 = 10^(4 * cur[[i]]$FL1.H / 1024) * 1
#     meanFL1[i]   = mean(linfl1)
#     varFL1[i]    = var(linfl1)
#     ncells[i]    = dim(cur[[i]])[1]
#     meanFSC[i]   = mean(cur[[i]]$FSC.H)
#     meanSSC[i]   = mean(cur[[i]]$SSC.H)
#     nb_on[i]   = data$anno$nb_on[i]
#     nb_off[i]   = data$anno$nb_off[i]
#     fraction_on[i]   = data$anno$fraction_on[i]
#   }
#   export <- cbind(data$anno[,c("date","time", "gal", "strain", "bg")], meanFL1, varFL1, nb_on, nb_off, fraction_on, ncells, meanFSC, meanSSC)
#   date=unique(data$anno$date)
#   write.table(export, file = paste("outputs/",date,"_gal3_priors.tab", sep = ""), quote = FALSE, row.names=FALSE)
#   return(export)
# }


#  export_data = function (data, type="norm" ){
#   meanFL1 = NULL
#   varFL1  = NULL
#   on_frac  = NULL
#   ncells  = NULL
#   meanFSC = NULL
#   meanSSC = NULL
#   on_threshold = NULL
#   for (i in 1:dim(data$anno)[1]){
#     cur = data$ngdat
#     if (type=="norm") {
#       cur = data$ngdat
#       on_threshold[i] = 10^(4 * data$anno$on_threshold[i] / 1024) * 1
#       on_frac[i]   = data$anno$ngon_frac[i]
#     } else if (type=="not_norm") {
#       cur = data$gdat
#       on_threshold[i] = 10^(4 * data$anno$on_threshold[i] / 1024) * 1
#       on_frac[i]   = data$anno$gon_frac[i]
#     }
#     linfl1 = 10^(4 * cur[[i]]$FL1.H / 1024) * 1
#     meanFL1[i]   = mean(linfl1)
#     varFL1[i]    = var(linfl1)
#     ncells[i]    = dim(cur[[i]])[1]
#     meanFSC[i]   = mean(cur[[i]]$FSC.H)
#     meanSSC[i]   = mean(cur[[i]]$SSC.H)
#   }
#   export <- cbind(data$anno[,c("date","time", "gal", "strain", "bg")], meanFL1, varFL1, on_frac, on_threshold, ncells, meanFSC, meanSSC)
#   date=unique(data$anno$date)
#   write.table(export, file = paste("outputs/",date,"_gal3_priors.tab", sep = ""), quote = FALSE, row.names=FALSE)
#   return(export)
# }
