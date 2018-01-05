#
# Author: Florent Chuffart, INSERM
# florent.chuffart@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Reads the facs data get the summary informations on data distributions.
#'
#' This function works on FL1H, FSC and SSC data acquired on 96-well plates. The data directory should be organized as follow: "date/plate/dat/your_fsc_data"
#' @param datadir The path to root directory of data
#' @param datadate The name of directory containing data, in date format: YYYMMDD, can be a vector of dates
#' @param annofile The experimenter annotation file, semi-colon separated, must at least contain the following informations :
#                 plate, well
#' @return A list of data.frames containing the data of ALL cells. data$alldat in a list of data, and data$anno is the summary of the informations.
#' @importFrom memoise memoise
#' @importFrom flowCore read.FCS
#' @export
read_facs_data <- function(datadir, datadate, annofile, anno=NULL){
  mread.fsc = memoise::memoise(flowCore::read.FCS)
  if (is.null(anno)) {
    anno = lapply(datadate, function(d) {
      anno = read.csv(paste(datadir, d, annofile, sep=""), sep=";", stringsAsFactors=FALSE)
      anno = anno[,1:length(colnames(anno))]
      date = gsub("/", "", d)
      anno <- cbind(anno, date = rep(date, times = dim(anno)[1]))
      return(anno)
    })
    anno = do.call(rbind, anno)
  }
  # for (i in 1:nrow(anno)){
  results = lapply(1:nrow(anno), function(i){
    a = anno[i,]
    plate = a$plate
    well = a$well
    date = a$date
    #time = a$time
    #gal = a$gal
   # strain = a$strain
    #if (nbcol_in_anno == 6) {rep = a$rep}
    # Reading fcs data
    path = paste(datadir, date, "/", plate, "/dat/", sep="")
    files = list.files(path)
    well_filename = files[match(well, substr(files, nchar(files) - 2, nchar(files)))]
    fcs_filename = paste(path, well_filename, sep = "")
    fcs_orig = mread.fsc(fcs_filename, transformation=FALSE)
    fcs = data.frame(flowCore::exprs(fcs_orig)[,1:3])
    return(fcs)
  })
  # quality control
  foo = lapply(results, function(l){
    list(
      median_fl1 = median(l$FL1.H),
      median_fsc = median(l$FSC.H),
      median_ssc = median(l$SSC.H),
      sd_fl1 = sd(l$FL1.H),
      sd_fsc = sd(l$FSC.H),
      sd_ssc = sd(l$SSC.H),
      nb_cells = length(l$FL1.H)
    )
  })
  foo = do.call(rbind, foo)
  foo = data.frame(lapply(data.frame(foo, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
  anno = cbind(anno, foo)
  return(list(alldat = results, anno = anno ))
}

