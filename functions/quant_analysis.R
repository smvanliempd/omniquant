quant.analysis <- function(project.dirs) {
  l <- length(project.dirs)
  i <- 1
  d.all <- sapply(project.dirs, function(project) {
    cat(paste0("Processing project ", i ," of ", l,"...\n",
               "Path: ",getwd(),project,
               "\n Data prepartion"))
    d1 <- data.prep(project)
    cat("...done!")
    cat(paste0("\n Deisotope labelled signals"))
    d2 <- deisotope(d1)
    cat("...done!")
    cat(paste0("\n MFC adjustments"))
    d3 <- mfc.adjust(d2)
    cat("...done!")
    cat(paste0("\n Calibration curve adjustments based on MFC values"))
    d4 <- mfc.adjust.curve(d3)
    cat("...done!")
    cat(paste0("\n QC adjustments"))
    d5 <- qc.adjust(d4,  ord = 2)
    cat("...done!")
    cat(paste0("\n Check calibration curves"))
    d6 <- curve.checks(d5)
    cat("...done!")
    cat(paste0("\n Quantifications"))
    d7 <- quantify(d6)
    cat("...done!")
    cat(paste0("\n Calculate biological concentrtions and normalized values"))
    d8 <- biol.conc(d7)
    cat("...done!\n___________\n\n\n")
    i <<- i + 1
    return(d8)
  }, USE.NAMES = T,simplify = F )
  return(d.all)
}