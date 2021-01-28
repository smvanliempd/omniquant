
adjustment.checks <- function(dat){
  
  # Get data
  project.path <- dat$project_path
  full.data    <- dat$data
  n.analytes <- length(dat$analytes$Analytes)
  
  # calculate RSDs for uncorrected and corrected signals in samples
  dc <- full.data[Sample.Class == "Sample", sapply(.SD, function(d) {
    m <- mean(d, na.rm = T)
    s <- sd(d, na.rm = T)
    rsd <- s/m
    list(rsd)
  } ) ,by = c("Analytes","Sample.ID") , .SDcols =  c("Area_deiso", "Signal_MFC", "Signal_MFC_QC") ]
  
  # melt
  dc <- melt(dc,id.vars = c("Sample.ID","Analytes"),value.name = "RSD_sample", variable.name = "Data_type")
  
  # calculate rsd over all samples per analyte
  dc[ , c("RSD_mean","RSD_sd","RSD_min","RSD_max") := {
    m <- mean(RSD_sample, na.rm = T)
    s <- sd(RSD_sample, na.rm = T)
    mn <- m - s
    if(mn < 0) mn <- 0
    mx <- m + s
    list(m,s,mn,mx)
  }, by = c("Analytes", "Data_type") ]
  
  pa <- ggplot(dc, aes(x = Data_type, y = RSD_mean ))+
    geom_point() +
    geom_linerange(aes(ymin = RSD_min, ymax = RSD_max)) +
    ylim(0,NA) +
    labs(x = "Adjustment type", y = "RSD (mean ± sd)",
         title = "RSDs over all samples per analyte") +
    facet_grid(~Analytes) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1))
  ggsave(paste0(project.path,"/graphs/adjust_checks.png"),pa,device = "png",dpi = 300, units = "in",
  height= 4, width = 1*n.analytes + 1)
  dat$adjustements <- pa

  
  }