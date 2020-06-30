# Calculate median fold-change values from markerlynx data.

mfc.adjust   <- function( dat ) { 
  
  # Read data and meta data
  data.path    <- dat$data_path
  project.path <- dat$project_path
  full.data    <- dat$data
  analyte.data <- dat$analytes
  analytes.include <- analyte.data[Apply_MFC_corr == "x", Analytes]
  sample.data <- dat$samples
  colnames(sample.data) <- make.names(colnames(sample.data))
  sample.data <- sample.data[Delete != "ALL" | is.na(Delete)]
  
  # Read and clean markerlynx data
  ml.file <- dir(data.path,pattern = "markerlynx_export")
  ml.data.raw <- data.table(read.csv(paste0(data.path,ml.file), sep = "\t",
                                     header = T,
                                     stringsAsFactors = F))
  colnames(ml.data.raw)[1] <- "File.Name"
  ml.data <- merge(sample.data[Sample.Class %in% c("QC","Sample","Blank"), .(File.Name,Sample.Class)], ml.data.raw, by = "File.Name", all.x = T )
  ml.data[ml.data == 0] <- NA
  ml.features   <- na.exclude(str_extract(colnames(ml.data),"[NXP]{1}[0-9]{1,}.[0-9]{2}_[0-9]{1,}.[0-9]{1,}"))
  
  # to long format
  ml.mfc <- melt(ml.data,measure.vars = ml.features,variable.name = "Feature",variable.factor = F)
  
  # if blank present delete features where one or more values in blank is over 500
  ml.del.cutoff <- 500
  ml.features.delete <- ml.mfc[Sample.Class == "Blank" & value > ml.del.cutoff, unique(Feature)]
  
  # calculate MFC
  ml.mfc <- ml.mfc[Sample.Class != "Blank" & !(Feature %in% ml.features.delete)]
  ml.ref <- ml.mfc[is.na(value) , .N, by = c("File.Name") ]
  ml.ref <- ml.ref[which.max(N)]
  ml.mfc[ , FC := value/value[File.Name == ml.ref$File.Name], by = Feature]
  
  ml.mfc <- ml.mfc[!is.na(FC)]
  ml.mfc <- merge(sample.data, ml.mfc, by = "File.Name")
  ml.mfc.red <- ml.mfc[,  .(N = .N, MFC = median(FC),y_max =  max(ml.mfc$FC)), by = c("File.Name","Injection.ID")]
  ml.mfc.red[ , MFC_c := MFC/median(MFC)]
  
  # plot FC values per s
  p <- ggplot() +
    geom_hline(yintercept = 1, col = "red")+
    geom_boxplot(data = ml.mfc, aes(x = Injection.ID, y = FC/median(ml.mfc.red$MFC)), alpha = 0.9)+
    geom_text(data = ml.mfc.red, aes(x = Injection.ID, y = y_max, label =N),angle = 90,hjust = 0,vjust = 0.5,size =3,nudge_y = 0.1)+
    scale_y_log10()+
    labs(y = "FC (scaled on median of MFC)")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.title.x = element_blank())
  n_samples <- nrow(sample.data)
  ggsave(paste0(project.path,"/graphs/fc_plots.png"),device = "png",dpi = 300, units = "in",
         height= 8, width = 0.18*n_samples + 1)
  
  # merge data with mfc values and apply correction
  full.data   <- merge(full.data, ml.mfc.red[,.(File.Name,MFC_c)], by = "File.Name", all.x = T)
  full.data[ , Signal_MFC := ifelse(Analytes %in% analytes.include, Area_deiso/MFC_c, Area_deiso)]
  mfc.data <- data.table(merge(sample.data[ , .(File.Name,Sample.Class,Sample.ID,Experiment,Injection.Number,Incubation.Time) ],
                               ml.mfc.red[,.(File.Name,MFC,MFC_c)], by = "File.Name"))
  
  # out
  dat$data <-  full.data
  dat$mfc  <-  list(data = mfc.data, plot = p ,included = analytes.include )
  
  return(dat)
  
}
