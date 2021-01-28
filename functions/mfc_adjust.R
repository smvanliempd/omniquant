# Calculate median fold-change values from markerlynx data.

mfc.adjust   <- function( dat ) { 
  
  # Read data and meta data
  data.path        <- dat$data_path
  project.path     <- dat$project_path
  full.data        <- dat$data
  analyte.data     <- dat$analytes
  analytes.include <- analyte.data[Apply_MFC_corr == "x", Analytes]
  sample.data      <- dat$samples[Delete != "ALL" | is.na(Delete)]
  colnames(sample.data) <- make.names(colnames(sample.data))
  
  # Read and clean markerlynx data
  ml.file <- dir(data.path, pattern = "markerlynx_export")
  ml.data.raw <- data.table(read.csv(paste0(data.path,ml.file), sep = "\t",
                                     header = T,
                                     stringsAsFactors = F))
  colnames(ml.data.raw)[1] <- "File.Name"
  ml.data <- merge(x = sample.data[Sample.Class %in% c("QC","Sample","Blank"), .(File.Name,
                                                                                 Sample.Class,
                                                                                 Sample.Group,
                                                                                 Sample.ID,
                                                                                 Injection.ID,
                                                                                 Injection.Number,
                                                                                 Injection.Replicate)], 
                   y = ml.data.raw, by = "File.Name", 
                   all.y = T )
  # ml.data[ml.data == 0] <- NA
  ml.features   <- na.exclude(str_extract(colnames(ml.data),"[NXP]{1}[0-9]{1,}.[0-9]{2}_[0-9]{1,}.[0-9]{1,}"))
  
  # to long format
  ml.mfc <- melt(ml.data,measure.vars = ml.features,variable.name = "Feature",variable.factor = F)
  ml.mfc[value == 0, value := NA ]
  
  # if blank present delete features where one or more values in blank is over 500
  ml.del.cutoff <- 500
  ml.features.delete <- ml.mfc[Sample.Class == "Blank" & value > ml.del.cutoff, unique(Feature)]
  
  ## calculate MFC
  # Delete blank samples and features present in the blank
  ml.mfc <- ml.mfc[Sample.Class != "Blank" & !(Feature %in% ml.features.delete)]
  
  # Get reference sample. Selection is based on highest number of features 
  ml.ref <- ml.mfc[is.na(value) , .N, by = c("File.Name") ]
  ml.ref <- ml.ref[which.max(N)]
  
  # Calculate fold-change values
  ml.mfc[ , FC := value/value[File.Name == ml.ref$File.Name], by = Feature]
  # ml.mfc <- ml.mfc[!is.na(FC)]
  
  # reduce mfc data
  ml.mfc.red <- ml.mfc[,  .(N = .N, 
                            MFC = median(FC, na.rm = T),
                            y_max =  max(ml.mfc$FC, na.rm = T),
                            Inj_nr = unique(as.integer(Injection.Number)),
                            Inj_rep = unique(as.integer(Injection.Replicate))), 
                       by = c("File.Name","Sample.Class","Sample.Group","Sample.ID","Injection.ID")]
  ml.mfc.red[ , MFC_c := MFC/median(MFC)]
  
  # merge data with mfc values and apply correction
  full.data   <- merge(full.data, ml.mfc.red[,.(File.Name,MFC_c)], by = "File.Name", all.x = T)
  full.data[ , Signal_MFC := ifelse(Analytes %in% analytes.include, Area_deiso/MFC_c, Area_deiso)]
  mfc.data <- data.table(merge(sample.data[ , .(File.Name,File.Text,Sample.Class,Sample.ID,Experiment,Injection.Number,Incubation.Time) ],
                               ml.mfc.red[,.(File.Name,MFC,MFC_c)], by = "File.Name"))
  
  # plot FC values per sample/injection number
  n_samples <- nrow(sample.data)
  p1 <- ggplot() +
    geom_hline(yintercept = 1, col = "red")+
    geom_boxplot(data = ml.mfc, aes(x = Injection.ID, y = FC/median(ml.mfc.red$MFC)), alpha = 0.9)+
    geom_text(data = ml.mfc.red, aes(x = Injection.ID, y = y_max, label =N),angle = 90,hjust = 0,vjust = 0.5,size =3,nudge_y = 0.1)+
    scale_y_log10()+
    labs(y = "FC (scaled on median of MFC)")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.title.x = element_blank())
  ggsave(paste0(project.path,"/graphs/fc_plots.png"),p1,device = "png",dpi = 300, units = "in",
         height= 8, width = 0.18*n_samples + 1)
  
  p2 <- ggplot(ml.mfc.red, aes(x = Inj_nr, y = MFC_c, col = Sample.Class )) +
    geom_smooth(method = "lm", formula = 'y ~ x',col = "black", lty = 2, size = .5)+
    geom_point() +
    scale_color_manual(values = c("red","black"))+
    scale_x_continuous(breaks = seq(0,max(ml.mfc.red$Inj_nr), by = 10)) +
    ylim(0,NA)+
    theme_bw()
  ggsave(paste0(project.path,"/graphs/mfc_vs_inj-nr.png"),p2,device = "png",dpi = 300, units = "in",
         height= 5, width = 0.05*n_samples + 1)
  
  
  p3 <- ggplot(ml.mfc.red[Sample.Class == "Sample"], aes(x = Sample.ID, y = MFC_c, col = factor(Inj_rep))) +
    geom_hline(yintercept = 1, lty = 2)+
    geom_hline(yintercept = 0)+
    geom_point() +
    scale_color_manual(name = "Injection\nnumber",values = c("blue","orangered","darkgreen"))+
    ylim(0,NA)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0(project.path,"/graphs/mfc_vs_sample-id.png"),p3,device = "png",dpi = 300, units = "in",
         height= 5, width = 0.05*n_samples + 1)
  
  # out
  dat$data <- full.data
  dat$mfc  <- list(raw = ml.mfc ,data = mfc.data, plot_fc = p1 ,plot_injnr = p2, plot_id = p3,included = analytes.include )
  
  return(dat)
  
}
