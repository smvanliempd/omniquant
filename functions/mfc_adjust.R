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
  ml.file <- dir(data.path, pattern = "markerlynx_export")
  
  # Read and clean markerlynx data
  if (length(ml.file) == 1) {
    
    # get raw markerlynx data
    ml.data.raw <- data.table(read.csv(paste0(data.path,ml.file), sep = "\t",
                                       header = T,
                                       stringsAsFactors = F))
    colnames(ml.data.raw)[1] <- "File.Name"
    ml.data.raw <- ml.data.raw[!(File.Name %in% c("Retention Time","Mass")) ]
    
    # bin raw with sample meta data
    id.cols <- c("File.Name","Sample.Class","Sample.Group","Sample.ID","Injection.ID","Injection.Number","Injection.Replicate")
    ml.data <- merge(x = sample.data[Sample.Class %in% c("QC","Sample","Blank"), ..id.cols], 
                     y = ml.data.raw,
                     by = "File.Name", 
                     all.y = TRUE)
    ml.features   <- na.exclude(str_extract(colnames(ml.data),"[NXP]{1}\\d{1,}.\\d{2}_\\d{1,}.\\d{1,}"))
    
    # get mz and RTs from feature names
    ft.meta <- data.table(Feature = ml.features,
                          RT = as.numeric(str_extract(ml.features,"\\d+\\.\\d{2}")),
                          mz = as.numeric(str_extract(ml.features,"\\d+\\.\\d{4}")),
                          mean_signal = ml.data[ , sapply(.SD, mean,simplify = TRUE), .SDcols = ml.features,],
                          p_null = ml.data[ , sapply(.SD, function(x) sum(x==0)/length(x),simplify = TRUE), .SDcols = ml.features,]
    )
    
    # chuck out features with more than x% missing
    ft.meta <- ft.meta[p_null < 0.3]
    
    # filter ringing signals
    ft.meta <- ft.meta[
      order(RT)
    ][
      # make RT groups with a 0.013' window (arbitrary, should be parameterized)
      , RT_group := as.integer(factor(floor( (RT - RT[1])/0.013)))
    ][
      order(mz)
    ][
      , c("mz_group","mz_rt_group") :={
        # make mz groups with a 1 Da window (arbitrary, should be parameterized)
        mg  <- as.integer(factor(floor(mz-mz[1])))
        
        # combine mz/Rt groups to unique groups
        mtg <- as.integer(as.factor(paste0(mg,"_",RT_group)))
        list(mg, mtg)
      }
    ][
      # get proportions of max signals for features in mz/RT groups
      , p_max := mean_signal/max(mean_signal), by = mz_rt_group
    ][
      order(-p_max)
    ][
      # get number of features in mz/RT group (not used in calcs, just for analysis)
      , n := .N, by = mz_rt_group 
    ][
      # filter ringing features. p_max lower than 5% of max and not more then 0.6 Da from max feature
      , c("mz_diff","flag_ring") := {
        md = mz - mz[1]
        fr =ifelse(p_max < 0.05 & md < 0.6 , 1, 0)
        list(md,fr)
      }, by = mz_rt_group
    ]
    ml.features.clean <- ft.meta[flag_ring == 0, Feature]
    ml.data <- ml.data[ , .SD, .SDcols = c(id.cols, ml.features.clean)]
    
    # to long format
    ml.mfc <- melt(ml.data,
                   measure.vars = ml.features.clean,
                   variable.name = "Feature",
                   variable.factor = FALSE)
    ml.mfc[value == 0, value := NA ]
    
    if (FALSE) { # ml.del.cutoff should be parameterized for instrument type 
                 #  (eg 500 G2, 5000 G2S), for the moment do not use.
      
      # if blank present delete features where one or more values in blank is over 500
      ml.del.cutoff <- 500
      ml.features.delete <- ml.mfc[Sample.Class == "Blank" & value > ml.del.cutoff, unique(Feature)]
      
      # Delete blank samples and features present in the blank
      ml.mfc <- ml.mfc[Sample.Class != "Blank" & !(Feature %in% ml.features.delete)]
    } else {
      ml.mfc <- ml.mfc[Sample.Class != "Blank"]
    }
    
    ## calculate MFC
    # Get reference sample. Selection is based on highest number of features 
    ml.ref <- ml.mfc[!is.na(value) , .N, by = c("File.Name") ]
    ml.ref <- ml.ref[which.max(N)]
    
    # Calculate fold-change values
    ml.mfc[ , FC := value/value[File.Name == ml.ref$File.Name], by = Feature]
    # ml.mfc <- ml.mfc[!is.na(FC)]
    
    # reduce mfc data
    ml.mfc.red <- ml.mfc[!is.na(FC),  .(N = .N, 
                                        MFC = median(FC, na.rm = T),
                                        y_max =  max(ml.mfc$FC, na.rm = T),
                                        Inj_nr = unique(as.integer(Injection.Number)),
                                        Inj_rep = unique(as.integer(Injection.Replicate))), 
                         by = c("File.Name",
                                "Sample.Class",
                                "Sample.Group",
                                "Sample.ID",
                                "Injection.ID")]
    ml.mfc.red[ , MFC_c := MFC/median(MFC)]
    
    # merge raw data with MFC values
    ml.mfc <- merge(ml.mfc, ml.mfc.red[,.(File.Name,MFC_c)], by = "File.Name")
    
    # merge data with mfc values and apply correction
    full.data   <- merge(full.data, ml.mfc.red[,.(File.Name,MFC_c)], by = "File.Name", all.x = T)
    full.data[ , Signal_MFC := ifelse(Analytes %in% analytes.include, Signal_deiso/MFC_c, Signal_deiso)]
    mfc.data <- data.table(merge(sample.data[ , .(File.Name,
                                                  File.Text,
                                                  Sample.Class,
                                                  Sample.ID,
                                                  # Experiment,
                                                  Injection.Number,
                                                  Cell.Count,
                                                  Weight,
                                                  Protein.Content) ],
                                 ml.mfc.red[,.(File.Name,MFC,MFC_c)], by = "File.Name"))
    
    # get number of samples to adjust graph widths accordingly
    n_samples <- nrow(sample.data)
    
    # plot FC values per sample/injection number
    p1 <- ggplot() +
      geom_hline(yintercept = 1, col = "red")+
      geom_boxplot(data = ml.mfc[!is.na(FC)],
                   aes(x = Injection.ID,
                       y = FC/median(ml.mfc.red$MFC)), 
                   alpha = 0.9)+
      geom_text(data = ml.mfc.red, 
                aes(x = Injection.ID, 
                    y = y_max, 
                    label = N),
                angle = 90,
                hjust = 0,
                vjust = 0.5,
                size = 3,
                nudge_y = 0.1)+
      scale_y_log10()+
      labs(y = "FC (scaled on median of MFC)")+
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
            axis.title.x = element_blank())
    ggsave(paste0(project.path,"/graphs/fc_plots.png"),
           plot = p1,
           dpi = 300, 
           height= 8, 
           width = 0.18 * n_samples + 1)
    
    p2 <- ggplot(ml.mfc.red, 
                 aes(x = Inj_nr,
                     y = MFC_c,
                     col = Sample.Class )) +
      geom_smooth(method = "lm", 
                  formula = 'y ~ x',
                  col = "black", 
                  lty = 2, 
                  linewidth = .5)+
      geom_point() +
      scale_color_manual(values = c("red","black"))+
      scale_x_continuous(breaks = seq(0,max(ml.mfc.red$Inj_nr), by = 10)) +
      ylim(0,NA)+
      theme_bw()
    ggsave(paste0(project.path,"/graphs/mfc_vs_inj-nr.png"),
           plot = p2,
           dpi = 300,
           height= 5, 
           width = 0.04*n_samples + 2)
    
    
    p3 <- ggplot(ml.mfc.red[Sample.Class == "Sample"], 
                 aes(x = Sample.ID, 
                     y = MFC_c,
                     col = factor(Inj_rep))) +
      geom_hline(yintercept = 1, lty = 2)+
      geom_hline(yintercept = 0)+
      geom_point(alpha = 0.8) +
      scale_color_manual(name = "Injection\nnumber",
                         values = c("blue","orangered","darkgreen"))+
      ylim(0,NA)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
            axis.title.x = element_blank())
    ggsave(paste0(project.path,"/graphs/mfc_vs_sample-id.png"),
          plot = p3,
          dpi = 300, 
          height= 5, 
          width = 0.04 * n_samples + 2)
    
    p4 <- ggplot(ml.mfc.red[Sample.Class == "Sample"], 
                 aes(
                   x = Inj_nr,
                   y = MFC_c,
                   col = factor(Inj_rep)
                 ) ) +
      geom_point() +
      geom_line(aes(group = Sample.ID), 
                col = "grey50")+
      scale_color_manual(values = c("blue","orangered"))+
      scale_x_continuous(breaks = seq(0,max(ml.mfc.red$Inj_nr), by = 10)) +
      ylim(0,NA)+
      theme_bw() 
    ggsave(paste0(project.path,"/graphs/mfc_vs_inj-nr_sample-id.png"),
           plot = p4,
           dpi = 300,
           height= 5, 
           width = 0.04 * n_samples + 2)
    
    
    # out
    dat$data <- full.data
    dat$mfc  <- list(done = TRUE,
                     raw = ml.mfc,
                     data = mfc.data,
                     features = ft.meta,
                     plot_fc = p1,
                     plot_injnr = p2, 
                     plot_id = p3,
                     plot_injnr_smplid = p4,
                     included = analytes.include )
    
  } else {
    
    dat$data[ , MFC_c := 1]
    dat$data[ , Signal_MFC := dat$data$Signal_deiso]
    dat$mfc$done <- FALSE
    
  }
  
  return(dat)
  
}

