# QC corrections of MFC adjusted data 
qc.adjust <- function( dat, ord = 3 ) {
  
  # get data
  full.data    <- dat$data
  qc.include   <- dat$analytes[Apply_QC_corr == "x" , Analytes ]
  project.path <- dat$project_path
  n.analytes   <- sum(nrow(dat$analytes))
  
  # Scale MFC corrected QC values on maximum values per Analyte
  full.data[Sample.Class == "QC"  , QC_signal_s := {
    max_s <- max(Signal_MFC)
    if(max_s > 0 ) {
      ss <- Signal_MFC/max_s
    } else {
      ss <- NA
    }
    list(ss)
  } , by = Analytes]
  
  # Make QC models based on the non-labeled analyte
  if (length(qc.include) > 0 ) { #check if there are any analytes that need QC correction
    full.data[Analytes %in% qc.include & Analytes == Metabolite, QC_model := { 
      
      # define qc data
      d   <- data.frame(signal = QC_signal_s, Inj = Injection.Number)
      
      # linear model with order depending on the number of QCs available
      o <- 1:ord
      f <- formula(paste0("signal ~ ",paste( paste0("I(Inj^", o, ")"), collapse  = "+") ) )
      m.lin  <- lm(f, data = d)
      c.lin  <- coef(m.lin)
      
      # calculate correction factor
      S.mod <- sapply( c(o, ord + 1)  , function(i) c.lin[i] * Injection.Number^(i - 1), simplify = T  )
      S.mod <- rowSums(S.mod)
      
      list(S.mod)
      
    }, by = Analytes ]
    
    # Bind QC-model of unlabeled metabolite with corresponding labeled one
    qc.mods <- full.data[Analytes %in% qc.include & Analytes == Metabolite,.(File.Name, Metabolite, QC_model)]
    full.data[ ,QC_model := NULL]
    full.data <- merge(full.data, qc.mods, by = c("File.Name", "Metabolite"), all.x = T)
    
  } else {
    
    full.data[Sample.Class == "QC"  , QC_model := as.double(NA)]
    
  }
  
  
  
  full.data[is.na(QC_model), QC_model := 1]
  
  # Adjust for QC drift (QC and Samples only)
  full.data[ Sample.Class %in% c("QC","Sample"), Signal_MFC_QC := Signal_MFC/QC_model] #
  
  # plot QC- drift models 
  qc.plots <- ggplot(full.data[Sample.Class %in% c("QC","Sample")] ) +
    geom_point(aes(x   = Injection.Number,
                   y   = QC_signal_s), na.rm = T) +
    geom_line(aes(x = Injection.Number,
                  y = QC_model), lty = 2, na.rm = T ) +
    geom_hline(yintercept = 0.85, lty = 2, size = 0.3) +
    scale_color_manual(name = "QC groups",values = c("#0029ff", "#ff5700"), na.translate = F ) +
    labs(title = "QC trents",
         subtitle = paste0("If necessary, QC adjustement done with \n",  ord,
                           ifelse(ord == 1, "st",
                                  ifelse(ord == 2, "nd", "rd") ),"-order polynomal."),
         x = "Injection Number",
         y = "QC-signals\n(scaled, max = 1; MFC adj)") +
    ylim(0,NA) +
    facet_wrap( ~ Analytes, ncol = 6) +
    theme_bw()
  ggsave(paste0(project.path, "/graphs/qc_plots.png"),
         qc.plots,
         dpi = 300,
         width = ifelse(n.analytes < 6 , n.analytes * 6/3 , 8),
         height = ceiling(n.analytes / 6) * 2)
  
  # Out
  dat$data <- full.data
  
  return(dat)
  
}
