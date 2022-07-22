# QC corrections of MFC adjusted data 
qc.adjust <- function( dat, ord ) {
  
  # get data
  full.data      <- dat$data
  qc.met.include <- unique(dat$analytes[Apply_QC_corr == "x" , Metabolite ])
  project.path   <- dat$project_path
  n.analytes     <- sum(nrow(dat$analytes))
  
  # Scale MFC corrected QC values on first QC value
  full.data[Sample.Class == "QC",
            QC_signal_s := {
              max_s <- Signal_MFC[Injection.Replicate == 1]
              if(max_s > 0 ) {
                ss <- Signal_MFC/max_s
              } else {
                ss <- as.double(NA)
              }
              list(ss)
            } , by = Analytes]
  
  # Make QC models based on the non-labeled analyte
  #  if there are any analytes that need QC correction
  if (length(qc.met.include) > 0 ) { 
    
    # get analyte with highest signal per metabolite
    d.qc.high <- full.data[
      Sample.Class == "QC" & Metabolite %in% qc.met.include         # get QC data for included metabolites
      ][
        , .(med_sig = median(Signal)), by = .(Metabolite, Analytes) # calculate median of QC signals per analyte
      ][ 
        , qc_rank := rank(-med_sig,), by = Metabolite               # rank on highest signal
      ][
      qc_rank == 1, Analytes                                        # select analytes to use for QC correction
      ]
    
    # make models
    qc.mods <- full.data[Analytes %in% d.qc.high,  { 
      
      # get qc data
      d <- .SD[Sample.Class == "QC", .(signal = QC_signal_s, Inj = Injection.Number)]

      # linear model with order depending on the number of QCs available
      o <- 1:ord
      f <- formula(paste0("signal ~ ", paste(paste0("I(Inj^", o, ")"), collapse = "+") ) )
      m.lin  <- lm(f, data = d)
      c.lin  <- coef(m.lin)
      
      # calculate correction factor
      S.mod <- sapply( c(o, ord + 1), function(i) c.lin[i] * Injection.Number^(i - 1), simplify = T  )
      S.mod <- rowSums(S.mod)

      list(QC_model  = S.mod, 
           File.Name = File.Name)

    }, by = .(Metabolite,Analytes) ]
    
    # merge qc models
    full.data <- merge(full.data, qc.mods[,.(File.Name, Metabolite, QC_model)],
                       by = c("Metabolite", "File.Name"),
                       all = T)
    
  } else {
    
    full.data[Sample.Class == "QC"  , QC_model := as.double(NA)]
    
  }
  full.data[is.na(QC_model), QC_model := 1]
  
  # Adjust for QC drift (QC and Samples only)
  full.data[Sample.Class %in% c("QC","Sample"), Signal_MFC_QC := Signal_MFC/QC_model] #
  
  # plot QC- drift models
  qc.dat <- full.data[Sample.Class %in% c("QC","Sample"),
                      .(Injection.Number, QC_signal_s, QC_model, Analytes)]
  qc.dat[ , Analytes_qc := ifelse(Analytes %in% d.qc.high, paste0(Analytes," (*)"), Analytes)]
  qc.plots <- ggplot(qc.dat ) +
    geom_point(
      aes(
        x   = Injection.Number,
        y   = QC_signal_s
      ), 
      na.rm = T
    ) +
    geom_line(
      aes(
        x = Injection.Number,
        y = QC_model
      ), 
      lty = 2,
      na.rm = T 
    ) +
    geom_hline(
      yintercept = 0.85, 
      lty = 2, 
      size = 0.3
    ) +
    labs(
      title = "QC trends",
      subtitle = paste0("If necessary, QC adjustement done with\n",  ord,
                       "-order polynomal. Reference analytes\nare indicated witn (*)"),
      x = "Injection Number",
      y = "QC-signals\n(MFC adjusted, scaled, 1st = 1)"
    ) +
    ylim(0,2) +
    facet_wrap( ~ Analytes_qc, ncol = 6) +
    theme_bw()
  ggsave(paste0(project.path, "/graphs/qc_plots.png"),
         qc.plots,
         dpi = 300,
         width = ifelse(n.analytes < 6 , n.analytes * 6/3 + 2 , 10),
         height = ceiling(n.analytes / 6) * 2 + 2)
  
  # Out
  dat$data <- full.data
  
  return(dat)
  
}
