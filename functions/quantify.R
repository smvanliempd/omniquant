quantify     <- function( dat ) {
  
  # dat <- readRDS("d6.rds")
  # Declare data input
  project.path  <- dat$project_path
  full.data     <- dat$data
  calibrants    <- dat$calibrants
  analyte.data  <- dat$analytes[ !is.na(Calibrants), .(Analytes,Calibrants) ]
  quant.targets <- analyte.data$Analytes
  
  if (length(calibrants) > 0){
    # replace adjusted curve values for calibration targets with those of the calibrants
    full.data <- merge(full.data, analyte.data, by = "Analytes", all.x = T)
    cal.data  <- full.data[Sample.Class == "Curve" & Analytes %in% calibrants, .(File.Name,Analytes,Signal_adj_cal) ]
    setnames(cal.data, "Signal_adj_cal","Dummy")
    full.data <- merge(full.data, cal.data , by.x = c("File.Name", "Calibrants"), by.y = c("File.Name", "Analytes"), all.x = T)
    full.data[ , Signal_adj_cal := Dummy]
    full.data[ , Dummy := NULL]
    
    # Define curve baseline
    full.data[Sample.Class == "Curve" & Analytes %in% quant.targets ,
              Baseline := .SD[ Curve.Concentrations == 0, mean(Signal_adj_cal,na.rm = T)  ] , 
              by = c("Analytes", "Injection.Replicate")  ]
    
    # Subtract baseline from curves and delete signals below 3 sd's from the baseline
    full.data[ , c("Curve_adj" , "Curve_min_adj") := {
      ai   <- Signal_adj_cal - Baseline
      mini <- 3 * Baseline
      ai   <- ifelse(ai <= mini, NA, ai)
      list(ai, mini )
    }   ]
    
    ## Focused calibration
    # Calculate median signal for each metabolite (was Analyte )
    full.data[  ,Signal_med := .SD[Sample.Class == "Sample", median(Signal_MFC_QC, na.rm = T) ], by = Metabolite]
    full.data[ Sample.Class == "Curve", idx_focus  := ifelse(Curve_adj > Signal_med | is.na(Curve_adj) | Delete %in% c("curve"), F, T )  ]
    curve.tol <- 0.15
    full.data[Analytes %in% quant.targets & !is.na(Curve_adj) , 
              c("quant_flag",
                "curve_mod", 
                "dev",
                "conc_foc",
                "conc_blq", 
                "included", 
                "alpha", 
                "beta") := {
                  
                  # prepare data
                  conc  <- Curve.Concentrations
                  sig   <- Curve_adj
                  repl  <- Injection.Replicate
                  halt  <- list(quant_flag = F,
                                curve_mod = as.double(NA),
                                dev = as.double(NA),
                                conc_focus = as.double(NA),
                                conc_blq =  as.double(NA),
                                included = "no",
                                alpha = as.double(NA),
                                beta = as.double(NA))
                  
                  if ( length(unique(conc)) >= 3 ) { # check if there are enough unique concentrations (more than 3)
                    
                    # determine focus concentration, based on median signal in the samples
                    conc_focus <- max(ifelse(idx_focus , conc, min(conc) ) )
                    row_focus  <- which(conc == conc_focus)
                    
                    # make all possible 4-point models
                    l_conc   <- length(conc)
                    d_combn  <- combn(1:l_conc, ifelse(l_conc >= 4, 4, l_conc))
                    d_coeffs <- apply(d_combn, 2, function(i) {
                      mod   <- lm(log10(sig[i]) ~ log10(conc[i]))
                      mod$coefficients
                    } )
                    
                    # Calculate concentrations and deviation from theoretical
                    d_dev <- sapply(1:ncol(d_coeffs), function(i) {
                      conc_calc <- ( sig/10^d_coeffs[1,i] )^(1/d_coeffs[2,i])
                      (conc_calc - conc)/conc
                    } )
                    
                    # Excluding all devs > 0.15
                    d_qual  <- d_dev
                    d_qual[abs(d_dev) > curve.tol ] <- NA
                    
                    # Check if at least one focus concentration deviation is less than 15%
                    d_focus  <- !is.na(d_qual[row_focus,  ] )
                    d_focus  <- if(is.vector(d_focus) ) { sum(d_focus) } else {colSums(d_focus,na.rm = T)}
                    d_focus  <- ifelse(d_focus == 0 , F, T)
                    
                    if ( sum(d_focus  ) >= 1 ) {
                      
                      # delete models where focus concentration deviation is over 15%
                      d_coeffs <- matrix( d_coeffs[,d_focus], nrow = nrow(d_coeffs) )
                      d_dev    <- matrix(    d_dev[,d_focus], nrow = nrow(d_dev)    )
                      d_qual   <- matrix(   d_qual[,d_focus], nrow = nrow(d_qual)   )
                      
                      # Check which model includes the most unique concentrations after deleting dev<0.15 ( >0.15? )
                      d_conc  <- apply(d_qual, 2, function(v) list( unique( conc[ !is.na(v) ] ) ) )
                      d_conc  <- sapply(d_conc,   function(v) length( unlist(v) ) , simplify = T )
                      
                      if ( max(d_conc) >= 3 ) { # check if there are enough unique concentrations left after deleting dev<0.15 ( >0.15? )
                        
                        col_max <- which( d_conc ==  max(d_conc) )
                        
                        # check/select which model with maximum unique concentrations include concentrations that are closest to focal conc
                        d_dist   <- sapply(col_max, function(c) sqrt(na.omit(sum((d_combn[,c] - conc_focus)^2))))
                        col_max  <- col_max[which.min(d_dist)]
                        coef_sel <- d_coeffs[,col_max]
                        
                        # calculate modeled signal and return relevant data
                        alpha     <- coef_sel[1]
                        beta      <- coef_sel[2]
                        curve_mod <- 10^alpha * Curve.Concentrations^beta
                        conc_blq  <- ( Curve_min_adj/(10^alpha) )^(1/beta)
                        dev       <- d_dev[,col_max]
                        included  <- ifelse(!is.na(d_qual[, col_max]), "yes", "no")
                        
                        
                        list(quant_flag = T, curve_mod, dev, conc_focus , conc_blq ,included, alpha, beta)
                        
                      } else {
                        
                        halt
                        
                      }
                      
                    } else {
                      
                      halt
                      
                    }
                    
                  } else {
                    
                    halt
                    
                  }
                  
                }, by = Analytes ] 
    
    # update final quantification targets based on the validity of the calibration curve; bad curve --> not included
    n.quant.targs <- length(quant.targets)
    
    # Extrapolate curves
    full.data[Analytes %in% quant.targets & Curve.Concentrations != 0, curve_mod  := {
      
      alpha     <- mean(alpha, na.rm = T)
      beta      <- mean(beta, na.rm = T)
      curve_mod <- 10^alpha * Curve.Concentrations^beta
      
      list(curve_mod)
    }, by = Analytes]
    
    # Determine sample signal extremes
    full.data[ , Signal_min := .SD[Sample.Class == "Sample" & Signal_MFC_QC >0,  min(Signal_MFC_QC, na.rm = T)], by = Metabolite ]#Analytes 
    full.data[ , Signal_max := .SD[Sample.Class == "Sample",  max(Signal_MFC_QC, na.rm = T)], by = Metabolite ] #Analytes 
    
    # Make curve plots for ALL INITIAL calibrants (also if curve is not possible)
    conc <- unique(full.data[Sample.Class == "Curve" & Curve.Concentrations != 0, Curve.Concentrations])
    max.conc <-full.data[,max(Curve.Concentrations, na.rm = T)]
    p1 <- ggplot(full.data[Sample.Class == "Curve" & Analytes %in% calibrants & Curve.Concentrations != 0, ], 
                 aes(x = Curve.Concentrations )) +
      geom_point(aes(y = Curve_adj, col = factor(Injection.Replicate) , alpha = included), na.rm = T) +
      geom_line(aes(y = curve_mod)) +
      geom_ribbon( aes(ymin = Signal_min, ymax = Signal_max  , fill = "Sample space") , alpha = 0.1, show.legend = T  ) +
      geom_line( aes( y = Signal_med ), col = "#0029ff" ,alpha = 0.3 ) +
      scale_fill_manual(name = NULL, values = "#0029ff") +
      scale_color_manual(name = NULL, values = c("#0029ff", "#ff5700"), labels = c("curve 1", "curve 2") ) +
      scale_alpha_manual(name = "Included", values = c(no = 0.2, yes = 0.8), labels = c("no", "yes"), na.translate = F  ) +
      scale_x_log10(breaks = conc, limits = c(0.01,max(max.conc)),
                    labels =format(conc, drop0trailing = T ) ) +
      scale_y_log10() +
      labs(title = "Calibration curves",
           subtitle = paste0("Deviation of included points is less than ±",round(100*curve.tol), "%"),
           x = expression(paste("Concentration (", mu, "M)" ) ),
           y = "Adjusted signal") +
      facet_wrap( ~ Analytes, ncol = 4, scales = "free_y" ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0(project.path, "/graphs/curves.png"), p1, dpi = 300,
           width  = 10, 
           height = 6 )
    
    # delete data for calibrant-only compounds e.g Methionine-SL
    full.data <- full.data[ !(Analytes %in% dat$analytes[Calibrant_Only == "x", Analytes]) ]
    
    # Calculate raw concentrations
    full.data[Analytes %in% quant.targets , c("Conc_raw", "Conc_raw_adj", "alpha", "beta")  := {
      alpha  <- mean(alpha, na.rm = T)
      beta   <- mean(beta, na.rm = T)
      conc <- ( Signal_MFC_QC/(10^alpha) )^(1/beta)
      conc_adj <- conc # * dil.homog
      list(conc, conc_adj, alpha, beta)
    }, by = Analytes]
    
    full.data[Analytes %in% quant.targets, Conc_raw_repl := {
      m <- mean(Conc_raw_adj, na.rm = T)
      m[ is.nan(m) ] <- NA
      list(m)
    }, by = c("Analytes","Sample.ID") ]
    
    # Set quantification labels based on QC/MFC-only corrected signals
    full.data[Analytes %in% quant.targets , c("Curve_min_adj","Curve_low_adj", "Curve_hgh_adj") := {
      sig_min  <- mean(Curve_min_adj, na.rm = T)
      conc_low <- min(.SD[included == "yes", Curve.Concentrations ])
      sig_low  <- 10^alpha * conc_low^beta
      conc_hgh <- max(.SD[included == "yes", Curve.Concentrations ])
      sig_hgh  <- 10^alpha * conc_hgh^beta
      list(sig_min ,sig_low,sig_hgh)
    }, by = Analytes ]
    
    full.data[Analytes %in% quant.targets , Quant_lab := {
      blq  <- mean(Curve_min_adj, na.rm = T ) # below LQ
      epol <- mean(Curve_low_adj, na.rm = T ) # extrapolated
      alq  <- mean(Curve_hgh_adj, na.rm = T ) # above LQ
      lab  <- as.character(ifelse(Signal_MFC_QC < blq,"BLQ"  , 
                                  ifelse(Signal_MFC_QC < epol ,"ExPol low" , 
                                         ifelse(Signal_MFC_QC > alq  , "ExPol high" , "OK" ) ) ) )
      list(lab)
    }, by = Analytes]
    
    full.data[Analytes %in% quant.targets, 
              Labs_repl :=  ifelse("BLQ" %in% Quant_lab, "BLQ",
                                   ifelse("ExPol high" %in% Quant_lab, "ExPol high", 
                                          ifelse("ExPol low" %in% Quant_lab, "ExPol low", "OK"))), by = c("Analytes", "Sample.ID") ]
    
    dat$curves$plot_curves <- p1
  }
  
  # Data output
  dat$data <- full.data
  
  
  return(dat)
  
}
