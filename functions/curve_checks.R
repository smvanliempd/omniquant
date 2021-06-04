# Heuristic check if there is systematic bias between the 
#  first and the second calibration curve after MFC correction.

curve.checks <- function( dat ) {
  
  # Get data
  project.path <- dat$project_path
  full.data    <- dat$data
  calibrants   <- dat$calibrants
  n.calibrants <- length(calibrants) 
  
  if (n.calibrants > 0 ){
    
    # Note that it does NOT use QC adjusted values for the curve. 
    #  This is OK because extrapolation from a polynomial of order 2 or higher 
    #  can go terribly wrong.
    logical.test <- expression(Analytes %in% calibrants & Sample.Class == "Curve" )
    full.data[ eval(logical.test) , Signal_adj_cal := Signal_MFC]
    
    # Calculate signed relative standard deviation (snRSD) between first and second curves
    full.data[  eval(logical.test),
                snRSD := sign(Signal_adj_cal[Injection.Replicate == 1] - Signal_adj_cal[Injection.Replicate == 2]) * sd(Signal_adj_cal)/mean(Signal_adj_cal)
                ,  by = c("Analytes", "Sample.ID") ]
    
    # less than 10% of the density of the snRSD distribution should be
    #  below (curve 1  < curve 2) or above (curve 1 > curve 2) zero (arbitrary)
    #  in order to not assign bias.
    pcurve.tol <- 0.1 #0.35
    logical.test <- expression(Analytes %in% calibrants & Sample.Class == "Curve" & Curve.Concentrations > 0.1)
    cb.data <- full.data[ eval(logical.test) ,{
      cbm <- mean(snRSD, na.rm = T)
      cbs <- sd(snRSD, na.rm = T)
      cbp <- if (cbm < 0) { pnorm(cbm, 0, cbs )} else { 1 - pnorm(cbm, 0, cbs ) }
      cbf <- if (cbp < pcurve.tol) {"Biased"} else {"OK"}
      cbmax <- max(snRSD)
      list(cbm,cbs,cbp,cbf,cbmax)
    }
    ,  by = Analytes ]
    
    # Select calibrants of which the first curve is biased compared to second curve
    cb.analytes <- cb.data[cbf == "Biased", Analytes]
    
    # Delete 2nd curve when curves are biased after MFC adjustment
    full.data[ Analytes %in% cb.analytes & Injection.Replicate == 2, Signal_adj_cal := NA  ]
    
    # Plot curve bias
    p1 <- ggplot() + 
      geom_point(data = full.data[eval(logical.test)],
                 aes(x = factor(Curve.Concentrations), y = snRSD) ) +
      geom_hline(yintercept = 0, lty = 2, col = "red") +
      geom_text(data = cb.data,
                aes(x = 1, y= 1.2 * max(cbmax), label = cbf, col = cbf), 
                show.legend = F, hjust = 0, size = 3 ) +
      scale_color_manual(values = c( "Biased" = "red","OK" = "green" ) ) +
      labs(x = "Curve concentrations", y ="signed RSD", 
           title = "Determine if there is a systematic bias between curves",
           subtitle = "Positive values means curve 1 is higher.") +
      facet_wrap( ~ Analytes, ncol = 6)  +
      theme_bw() +
      theme(axis.text.x  = element_text(angle = 90, hjust = 1,vjust = 0.5) )
    
    ggsave(paste0(project.path, "/graphs/curve_bias.png"), p1, dpi = 300,
           width  = ifelse(n.calibrants <= 6, 1.25 * n.calibrants,  7.5), 
           height = ceiling(n.calibrants/6) * 2)
    
    # Data output
    dat$curves <- list(biased_calibrants = cb.analytes,
                       data_curve_bias   = cb.data,
                       plot_curve_bias  = p1)
  }
  
  dat$data <- full.data
  
  return(dat)
  
}
