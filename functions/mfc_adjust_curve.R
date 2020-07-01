# adjust calibration curves based on linear drift of MFC values 

mfc.adjust.curve <- function( dat ) {
  
  # Get data
  full.data <- dat$data
  mfc.data  <- dat$mfc$data
  analytes.include <- dat$mfc$included
  calibrants <- na.omit(unique(dat$analytes[, Calibrants]))
  
  # Make 1st order linear model from MFC values
  mod <- lm(MFC_c ~ Injection.Number, data = mfc.data)
  coefs <- mod$coefficients
  
  # Extrapolate model to adjust curves
  full.data[, MFC_mod := ifelse(Analytes %in% calibrants & 
                                  Analytes %in% analytes.include & 
                                  Sample.Class == "Curve",
                                coefs[1] + coefs[2] * Injection.Number, 1)]
  full.data[Analytes %in% analytes.include & Sample.Class == "Curve" , Signal_MFC := Area_deiso/MFC_mod]
  
  # Out
  dat$data <- full.data
  dat$calibrants <- calibrants
  
  return(dat)
  
}
