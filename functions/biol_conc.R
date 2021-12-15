# function to calculate biological concentrations from sample concentrations
#  and to 

biol.conc    <- function( dat ) {
  
  # Get data
  n.calibrants <- length(dat$calibrants)
  full.dat <- dat$data
  assay.data <- dat$assay_pars
  vol_hom <- assay.data[Metric == "Homogenization_volume", Value] # µL
  vol_sv  <- assay.data[Metric == "Speedvac_volume", Value]       # µL
  vol_res <- assay.data[Metric == "Resuspention_volume", Value]   # µL
  cell_nr <- assay.data[Metric == "Cell_count_per_well", Value]
  
  #  calculate loss due to sample preparation
  loss_fact <- vol_hom/vol_sv
  
  if (n.calibrants > 0) {
    
    # calculate amount from original sample
    full.dat[ Sample.Class == "Sample", Sample_amount := (loss_fact  * vol_res * Conc_raw_adj)] # pmol
    
    # calculate amount per tissue weight depending on whether values are MFC adjusted or not.
    mets.mfc <- dat$mfc$included
    full.dat[ Sample.Class == "Sample", Amount_mg_tissue := ifelse(Analytes %in% mets.mfc, 
                                                                   Sample_amount/median(Weight),
                                                                   Sample_amount/Weight)]
    
    # calculate concentrations per amount of cells
    full.dat[ Sample.Class == "Sample" , Conc_biol := (loss_fact  * vol_res * Conc_raw_adj)/cell_nr] #nmol/million cells
  }
  
  # correct signals for different resuspension and speedvac volumes in each assay
  full.dat[ Sample.Class == "Sample" , Signal_assay_norm := (loss_fact  * vol_res * Signal_MFC_QC)/100 ]
  
  # scale signals on median of unlabeled metabolite per metabolite
  full.dat[ Sample.Class == "Sample" , Signal_median_norm := Signal_assay_norm/.SD[ Analytes == Metabolite , median(Signal_assay_norm, na.rm = T)  ] , by = Metabolite]
  
  # scale signals on tissue weight
  full.dat[ Sample.Class == "Sample" , Signal_weight_norm := Signal_deiso/Weight]
  
  # sort
  setkeyv(x = full.dat, cols= c("Metabolite", "Analytes"))
  
  # out
  dat$data <- full.dat
  return(dat)
  
}

