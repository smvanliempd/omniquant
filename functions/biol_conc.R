# function to calculate biological concentrations from sample concentrations
#  and to 

biol.conc    <- function( dat ) {
  
  # Get data
  full.dat <- dat$data
  assay.data <- dat$assay_pars
  vol_hom <- assay.data[Metric == "Homogenization_volume", Value] # µL
  vol_sv  <- assay.data[Metric == "Speedvac_volume", Value]       # µL
  vol_res <- assay.data[Metric == "Resuspention_volume", Value]   # µL
  cell_nr <- assay.data[Metric == "Cell_count_per_well", Value]   # million
  
  #  calculate loss due to sample preparation
  loss_fact <- vol_hom/vol_sv
  
  # calculate concentrations per amount of cells
  full.dat[ Sample.Class == "Sample" , Conc_biol := (loss_fact  * vol_res * Conc_raw_adj)/cell_nr] #nmol/million cells
  
  # correct signals for different resuspension and speedvac volumes in each assay
  full.dat[ Sample.Class == "Sample" , Signal_assay_norm := (loss_fact  * vol_res * Signal_MFC_QC)/100 ]
  
  # scale signals on median of unlabeled metabolite per metabolite
  full.dat[ Sample.Class == "Sample" , Signal_median_norm := Signal_assay_norm/.SD[ Analytes == Metabolite , median(Signal_assay_norm, na.rm = T)  ] , by = Metabolite]
  
  # sort
  setkeyv(x = full.dat, cols= c("Metabolite", "Analytes"))
  
  # out
  dat$data <- full.dat
  return(dat)
  
}

