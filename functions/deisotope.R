# Function for subtracting natural isotope signals from artificially introduced isotopes.
# It needs the output from the data.prep() function.

deisotope    <- function( dat ){
  
  # Get all data
  full.dat <- dat$data
  
  # Get 13C0-data (natural abundance)
  dat.null <- full.dat[Analytes == Metabolite, .(File.Name,Metabolite, Signal_sub = Signal )]
  
  # merge null data with full data
  full.dat <- merge(full.dat , dat.null, by = c("File.Name","Metabolite") , all.y = T)
  
  # subtract natural 13Cn (n>0) from labeled if not NA, otherwise take original Signal
  full.dat[Analytes != Metabolite , Signal_deiso := as.double(ifelse(is.na(Signal_sub), Signal, Signal - Isotope_Corr * Signal_sub ))]
  
  # if subtractions lead to negative numbers set Signal to 0
  # full.dat[Signal_deiso < 0, Signal_deiso:= 0]
  full.dat[Signal_deiso < 0, Signal_deiso:= ifelse(Sample.Class == "Sample", NA, 0) ]
  
  # for 13C0 set Signal_deiso to the original Signal
  full.dat[Analytes == Metabolite ,Signal_deiso := Signal  ]
  
  # out
  dat$data <- full.dat
  
  return(dat)
  
}
