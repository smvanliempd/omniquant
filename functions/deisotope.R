# Function for subtracting natural isotope signals from artificially introduced isotopes.
# It needs the output from the data.prep() function.

deisotope    <- function( dat ){
  
  # Get all data
  full.dat <- dat$data
  
  # Get 13C0-data (natural abundance)
  dat.null <- full.dat[Analytes == Metabolite, .(File.Name,Metabolite, Area_sub = Area )]
  
  # merge null data with full data
  full.dat <- merge(full.dat , dat.null, by = c("File.Name","Metabolite") , all.y = T)
  
  # subtract natural 13Cn (n>0) from labeled if not NA, otherwise take original area
  full.dat[Analytes != Metabolite , Area_deiso := as.double(ifelse(is.na(Area_sub), Area, Area - Isotope_Corr * Area_sub ))]
  
  # if subtractions lead to negative numbers set area to 0
  full.dat[Area_deiso < 0, Area_deiso:= 0]
  
  # for 13C0 set Area_deiso to the original area
  full.dat[Analytes == Metabolite ,Area_deiso := Area  ]
  
  # out
  dat$data <- full.dat
  
  return(dat)
  
}
