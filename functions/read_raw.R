# Function to read raw exported TargetLynx files.
# This function can be used as stand-alone without meta data (sample.data) or
#  within the data.prep() function.

read.raw <- function(data.path, raw.files, sample.data) {
  
  # check if sample meta data is available
  no.meta <- missing(sample.data)
  
  # Cycle trough raw files
  dat <- sapply( raw.files , function(f) { 
    
    # Get raw text file
    ql.data.raw <- read.csv(paste0(data.path,f), sep = "\t", header = FALSE, stringsAsFactors = F)
    
    # Extract metabolite names, metabolite number
    target.string <-  "Compound\\s[0-9]{1,}:\\s{1,}"
    analytes <- str_replace(ql.data.raw$V1[grepl(target.string,ql.data.raw$V1)],target.string,"")
    
    ## Tidy up TargetLynx data
    # coerce sample numbers to integers and text fields to NA 
    row.select <- suppressWarnings(as.integer(ql.data.raw$V1))
    
    # get the maximum sample numbers to get the total number of samples
    n.samples <- max(row.select,na.rm = T)
    
    # Select rows
    row.select <- !is.na(row.select)
    ql.data <- ql.data.raw[row.select,-1]
    
    # Set proper column names
    col.names <- c(ql.data.raw[grep("Name",ql.data.raw$V3)[1],-1])
    colnames(ql.data)      <- make.names(col.names)
    colnames(ql.data)[1:2] <- c("Order","File.Name")
    
    # Set data types
    ql.data$Order              <- as.integer(ql.data$Order)
    ql.data$RT                 <- as.numeric(ql.data$RT)
    ql.data$Area               <- as.numeric(ql.data$Area)
    ql.data$Peak.Start.Height  <- as.numeric(ql.data$Peak.Start.Height)
    ql.data$Height             <- as.numeric(ql.data$Height)
    ql.data$S.N                <- as.numeric(ql.data$S.N)
    ql.data$Analyte            <- as.vector(sapply(analytes, function(m)  rep(m,n.samples) ) )
    
    # Transform to data.table
    ql.data <- data.table(ql.data)
    
    ## If sample meta data is available merge it with raw data
    if(no.meta == T) {
      
      full.data <- ql.data
      
    } else {
      
      # Check if all metadata file names are present in the quanlynx files
      files.md <- sample.data$File.Name
      files.ql <- unique(ql.data$File.Name)
      if(sum(files.md %in% files.ql, na.rm = T) < length(files.md) ) {
        stop( paste0("\n\n\n\n\nSome samples are not integrated and are missing\nfrom the QuanLynx results.\n\n\n\n\n\n") )
      }
      
      # include samples included in the meta data only 
      full.data <- merge(sample.data, ql.data, by = "File.Name", all.x = T) 
      full.data <- full.data[ order( Analyte, File.Name ) ]
      
    }
    
    return(full.data)
    
  },USE.NAMES = T , simplify = F)
  
  return(dat)
  
}
