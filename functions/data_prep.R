# Function to load and tidy raw LCMS data

data.prep    <- function( project.path ) {
  
  # get raw data and meta data file names
  data.path <- paste0(project.path,"/data/")
  raw.files <- dir(data.path,pattern = "SET_") 
  meta.file <- dir(data.path,pattern = "^meta_data")
  
  # Read meta data
  # Assay data
  assay.data <- read_excel(paste0(data.path,meta.file) , sheet = "Assay_conditions" )
  assay.data <- data.table(assay.data)
  
  # Sample data
  sample.data  <- read_excel(paste0(data.path,meta.file) , sheet = "Samples", )
  c.names <- make.names(colnames(sample.data))
  colnames(sample.data) <- c.names
  sample.data <- data.table(sample.data)
  
  #Delete superfluous columns in sample data
  supcols <- c("Inlet.File","MS.File")
  if (all(supcols %in% colnames(sample.data))) {
    sample.data[,get("supcols") := NULL]
  }
  
  # Analyte data
  analyte.data <- data.table(read_excel(paste0(data.path,meta.file) , sheet = "Compounds" ))
  
  # Parse all raw data to list of data.tables
  dat <- read.raw(data.path = data.path,
                  raw.files = raw.files,
                  sample.data = sample.data)
  
  # Bind all integration data sets 
  dat <- do.call(rbind, dat)
  
  # Check inconsistencies between integrations and meta data
  # If names in TargetLynx export and sample.data differ, throw error.
  analytes.md <- na.exclude(analyte.data$QL_names)
  analytes    <- unique(dat$Analyte)
  if( !setequal(analytes,analytes.md ) ) {
    stop(
      paste0("\nThese analytes are not in the MASS SPEC data:\n", paste(analytes.md[!(analytes.md %in% analytes)],collapse = "\n"),
             "\n\nThese analytes are not in the META DATA data:\n",paste(analytes[!(analytes %in% analytes.md)],collapse = "\n")
      ) 
    )
  }
  
  # Merge QL data with compound data and set Analyte levels based on meta data order
  dat <- merge(analyte.data[ , .(QL_names,
                                 Analytes,
                                 Label_Comp,
                                 Cycle,
                                 Not_Use,
                                 Label_enrichment,
                                 Label_ID,
                                 Metabolite,
                                 Metabolite_Class,
                                 Isotope_Corr)] , dat  , by.x = "QL_names", by.y = "Analyte", all.y = T) #
  
  # delete analytes
  dat <- dat[is.na(Not_Use)]
  dat$Not_Use <- NULL
  analyte.data <- analyte.data[is.na(Not_Use)]
  analyte.data$Not_Use <- NULL
  
  # Delete samples for all metabolites
  dat <- dat[Delete != "ALL" | is.na(Delete)]
  
  # Delete samples for specifically defined Analytes
  dat[!is.na(Delete)  , Delete := ifelse(Analytes %in% c(str_extract_all(Delete, "[A-Za-z0-9\\-\\_\\s\\(\\)\\:]{1,}",simplify = T) ), "x", NA ) , by = c("File.Name") ]
  dat <- dat[Delete != "x" | is.na(Delete)]
  
  # get polarity and assay from file location
  dat$Polarity <- str_extract(project.path , "POS|NEG")
  dat$Assay <- project.path
  
  # out
  dat <- list(data         = dat,
              data_path    = data.path,
              analytes     = analyte.data,
              samples      = sample.data,
              assay_pars   = assay.data,
              project_path = project.path,
              meta_file    = meta.file)
  
  return(dat)
  
}
