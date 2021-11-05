#########################
## STAT 8053 - Model-X ##
##   Cleaning script   ##
##     Linh Nguyen     ##
##     2021-11-04      ##
#########################

# Meta ----

libraries <- c("tidyverse") # general wrangling

lapply(libraries, library, character = TRUE)

# Phenotypes ----

# read in file and order by ID
phenotypes <- readRDS(file = "phenotypes.RData")

phenotypes <- phenotypes[order(phenotypes$user_id), ]

# delete quotation mark and space in height
phenotypes$Height <- gsub("[[:punct:]]", "", phenotypes$Height)
phenotypes$Height <- gsub(" ", "", phenotypes$Height)

# delete cm in height
phenotypes$Height <- gsub("cm", "", phenotypes$Height)

# clean height values
height_dict <- rbind(
  c("Tall180",                              180),
  c("Average173",                           173),
  c("510",                                  177.8),
  c("51or155",                              155),
  c("58",                                   172.72),
  c("59",                                   175.26),
  c("Average165x180",                       172.5),
  c("Average",                              NA),
  c("56",                                   167.64),
  c("Tall",                                 NA),
  c("62",                                   187.96),
  c("59176",                                176),
  c("55",                                   165.1),
  c("GG",                                   NA),
  c("61",                                   185.42),
  c("54",                                   162.56),
  c("52",                                   157.48),
  c("6",                                    182.88),
  c("53",                                   160.02),
  c("57",                                   170.18),
  c("63",                                   190.5),
  c("511",                                  180.34),
  c("61185",                                185),
  c("585",                                  173.99),
  c("ttcc",                                 NA),
  c("68",                                   203.2),
  c("64",                                   193.04),
  c("565",                                  168.91),
  c("48",                                   142.24),
  c("50",                                   152.40),
  c("60183",                                183),
  c("65",                                   195.58),
  c("50'",                                  152.40),
  c("rs6060371",                            NA),
  c("411",                                  149.86),
  c("66",                                   198.12),
  c("62",                                   187.96),
  c("5912",                                 176.53),
  c("Wellabovetheaverageforwomen",          NA),
  c("5212",                                 158.75),
  c("1m87",                                 187),
  c("5434tallerthanallfemalesinmyfamily",   161.925),
  c("515",                                  156.21),
  c("535",                                  161.29),
  c("6ft0in",                               182.88),
  c("5275atmax",                            159.385),
  c("200",                                  200),
  c("5",                                    152.4),
  c("5â€™6",                                167.64),
  c("64",                                   193.04),
  c("1778",                                 177.8),
  c("17272",                                172.72),
  c("169316",                               169.316),
  c("16256",                                162.56),
  c("163195",                               163.195)
)

for(i in seq(nrow(height_dict))) {
  phenotypes <- phenotypes %>% 
    mutate(Height = ifelse(Height == height_dict[i, 1],
                           # replace with correct value
                           height_dict[i, 2], 
                           # or else return original value
                           Height))
}

phenotypes$Height <- as.numeric(phenotypes$Height)

# keep only first genotype file per ID
phenotypes <- phenotypes[!duplicated(phenotypes$user_id), ]
rownames(phenotypes) <- NULL

# clean genotype_filename
phenotypes$genotype_filename <- 
  # delete everything before first .
  gsub(pattern = "^.*?\\.", "", phenotypes$genotype_filename)


# split everything after . to another column called file
phenotypes <- phenotypes %>% 
  separate(genotype_filename, c("genotype", "file"), sep = "[.]")

# add filename variable to conform to actual .txt file names
# to zip automatically later
phenotypes[is.na(phenotypes)] <- "unknown" # temporary to match file names
phenotypes <- phenotypes %>% 
  mutate(filename = paste0("user", user_id, "_file", file,
                           "_yearofbirth_", date_of_birth,
                           "_sex_", chrom_sex, 
                           ".", genotype, ".txt"))
phenotypes[phenotypes == "unknown"] <- NA # change back to NA

# extract file name list to unzip and delete unnecessary columns
genotype_filename <- pull(phenotypes, filename)
phenotypes <- phenotypes %>% select(user_id, date_of_birth, chrom_sex, Height)

# Genotypes ----

# unzip all relevant genotype files
# unzip(zipfile = "opensnp_datadump.current.zip", files = genotype_filename)

# read in all files
genotypes <- lapply(genotype_filename[1:648], function(x) {
  tmp <- try(read.table(paste(x), sep = "\t", header = FALSE))
  if (!inherits(tmp, 'try-error')) tmp
})