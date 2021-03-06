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

# filter to keep only 23andMe genotype files
phenotypes <- phenotypes %>% filter(genotype == "23andme")

# keep only first genotype file per ID
phenotypes <- phenotypes[!duplicated(phenotypes$user_id), ]
rownames(phenotypes) <- NULL

# variable types
# Height: num (cm)
# date_of_birth: num (year)
# chrom_sex: factor (XX, XY)
phenotypes$Height <- as.numeric(phenotypes$Height)
phenotypes$date_of_birth <- as.numeric(phenotypes$date_of_birth)
phenotypes$chrom_sex[phenotypes$chrom_sex == "other"] <- NA
phenotypes$chrom_sex <- as.factor(phenotypes$chrom_sex)

# extract file name list to unzip and delete unnecessary columns
genotype_filename <- phenotypes %>% select(user_id, filename)
phenotypes <- phenotypes %>% select(user_id, date_of_birth, chrom_sex, Height)

# sample 250 genotype files
set.seed(8053)
genotype <- sample(nrow(genotype_filename), size = 250, replace = FALSE)
genotype <- genotype_filename[genotype, ]

# save phenotypes objects to include only 250 selected genotypes
phenotypes_clean <- phenotypes %>% filter(user_id %in% genotype[, "user_id"]) 
#saveRDS(object = phenotypes_clean, file = "phenotypes_clean.RData")

# clean environment
rm(height_dict, libraries, phenotypes, phenotypes_clean, genotype_filename, i)
gc()

# Genotypes ----

# unzip all selected genotype files - 250 selected
# setwd() appropriately
# unzip(zipfile = "opensnp_datadump.current.zip", files = genotype[, "filename"])

# read in all files
# keeping 2 variables: (1) SNP ID and (2) genotype variant
genotypes <- cbind(user_id = genotype[1, "user_id"],
                   read.table(genotype[1, "filename"], 
                              sep = "\t", header = FALSE) %>%
                     # remove mitochrondrial DNA and sex chromosome
                     filter(V2!= "MT" & V2 != "X" & V2 != "Y") %>%
                     # keep only known SNPs in standard databases
                     filter(grepl("rs", V1)) %>%
                     # select SNP id and genotype
                     select(V1, V4) %>%
                     # long to wide format
                     spread(., key = V1, value = V4))

for(file in 2:250) {
  try <- try(bind_rows(genotypes, 
                           cbind(user_id = genotype[file, "user_id"],
                                 read.table(genotype[file, "filename"], 
                                            sep = "\t", header = FALSE) %>%
                                   # remove mitochrondrial DNA and sex chromosome
                                   filter(V2!= "MT" & V2 != "X" & V2 != "Y") %>%
                                   # keep only known SNPs in standard databases
                                   filter(grepl("rs", V1)) %>%
                                   # select SNP id and genotype
                                   select(V1, V4) %>%
                                   # long to wide format
                                   spread(., key = V1, value = V4))))
  
  if (class(try) != "try-error") {
    genotypes <- try
    cat(file, "genotype files have been loaded and combined. \n")
  }
}

# keep only non-missing columns
genotypes_clean <- genotypes[, colSums(is.na(genotypes)) == 0]

# mark -- as NA
genotypes_clean[genotypes_clean == "--"] <- NA
genotypes_clean <- genotypes_clean[, colSums(is.na(genotypes_clean)) == 0]

# remove observation ID 6965 and 2925 (genotype only has one allele)
genotypes_clean <- genotypes_clean %>% filter(user_id != 6965 & user_id != 2925)

# remove observation ID 9262 and 77 due to questionable genotype data
genotypes_clean <- genotypes_clean %>% filter(user_id != 9262 & user_id != 77)

# recode some alleles for consistency
genotypes_clean[genotypes_clean == "TC"] <- "CT"
genotypes_clean[genotypes_clean == "TG"] <- "GT"

# export RData
# saveRDS(object = genotypes_clean, file = "genotypes_clean.RData")

### save at most 5% missing ###
genotypes_5 <- genotypes[, colSums(is.na(genotypes)) < 11]

# mark -- as NA
genotypes_5[genotypes_5 == "--"] <- NA

# remove observation ID 6965 and 2925 (genotype only has one allele)
genotypes_5 <- genotypes_5 %>% filter(user_id != 6965 & user_id != 2925)

# remove observation ID 9262 and 77 due to questionable genotype data
genotypes_5 <- genotypes_5 %>% filter(user_id != 9262 & user_id != 77)

# recode some alleles for consistency
genotypes_5[genotypes_5 == "TC"] <- "CT"
genotypes_5[genotypes_5 == "TG"] <- "GT"


# Full data ----

# merge genotype and phenotype 
data <- merge(phenotypes_clean, genotypes_clean)

# select only users with phenotype data
data <- data %>% filter(!is.na(Height))

# export RData
# saveRDS(object = data, file = "data_clean.RData")

### full data with at most 10% missing ###
phenotypes_clean <- readRDS(file = "phenotypes_clean.RData")
data_5 <- merge(phenotypes_clean, genotypes_5)

# select only users with phenotypes data
data_5 <- data_5 %>% filter(!is.na(Height))
                               
# export RData
saveRDS(object = data_5, file = "data_5.RData")

