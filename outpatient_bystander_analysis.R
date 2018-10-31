# 1a. SETUP - WORKSPACE ###################################################################

# Load required packages
library(gdata)
library(ggplot2)
library(dplyr)

# 1b. SETUP - DEFINE VARIABLES OF INTEREST ################################################

# Define species of interest
species <- c("Streptococcus_pneumoniae", "Haemophilus_influenzae", "Moraxella_catarrhalis", 
             "Staphylococcus_aureus", "Escherichia_coli", "Pseudomonas_aeruginosa", 
             "Klebsiella_pneumoniae", "Streptococcus_agalactiae", "Streptococcus_pyogenes")

# Define conditions of interest
conditions <- c("acuteSinusitis", "chronicSinusitis", "suppOtitisMedia", "pharyngitis.nostrep", "strep", "uti",
                "bronchitis", "cellulitis", "pneumonia", "viralPneumonia", "viralURI",
                "influenza", "miscBacInf", "allergicRhinitis", "nonSuppOtitisMedia",
                "pneumoniaSP", "pneumoniaKP", "pneumoniaPA", "pneumoniaHI", "pneumoniaSgpA",
                "pneumoniaSgpB", "pneumoniaSA", "pneumoniaEC", "acne")

# Define antibiotics of interest
# "inclab" = any of the explicitly included antibiotics
# "anyab" = any antibiotic, as identified by MULTUM Lexicon system in NAMCS/NHAMCS
antibiotics <- c("AMOXICILLIN", "AMOXICILLIN-CLAVULANATE", "PENICILLIN",
                 "AZITHROMYCIN", "CLINDAMYCIN", "CLARITHROMYCIN",
                 "CIPROFLOXACIN", "LEVOFLOXACIN", "MOXIFLOXACIN",
                 "CEFTRIAXONE", "CEPHALEXIN", "CEFDINIR",
                 "SULFAMETHOXAZOLE-TRIMETHOPRIM", "NITROFURANTOIN", "DOXYCYCLINE", 
                 "TETRACYCLINE", "MINOCYCLINE", "inclab", "anyab")

# Define drug classes
# `drugclasses` vector matches respective class to antibiotics as defined in `antibiotics` vector
drugclasses <- c("PENICILLINS", "PENICILLINS", "PENICILLINS",
                 "MACROLIDES", "MACROLIDES", "MACROLIDES",
                 "FLUOROQUINOLONES", "FLUOROQUINOLONES", "FLUOROQUINOLONES",
                 "CEPHALOSPORINS", "CEPHALOSPORINS", "CEPHALOSPORINS", 
                 "SULFAMETHOXAZOLE-TRIMETHOPRIM", "NITROFURANTOIN", "TETRACYCLINES", 
                 "TETRACYCLINES", "TETRACYCLINES", "inclab", "anyab")
drugclasses.key <- as.data.frame(cbind(antibiotic=antibiotics, drugclass=drugclasses))
drugclasses.unique <- c("PENICILLINS", "MACROLIDES", "FLUOROQUINOLONES",
                        "CEPHALOSPORINS", "SULFAMETHOXAZOLE-TRIMETHOPRIM", "NITROFURANTOIN", 
                        "TETRACYCLINES", "inclab", "anyab")

# Age groups of interest
groups <- c("kids0", "kids1to5", "adults")

# 1c. SETUP - DATA FILES ##################################################################

# Load in required files
load("outpatient_bystander_analysis.Rdata")

# 2a. PRE-ANALYSIS - MICROBIOME DATASETS #########################################################
# Goal of processing steps:
# 1. For data not from HMP (kids0 and kids1to5), need to consolidate from different studies to get overall prevalence for each bug
# 2. Create one master microbiome file (with data from all 3 age groups)

### Pre-processing steps
# Create list of body sites in all separate microbiome studies
microbiomeprep.bodysites <- unique(c(microbiome.kids0$Body_site, microbiome.kids1to5$Body_site))

# Change sample size and prevalence columns to numeric 
microbiome.kids0[,5:ncol(microbiome.kids0)] <- apply(microbiome.kids0[,5:ncol(microbiome.kids0)], 2, function(x) as.numeric(as.character(x)))
microbiome.kids1to5[,5:ncol(microbiome.kids1to5)] <- apply(microbiome.kids1to5[,5:ncol(microbiome.kids1to5)], 2, function(x) as.numeric(as.character(x)))

### Calculate overall prevalence numbers for each bug in kids age groups
microbiomeprep.func <- function(kidsdf){
  microbiomeprep.species <- names(kidsdf)[6:ncol(kidsdf)]
  names(kidsdf)[6:ncol(kidsdf)] <- paste0("prev.",names(kidsdf)[6:ncol(kidsdf)])
  
  # Reformat into long dataframe with row for each prevalence from each study at each body site and age
  kidsdf.long <- reshape(kidsdf, direction="long", varying=list(6:ncol(kidsdf)), sep=".", times=microbiomeprep.species,
                         v.names="prev", timevar="species")
  kidsdf.long[,c("n","prev")] <- apply(kidsdf.long[,c("n","prev")], 2, function(x) as.numeric(as.character(x)))
  row.names(kidsdf.long) <- NULL
  
  # Take weighted average of prevalence for each study, species, and body site (over different intermediate ages)
  kidsdf.samestudy <- kidsdf.long %>% group_by(Study,Body_site,species)%>%
                      summarise(total_n = sum(n), wtprev_samestudy = weighted.mean(prev,n,na.rm=TRUE))
  
  # Take weighted average of prevalence for each species and body site (over studies)
  kidsdf.summary <- kidsdf.samestudy %>% group_by(Body_site,species)%>%
                    summarise(wtprev_overstudy = weighted.mean(wtprev_samestudy,total_n,na.rm=TRUE))
  
  kidsdf.summary$wtprev_overstudy[is.na(kidsdf.summary$wtprev_overstudy)] <- 0
  
  # Current df should be three columns: Body_site, species, wtprev_overstudy (weighted prevalence over all studies)
  # Assume independence between body sites to calculate overall prevalences
  p.notpresent <- 1-subset(kidsdf.summary, Body_site==microbiomeprep.bodysites[1])$wtprev_overstudy
  for (i in 2:length(microbiomeprep.bodysites)){
    p.notpresent <- p.notpresent*(1-subset(kidsdf.summary, Body_site==microbiomeprep.bodysites[i])$wtprev_overstudy)
  }
  wtprev_kids <- 1-p.notpresent
  wtprev_kidsdf <- as.data.frame(cbind(microbiomeprep.species, wtprev_kids))
  
  return(wtprev_kidsdf)
}

wtprev_kids0.df <- microbiomeprep.func(microbiome.kids0)
names(wtprev_kids0.df)[2] <- "wtprev_kids0"
wtprev_kids1to5.df <- microbiomeprep.func(microbiome.kids1to5)
names(wtprev_kids1to5.df)[2] <- "wtprev_kids1to5"

### Merge into one microbiome file
microbiome_temp <- merge(microbiome.adults, wtprev_kids0.df, by.x="Species.strain.Key", by.y="microbiomeprep.species")
microbiome <- merge(microbiome_temp, wtprev_kids1to5.df, by.x="Species.strain.Key", by.y="microbiomeprep.species")
microbiome[,2:ncol(microbiome)] <- apply(microbiome[,2:ncol(microbiome)], 2, function(x) as.numeric(as.character(x)))

# No data available for P. aeruginosa and S. agalactiae for kids 1 to 5, so use data for kids under 1
microbiome$wtprev_kids1to5[which(microbiome$Species.strain.Key=="Pseudomonas_aeruginosa")] <- microbiome$wtprev_kids0[which(microbiome$Species.strain.Key=="Pseudomonas_aeruginosa")]
microbiome$wtprev_kids1to5[which(microbiome$Species.strain.Key=="Streptococcus_agalactiae")] <- microbiome$wtprev_kids0[which(microbiome$Species.strain.Key=="Streptococcus_agalactiae")]

# 2b. PRE-ANALYSIS - NAMCS/NHAMCS PREP #########################################################
library(foreign)
library(readstata13)

# i)    Read in NAMCS/NHAMCS -----------------------------------------------------------------------------

# Download NAMCS and NHAMCS 2010-2011 stata files: 
# NAMCS: http://www.nber.org/data/national-ambulatory-medical-care-survey.html
# NHAMCS: http://www.nber.org/data/national-hospital-ambulatory-medical-care-survey.html
namcs2010 <- read.dta13("namcs2010.dta")
namcs2011 <- read.dta13("namcs2011.dta")
nhamcsopd2010 <- read.dta13("nhamcsopd2010.dta")
nhamcsopd2011 <- read.dta13("nhamcsopd2011.dta")
nhamcsed2010 <- read.dta13("nhamcsed2010.dta")
nhamcsed2011 <- read.dta13("nhamcsed2011.dta")
totalVisits.noexcl <- sum(nrow(namcs2010), nrow(namcs2011), nrow(nhamcsopd2010), nrow(nhamcsopd2011), nrow(nhamcsed2010), nrow(nhamcsed2011))

# Mark visits which did not result in hospital or observation visit (1 if no/included, 0 otherwise)
# If doing this in R (or other case-sensitive language), note that in 2010 case answers are NOT capitalized, but ARE capitalized in 2011
namcs2010$nonhospadm <- namcs2010$eradmhos=="no"
namcs2011$nonhospadm <- namcs2011$eradmhos=="No"
nhamcsopd2010$nonhospadm <- nhamcsopd2010$eradmhos=="no"
nhamcsopd2011$nonhospadm <- nhamcsopd2011$eradmhos=="No"
nhamcsed2010$nonhospadm <- nhamcsed2010$admithos=="no" & nhamcsed2010$obshos=="no" & nhamcsed2010$obsdis=="no"
nhamcsed2011$nonhospadm <- nhamcsed2011$admithos=="No" & nhamcsed2011$obshos=="No" & nhamcsed2011$obsdis=="No"

# Select variables to keep and create restricted datasets
varstokeep <- c("drugid1", "drugid2", "drugid3", "drugid4", "drugid5", "drugid6", "drugid7", "drugid8",
               "rx1v2c1", "rx1v2c2", "rx1v2c3", "rx1v2c4","rx2v2c1", "rx2v2c2", "rx2v2c3", "rx2v2c4",
               "rx3v2c1", "rx3v2c2", "rx3v2c3", "rx3v2c4","rx4v2c1", "rx4v2c2", "rx4v2c3", "rx4v2c4",
               "rx5v2c1", "rx5v2c2", "rx5v2c3", "rx5v2c4","rx6v2c1", "rx6v2c2", "rx6v2c3", "rx6v2c4",
               "rx7v2c1", "rx7v2c2", "rx7v2c3", "rx7v2c4","rx8v2c1", "rx8v2c2", "rx8v2c3", "rx8v2c4",
               "med1", "med2", "med3", "med4", "med5", "med6", "med7", "med8",
               "diag1", "diag2", "diag3", "patwt", "cstratm", "cpsum", "age", "nonhospadm", "sex",
               "rx1v3c1", "rx1v3c2", "rx1v3c3", "rx1v3c4","rx2v3c1", "rx2v3c2", "rx2v3c3", "rx2v3c4",
               "rx3v3c1", "rx3v3c2", "rx3v3c3", "rx3v3c4","rx4v3c1", "rx4v3c2", "rx4v3c3", "rx4v3c4",
               "rx5v3c1", "rx5v3c2", "rx5v3c3", "rx5v3c4","rx6v3c1", "rx6v3c2", "rx6v3c3", "rx6v3c4",
               "rx7v3c1", "rx7v3c2", "rx7v3c3", "rx7v3c4","rx8v3c1", "rx8v3c2", "rx8v3c3", "rx8v3c4")
namcs2010.res <- namcs2010[varstokeep]
namcs2011.res <- namcs2011[varstokeep]
nhamcsopd2010.res <- nhamcsopd2010[varstokeep]
nhamcsopd2011.res <- nhamcsopd2011[varstokeep]
nhamcsed2010.res <- nhamcsed2010[varstokeep]
nhamcsed2011.res <- nhamcsed2011[varstokeep]

# ii)   Pull in drug and diagnosis names ----------------------------------------------------------------------
# Use MULTUM Lexicon classification system (in-house system does not have classes)

# Pull in drug identifier codebooks
# Take unique - some have duplicates
DRUGID.10 <- read.csv("DRUGID.10.csv", header=FALSE, col.names=c("Code","Drug"), colClasses="character")
DRUGID.10 <- unique(DRUGID.10)
DRUGID.11 <- read.csv("DRUGID.11.csv", header=FALSE, col.names=c("Code","Drug"), colClasses="character")
DRUGID.11 <- unique(DRUGID.11)

DRGLV2.10 <- read.csv("DRGLV2.10.csv", header=FALSE, col.names=c("Code","Category"), colClasses="character")
DRGLV2.10 <- unique(DRGLV2.10)

DRGLV2.A11 <- read.csv("DRGLV2.A11.csv", header=FALSE, colClasses="character")
DRGLV2.A11 <- DRGLV2.A11[1:2]
names(DRGLV2.A11) <- c("Code", "Category")
DRGLV2.A11 <- unique(DRGLV2.A11)

DRGLV2.H11 <- read.csv("DRGLV2.H11.csv", header=FALSE, colClasses="character")
DRGLV2.H11 <- DRGLV2.H11[1:2]
names(DRGLV2.H11) <- c("Code", "Category")
DRGLV2.H11 <- unique(DRGLV2.H11)

DRGLV2.HED11 <- read.csv("DRGLV2.HED11.csv", header=FALSE, colClasses="character")
DRGLV2.HED11 <- DRGLV2.HED11[1:2]
names(DRGLV2.HED11) <- c("Code", "Category")
DRGLV2.HED11 <- unique(DRGLV2.HED11)

# Read in diagnosis identifier codebooks (no duplicates found in any of the diag codebooks)
DIAG.A10 <- read.csv("DIAG.A10.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")
DIAG.A11 <- read.csv("DIAG.A11.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")

DIAG.H10 <- read.csv("DIAG.H10.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")
DIAG.H11 <- read.csv("DIAG.H11.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")

DIAG.HED10 <- read.csv("DIAG.HED10.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")
DIAG.HED11 <- read.csv("DIAG.HED11.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")

# Function to pull in drug names, MULTUM category names, and full diagnoses with appropriate codebook
decodeDrugDiag <- function(dataset, df, drugCodeBook, catCodeBook, diagCodeBook){
  df.temp <- df
  
  # Merge in drug names
  drugid1.1st <- which(names(namcs2010.res)=="drugid1")
  for(i in drugid1.1st:(drugid1.1st+7)){
    df.temp <- merge(df.temp, drugCodeBook, by.x=c(names(df[i])), by.y=c("Code"), all.x=TRUE, all.y=FALSE)
    names(df.temp)[ncol(df.temp)] <- c(paste(names(df)[i],".name", sep=""))
  }
  
  # Merge in category names
  # Level 2 for 8 medications, each one can have up to 4 associated categories
  rx1v2c1.1st <- which(names(namcs2010.res)=="rx1v2c1")
  for(i in rx1v2c1.1st:(rx1v2c1.1st + 8*4 - 1)){
    df.temp <- merge(df.temp, catCodeBook, by.x=c(names(df[i])), by.y=c("Code"), all.x=TRUE, all.y=FALSE)
    names(df.temp)[ncol(df.temp)] <- c(paste(names(df)[i],".name", sep=""))
  }
  
  # Merge in diagnosis names
  diag1.1st <- which(names(namcs2010.res)=="diag1")
  for (i in diag1.1st:(diag1.1st+2)){
    df.temp <- merge(df.temp, diagCodeBook, by.x=c(names(df[i])), by.y=c("Code"), all.x=TRUE, all.y=FALSE)
    names(df.temp)[ncol(df.temp)] <- c(paste(names(df)[i],".name", sep=""))
  }
  
  df.temp <- cbind(dataset=dataset, df.temp)
  
  return(df.temp)
}

namcs2010.drugdecode <- decodeDrugDiag("namcs2010", namcs2010.res, DRUGID.10, DRGLV2.10, DIAG.A10)
namcs2011.drugdecode <- decodeDrugDiag("namcs2011", namcs2011.res, DRUGID.11, DRGLV2.A11, DIAG.A11)
nhamcsopd2010.drugdecode <- decodeDrugDiag("nhamcsopd2010", nhamcsopd2010.res, DRUGID.10, DRGLV2.10, DIAG.H10)
nhamcsopd2011.drugdecode <- decodeDrugDiag("nhamcsopd2011", nhamcsopd2011.res, DRUGID.11, DRGLV2.H11, DIAG.H11)
nhamcsed2010.drugdecode <- decodeDrugDiag("nhamcsed2010", nhamcsed2010.res, DRUGID.10, DRGLV2.10, DIAG.HED10)
nhamcsed2011.drugdecode <- decodeDrugDiag("nhamcsed2011", nhamcsed2011.res, DRUGID.11, DRGLV2.HED11, DIAG.HED11)

# iii)  Add antibiotic flags ('antibiotic<#>_flag') ----------------------------------------------------------------------------
# Create flag for each med (1-8) where 1=med is antibiotic, 0=otherwise (based on specification within anti-infectives, see link below)
# Lexicon: https://wwwn.cdc.gov/nchs/nhanes/1999-2000/RXQ_DRUG.htm#Appendix_3:_Multum_Lexicon_Therapeutic_Classification_Scheme_
# Additional NAMCS/NHAMCS resource: https://www.cdc.gov/nchs/ppt/nchs2010/13_schappert.pdf
# 2018.09.25 Added parentheses around grepl statement, which affects summary stats at beginning of results and McAdams (but I think did not use any of these restrictions for that analysis)
antibioticsflag <- function(df){
  rx1v2c1.name.1st <- which(names(df)=="rx1v2c1.name")
  
  # Create empty dataframe to store flags
  antibioticsflag.df <- data.frame(matrix(vector(),nrow(df),8))
  
  # j represents the antibiotic # (can record up to 8 medications)
  for (j in 0:7){
    temp <- rep(0,nrow(df))
    
    # For each med (1-8), look through all level 2 categories (up to 4 may exist)
    # Add if any of them is an "Anti-infective" but NOT one of the other categories listed
    for (i in c((rx1v2c1.name.1st+4*j):(rx1v2c1.name.1st+4*j+3))){
      temp <- temp + (grepl("Anti-infectives", df[,i]) & !grepl("amebicides", df[,i]) & !grepl("anthelmintics", df[,i]) & !grepl("antifungals", df[,i]) & !grepl("antimalarial agents", df[,i]) & !grepl("antiviral agents", df[,i]))
    }
    
    antibioticsflag.df[,(j+1)] <- as.numeric(temp>0)
  }
  
  # Update names of antibiotic flag columns
  for (i in 1:8){
    names(antibioticsflag.df)[i] <- paste("antibioticflag_",i,sep="")
  }
  
  df <- cbind(df, antibioticsflag.df)
  return(df)
}

namcs2010.drugdecode.ab <- antibioticsflag(namcs2010.drugdecode)
namcs2011.drugdecode.ab <- antibioticsflag(namcs2011.drugdecode)
nhamcsopd2010.drugdecode.ab <- antibioticsflag(nhamcsopd2010.drugdecode)
nhamcsopd2011.drugdecode.ab <- antibioticsflag(nhamcsopd2011.drugdecode)
nhamcsed2010.drugdecode.ab <- antibioticsflag(nhamcsed2010.drugdecode)
nhamcsed2011.drugdecode.ab <- antibioticsflag(nhamcsed2011.drugdecode)

# iv)   Create flag for inclusion criteria ('inclusion') ----------------------------------------------------------------------------

# Create overall inclusion flag (1 to include, 0 otherwise)
# Removed criteria relating to oral, parenteral and topical antibiotics from Fleming-Dutra paper (code in v17)
# Assumed that `nonhospadm` flag removed all inpatient visits

namcs2010.drugdecode.ab$inclusion <- namcs2010.drugdecode.ab$nonhospadm
namcs2011.drugdecode.ab$inclusion <- namcs2011.drugdecode.ab$nonhospadm
nhamcsopd2010.drugdecode.ab$inclusion <- nhamcsopd2010.drugdecode.ab$nonhospadm
nhamcsopd2011.drugdecode.ab$inclusion <- nhamcsopd2011.drugdecode.ab$nonhospadm
nhamcsed2010.drugdecode.ab$inclusion <- nhamcsed2010.drugdecode.ab$nonhospadm
nhamcsed2011.drugdecode.ab$inclusion <- nhamcsed2011.drugdecode.ab$nonhospadm

# Count total visits
totalVisits.incl <- sum(namcs2010.drugdecode.ab$inclusion, namcs2011.drugdecode.ab$inclusion,
                    nhamcsopd2010.drugdecode.ab$inclusion, nhamcsopd2011.drugdecode.ab$inclusion,
                    nhamcsed2010.drugdecode.ab$inclusion, nhamcsed2011.drugdecode.ab$inclusion)

# v)    Add flags for antibiotics of interest ('<ANTIBIOTIC>_flag')----------------------------------------------------------------------------
# Create flag for each antibiotic where 1=antibiotic a was prescribed at that visits, 0=otherwise

# length(antibiotics)-2 because do not create flags for "inclabnum" and "anyab" here
specificABflag <- function(df){
  for(i in 1:(length(antibiotics)-2)){
    temp <- df$drugid1.name==antibiotics[i] | df$drugid2.name==antibiotics[i] |
            df$drugid3.name==antibiotics[i] | df$drugid4.name==antibiotics[i] |
            df$drugid5.name==antibiotics[i] | df$drugid6.name==antibiotics[i] |
            df$drugid7.name==antibiotics[i] | df$drugid8.name==antibiotics[i]
    df[,ncol(df)+1] <- temp
    names(df)[ncol(df)] <- c(paste(antibiotics[i], "flag", sep="_"))
  } 
  return(df)
}

namcs2010.drugdecode.ab <- specificABflag(namcs2010.drugdecode.ab)
namcs2011.drugdecode.ab <- specificABflag(namcs2011.drugdecode.ab)
nhamcsopd2010.drugdecode.ab <- specificABflag(nhamcsopd2010.drugdecode.ab)
nhamcsopd2011.drugdecode.ab <- specificABflag(nhamcsopd2011.drugdecode.ab)
nhamcsed2010.drugdecode.ab <- specificABflag(nhamcsed2010.drugdecode.ab)
nhamcsed2011.drugdecode.ab <- specificABflag(nhamcsed2011.drugdecode.ab)

# vi)   Add counts for drug classes of interest ('<CLASS>_count') ---------------------------------------------------------
# Count of number of drugs of that drug class prescribed at each visit
specificABclasscount <- function(df){
  rx1v2c1.name.1st <- which(names(df)=="rx1v2c1.name")
  
  # Identifiers (based on level 2 of MULTUM Lexicon classification system) for each drug class within anti-infectives
  # Order of identifiers here matches with `drugclasses.unique` -- MULTUM identifier is irrelevant for ones with "x"
  multum.names <- c("penicillins", "macrolide", "quinolones", "cephalosporins", "x", "x", "tetracyclines", "x", "x")
  
  for(i in c(1,2:4,7)){
    df[,paste0(drugclasses.unique[i], "_count")] <- rowSums(apply(df[, c(rx1v2c1.name.1st:(rx1v2c1.name.1st + 8*4 - 1))], 2, function(x) grepl(multum.names[i], x, ignore.case=TRUE)))
  }
  
  # Macrolides class (in this analysis) additionally includes lincomycin derivatives
  df[,paste0(drugclasses.unique[2], "_count")] <- df[,paste0(drugclasses.unique[2], "_count")] +
                                                rowSums(apply(df[, c(rx1v2c1.name.1st:(rx1v2c1.name.1st + 8*4 - 1))], 2, function(x) grepl("lincomycin derivatives", x, ignore.case=TRUE)))
  
  # For nitrofurantoin and TMP/SMX, they are only drugs in their class
  df$NITROFURANTOIN_count <- df$NITROFURANTOIN_flag
  df$`SULFAMETHOXAZOLE-TRIMETHOPRIM_count` <- df$`SULFAMETHOXAZOLE-TRIMETHOPRIM_flag`
  
  return(df)
}

namcs2010.drugdecode.ab <- specificABclasscount(namcs2010.drugdecode.ab)
namcs2011.drugdecode.ab <- specificABclasscount(namcs2011.drugdecode.ab)
nhamcsopd2010.drugdecode.ab <- specificABclasscount(nhamcsopd2010.drugdecode.ab)
nhamcsopd2011.drugdecode.ab <- specificABclasscount(nhamcsopd2011.drugdecode.ab)
nhamcsed2010.drugdecode.ab <- specificABclasscount(nhamcsed2010.drugdecode.ab)
nhamcsed2011.drugdecode.ab <- specificABclasscount(nhamcsed2011.drugdecode.ab)

# vii)  Add counts for overall antibiotics ('inclab_count' and 'anyab_count') ---------------------------------------------------------
# 'inclab_count' = count of antibiotic in any of the explicitly included classes prescribed at that visit
# 'anyab_count' = count of any antibiotic prescribed at that visit, as identified by MULTUM Lexicon system in NAMCS/NHAMCS
# 2018.06.25 Change to now include any of the drug classes of interest (since shown in Figure 1)
# 2018.06.28 Changed antibiotics.unique to drugclasses.unique
abcounts <- function(df){
  inclab.temp <- rep(0, nrow(df))
  #drugclasses.unique-2 because does not include "inclab" or "anyab"
  for (i in 1:(length(drugclasses.unique)-2)){
    inclab.temp <- inclab.temp + df[,paste0(drugclasses.unique[i], "_count")]
  }
  df$inclab_count <- inclab.temp
  
  anyab.temp <- df$antibioticflag_1 + df$antibioticflag_2 +
                df$antibioticflag_3 + df$antibioticflag_4 +
                df$antibioticflag_5 + df$antibioticflag_6 + 
                df$antibioticflag_7 + df$antibioticflag_8
  df$anyab_count <- anyab.temp
  return(df)
}

namcs2010.drugdecode.ab <- abcounts(namcs2010.drugdecode.ab)
namcs2011.drugdecode.ab <- abcounts(namcs2011.drugdecode.ab)
nhamcsopd2010.drugdecode.ab <- abcounts(nhamcsopd2010.drugdecode.ab)
nhamcsopd2011.drugdecode.ab <- abcounts(nhamcsopd2011.drugdecode.ab)
nhamcsed2010.drugdecode.ab <- abcounts(nhamcsed2010.drugdecode.ab)
nhamcsed2011.drugdecode.ab <- abcounts(nhamcsed2011.drugdecode.ab)

# viii) Add flags for conditions of interest ('<CONDITION>_flag') ----------------------------------------------------------------------------
# Add flag (1=diagnosed during that visit, 0=otherwise) based on codes in condition-code key
# Condition-code key contains diagnosis codes to include in that condition and diagnosis codes to exclude
# E.g. Fleming-Dutra exclude diagnoses of bronchitis that have a co-diagnosis of COPD

#'creategreplstring' function creates the string of code needed to appropriately include/exclude codes
creategreplstring <- function(df, inclcode, exclcode){
  m <- length(inclcode)
  n <- length(exclcode)

  # 'grepl' marks as TRUE if string (e.g. df$diag1) includes the given string segment (inclcode[1], in the very first part of this line)
  greplstr <- paste0("grepl('", inclcode[1], "', df$diag1) | grepl('", inclcode[1], "', df$diag2) | grepl('", inclcode[1], "', df$diag3) ")

  if (m>1){
    for (i in 2:m){
      greplstr <- paste0(greplstr, "| grepl('", inclcode[i], "', df$diag1) | grepl('", inclcode[i], "', df$diag2) | grepl('", inclcode[i], "', df$diag3)")
    }
  }

  if(n>0){
    greplstr <- paste0("(", greplstr, ")", " & !(grepl('", exclcode[1], "', df$diag1) | grepl('", exclcode[1], "', df$diag2) | grepl('", exclcode[1], "', df$diag3)")
    if (n>1){
      for (i in 2:n){
        greplstr <- paste0(greplstr, " | grepl('", exclcode[i], "', df$diag1) | grepl('", exclcode[i], "', df$diag2) | grepl('", exclcode[i], "', df$diag3)")
      }
    }
    greplstr <- paste0(greplstr, ")")
  }

  return(eval(parse(text=greplstr)))
}

specificconditionflag <- function(df){
  for (i in 1:length(conditions)){
    incl <- subset(conditionCodes, conditionCodes$Condition==conditions[i] & conditionCodes$Exclude==0)$Code
    excl <- subset(conditionCodes, conditionCodes$Condition==conditions[i] & conditionCodes$Exclude==1)$Code
    df[,paste0(conditions[i], "_flag")] <- creategreplstring(df, incl, excl)
  }
  return(df)
}

namcs2010.drugdecode.ab.cond <- specificconditionflag(namcs2010.drugdecode.ab)
namcs2011.drugdecode.ab.cond <- specificconditionflag(namcs2011.drugdecode.ab)
nhamcsopd2010.drugdecode.ab.cond <- specificconditionflag(nhamcsopd2010.drugdecode.ab)
nhamcsopd2011.drugdecode.ab.cond <- specificconditionflag(nhamcsopd2011.drugdecode.ab)
nhamcsed2010.drugdecode.ab.cond <- specificconditionflag(nhamcsed2010.drugdecode.ab)
nhamcsed2011.drugdecode.ab.cond <- specificconditionflag(nhamcsed2011.drugdecode.ab)

# viii.b) Create TIERED diagosis (from Fleming-Dutra) - COMMENT OUT IF DOING BASELINE ANALYSIS ----------------------------------------------------------------------------
# # 'creategreplstring' function creates the string of code needed to appropriately include/exclude codes
# creategreplstring <- function(df, diagnum, inclcode, exclcode){
#   m <- length(inclcode)
#   n <- length(exclcode)
# 
#   # 'grepl' marks as TRUE if string (e.g. df$diag1) includes the given string segment (inclcode[1], in the very first part of this line)
#   greplstr <- paste0("grepl('", inclcode[1], "', df$diag", diagnum, ") ")
# 
#   if (m>1){
#     for (i in 2:m){
#       greplstr <- paste0(greplstr, "| grepl('", inclcode[i], "', df$diag", diagnum, ")")
#     }
#   }
# 
#   if(n>0){
#     greplstr <- paste0("(", greplstr, ")", " & !(grepl('", exclcode[1], "', df$diag1) | grepl('", exclcode[1], "', df$diag2) | grepl('", exclcode[1], "', df$diag3)")
#     if (n>1){
#       for (i in 2:n){
#         greplstr <- paste0(greplstr, " | grepl('", exclcode[i], "', df$diag1) | grepl('", exclcode[i], "', df$diag2) | grepl('", exclcode[i], "', df$diag3)")
#       }
#     }
#     greplstr <- paste0(greplstr, ")")
#   }
# 
#   return(eval(parse(text=greplstr)))
# }
# 
# specificconditionflag <- function(df){
#   for (i in 1:length(conditions)){
#     # for diagnosis numbers
#     for (j in 1:3){
#       incl <- subset(conditionCodes, conditionCodes$Condition==conditions[i] & conditionCodes$Exclude==0)$Code
#       excl <- subset(conditionCodes, conditionCodes$Condition==conditions[i] & conditionCodes$Exclude==1)$Code
#       temp <- creategreplstring(df, j, incl, excl)
#       df[which(temp==TRUE), paste0("diag",j,".condition")] <- conditions[i]
#     }
#   }
#   return(df)
# }
# 
# namcs2010.drugdecode.ab.cond <- specificconditionflag(namcs2010.drugdecode.ab)
# namcs2011.drugdecode.ab.cond <- specificconditionflag(namcs2011.drugdecode.ab)
# nhamcsopd2010.drugdecode.ab.cond <- specificconditionflag(nhamcsopd2010.drugdecode.ab)
# nhamcsopd2011.drugdecode.ab.cond <- specificconditionflag(nhamcsopd2011.drugdecode.ab)
# nhamcsed2010.drugdecode.ab.cond <- specificconditionflag(nhamcsed2010.drugdecode.ab)
# nhamcsed2011.drugdecode.ab.cond <- specificconditionflag(nhamcsed2011.drugdecode.ab)
# 
# assigntier <- function(df){
#   for(j in 1:3){
#     tier1.conditions <- c("pneumonia", "pneumoniaSP", "pneumoniaKP", "pneumoniaPA", "pneumoniaHI", "pneumoniaSgpA", "pneumoniaSgpB", "pneumoniaSA", "pneumoniaEC", "uti")
#     tier2.conditions <- c("acne", "pharyngitis.nostrep", "strep", "acuteSinusitis", "chronicSinusitis", "suppOtitisMedia", "cellulitis")
#     tier3.conditions <- c("bronchitis", "viralPneumonia", "viralURI", "influenza", "miscBacInf", "allergicRhinitis", "nonSuppOtitisMedia")
# 
#     df[df[,paste0("diag",j,".condition")] %in% tier1.conditions, paste0("diag",j,".tier")] <- 1
#     df[df[,paste0("diag",j,".condition")] %in% tier2.conditions, paste0("diag",j,".tier")] <- 2
#     df[df[,paste0("diag",j,".condition")] %in% tier3.conditions, paste0("diag",j,".tier")] <- 3
#   }
#   df$mintier <- unlist(lapply(1:nrow(df), function(x) min(df$diag1.tier[x], df$diag2.tier[x], df$diag3.tier[x], na.rm=TRUE)))
#   df$mintier[which(df$mintier==Inf)] <- NA
#   return(df)
# }
# 
# namcs2010.drugdecode.ab.cond <- assigntier(namcs2010.drugdecode.ab.cond)
# namcs2011.drugdecode.ab.cond <- assigntier(namcs2011.drugdecode.ab.cond)
# nhamcsopd2010.drugdecode.ab.cond <- assigntier(nhamcsopd2010.drugdecode.ab.cond)
# nhamcsopd2011.drugdecode.ab.cond <- assigntier(nhamcsopd2011.drugdecode.ab.cond)
# nhamcsed2010.drugdecode.ab.cond <- assigntier(nhamcsed2010.drugdecode.ab.cond)
# nhamcsed2011.drugdecode.ab.cond <- assigntier(nhamcsed2011.drugdecode.ab.cond)
# 
# pullcondition <- function(df){
#   # Add 1 to incorporate the 'NA' (if not found, `which` returns integer(0))
#   df$finalcondition <- unlist(lapply(1:nrow(df), function(x) c(NA, df$diag1.condition[x], df$diag2.condition[x], df$diag3.condition[x])[min(which(c(df$diag1.tier[x], df$diag2.tier[x], df$diag3.tier[x])==df$mintier[x]))+1]))
#   return(df)
# }
# 
# namcs2010.drugdecode.ab.cond <- pullcondition(namcs2010.drugdecode.ab.cond)
# namcs2011.drugdecode.ab.cond <- pullcondition(namcs2011.drugdecode.ab.cond)
# nhamcsopd2010.drugdecode.ab.cond <- pullcondition(nhamcsopd2010.drugdecode.ab.cond)
# nhamcsopd2011.drugdecode.ab.cond <- pullcondition(nhamcsopd2011.drugdecode.ab.cond)
# nhamcsed2010.drugdecode.ab.cond <- pullcondition(nhamcsed2010.drugdecode.ab.cond)
# nhamcsed2011.drugdecode.ab.cond <- pullcondition(nhamcsed2011.drugdecode.ab.cond)
# 
# finalconditionflag <- function(df){
#   for(i in 1:length(conditions)){
#     df[,paste0(conditions[i], "_flag")] <- (df$finalcondition==conditions[i])
#     df[,paste0(conditions[i], "_flag")][is.na(df[,paste0(conditions[i], "_flag")])] <- 0
#     df[,paste0(conditions[i], "_flag")][df[,paste0(conditions[i], "_flag")]==FALSE] <- 0
#     df[,paste0(conditions[i], "_flag")][df[,paste0(conditions[i], "_flag")]==TRUE] <- 1
#   }
#   return(df)
# }
# 
# namcs2010.drugdecode.ab.cond <- finalconditionflag(namcs2010.drugdecode.ab.cond)
# namcs2011.drugdecode.ab.cond <- finalconditionflag(namcs2011.drugdecode.ab.cond)
# nhamcsopd2010.drugdecode.ab.cond <- finalconditionflag(nhamcsopd2010.drugdecode.ab.cond)
# nhamcsopd2011.drugdecode.ab.cond <- finalconditionflag(nhamcsopd2011.drugdecode.ab.cond)
# nhamcsed2010.drugdecode.ab.cond <- finalconditionflag(nhamcsed2010.drugdecode.ab.cond)
# nhamcsed2011.drugdecode.ab.cond <- finalconditionflag(nhamcsed2011.drugdecode.ab.cond)

# ix)   Add flag and count for any condition of interest ('inclcond_flag', 'inclcond_count', 'other_flag') ----------------------------------------------------------------------------
# 'inclcond_flag'=1 if any of our conditions of interest was diagnosed during that visit, 0 otherwise
# 'inclcond_count' is the count of conditions of interest diagnosed during that visit
# 'other_flag'=1 if none of the conditions of interest was diagnosed during that visit (inverse of 'inclcond_flag')
inclcondfunc <- function(df){
  temp <- rep(0, nrow(df))
  for (i in 1:length(conditions)){
    temp <- temp + df[,paste0(conditions[i], "_flag")]
  }
  df$inclcond_flag <- as.numeric(temp>=1)
  df$inclcond_count <- temp
  df$other_flag <- as.numeric(df$inclcond_flag==0)
  return(df)
}

namcs2010.drugdecode.ab.cond <- inclcondfunc(namcs2010.drugdecode.ab.cond)
namcs2011.drugdecode.ab.cond <-inclcondfunc(namcs2011.drugdecode.ab.cond)
nhamcsopd2010.drugdecode.ab.cond <- inclcondfunc(nhamcsopd2010.drugdecode.ab.cond)
nhamcsopd2011.drugdecode.ab.cond <- inclcondfunc(nhamcsopd2011.drugdecode.ab.cond)
nhamcsed2010.drugdecode.ab.cond <- inclcondfunc(nhamcsed2010.drugdecode.ab.cond)
nhamcsed2011.drugdecode.ab.cond <- inclcondfunc(nhamcsed2011.drugdecode.ab.cond)

# x)    Add flags for age groups of interest ----------------------------------------------------------------------------

agegroupsfunc <- function(df){
  df$agegroup[df$age==0] <- "kids0"
  df$agegroup[df$age>0 & df$age<=5] <- "kids1to5"
  df$agegroup[df$age>5] <- "adults"
  return(df)
}

namcs2010.drugdecode.ab.cond.grp <- agegroupsfunc(namcs2010.drugdecode.ab.cond)
namcs2011.drugdecode.ab.cond.grp <- agegroupsfunc(namcs2011.drugdecode.ab.cond)
nhamcsopd2010.drugdecode.ab.cond.grp <- agegroupsfunc(nhamcsopd2010.drugdecode.ab.cond)
nhamcsopd2011.drugdecode.ab.cond.grp <- agegroupsfunc(nhamcsopd2011.drugdecode.ab.cond)
nhamcsed2010.drugdecode.ab.cond.grp <- agegroupsfunc(nhamcsed2010.drugdecode.ab.cond)
nhamcsed2011.drugdecode.ab.cond.grp <- agegroupsfunc(nhamcsed2011.drugdecode.ab.cond)

# 2c. PRE-ANALYSIS - BOOTSTRAP SETUP --------------------------------------------------------------------------------
# Calculate means and CI for proportions of diagnoses with given medication using 'survey' package
# R documentation: https://cran.r-project.org/web/packages/survey/survey.pdf
# For multistage sampling:
# id = formula with the cluster identifiers at each stage
# strata = if subsequent stages are stratified, should be a formula with stratum identifiers at each stage
# fpc = population size for each level of sampling; if not specified, then only does one-level with replacement
# nest = TRUE --> PSUs reuse the same identifiers across strata
library(survey)
library(srvyr)

# According to this documentation: https://www.cdc.gov/nchs/ppt/nchs2015/Hing_Monday_GlenEcho_C2.pdf
# NHAMCS researchers should use combined ED and OPD files when computing variances for emergency and/or outpatient department estimates. 
# Including both files takes into account NHAMCS’s complex sample design.
# Thus from now on, we combine both NHAMCS files
# The "names..." part ensures that columns are in the same order
nhamcs2010.drugdecode.ab.cond.grp <- rbind(nhamcsopd2010.drugdecode.ab.cond.grp, nhamcsed2010.drugdecode.ab.cond.grp[,names(nhamcsopd2010.drugdecode.ab.cond.grp)])
nhamcs2011.drugdecode.ab.cond.grp <- rbind(nhamcsopd2011.drugdecode.ab.cond.grp, nhamcsed2011.drugdecode.ab.cond.grp[,names(nhamcsopd2011.drugdecode.ab.cond.grp)])

# `surveyCalcs` function creates summaries of NAMCS/NHAMCS with weighted visits by condition and antibiotic and corresponding SE
# options(survey.lonely.psu=) accounts for strata that have a single PSU
# Based on R documentation: http://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
# With options(survey.lonely.psu="adjust") the data for the single-PSU stratum are centered at the sample grand mean rather than the stratum mean. 
# This is conservative.
# Includes "other" as a condition - visits where none of our included conditions were diagnosed
options(survey.lonely.psu = "adjust")
surveyCalcs <- function(dataset, class){
  # Create empty dataframe for results
  results <- data.frame(dataset=character(),
                        condition=character(),
                        agegroup=character(),
                        antibiotic=character(),
                        mean=double(),
                        SE=double(),
                        stringsAsFactors = FALSE)
  
  # Set up survey design
  design <- svydesign(ids=~cpsum, strata=~cstratm, weights=~patwt, data=get(paste0(dataset, ".drugdecode.ab.cond.grp")), nest=TRUE)
  
  # `n` to be appropriate list (antibiotics or drug classes)
  if(class=="no"){n<-antibiotics[1:(length(antibiotics)-2)]}else{n<-drugclasses.unique}
  
  # Create special conditions list (just for this function) with 'other' included as a condition
  surveyCalcs.conditions <- c(conditions, "other")
  
  for (i in 1:(length(n))){
    for (j in 1:length(surveyCalcs.conditions)){
      # Set flags for this iteration
      if(class=="no"){tempabflag <- paste0(antibiotics[i], "_flag")}else{tempabflag <- paste0(drugclasses.unique[i], "_count")}
      tempcondflag <- paste0(surveyCalcs.conditions[j], "_flag")
      
      # Create survey summary: select for included counts (weighted by patwt) by age group, condition, and antibiotic
      # Add additional columns for labelling dataset, antibiotic and condition
      tempsurvey <- as.data.frame(svyby(~get(tempabflag), ~inclusion+agegroup+get(tempcondflag), design, svytotal))
      tempsurvey.res <- subset(tempsurvey, inclusion==1 & tempsurvey[,which(names(tempsurvey)=="get(tempcondflag)")]==TRUE)
      rownames(tempsurvey.res) <- NULL
      # Fix (in the event that there are NO visits with that condition in any age group)
      # Otherwise throws error because dataframe is completely empty
      if(nrow(tempsurvey.res)==0){tempsurvey.res[1,] <- rep(0,ncol(tempsurvey.res))}
      tempsurvey.res <- merge(data.frame(agegroup=groups), tempsurvey.res, by.x="agegroup", by.y="agegroup", all.x=TRUE)
      tempsurvey.res$dataset <- dataset
      if(class=="no"){tempsurvey.res$antibiotic <- antibiotics[i]}else{tempsurvey.res$antibiotic <- drugclasses.unique[i]}
      tempsurvey.res$condition <- surveyCalcs.conditions[j]
      tempsurvey.res[is.na(tempsurvey.res)] <- 0
      
      # Results have different column names depending on whether column is binary or count
      # Keep columns of interest
      if((class=="no") | (class=="yes" & (drugclasses.unique[i]=="NITROFURANTOIN"|drugclasses.unique[i]=="SULFAMETHOXAZOLE-TRIMETHOPRIM"))){
        tempresults <- tempsurvey.res[, c("dataset", "condition", "agegroup", "antibiotic","get(tempabflag)TRUE", "se.get(tempabflag)TRUE")]
      } else {
        tempresults <- tempsurvey.res[, c("dataset", "condition", "agegroup", "antibiotic", "get(tempabflag)", "se")]
      }
      
      if(class=="no"){names(tempresults) <- c("dataset", "condition", "agegroup", "antibiotic", "wtVisits", "wtVisits.se")}
      else{names(tempresults) <- c("dataset", "condition", "agegroup", "drugclass", "wtVisits", "wtVisits.se")}
      
      results <- rbind(results, tempresults)
    }
  }
  
  return(results)
}

namcs2010.summary <- surveyCalcs("namcs2010", "no")
namcs2011.summary <- surveyCalcs("namcs2011", "no")
nhamcs2010.summary <- surveyCalcs("nhamcs2010", "no")
nhamcs2011.summary <- surveyCalcs("nhamcs2011", "no")

namcs2010.summary.byclass <- surveyCalcs("namcs2010", "yes")
namcs2011.summary.byclass <- surveyCalcs("namcs2011", "yes")
nhamcs2010.summary.byclass <- surveyCalcs("nhamcs2010", "yes")
nhamcs2011.summary.byclass <- surveyCalcs("nhamcs2011", "yes")

# i)    Count visits by antibiotic and total visits ----------------------------------------------------------
# Identify column index numbers in summary datasets
condition.col <- which(names(namcs2010.summary)=="condition")
antibiotic.col <- which(names(namcs2010.summary)=="antibiotic")
agegroup.col <- which(names(namcs2010.summary)=="agegroup")

# "visits" = raw number of visits with given condition and given antibiotic prescribed
namcs2010.summary$visits <- apply(namcs2010.summary, 1, function(x) nrow(subset(namcs2010.drugdecode.ab.cond.grp, agegroup==x[agegroup.col] & 
                                                                                get(paste0(x[condition.col], "_flag"))==1 & get(paste0(x[antibiotic.col], "_flag"))==1)))
namcs2011.summary$visits <- apply(namcs2011.summary, 1, function(x) nrow(subset(namcs2011.drugdecode.ab.cond.grp, agegroup==x[agegroup.col] & 
                                                                                get(paste0(x[condition.col], "_flag"))==1 & get(paste0(x[antibiotic.col], "_flag"))==1)))
nhamcs2010.summary$visits <- apply(nhamcs2010.summary, 1, function(x) nrow(subset(nhamcs2010.drugdecode.ab.cond.grp, agegroup==x[agegroup.col] & 
                                                                                  get(paste0(x[condition.col], "_flag"))==1 & get(paste0(x[antibiotic.col], "_flag"))==1)))
nhamcs2011.summary$visits <- apply(nhamcs2011.summary, 1, function(x) nrow(subset(nhamcs2011.drugdecode.ab.cond.grp, agegroup==x[agegroup.col] & 
                                                                                  get(paste0(x[condition.col], "_flag"))==1 & get(paste0(x[antibiotic.col], "_flag"))==1)))

# "totalVisits" = raw number of visits with given condition
namcs2010.summary$totalVisits <- apply(namcs2010.summary, 1, function(x) nrow(subset(namcs2010.drugdecode.ab.cond.grp, 
                                                                                     agegroup==x[agegroup.col] & get(paste0(x[condition.col], "_flag"))==1)))
namcs2011.summary$totalVisits <- apply(namcs2011.summary, 1, function(x) nrow(subset(namcs2011.drugdecode.ab.cond.grp, 
                                                                                     agegroup==x[agegroup.col] & get(paste0(x[condition.col], "_flag"))==1)))
nhamcs2010.summary$totalVisits <- apply(nhamcs2010.summary, 1, function(x) nrow(subset(nhamcs2010.drugdecode.ab.cond.grp, 
                                                                                       agegroup==x[agegroup.col] & get(paste0(x[condition.col], "_flag"))==1)))
nhamcs2011.summary$totalVisits <- apply(nhamcs2011.summary, 1, function(x) nrow(subset(nhamcs2011.drugdecode.ab.cond.grp, 
                                                                                       agegroup==x[agegroup.col] & get(paste0(x[condition.col], "_flag"))==1)))

# ii)   Merge NAMCS/NHAMCS dataframes ----------------------------------------------------
# Merge full NAMCS/NHAMCS datasets
fullNAMCS <- rbind(namcs2010.drugdecode.ab.cond.grp,
                   namcs2011.drugdecode.ab.cond.grp[names(namcs2010.drugdecode.ab.cond.grp)],
                   nhamcs2010.drugdecode.ab.cond.grp[names(namcs2010.drugdecode.ab.cond.grp)],
                   nhamcs2011.drugdecode.ab.cond.grp[names(namcs2010.drugdecode.ab.cond.grp)])
fullNAMCS.incl <- subset(fullNAMCS, inclusion==1)

# Merge NAMCS/NHAMCS summaries
NAMCS.summary <- rbind(namcs2010.summary, 
                       namcs2011.summary[names(namcs2010.summary)],
                       nhamcs2010.summary[names(namcs2010.summary)], 
                       nhamcs2011.summary[names(namcs2010.summary)])
NAMCS.summary[is.na(NAMCS.summary)] <- 0
# Pull in drug classes
NAMCS.summary <- merge(NAMCS.summary, drugclasses.key, by.x="antibiotic")

# Merge NAMCS/NHAMCS summaries by class
NAMCS.summary.byclass <- rbind(namcs2010.summary.byclass, 
                               namcs2011.summary.byclass[names(namcs2010.summary.byclass)],
                               nhamcs2010.summary.byclass[names(namcs2010.summary.byclass)], 
                               nhamcs2011.summary.byclass[names(namcs2010.summary.byclass)])
NAMCS.summary.byclass[is.na(NAMCS.summary.byclass)] <- 0

# 3. BYSTANDER ANALYSIS ###################################################################
# Weighted visits are used for all bystander analyses
# Setup bystander analysis --------------------------------------------------------------
# Set up dataframes for mapply
# Remove last two ("inclab" and "anyab")
antibiotic.species <- expand.grid(antibiotics[1:(length(antibiotics)-2)], species)
names(antibiotic.species) <- c("antibiotic", "species")
antibiotic.species[,c(1:2)] <- sapply(antibiotic.species[,c(1:2)], as.character)

drugclass.species <- expand.grid(drugclasses.unique, species)
names(drugclass.species) <- c("drugclass", "species")
drugclass.species[,c(1:2)] <- sapply(drugclass.species[,c(1:2)], as.character)

# Bystander prop by AB and species -------------------------------------------------------------------
# Preliminary setup for bootstrapping
etiologies[is.na(etiologies)] <- 0
etiologies.bs <- etiologies
microbiome.bs <- microbiome
NAMCS.summary$wtVisits_curr <- NAMCS.summary$wtVisits

# Calculate proportion of bystander exposures for a given antibiotic and species, summed across all conditions and groups
bystanderCalc <- function(species_temp, antibiotic_temp){
  # Initialize T (denominator) and N (numerator) at 0
  T.sum <- 0
  N.sum <- 0
  
  # For each condition, pull t.weighted, e and p for this antibiotic, group, and species
  cond.func <- function(cond.curr, group.curr){
    t.weighted <- sum(subset(NAMCS.summary, antibiotic==antibiotic_temp & condition==cond.curr & agegroup==group.curr)$wtVisits_curr)
    e <- etiologies.bs[which(etiologies.bs[,1]==species_temp), which(names(etiologies.bs)==paste0(cond.curr, "_", group.curr))]
    p <- e + (1-e)*microbiome.bs[which(microbiome.bs$Species.strain.Key==species_temp), which(names(microbiome.bs)==paste0("wtprev_", group.curr))]  
    
    return(c(t.weighted*p, t.weighted*e))
  }
  
  group.func <- function(group.curr){
    cond.result <- as.data.frame(lapply(conditions, function(x) cond.func(x,group.curr)))
    T.sum <- rowSums(cond.result)[1]
    N.sum <- rowSums(cond.result)[2]
    
    # Add "other" antibiotic use into the denominator (NOT for one of our prespecified conditions)
    p.sg <- microbiome.bs[which(microbiome.bs$Species.strain.Key==species_temp), paste0("wtprev_", group.curr)]
    T.sum <- T.sum + sum(subset(NAMCS.summary, condition=="other" & agegroup==group.curr & antibiotic==antibiotic_temp)$wtVisits_curr)*p.sg
    
    return(c(T.sum, N.sum))
  }
  
  group.result <- as.data.frame(lapply(groups, function(x) group.func(x)))
  T.sum <- rowSums(group.result)[1]
  N.sum <- rowSums(group.result)[2]
  
  bystander <- 1-(N.sum/T.sum)
  
  return(c(species_temp, antibiotic_temp, bystander, N.sum, T.sum))
}

bystander.cbar.df <- as.data.frame(t(mapply(bystanderCalc, species_temp=antibiotic.species$species, antibiotic_temp=antibiotic.species$antibiotic)))
rownames(bystander.cbar.df) <- NULL
colnames(bystander.cbar.df) <- c("Species", "Antibiotic", "Bystander_prop", "N.sum", "T.sum")
bystander.cbar.df[,3:ncol(bystander.cbar.df)] <- sapply(bystander.cbar.df[,3:ncol(bystander.cbar.df)], as.character)
bystander.cbar.df[,3:ncol(bystander.cbar.df)] <- sapply(bystander.cbar.df[,3:ncol(bystander.cbar.df)], as.numeric)
#write.table(bystander.cbar.df, "/Users/ctedijanto/Documents/03 Research/02 Bystander selection/01 Data/2018.09.24 bystander.cbar.df.txt", sep="\t")

# Bystander prop by AB class and species -------------------------------------------------------------
etiologies[is.na(etiologies)] <- 0
etiologies.bs <- etiologies
microbiome.bs <- microbiome
NAMCS.summary.byclass$wtVisits_curr <- NAMCS.summary.byclass$wtVisits 

bystanderCalc.byclass <- function(species_temp, drugclass_temp){
  T.sum <- 0
  N.sum <- 0
  
  cond.func <- function(cond.curr, group.curr){
    t.weighted <- sum(subset(NAMCS.summary.byclass, drugclass==drugclass_temp & condition==cond.curr & agegroup==group.curr)$wtVisits_curr)
    e <- etiologies.bs[which(etiologies.bs[,1]==species_temp), which(names(etiologies.bs)==paste0(cond.curr, "_", group.curr))]
    p <- e + (1-e)*microbiome.bs[which(microbiome.bs$Species.strain.Key==species_temp), which(names(microbiome.bs)==paste0("wtprev_", group.curr))]   
   
    return(c(t.weighted*p, t.weighted*e))
  }
  
  group.func <- function(group.curr){
    cond.result <- as.data.frame(lapply(conditions, function(x) cond.func(x,group.curr)))
    T.sum <- rowSums(cond.result)[1]
    N.sum <- rowSums(cond.result)[2]
    
    # Add "other" antibiotic use into the denominator (NOT for one of our prespecified conditions)
    p.sg <- microbiome.bs[which(microbiome.bs$Species.strain.Key==species_temp), paste0("wtprev_", group.curr)]
    T.sum <- T.sum + sum(subset(NAMCS.summary.byclass, condition=="other" & agegroup==group.curr & drugclass==drugclass_temp)$wtVisits_curr)*p.sg
    
    return(c(T.sum, N.sum))
  }
  
  group.result <- as.data.frame(lapply(groups, function(x) group.func(x)))
  T.sum <- rowSums(group.result)[1]
  N.sum <- rowSums(group.result)[2]
  
  bystander <- 1-(N.sum/T.sum)
  
  return(c(species_temp, drugclass_temp, bystander, N.sum, T.sum))
}

bystanderbyclass.cbar.FD.df <- as.data.frame(t(mapply(bystanderCalc.byclass, species_temp=drugclass.species$species, drugclass_temp=drugclass.species$drugclass)))
rownames(bystanderbyclass.cbar.FD.df) <- NULL
colnames(bystanderbyclass.cbar.FD.df) <- c("Species", "drugclass", "bystanderbyclass_prop", "N.sum", "T.sum")
bystanderbyclass.cbar.FD.df[,3:ncol(bystanderbyclass.cbar.FD.df)] <- sapply(bystanderbyclass.cbar.FD.df[,3:ncol(bystanderbyclass.cbar.FD.df)], as.character)
bystanderbyclass.cbar.FD.df[,3:ncol(bystanderbyclass.cbar.FD.df)] <- sapply(bystanderbyclass.cbar.FD.df[,3:ncol(bystanderbyclass.cbar.FD.df)], as.numeric)
write.table(bystanderbyclass.cbar.FD.df, "/Users/ctedijanto/Documents/03 Research/02 Bystander selection/01 Data/2018.09.27 bystanderbyclass.cbar.FD.df.txt", sep="\t")

# Produce wide table of bystander proportions ----------------------------------------------------------
library(tidyr)
temp <- bystander.cbar.df
temp$Bystander_prop <- round(temp$Bystander_prop,3)
View(spread(temp[,c(1:3)], Antibiotic, Bystander_prop))

temp <- bystanderbyclass.cbar.df
temp$bystanderbyclass_prop <- round(temp$bystanderbyclass_prop,3)
View(spread(temp[,c(1:3)], drugclass, bystanderbyclass_prop))

# Total exposures by AB, species and condition -------------------------------------------------------
totalExp.byCond <- function(){
  results <- data.frame(antibiotic=character(),
                        species=character(),
                        condition=character(),
                        exposures=double(),
                        stringsAsFactors = FALSE)
  
  # Add "other" to conditions (to account for conditions that are not explicitly included)
  totalExp.byCond.conditions <- c(conditions, "other")
  
  for (s in 1:length(species)){
    for (a in 1:(length(antibiotics)-2)){
      for (i in 1:(length(totalExp.byCond.conditions))){
        exp.sum <- 0
        for (j in 1:length(groups)){
          t.weighted <- sum(subset(NAMCS.summary, antibiotic==antibiotics[a] & condition==totalExp.byCond.conditions[i] & agegroup==groups[j])$wtVisits)
          if(totalExp.byCond.conditions[i]=="other"){
            e <- 0
          } else {
            e <- etiologies[which(etiologies[,1]==species[s]), which(names(etiologies)==paste0(totalExp.byCond.conditions[i], "_", groups[j]))]
          }
          
          p <- e + (1-e)*microbiome[which(microbiome$Species.strain.Key==species[s]), which(names(microbiome)==paste0("wtprev_", groups[j]))]  
          
          exp.sum <- exp.sum +  t.weighted*p
        }
        results[nrow(results)+1,1] <- antibiotics[a]
        results[nrow(results),2] <- species[s]
        results[nrow(results),3] <- totalExp.byCond.conditions[i]
        results[nrow(results),4] <- exp.sum
      }
    }
  }
  return(results)
}

etiologies[is.na(etiologies)] <- 0
totalExp.byCond.df <- totalExp.byCond()
totalExp.byCond.df$totalExp <- NULL
for (i in 1:nrow(totalExp.byCond.df)){
  totalExp.byCond.df$totalExp[i] <- sum(subset(totalExp.byCond.df, species==totalExp.byCond.df$species[i] & antibiotic==totalExp.byCond.df$antibiotic[i])$exposures)
}
totalExp.byCond.df$exp.prop <- totalExp.byCond.df$exposures/totalExp.byCond.df$totalExp
#write.table(totalExp.byCond.df, "/Users/ctedijanto/Documents/03 Research/02 Bystander selection/01 Data/2018.07.05 totalExp.txt", sep="\t")

# 4. BOOTSTRAP #############################################################
# Resample NAMCS/NHAMCS ---------------------------------------------------
# SE for natural log of weighted visits based on delta method
NAMCS.summary$wtVisits.ln <- log(NAMCS.summary$wtVisits)
NAMCS.summary$wtVisits.se.ln <- NAMCS.summary$wtVisits.se*(1/NAMCS.summary$wtVisits)

NAMCS.summary.byclass$wtVisits.ln <- log(NAMCS.summary.byclass$wtVisits)
NAMCS.summary.byclass$wtVisits.se.ln <- NAMCS.summary.byclass$wtVisits.se*(1/NAMCS.summary.byclass$wtVisits)

# Resample microbiome -----------------------------------------------------
# In the microbiome data, values of 0 maybe be given slight value due to use of Jeffreys (uninformative) prior
# Info on Jeffreys prior: http://www.stat.cmu.edu/~larry/=sml/Bayes.pdf
bootstrap_microbiome_kids <- function(kidsdf_temp){
  # For each bug reported from each of the studies (and age groups), resample using Jeffreys prior
  for(i in 6:ncol(kidsdf_temp)){
   p <- kidsdf_temp[,i]
   n <- kidsdf_temp$n
   alpha <- p*n
   beta <- n-alpha
   kidsdf_temp[,i] <- unlist(lapply(c(1:nrow(kidsdf_temp)), function(x) rbeta(1, alpha[x]+0.5, beta[x]+0.5)))
  }
  return(kidsdf_temp)
}

bootstrap_microbiome <- function(microbiome_df, kids0_temp, kids1to5_temp){
  
  # Bootstrap for kids
  kids0_temp <- bootstrap_microbiome_kids(kids0_temp)
  kids1to5_temp <- bootstrap_microbiome_kids(kids1to5_temp)
  wtprev_kids0_bs <- microbiomeprep.func(kids0_temp)
  names(wtprev_kids0_bs)[2] <- "wtprev_kids0"
  wtprev_kids1to5_bs <- microbiomeprep.func(kids1to5_temp)
  names(wtprev_kids1to5_bs)[2] <- "wtprev_kids1to5"
  
  # Bootstrap adult values
  # For each bug reported from each of the studies (and age groups), resample using Jeffreys prior
  for(i in 1:3){
    p <- microbiome_df[, paste0("prev_", i)]
    n <- microbiome_df[, paste0("n_", i)]
    alpha <- p*n
    beta <- n-alpha
    microbiome_df[,ncol(microbiome_df)+1] <- unlist(lapply(c(1:nrow(microbiome_df)), function(x) rbeta(1, alpha[x]+0.5, beta[x]+0.5)))
    names(microbiome_df)[ncol(microbiome_df)] <- paste("prev", i, "bs", sep="_")
  }
  microbiome_df$wtprev_adults <- microbiome_df$prev_1_bs*(microbiome_df$n_1/microbiome_df$n_123) +
                                 microbiome_df$prev_2_bs*(microbiome_df$n_2/microbiome_df$n_123) +
                                 microbiome_df$prev_3_bs*(microbiome_df$n_3/microbiome_df$n_123)
  
  # Merge microbiome datasets
  microbiome_temp <- merge(microbiome_df, wtprev_kids0_bs, by.x="Species.strain.Key", by.y="microbiomeprep.species")
  microbiome_df <- merge(microbiome_temp, wtprev_kids1to5_bs, by.x="Species.strain.Key", by.y="microbiomeprep.species")
  microbiome_df[,2:ncol(microbiome_df)] <- apply(microbiome_df[,2:ncol(microbiome_df)], 2, function(x) as.numeric(as.character(x)))
  
  # Imputing values
  microbiome_df$wtprev_kids1to5[which(microbiome_df$Species.strain.Key=="Pseudomonas_aeruginosa")] <- microbiome_df$wtprev_kids0[which(microbiome_df$Species.strain.Key=="Pseudomonas_aeruginosa")]
  microbiome_df$wtprev_kids1to5[which(microbiome_df$Species.strain.Key=="Streptococcus_agalactiae")] <- microbiome_df$wtprev_kids0[which(microbiome_df$Species.strain.Key=="Streptococcus_agalactiae")]
  
  return(microbiome_df)
}

# Resample etiologies -----------------------------------------------------
# For etiologies, we assume that values of 0 or 1 remain fixed (e.g. strep or viral URI) 
bootstrap_etiologies <- function(etiologies_df){
  for (i in 1:length(conditions)){
    # All resampled using Jeffreys prior (assume binary because that most studies note that the pathogen was either present or not during infection)
    adults.p <- etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "adults", sep="_"))]
    adults.n <- etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "adults", "n", sep="_"))]
    adults.alpha <- adults.p*adults.n
    adults.beta <- adults.n-adults.alpha
    
    # Conditions where etiological data was available for adults and kids -- resampled separately for both groups
    if(conditions[i] == "suppOtitisMedia" | conditions[i] == "pneumonia" | conditions[i] == "uti"){
      etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "adults", sep="_"))] <- unlist(lapply(c(1:nrow(etiologies_df)), function(x) rbeta(1, adults.alpha[x]+0.5, adults.beta[x]+0.5)))
      kids.p <- etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "kids1to5", sep="_"))]
      kids.n <- etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "kids1to5", "n", sep="_"))]
      kids.alpha <- kids.p*kids.n
      kids.beta <- kids.n-kids.alpha
      
      etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "kids1to5", sep="_"))] <- unlist(lapply(c(1:nrow(etiologies_df)), function(x) rbeta(1, kids.alpha[x]+0.5, kids.beta[x]+0.5)))
      etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "kids0", sep="_"))] <- etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "kids1to5", sep="_"))]
      
    # Conditions where etiological data is only available for adults -- resampled for adults and then applied to both children groups
    } else if(conditions[i] == "acuteSinusitis" | conditions[i] == "chronicSinusitis" | conditions[i] == "pharyngitis.nostrep" |
              conditions[i] == "cellulitis" |  conditions[i] == "bronchitis" | conditions[i] == "acne"){
      etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "adults", sep="_"))] <- unlist(lapply(c(1:nrow(etiologies_df)), function(x) rbeta(1, adults.alpha[x]+0.5, adults.beta[x]+0.5)))
      etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "kids1to5", sep="_"))] <- etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "adults", sep="_"))]
      etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "kids0", sep="_"))] <- etiologies_df[, which(names(etiologies_df)==paste(conditions[i], "adults", sep="_"))]
    } 
    
    # For rest of conditions, etiology is just based on assumption (e.g. 0 for all bacteria for viral URI) -- thus does not need to be resampled
  }
  return (etiologies_df)
}

# Run bootstrap! ----------------------------------------------------------
etiologies[etiologies==0] <- NA
bystanderbyclass.cbar.FD.df_bs <- bystanderbyclass.cbar.FD.df
bs_reps <- 1000
for(i in 1:bs_reps){
  NAMCS.summary.byclass$wtVisits_curr <- apply(NAMCS.summary.byclass, 1, 
                                               function(x) exp(rnorm(1, mean=as.numeric(x[which(names(NAMCS.summary.byclass)=="wtVisits.ln")]), sd=as.numeric(x[which(names(NAMCS.summary.byclass)=="wtVisits.se.ln")]))))
  NAMCS.summary.byclass[is.na(NAMCS.summary.byclass)] <- 0
  NAMCS.summary$wtVisits_curr <- apply(NAMCS.summary, 1, 
                                       function(x) exp(rnorm(1, mean=as.numeric(x[which(names(NAMCS.summary)=="wtVisits.ln")]), sd=as.numeric(x[which(names(NAMCS.summary)=="wtVisits.se.ln")]))))
  NAMCS.summary[is.na(NAMCS.summary)] <- 0
  etiologies.bs <- bootstrap_etiologies(etiologies)
  microbiome.bs <- bootstrap_microbiome(microbiome.adults, microbiome.kids0, microbiome.kids1to5)
  etiologies.bs[is.na(etiologies.bs)] <- 0
  temp <- as.data.frame(t(as.data.frame(mapply(bystanderCalc.byclass, species_temp=drugclass.species$species, drugclass_temp=drugclass.species$drugclass, SIMPLIFY = FALSE))))
  rownames(temp) <- NULL
  colnames(temp) <- c("Species", "drugclass", "bystanderbyclass_prop", "N.sum", "T.sum")
  temp[,3] <- as.numeric(as.character(temp[,3]))
  bystanderbyclass.cbar.FD.df_bs[,paste0("bs_",i)] <- temp[,3]
  print(i)
}  
bystanderbyclass.cbar.FD.df_bs$q2.5 <- apply(bystanderbyclass.cbar.FD.df_bs[,6:1005], 1, function(x) quantile(x, 0.025))
bystanderbyclass.cbar.FD.df_bs$q97.5 <- apply(bystanderbyclass.cbar.FD.df_bs[,6:1005], 1, function(x) quantile(x, 0.975)) 
write.table(bystanderbyclass.cbar.FD.df_bs, "/Users/ctedijanto/Documents/03 Research/02 Bystander selection/01 Data/2018.09.27 bystanderbyclass.cbar.FD.df_bs.txt", sep="\t")

etiologies[etiologies==0] <- NA
bystander.cbar.df_bs <- bystander.cbar.df
bs_reps <- 1000
for(i in 1:bs_reps){
  NAMCS.summary$wtVisits_curr <- apply(NAMCS.summary, 1, 
                                              function(x) exp(rnorm(1, mean=as.numeric(x[which(names(NAMCS.summary)=="wtVisits.ln")]), sd=as.numeric(x[which(names(NAMCS.summary)=="wtVisits.se.ln")]))))
  NAMCS.summary[is.na(NAMCS.summary)] <- 0
  etiologies.bs <- bootstrap_etiologies(etiologies)
  microbiome.bs <- bootstrap_microbiome(microbiome.adults, microbiome.kids0, microbiome.kids1to5)
  etiologies.bs[is.na(etiologies.bs)] <- 0
  temp <- as.data.frame(t(as.data.frame(mapply(bystanderCalc, species_temp=antibiotic.species$species, antibiotic_temp=antibiotic.species$antibiotic, SIMPLIFY = FALSE))))
  rownames(temp) <- NULL
  colnames(temp) <- c("Species", "antibiotic", "bystander_prop", "N.sum", "T.sum")
  temp[,3] <- as.numeric(as.character(temp[,3]))
  bystander.cbar.df_bs[,paste0("bs_",i)] <- temp[,3]
  print(i)
}
bystander.cbar.df_bs$q2.5 <- apply(bystander.cbar.df_bs[,6:1005], 1, function(x) quantile(x, 0.025))
bystander.cbar.df_bs$q97.5 <- apply(bystander.cbar.df_bs[,6:1005], 1, function(x) quantile(x, 0.975)) 
write.table(bystander.cbar.df_bs, "/Users/ctedijanto/Documents/03 Research/02 Bystander selection/01 Data/2018.09.24 bystander.cbar.df_bs.txt", sep="\t")
