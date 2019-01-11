# A1. SETUP - Install and load required packages ############################

if (!require(gdata)) install.packages('gdata')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(dplyr)) install.packages('dplyr')
if (!require(foreign)) install.packages('foreign')
if (!require(readstata13)) install.packages('readstata13')
if (!require(survey)) install.packages('survey')
if (!require(srvyr)) install.packages('srvyr')

library(gdata)
library(ggplot2)
library(dplyr)
library(foreign)
library(readstata13)
library(survey)
library(srvyr)

# A2. SETUP - Load required data ############################################

load("outpatient_bystander_analysis.Rdata")

# Load NAMCS and NHAMCS 2010-2011 Stata files (available at links below)
# NAMCS: http://www.nber.org/data/national-ambulatory-medical-care-survey.html
# NHAMCS: http://www.nber.org/data/national-hospital-ambulatory-medical-care-survey.html
namcs2010 <- read.dta13("namcs2010.dta")
namcs2011 <- read.dta13("namcs2011.dta")
nhamcsopd2010 <- read.dta13("nhamcsopd2010.dta")
nhamcsopd2011 <- read.dta13("nhamcsopd2011.dta")
nhamcsed2010 <- read.dta13("nhamcsed2010.dta")
nhamcsed2011 <- read.dta13("nhamcsed2011.dta")

# A3. SETUP - Define variables of interest ##################################

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
                 "SULFAMETHOXAZOLE-TRIMETHOPRIM", "NITROFURANTOIN", 
                 "DOXYCYCLINE", "TETRACYCLINE", "MINOCYCLINE",
                 "inclab", "anyab")

# Define drug classes of interest
# `drugclasses` vector matches respective class to antibiotics as defined in `antibiotics` vector
drugclasses <- c("PENICILLINS", "PENICILLINS", "PENICILLINS",
                 "MACROLIDES", "MACROLIDES", "MACROLIDES",
                 "QUINOLONES", "QUINOLONES", "QUINOLONES",
                 "CEPHALOSPORINS", "CEPHALOSPORINS", "CEPHALOSPORINS", 
                 "SULFAMETHOXAZOLE-TRIMETHOPRIM", "NITROFURANTOIN", 
                 "TETRACYCLINES", "TETRACYCLINES", "TETRACYCLINES", 
                 "inclab", "anyab")
drugclasses.key <- as.data.frame(cbind(antibiotic=antibiotics, drugclass=drugclasses))
drugclasses.unique <- unique(drugclasses)

# Define age groups of interest
groups <- c("kids0", "kids1to5", "adults")

# B1. DATA PREP - Microbiome #########################################################
# Goal of microbiome processing steps:
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
  
  # Take weighted average of prevalence for each study, species, and body site (over different age ranges)
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

# B2. DATA PREP - NAMCS/NHAMCS #########################################################
# i) Inclusion/exclusion -----------------------------------------------------------------------------
# To reduce size of dataframes, only keep necessary variables

# Check: sum raw number of visits across all datasets
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
# 'rx' variables are Level 2 Multum Lexicon classifications (each drug is allowed up to 4)
varstokeep <- c("drugid1", "drugid2", "drugid3", "drugid4", "drugid5", "drugid6", "drugid7", "drugid8",
               "rx1v2c1", "rx1v2c2", "rx1v2c3", "rx1v2c4","rx2v2c1", "rx2v2c2", "rx2v2c3", "rx2v2c4",
               "rx3v2c1", "rx3v2c2", "rx3v2c3", "rx3v2c4","rx4v2c1", "rx4v2c2", "rx4v2c3", "rx4v2c4",
               "rx5v2c1", "rx5v2c2", "rx5v2c3", "rx5v2c4","rx6v2c1", "rx6v2c2", "rx6v2c3", "rx6v2c4",
               "rx7v2c1", "rx7v2c2", "rx7v2c3", "rx7v2c4","rx8v2c1", "rx8v2c2", "rx8v2c3", "rx8v2c4",
               "med1", "med2", "med3", "med4", "med5", "med6", "med7", "med8",
               "diag1", "diag2", "diag3", "patwt", "cstratm", "cpsum", "age", "nonhospadm", "sex")
namcs2010.res <- namcs2010[varstokeep]
namcs2011.res <- namcs2011[varstokeep]
nhamcsopd2010.res <- nhamcsopd2010[varstokeep]
nhamcsopd2011.res <- nhamcsopd2011[varstokeep]
nhamcsed2010.res <- nhamcsed2010[varstokeep]
nhamcsed2011.res <- nhamcsed2011[varstokeep]

# ii) Pull in drug and diagnosis names ----------------------------------------------------------------------
# Pull in drug identifier codebooks (from SAS format commands) - be sure to update for the correct year if using updated data
# Some drug identifier codebooks may contain duplicates (cleaned in Excel)
DRUGID.10 <- read.csv("DRUGID.10.csv", header=FALSE, col.names=c("Code","Drug"), colClasses="character")
DRUGID.11 <- read.csv("DRUGID.11.csv", header=FALSE, col.names=c("Code","Drug"), colClasses="character")

# Pull in drug Level 2 codebooks
DRGLV2.10 <- read.csv("DRGLV2.10.csv", header=FALSE, col.names=c("Code","Category"), colClasses="character")

# 2011 has different drug Level 2 codebooks for NAMCS, NHAMCS OPD, and NHAMCS ED
DRGLV2.A11 <- read.csv("DRGLV2.A11.csv", header=FALSE, col.names=c("Code", "Category"), colClasses="character")
DRGLV2.H11 <- read.csv("DRGLV2.H11.csv", header=FALSE, col.names=c("Code", "Category"), colClasses="character")
DRGLV2.HED11 <- read.csv("DRGLV2.HED11.csv", header=FALSE, col.names=c("Code", "Category"), colClasses="character")

# Read in diagnosis identifier codebooks
DIAG.A10 <- read.csv("DIAG.A10.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")
DIAG.A11 <- read.csv("DIAG.A11.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")

DIAG.H10 <- read.csv("DIAG.H10.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")
DIAG.H11 <- read.csv("DIAG.H11.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")

DIAG.HED10 <- read.csv("DIAG.HED10.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")
DIAG.HED11 <- read.csv("DIAG.HED11.csv", header=FALSE, col.names=c("Code", "Diag"), colClasses="character")

# Function to pull in drug names, Multum Level 2 category names, and full diagnoses with appropriate codebook
decodeDrugDiag <- function(dataset, df, drugCodeBook, catCodeBook, diagCodeBook){
  df.temp <- df
  
  # Merge in drug names
  # Up to 8 drugs for each visit
  drugid1.1st <- which(names(df)=="drugid1")
  for(i in drugid1.1st:(drugid1.1st+7)){
    df.temp <- merge(df.temp, drugCodeBook, by.x=c(names(df[i])), by.y=c("Code"), all.x=TRUE, all.y=FALSE)
    df.temp[is.na(df.temp)] <- ""
    names(df.temp)[ncol(df.temp)] <- c(paste(names(df)[i],".name", sep=""))
  }
  
  # Merge in category names
  # Level 2 for 8 medications, each one can have up to 4 associated categories
  rx1v2c1.1st <- which(names(df)=="rx1v2c1")
  for(i in rx1v2c1.1st:(rx1v2c1.1st + 8*4 - 1)){
    df.temp <- merge(df.temp, catCodeBook, by.x=c(names(df[i])), by.y=c("Code"), all.x=TRUE, all.y=FALSE)
    df.temp[is.na(df.temp)] <- ""
    names(df.temp)[ncol(df.temp)] <- c(paste(names(df)[i],".name", sep=""))
  }
  
  # Merge in diagnosis names
  # Up to 3 diagnoses for each visit
  diag1.1st <- which(names(df)=="diag1")
  for (i in diag1.1st:(diag1.1st+2)){
    df.temp <- merge(df.temp, diagCodeBook, by.x=c(names(df[i])), by.y=c("Code"), all.x=TRUE, all.y=FALSE)
    df.temp[is.na(df.temp)] <- ""
    names(df.temp)[ncol(df.temp)] <- c(paste(names(df)[i],".name", sep=""))
  }
  
  # Label with name of dataset
  df.temp <- cbind(dataset=dataset, df.temp)
  
  return(df.temp)
}

namcs2010.drugdecode <- decodeDrugDiag("namcs2010", namcs2010.res, DRUGID.10, DRGLV2.10, DIAG.A10)
namcs2011.drugdecode <- decodeDrugDiag("namcs2011", namcs2011.res, DRUGID.11, DRGLV2.A11, DIAG.A11)
nhamcsopd2010.drugdecode <- decodeDrugDiag("nhamcsopd2010", nhamcsopd2010.res, DRUGID.10, DRGLV2.10, DIAG.H10)
nhamcsopd2011.drugdecode <- decodeDrugDiag("nhamcsopd2011", nhamcsopd2011.res, DRUGID.11, DRGLV2.H11, DIAG.H11)
nhamcsed2010.drugdecode <- decodeDrugDiag("nhamcsed2010", nhamcsed2010.res, DRUGID.10, DRGLV2.10, DIAG.HED10)
nhamcsed2011.drugdecode <- decodeDrugDiag("nhamcsed2011", nhamcsed2011.res, DRUGID.11, DRGLV2.HED11, DIAG.HED11)

# iii) Add antibiotic flags ('antibiotic<med#>_flag') ----------------------------------------------------------------------------
# Create flag for each medication (1-8) where 1=med is antibiotic, 0=otherwise (based on specification within anti-infectives, see NAMCS documentation below)
# Additional NAMCS/NHAMCS resource: https://www.cdc.gov/nchs/ppt/nchs2010/13_schappert.pdf
# To be used for "anyab" calculations -- not used directly in bystander paper

antibioticsflag <- function(df){
  rx1v2c1.name.1st <- which(names(df)=="rx1v2c1.name")
  
  # Create empty dataframe to store flags (same number of rows as df, 8 columns)
  antibioticsflag.df <- data.frame(matrix(vector(), nrow(df), 8))
  
  # j represents the antibiotic # (can record up to 8 medications)
  for (j in 0:7){
    temp <- rep(0, nrow(df))
    
    # For each med (1-8), look through all level 2 categories (up to 4 may exist)
    # Add if any of them is an "Anti-infective" but NOT one of the other categories listed
    for (i in c((rx1v2c1.name.1st+4*j):(rx1v2c1.name.1st+4*j+3))){
      temp <- temp + (grepl("Anti-infectives", df[,i]) & 
                        !grepl("amebicides", df[,i]) & 
                        !grepl("anthelmintics", df[,i]) & 
                        !grepl("antifungals", df[,i]) & 
                        !grepl("antimalarial agents", df[,i]) & 
                        !grepl("antiviral agents", df[,i]))
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

# iv) Create flag for inclusion criteria ('inclusion') ----------------------------------------------------------------------------
# Create overall inclusion flag (1 to include, 0 otherwise)
# Use `nonhospadm` flag (previously created) to remove visits that resulted in inpatient admission

namcs2010.drugdecode.ab$inclusion <- namcs2010.drugdecode.ab$nonhospadm
namcs2011.drugdecode.ab$inclusion <- namcs2011.drugdecode.ab$nonhospadm
nhamcsopd2010.drugdecode.ab$inclusion <- nhamcsopd2010.drugdecode.ab$nonhospadm
nhamcsopd2011.drugdecode.ab$inclusion <- nhamcsopd2011.drugdecode.ab$nonhospadm
nhamcsed2010.drugdecode.ab$inclusion <- nhamcsed2010.drugdecode.ab$nonhospadm
nhamcsed2011.drugdecode.ab$inclusion <- nhamcsed2011.drugdecode.ab$nonhospadm

# Count total included visits
totalVisits.incl <- sum(namcs2010.drugdecode.ab$inclusion, namcs2011.drugdecode.ab$inclusion,
                    nhamcsopd2010.drugdecode.ab$inclusion, nhamcsopd2011.drugdecode.ab$inclusion,
                    nhamcsed2010.drugdecode.ab$inclusion, nhamcsed2011.drugdecode.ab$inclusion)

# v) Add flags for antibiotics of interest ('<ANTIBIOTIC>_flag')----------------------------------------------------------------------------
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

namcs2010.drugdecode.ab2 <- specificABflag(namcs2010.drugdecode.ab)
namcs2011.drugdecode.ab2 <- specificABflag(namcs2011.drugdecode.ab)
nhamcsopd2010.drugdecode.ab2 <- specificABflag(nhamcsopd2010.drugdecode.ab)
nhamcsopd2011.drugdecode.ab2 <- specificABflag(nhamcsopd2011.drugdecode.ab)
nhamcsed2010.drugdecode.ab2 <- specificABflag(nhamcsed2010.drugdecode.ab)
nhamcsed2011.drugdecode.ab2 <- specificABflag(nhamcsed2011.drugdecode.ab)

# vi) Add counts for drug classes of interest ('<CLASS>_count') ---------------------------------------------------------
# Count of number of drugs of that drug class prescribed at each visit
specificABclasscount <- function(df){
  rx1v2c1.name.1st <- which(names(df)=="rx1v2c1.name")
  
  # Identifiers (based on level 2 of MULTUM Lexicon classification system) for each drug class within anti-infectives
  # Order of identifiers here matches with `drugclasses.unique` -- MULTUM identifier is irrelevant for ones with "x"
  multum.names <- c("penicillins", "macrolide derivatives", "quinolones", "cephalosporins", "x", "x", "tetracyclines", "x", "x")
  
  for(i in c(1:4,7)){
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

namcs2010.drugdecode.ab3 <- specificABclasscount(namcs2010.drugdecode.ab2)
namcs2011.drugdecode.ab3 <- specificABclasscount(namcs2011.drugdecode.ab2)
nhamcsopd2010.drugdecode.ab3 <- specificABclasscount(nhamcsopd2010.drugdecode.ab2)
nhamcsopd2011.drugdecode.ab3 <- specificABclasscount(nhamcsopd2011.drugdecode.ab2)
nhamcsed2010.drugdecode.ab3 <- specificABclasscount(nhamcsed2010.drugdecode.ab2)
nhamcsed2011.drugdecode.ab3 <- specificABclasscount(nhamcsed2011.drugdecode.ab2)

# vii) Add counts for overall antibiotics ('inclab_count' and 'anyab_count') ---------------------------------------------------------
# 'inclab_count' = count of antibiotic in any of the explicitly included classes prescribed at that visit
# 'anyab_count' = count of any antibiotic prescribed at that visit, as identified by Multum Lexicon system in NAMCS/NHAMCS
abcounts <- function(df){
  inclab.temp <- rep(0, nrow(df))
  #drugclasses.unique-2 because does not include "inclab" or "anyab"
  for (i in 1:(length(drugclasses.unique)-2)){
    inclab.temp <- inclab.temp + df[,paste0(drugclasses.unique[i], "_count")]
  }
  df$inclab_count <- inclab.temp
  
  df$anyab_count <- df$antibioticflag_1 + df$antibioticflag_2 +
                    df$antibioticflag_3 + df$antibioticflag_4 +
                    df$antibioticflag_5 + df$antibioticflag_6 + 
                    df$antibioticflag_7 + df$antibioticflag_8
  return(df)
}

namcs2010.drugdecode.ab4 <- abcounts(namcs2010.drugdecode.ab3)
namcs2011.drugdecode.ab4 <- abcounts(namcs2011.drugdecode.ab3)
nhamcsopd2010.drugdecode.ab4 <- abcounts(nhamcsopd2010.drugdecode.ab3)
nhamcsopd2011.drugdecode.ab4 <- abcounts(nhamcsopd2011.drugdecode.ab3)
nhamcsed2010.drugdecode.ab4 <- abcounts(nhamcsed2010.drugdecode.ab3)
nhamcsed2011.drugdecode.ab4 <- abcounts(nhamcsed2011.drugdecode.ab3)

# viii) Add flags for conditions of interest ('<CONDITION>_flag') ----------------------------------------------------------------------------
# Add flag (1=diagnosed during that visit, 0=otherwise) based on codes in condition-code key
# Condition-code key contains diagnosis codes to include in that condition and diagnosis codes to exclude
# E.g. Fleming-Dutra exclude diagnoses of bronchitis that have a co-diagnosis of COPD
#'creategreplstring' function creates the string of code (from condition-code key) needed to appropriately include/exclude codes
creategreplstring <- function(df, inclcode, exclcode){
  m <- length(inclcode)
  n <- length(exclcode)

  # 'grepl' marks as TRUE if string (e.g. df$diag1) includes the given string segment
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

namcs2010.drugdecode.ab.cond <- specificconditionflag(namcs2010.drugdecode.ab4)
namcs2011.drugdecode.ab.cond <- specificconditionflag(namcs2011.drugdecode.ab4)
nhamcsopd2010.drugdecode.ab.cond <- specificconditionflag(nhamcsopd2010.drugdecode.ab4)
nhamcsopd2011.drugdecode.ab.cond <- specificconditionflag(nhamcsopd2011.drugdecode.ab4)
nhamcsed2010.drugdecode.ab.cond <- specificconditionflag(nhamcsed2010.drugdecode.ab4)
nhamcsed2011.drugdecode.ab.cond <- specificconditionflag(nhamcsed2011.drugdecode.ab4)

# ix) ALTERNATE viii creating tiered diagnosis (from Fleming-Dutra et al. 2016) ----------------------------------------------------------------------------
# 
# # 'creategreplstring' function creates the string of code needed to appropriately include/exclude codes
# # In order to do tiered diagnosis analysis, have changed this to be diag# specific (instead of looking across any of the diagnoses per row)
# creategreplstring <- function(df, diagnum, inclcode, exclcode){
#   m <- length(inclcode)
#   n <- length(exclcode)
# 
#   # 'grepl' marks as TRUE if string (e.g. df$diag1) includes the given string segment
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
# # Alter to create a condition flag for each listed diagnosis (instead of for each visit)
# specificconditionflag <- function(df){
#   for (i in 1:length(conditions)){
#     # for diagnosis numbers
#     for (j in 1:3){
#       incl <- subset(conditionCodes, conditionCodes$Condition==conditions[i] & conditionCodes$Exclude==0)$Code
#       excl <- subset(conditionCodes, conditionCodes$Condition==conditions[i] & conditionCodes$Exclude==1)$Code
#       temp <- creategreplstring(df, j, incl, excl)
#       # Create new columns called diag1.condition, diag2.condition, diag3.condition
#       # When grepl function returns true, replace value in appropriate diag<#>.condition with name of condition
#       df[which(temp==TRUE), paste0("diag",j,".condition")] <- conditions[i]
#     }
#   }
#   return(df)
# }
# 
# namcs2010.drugdecode.ab.cond1 <- specificconditionflag(namcs2010.drugdecode.ab4)
# namcs2011.drugdecode.ab.cond1 <- specificconditionflag(namcs2011.drugdecode.ab4)
# nhamcsopd2010.drugdecode.ab.cond1 <- specificconditionflag(nhamcsopd2010.drugdecode.ab4)
# nhamcsopd2011.drugdecode.ab.cond1 <- specificconditionflag(nhamcsopd2011.drugdecode.ab4)
# nhamcsed2010.drugdecode.ab.cond1 <- specificconditionflag(nhamcsed2010.drugdecode.ab4)
# nhamcsed2011.drugdecode.ab.cond1 <- specificconditionflag(nhamcsed2011.drugdecode.ab4)
# 
# # Create new columns diag1.tier, diag2.tier, diag3.tier which contain the appropriate tier from Fleming-Dutra 2016
# # Antibiotics almost always appropriate for tier 1 and less so for 2 and 3
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
#   # For each row of the dataframe, find the minimum tier of all diagnoses conditions
#   df$mintier <- unlist(lapply(1:nrow(df), function(x) min(df$diag1.tier[x], df$diag2.tier[x], df$diag3.tier[x], na.rm=TRUE)))
#   # If `mintier`` is Inf, there were no tiers available
#   df$mintier[which(df$mintier==Inf)] <- NA
#   return(df)
# }
# 
# namcs2010.drugdecode.ab.cond2 <- assigntier(namcs2010.drugdecode.ab.cond1)
# namcs2011.drugdecode.ab.cond2 <- assigntier(namcs2011.drugdecode.ab.cond1)
# nhamcsopd2010.drugdecode.ab.cond2 <- assigntier(nhamcsopd2010.drugdecode.ab.cond1)
# nhamcsopd2011.drugdecode.ab.cond2 <- assigntier(nhamcsopd2011.drugdecode.ab.cond1)
# nhamcsed2010.drugdecode.ab.cond2 <- assigntier(nhamcsed2010.drugdecode.ab.cond1)
# nhamcsed2011.drugdecode.ab.cond2 <- assigntier(nhamcsed2011.drugdecode.ab.cond1)
# 
# # Create new column called `finalcondition` is a list of the conditions at the same visit with the lowest tier
# # If two conditions have the minimum tier, then the one that is listed first is used
# pullcondition <- function(df){
#   # Add 1 to incorporate the 'NA' (if not found, `which` returns integer(0))
#   df$finalcondition <- unlist(lapply(1:nrow(df),
#                               function(x) c(NA, df$diag1.condition[x], df$diag2.condition[x], df$diag3.condition[x])[min(which(c(df$diag1.tier[x], df$diag2.tier[x], df$diag3.tier[x])==df$mintier[x]))+1]))
#   return(df)
# }
# 
# namcs2010.drugdecode.ab.cond3 <- pullcondition(namcs2010.drugdecode.ab.cond2)
# namcs2011.drugdecode.ab.cond3 <- pullcondition(namcs2011.drugdecode.ab.cond2)
# nhamcsopd2010.drugdecode.ab.cond3 <- pullcondition(nhamcsopd2010.drugdecode.ab.cond2)
# nhamcsopd2011.drugdecode.ab.cond3 <- pullcondition(nhamcsopd2011.drugdecode.ab.cond2)
# nhamcsed2010.drugdecode.ab.cond3 <- pullcondition(nhamcsed2010.drugdecode.ab.cond2)
# nhamcsed2011.drugdecode.ab.cond3 <- pullcondition(nhamcsed2011.drugdecode.ab.cond2)
# 
# # Create `<CONDITION_flag>`` as in untiered version
# finalconditionflag <- function(df){
#   for(i in 1:length(conditions)){
#     temp <- (df$finalcondition==conditions[i])
#     temp[is.na(temp)] <- FALSE
#     df[,paste0(conditions[i], "_flag")] <- as.numeric(temp)
#   }
#   return(df)
# }
# 
# namcs2010.drugdecode.ab.cond <- finalconditionflag(namcs2010.drugdecode.ab.cond3)
# namcs2011.drugdecode.ab.cond <- finalconditionflag(namcs2011.drugdecode.ab.cond3)
# nhamcsopd2010.drugdecode.ab.cond <- finalconditionflag(nhamcsopd2010.drugdecode.ab.cond3)
# nhamcsopd2011.drugdecode.ab.cond <- finalconditionflag(nhamcsopd2011.drugdecode.ab.cond3)
# nhamcsed2010.drugdecode.ab.cond <- finalconditionflag(nhamcsed2010.drugdecode.ab.cond3)
# nhamcsed2011.drugdecode.ab.cond <- finalconditionflag(nhamcsed2011.drugdecode.ab.cond3)

# x) Add flag and count for any condition of interest ('inclcond_flag', 'inclcond_count', 'other_flag') ----------------------------------------------------------------------------
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

namcs2010.drugdecode.ab.cond.incl <- inclcondfunc(namcs2010.drugdecode.ab.cond)
namcs2011.drugdecode.ab.cond.incl <-inclcondfunc(namcs2011.drugdecode.ab.cond)
nhamcsopd2010.drugdecode.ab.cond.incl <- inclcondfunc(nhamcsopd2010.drugdecode.ab.cond)
nhamcsopd2011.drugdecode.ab.cond.incl <- inclcondfunc(nhamcsopd2011.drugdecode.ab.cond)
nhamcsed2010.drugdecode.ab.cond.incl <- inclcondfunc(nhamcsed2010.drugdecode.ab.cond)
nhamcsed2011.drugdecode.ab.cond.incl <- inclcondfunc(nhamcsed2011.drugdecode.ab.cond)

# xi) Add flags for age groups of interest ----------------------------------------------------------------------------

agegroupsfunc <- function(df){
  df$agegroup[df$age==0] <- "kids0"
  df$agegroup[df$age>0 & df$age<=5] <- "kids1to5"
  df$agegroup[df$age>5] <- "adults"
  return(df)
}

namcs2010.drugdecode.ab.cond.grp <- agegroupsfunc(namcs2010.drugdecode.ab.cond.incl)
namcs2011.drugdecode.ab.cond.grp <- agegroupsfunc(namcs2011.drugdecode.ab.cond.incl)
nhamcsopd2010.drugdecode.ab.cond.grp <- agegroupsfunc(nhamcsopd2010.drugdecode.ab.cond.incl)
nhamcsopd2011.drugdecode.ab.cond.grp <- agegroupsfunc(nhamcsopd2011.drugdecode.ab.cond.incl)
nhamcsed2010.drugdecode.ab.cond.grp <- agegroupsfunc(nhamcsed2010.drugdecode.ab.cond.incl)
nhamcsed2011.drugdecode.ab.cond.grp <- agegroupsfunc(nhamcsed2011.drugdecode.ab.cond.incl)

# B3. DATA PREP - Create NAMCS/NHAMCS '.summary' dataframes --------------------------------------------------------------------------------
# Calculate means and CI for proportions of diagnoses with given medication using 'survey' package
# R documentation: https://cran.r-project.org/web/packages/survey/survey.pdf
# For multistage sampling:
# id = formula with the cluster identifiers at each stage
# strata = if subsequent stages are stratified, should be a formula with stratum identifiers at each stage
# fpc = population size for each level of sampling; if not specified, then only does one-level with replacement
# nest = TRUE --> PSUs reuse the same identifiers across strata

# According to this documentation: https://www.cdc.gov/nchs/ppt/nchs2015/Hing_Monday_GlenEcho_C2.pdf
# NHAMCS researchers should use combined ED and OPD files when computing variances for emergency and/or outpatient department estimates. 
# Including both files takes into account NHAMCS’s complex sample design. Thus, in this step, we combine NHAMCS files
# "names..." ensures that columns are in the same order
nhamcs2010.drugdecode.ab.cond.grp <- rbind(nhamcsopd2010.drugdecode.ab.cond.grp, nhamcsed2010.drugdecode.ab.cond.grp[,names(nhamcsopd2010.drugdecode.ab.cond.grp)])
nhamcs2011.drugdecode.ab.cond.grp <- rbind(nhamcsopd2011.drugdecode.ab.cond.grp, nhamcsed2011.drugdecode.ab.cond.grp[,names(nhamcsopd2011.drugdecode.ab.cond.grp)])

# `surveyCalcs` function creates summaries of NAMCS/NHAMCS with weighted visits by condition and antibiotic and corresponding SE (for bootstrapping)
# options(survey.lonely.psu=) accounts for strata that have a single PSU
# Based on R documentation: http://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
# With options(survey.lonely.psu="adjust") the data for the single-PSU stratum are centered at the sample grand mean rather than the stratum mean. 
# This is conservative.
# Includes "other" as a condition - visits where none of our included conditions were diagnosed.
options(survey.lonely.psu = "adjust")
surveyCalcs <- function(dataset, class){
  # Create empty dataframe for results
  results <- data.frame(dataset=character(),
                        condition=character(),
                        agegroup=character(),
                        antibiotic=character(),
                        total=double(),
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
      # Select only those rows which meet inclusion criteria and which had the condition of interest
      tempsurvey.res <- subset(tempsurvey, inclusion==1 & tempsurvey[,which(names(tempsurvey)=="get(tempcondflag)")]==TRUE)
      rownames(tempsurvey.res) <- NULL
      
      # For the case where there are NO visits with that condition in any age group (completely empty dataframe):
      if(nrow(tempsurvey.res)==0){tempsurvey.res[1,] <- rep(0,ncol(tempsurvey.res))}
      
      # Now add all other important data columns
      tempsurvey.res <- merge(data.frame(agegroup=groups), tempsurvey.res, by.x="agegroup", by.y="agegroup", all.x=TRUE)
      tempsurvey.res$dataset <- dataset
      if(class=="no"){tempsurvey.res$antibiotic <- antibiotics[i]}else{tempsurvey.res$antibiotic <- drugclasses.unique[i]}
      tempsurvey.res$condition <- surveyCalcs.conditions[j]
      tempsurvey.res[is.na(tempsurvey.res)] <- 0
      
      # Results have different column names depending on whether column is binary (only contains 0 and 1) or count (more values than 0 and 1)
      # Keep columns of interest
      if((class=="no") | (class=="yes" & (drugclasses.unique[i]=="NITROFURANTOIN"|drugclasses.unique[i]=="SULFAMETHOXAZOLE-TRIMETHOPRIM"))){
        tempresults <- tempsurvey.res[, c("dataset", "condition", "agegroup", "antibiotic","get(tempabflag)TRUE", "se.get(tempabflag)TRUE")]
      } else {
        tempresults <- tempsurvey.res[, c("dataset", "condition", "agegroup", "antibiotic", "get(tempabflag)", "se")]
      }
      
      # Rename columns
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



# i) Merge NAMCS/NHAMCS dataframes ----------------------------------------------------
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

# C1. BYSTANDER ANALYSIS - Setup ###################################################################
# Weighted visits are used for all bystander analyses
# Set up dataframes for mapply
# Remove last two ("inclab" and "anyab")
antibiotic.species <- expand.grid(antibiotics[1:(length(antibiotics)-2)], species)
names(antibiotic.species) <- c("antibiotic", "species")
antibiotic.species[,c(1:2)] <- sapply(antibiotic.species[,c(1:2)], as.character)

drugclass.species <- expand.grid(drugclasses.unique, species)
names(drugclass.species) <- c("drugclass", "species")
drugclass.species[,c(1:2)] <- sapply(drugclass.species[,c(1:2)], as.character)

# C2. BYSTANDER ANALYSIS - Bystander prop by AB and species -------------------------------------------------------------------
# Calculate proportion of bystander exposures for a given antibiotic and species, summed across all conditions and groups
bystanderCalc <- function(species_temp, antibiotic_temp){
  # Initialize T (denominator) and N (numerator) at 0
  T.sum <- 0
  N.sum <- 0
  
  # For each condition, pull t.weighted, e and p for this antibiotic, group, and species
  cond.func <- function(cond.curr, group.curr){
    t.weighted <- sum(subset(NAMCS.summary, antibiotic==antibiotic_temp & condition==cond.curr & agegroup==group.curr)$wtVisits)
    e <- etiologies[which(etiologies[,1]==species_temp), which(names(etiologies)==paste0(cond.curr, "_", group.curr))]
    p <- e + (1-e)*microbiome[which(microbiome$Species.strain.Key==species_temp), which(names(microbiome)==paste0("wtprev_", group.curr))]  
    
    return(c(t.weighted*p, t.weighted*e))
  }
  
  # Calculate numerator and denominator for each group
  group.func <- function(group.curr){
    cond.result <- as.data.frame(lapply(conditions, function(x) cond.func(x,group.curr)))
    T.sum <- rowSums(cond.result)[1]
    N.sum <- rowSums(cond.result)[2]
    
    # Add "other" antibiotic use into the denominator (NOT for one of our prespecified conditions)
    p.sg <- microbiome[which(microbiome$Species.strain.Key==species_temp), paste0("wtprev_", group.curr)]
    T.sum <- T.sum + sum(subset(NAMCS.summary, condition=="other" & agegroup==group.curr & antibiotic==antibiotic_temp)$wtVisits)*p.sg
    
    return(c(T.sum, N.sum))
  }
  
  group.result <- as.data.frame(lapply(groups, function(x) group.func(x)))
  T.sum <- rowSums(group.result)[1]
  N.sum <- rowSums(group.result)[2]
  
  bystander <- 1-(N.sum/T.sum)
  
  return(c(species_temp, antibiotic_temp, bystander, N.sum, T.sum))
}

bystander.df <- as.data.frame(t(mapply(bystanderCalc, species_temp=antibiotic.species$species, antibiotic_temp=antibiotic.species$antibiotic)))
rownames(bystander.df) <- NULL
colnames(bystander.df) <- c("species", "antibiotic", "bystander_prop", "N.sum", "T.sum")
bystander.df[,3:ncol(bystander.df)] <- sapply(bystander.df[,3:ncol(bystander.df)], as.character)
bystander.df[,3:ncol(bystander.df)] <- sapply(bystander.df[,3:ncol(bystander.df)], as.numeric)

# C3. BYSTANDER ANALYSIS -  Bystander prop by AB class and species -------------------------------------------------------------
bystanderCalc.byclass <- function(species_temp, drugclass_temp){
  T.sum <- 0
  N.sum <- 0
  
  cond.func <- function(cond.curr, group.curr){
    t.weighted <- sum(subset(NAMCS.summary.byclass, drugclass==drugclass_temp & condition==cond.curr & agegroup==group.curr)$wtVisits)
    e <- etiologies[which(etiologies[,1]==species_temp), which(names(etiologies)==paste0(cond.curr, "_", group.curr))]
    p <- e + (1-e)*microbiome[which(microbiome$Species.strain.Key==species_temp), which(names(microbiome)==paste0("wtprev_", group.curr))]   
   
    return(c(t.weighted*p, t.weighted*e))
  }
  
  group.func <- function(group.curr){
    cond.result <- as.data.frame(lapply(conditions, function(x) cond.func(x,group.curr)))
    T.sum <- rowSums(cond.result)[1]
    N.sum <- rowSums(cond.result)[2]
    
    # Add "other" antibiotic use into the denominator (NOT for one of our prespecified conditions)
    p.sg <- microbiome[which(microbiome$Species.strain.Key==species_temp), paste0("wtprev_", group.curr)]
    T.sum <- T.sum + sum(subset(NAMCS.summary.byclass, condition=="other" & agegroup==group.curr & drugclass==drugclass_temp)$wtVisits)*p.sg
    
    return(c(T.sum, N.sum))
  }
  
  group.result <- as.data.frame(lapply(groups, function(x) group.func(x)))
  T.sum <- rowSums(group.result)[1]
  N.sum <- rowSums(group.result)[2]
  
  bystander <- 1-(N.sum/T.sum)
  
  return(c(species_temp, drugclass_temp, bystander, N.sum, T.sum))
}

bystanderbyclass.df <- as.data.frame(t(mapply(bystanderCalc.byclass, species_temp=drugclass.species$species, drugclass_temp=drugclass.species$drugclass)))
rownames(bystanderbyclass.df) <- NULL
colnames(bystanderbyclass.df) <- c("species", "drugclass", "bystanderbyclass_prop", "N.sum", "T.sum")
bystanderbyclass.df[,3:ncol(bystanderbyclass.df)] <- sapply(bystanderbyclass.df[,3:ncol(bystanderbyclass.df)], as.character)
bystanderbyclass.df[,3:ncol(bystanderbyclass.df)] <- sapply(bystanderbyclass.df[,3:ncol(bystanderbyclass.df)], as.numeric)
