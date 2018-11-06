### Bystander selection in the outpatient setting
*Author*: Christine Tedijanto  
*Email*: ctedijanto@g.harvard.edu  
*Last update*: November 6, 2018  
*Publication*: [Estimating the proportion of bystander selection for antibiotic resistance among potentially pathogenic bacterial flora](https://www.biorxiv.org/content/early/2018/09/28/288704) (*in press at PNAS*)

#### Description:
Antibiotic use creates a selective pressure for resistance. In this paper, we estimate the extent of bystander selection, which we define as selective pressures experienced by an organism due to antibiotics that were intended to treat another pathogen. We estimate a metric termed the "proportion of bystander exposures" using existing data from the National Ambulatory Medical Care Survey (NAMCS) and the National Hospital Ambulatory Medical Care Survey (NHAMCS), the Human Microbiome Project (HMP), and published carriage and etiology studies. This repository contains the analysis code and data needed to estimate the proportion of bystander exposures for 9 commonly carried bacterial species and 17 commonly used antibiotics in the United States.

#### Data:
* The `.Rdata` file contains the relevant Human Microbiome Project, ICD-9 code, and etiology data that are required for this project. The included dataframes are detailed below:
    1. `ICD9v28`: ICD9 diagnosis codes (v28) with long and short description from [NBER](http://www.nber.org/data/icd-9-cm-diagnosis-and-procedure-codes-and-titles.html).
    2. `conditionCodes`: ICD9 diagnosis codes to include or exclude in condition categories. Based on conditions in [Fleming-Dutra et al. 2016](https://jamanetwork.com/journals/jama/fullarticle/2518263). Codes are written using regular expressions ([regex](https://medium.com/factory-mind/regex-tutorial-a-simple-cheatsheet-by-examples-649dc1c3f285)) for use in `grepl` function.
    3. `etiologies`: Estimated etiologies by condition, species, and age group. See Supporting Information in paper for sources.
    4. `microbiome.adults`: Carriage prevalence by species and visit number among adults participating in the Human Microbiome Project.
    5. `microbiome.kids0`: Carriage prevalence by species, body site, and study among children under 1yo. See Supporting Information in paper for sources.
    6. `microbiome.kids1to5`: Carriage prevalence by species, body site, and study measured in children 1-5yo. See Supporting Information in paper for sources.  
<br>
* Data from the National Ambulatory Medical Care Survey (NAMCS) and the National Hospital Ambulatory Medical Care Survey (NHAMCS) are made publicly available by the CDC [here](https://www.cdc.gov/nchs/ahcd/datasets_documentation_related.htm). This analysis uses the Stata files from the 2010 and 2011 surveys.

#### Instructions:
1. Place `outpatient_bystander_analysis.R`, `outpatient_bystander_analysis.Rdata`, and downloaded NAMCS/NHAMCS Stata files in the same directory.



