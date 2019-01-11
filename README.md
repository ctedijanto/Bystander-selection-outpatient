### Bystander selection in the outpatient setting
*Author*: Christine Tedijanto  
*Email*: ctedijanto@g.harvard.edu  
*Last update*: January 10, 2019  
*Publication*: [Estimating the proportion of bystander selection for antibiotic resistance among potentially pathogenic bacterial flora](https://www.pnas.org/lookup/doi/10.1073/pnas.1810840115)

#### Description:
Antibiotic use creates a selective pressure for resistance. In this paper, we estimate the extent of bystander selection, which we define as selective pressures experienced by an organism due to antibiotics that were intended to treat another pathogen. We estimate a metric termed the "proportion of bystander exposures" using existing data from the National Ambulatory Medical Care Survey (NAMCS) and the National Hospital Ambulatory Medical Care Survey (NHAMCS), the Human Microbiome Project (HMP), and published carriage and etiology studies. This repository containss the analysis code and data needed to estimate the proportion of bystander exposures for 9 commonly carried bacterial species and 17 commonly used antibiotics in the United States.

#### Data:
* The `.Rdata` file contains the relevant Human Microbiome Project, ICD-9 code, and etiology data that are required for this project. The included dataframes are detailed below:  

    1. `ICD9v28`: ICD9 diagnosis codes (v28) with long and short description from [NBER](http://www.nber.org/data/icd-9-cm-diagnosis-and-procedure-codes-and-titles.html).
    2. `conditionCodes`: ICD9 diagnosis codes to include or exclude in condition categories. Based on conditions in [Fleming-Dutra et al. 2016](https://jamanetwork.com/journals/jama/fullarticle/2518263). Codes are written using regular expressions ([regex](https://medium.com/factory-mind/regex-tutorial-a-simple-cheatsheet-by-examples-649dc1c3f285)) for use in `grepl` function.
    3. `etiologies`: Estimated etiologies by condition, species, and age group. See Supporting Information in paper for sources.
    4. `microbiome.adults`: Carriage prevalence by species and visit number among adults participating in the Human Microbiome Project. Estimates for several species were based on additional carriage studies; see Supporting Information in paper for sources.
    5. `microbiome.kids0`: Carriage prevalence by species, body site, and study among children under 1yo. See Supporting Information in paper for sources.
    6. `microbiome.kids1to5`: Carriage prevalence by species, body site, and study measured in children 1-5yo. See Supporting Information in paper for sources.  
<br>
* Data from the National Ambulatory Medical Care Survey (NAMCS) and the National Hospital Ambulatory Medical Care Survey (NHAMCS) are made publicly available by the CDC [here](https://www.cdc.gov/nchs/ahcd/datasets_documentation_related.htm). This analysis uses the Stata files from the 2010 and 2011 surveys.

#### Instructions:
Place `outpatient_bystander_analysis.R`, `outpatient_bystander_analysis.Rdata`, and NAMCS/NHAMCS Stata files in the same directory. For baseline analysis, code chunk B2(ix) *ALTERNATE viii creating tiered diagnosis (from Fleming-Dutra et al. 2016)* should be commented out. Run entire file. Dataframes `bystander.df` and `bystanderbyclass.df` contain bystander proportions for all species across individual antibiotics and antibiotic classes, respectively.

In order the run the sensitivity analysis using the tiered diagnosis, comment out code chunk B2(viii) and comment in code chunk B2(ix). Rerun all code.

#### Checks:
The baseline analysis may take up to 30 minutes to run (majority of runtime due to code chunk B3).

##### 1. *Microbiome file*
Using the code and data included above, the age group-specific weighted prevalences for the first 5 rows of the `microbiome` dataframe should be as follows:

| Species.strain.Key    | wtprev_adults     | wtprev_kids0    | wtprev_kids1to5 |
| :-------------        | :----------:      | :-----------:   | :-----------:   |
| Escherichia_coli      | 0.66349810        | 0.94870000      | 1.00000000      |
| Haemophilus_influenzae| 0.68631179        | 1.00000000      | 0.95896469      |
| Klebsiella_pneumoniae | 0.07414449        | 0.39097000      | 0.15000000      |
| Moraxella_catarrhalis | 0.02281369        | 0.45485113      | 0.50790000      |
| Pseudomonas_aeruginosa| 0.01901141        | 0.01359456      | 0.01359456      |

##### 2. *NAMCS/NHAMCS summary file*
In the baseline analysis, using the code and data included above, the first 5 rows of the `NAMCS.summary` dataframe should be as follows:

| antibiotic  | dataset     | condition        | agegroup | wtVisits  | wtVisits.se | drugclass   |
| :---------- | :------     | :-------         | :------- | :------:  | :---------: | :---------- |
| AMOXICILLIN | namcs2010   | acuteSinusitis   | adults   | 2356812   | 565377.953  | PENICILLINS |
| AMOXICILLIN | namcs2010   | acuteSinusitis   | kids0    | 70355     | 70355.000   | PENICILLINS |
| AMOXICILLIN | namcs2010   | acuteSinusitis   | kids1to5 | 158253    | 67000.831   | PENICILLINS |
| AMOXICILLIN | namcs2010   | chronicSinusitis | adults   | 1769606   | 337174.213  | PENICILLINS |
| AMOXICILLIN | namcs2010   | chronicSinusitis | kids0    | 166430    | 106789.037  | PENICILLINS |

In the baseline analysis, using the code and data included above, the first 5 rows of the `NAMCS.summary.byclass` dataframe should be as follows:

| dataset     | condition        | agegroup | drugclass   | wtVisits  | wtVisits.se |
| :---------- | :-----------     | :------- | :-------    | :------:  | :---------: | 
| namcs2010   | acuteSinusitis   | adults   | PENICILLINS | 3448543   | 681417.82   |
| namcs2010   | acuteSinusitis   | kids0    | PENICILLINS | 70355     | 70355.00    |
| namcs2010   | acuteSinusitis   | kids1to5 | PENICILLINS | 158253    | 67000.83    |
| namcs2010   | chronicSinusitis | adults   | PENICILLINS | 3650426   | 562887.03   |
| namcs2010   | chronicSinusitis | kids0    | PENICILLINS | 166430    | 106789.04   |

##### 3. *Bystander output*
In the baseline analysis, using the code and data included above, the first 3 columns of the first 5 rows of the `bystander.df` dataframe should be as follows:

| species                   | drugclass                     | bystander_prop        | 
| :-------------            | :----------                   | :-----------:          | 
| Streptococcus_pneumoniae  | AMOXICILLIN                   | 0.8185041             | 
| Streptococcus_pneumoniae  | AMOXICILLIN-CLAVULANATE       | 0.8576846             | 
| Streptococcus_pneumoniae  | PENICILLIN                    | 0.9867820             | 
| Streptococcus_pneumoniae  | AZITHROMYCIN                  | 0.9160512             | 
| Streptococcus_pneumoniae  | CLINDAMYCIN                   | 0.9825318             | 





