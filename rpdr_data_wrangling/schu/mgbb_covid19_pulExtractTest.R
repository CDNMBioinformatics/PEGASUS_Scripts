# prelim analysis

library(data.table)
#library(doMC)
#library(ff)
library(magrittr)
library(readxl)
library(openxlsx)
source('./code/fxns.rpdr.drv.phen_1022.R')


filePathToPulTxt <- 'file/path/to/PartnersID_RandomNumbers_Pul.txt'

pul.rpt<-getTextReport(filePathToPulTxt, 
                       "Spiro|SPIROMETRY|Interpretation", "report_end", 
                       "[|]")

pul.id <- getPULIDs(filePathToPulTxt)

#name the elements of pul.rpt text list
names(pul.rpt) <- pul.id$EMPI

#drop header line that we get in first report of the text file
pul.rpt[[1]] <- pul.rpt[[1]][-1]

# number of unique subjects with any pft data in pul.txt file
unique(pul.id$EMPI) %>% length 



#report type/prefix strings:
strsplit(pul.id$Report_Number, '[[:punct:]]|[:0-9:]') %>% lapply(., '[[', 1) %>% unlist  %>% table
strsplit(pul.id$Report_Number, '[[:punct:]]|[:0-9:]') %>% lapply(., '[[', 1) %>% unlist %>% unique

# create ReportClass variable for sanity checking/debugging later
pul.id[ , ReportClass := strsplit(Report_Number, '[[:punct:]]|[:0-9:]') %>% lapply(., '[[', 1) %>% unlist]


# N.B.: some of the lung function text files don't have spirometry, or don't have pdfs to confirm accuracy of extraction
# DFCIBICSPUL, MGHPCISPUL, FHMT = no pdfs, but do have values {File 1}
# EPIC files from FH have no quantitative spirometry values, but do have qualitative ("normal")
# 
pul.id$MRN_Type %>% table
which(pul.id$MRN_Type=='FH')
pul.id[which(pul.id$MRN_Type=='FH'),]
pul.rpt[which(pul.id$MRN_Type=='FH')]

# Each subject has UP TO two records of each PFT they had (if there are two, 
# usually denoted as "FILE1" and "FILE2" with same prefix)
# Note that these functions pull out values from BOTH of the FILE1 and FILE2 reports
# Sometimes there will be more values from File 1 extractions, sometimes more from File 2 Extractions
# Still need to write a function to get the most complete set of data. 
# Note that postBD_Pred_XX will always be NA. Because I realized this should always just be the same as preBD_Pred_XX. 
# Function can be cleaned up to remove that column later.


fev1 <- dtPFTData(pul.rpt, measure='fev1', F)
fvc <- dtPFTData(pul.rpt, measure='fvc', F)
fvr <- dtPFTData(pul.rpt[1:10], measure='fvr', F)

