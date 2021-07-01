# Functions for deriving Phenotype Data from MGBB/RPDR (Based on COPD Project)
# Su Chu, Ph.D.
# su.chu@channing.harvvard.edu 

require(data.table)
require(stringr)
require(Hmisc)

#### 0. UTILITY FUNCTIONS ####
# data-wrangling
df.char2num=function(df){
        var.text=which(apply(df, 2, function(x){flag=grep('[:alpha:]',x); ifelse(length(flag)>1, 1,0)})>0)
        df[,-var.text]=apply(df[,-var.text], 2, as.numeric)
        return(df)
}

nunq <- function(x){ x %>% unique %>% length}

tabPercent <- function(x){return(rbind(table(x, exclude=NULL)/length(x)*100, table(x, exclude=F)))}

summary.sd <- function(x){return(c(summary(x), SD=sd(x)))}

expit <- function(x){
        y <- exp(x) / (1+exp(x))
        return(y)
}


#Formatting dates into R Date objects
fmtDate <- function(dt, colName='Date'){
        dt[ , Date_Fmt := strsplit(dt[, get(colName)], split=' ') %>% 
                    lapply(., '[[', 1) %>% 
                    unlist %>% 
                    as.Date(., '%m/%d/%Y')]
}

`%notin%` <- Negate(`%in%`)

#### 1. READING DATA ####
dataLoc=function(file, loc='mbp'){
        if(loc=='mbp') {return(paste("/Users/shchu/Dropbox (Partners HealthCare)/mgbb/covid19/COVID Data/RPDR_Files/Plasma_Files/07-15-2020/Base_RPDR_Files/SW029_20200715_153908_", file, sep=''))}
        #if(loc=='cdnm') {return(paste('C:/Users/rechu/Local-HMS-BWH/RPDR/raw/shc40_111317112943523114_', file, sep=''))}
}

# function to break lines of text by a specified delimiter, dlm
splitLines <- function(txtLine, dlm='\\|'){
        out <- strsplit(txtLine, split=dlm) %>% unlist
}

#function to read in data and extract lines with specific string patterns
processFile = function(filepath, pattern) {
        con = file(filepath, "r")
        l=list()
        i=1
        while ( TRUE ) {
                line = readLines(con, n = 1)
                if ( length(grep(paste(pattern), line))>0 ) { l[[i]] = line; i=i+1}
                if ( length(line) == 0 ) {
                        print('File processed.')
                        break
                }
        }
        
        close(con)
        return(l)
}

#function to use when fread doesn't work (slow)
get.dat=function(
        data=NA,         #input data file
        header=T,
        eol='\n',
        sep='\t',
        skl=0            #number of lines to skip
){
        f=file(data)
        row_list=sapply(scan(data, what=character(), sep=eol, skip=skl), strsplit, split=sep) 
        close(f)
        df=data.frame(do.call(rbind,row_list[2:length(row_list)]), stringsAsFactors=F)
        row.names(df) <- NULL; names(df) <- row_list[[1]]
        return(df)
}


#fn to extract text data in RPDR files with text reports (was initially built specifically for Pul.txt)
getTextReport = function(filepath, patternStart, patternEnd, idSplit, igc=T) {
        con = file(filepath, "r")
        l=list(); L=list()
        i=1; p=1
        recSwitch=FALSE
        nullSwitch=FALSE
        while ( TRUE ) {
                
                line = readLines(con, n = 1)
                
                if ( length(line) == 0 ) {
                        if(recSwitch){next}
                        if(nullSwitch){break}
                }
                
                #save patient information
                if ( length(grep(paste(idSplit), line, ignore.case=igc))>0 ) { l[[i]] = line; i=i+1 }
                
                #start with recording switch off / turn off at end of report
                if ( length(grep(paste(patternEnd), line, ignore.case=igc)) > 0 ) { 
                        recSwitch=FALSE
                        nullSwitch=TRUE
                        L[[p]]=unlist(l)
                        p=p+1
                        l=list()
                        # print(p)
                }
                
                #turn on recording switch if patternStart observed
                if ( length(grep(paste(patternStart), line, ignore.case=igc))>0 ) { recSwitch=TRUE; nullSwitch=FALSE }
                if(recSwitch){
                        l[[i]]=line
                        i=i+1
                }
                
        }
        
        close(con)
        return(L)
}

#### 2. SUMMARIZING DATA FUNCTIONS ####
# Fn to extract sets with <VAR> occurring > X TIMES/DATES (Returns subsetted data.table)
minOccur = function(sampleIDvec, dt, eventVarName, minCount){
        check.flag=rep(NA, length(sampleIDvec))
        for(i in 1:length(sampleIDvec)){
                check.flag[i]=ifelse(length(unique(dt[EMPI==sampleIDvec[i], paste(eventVarName)]))>minCount, 1, 0)
        }
        return(dt[EMPI %in% sampleIDvec[check.flag==1]])
}

library(RcppRoll) #for roll_sum function
# function to get the number of times an event occurs (e.g. in dx file, EMPI ID is equivalent to 
# COPD counts, as we restricted to only COPD ICDs in Diagnosis file)
# dt= data.table, varName='EMPI'
# eventCount = the number of consecutive dates to count the days between (e.g. if we want to assess time that elapsed between
#  3 events, we would denote 3 here)
# timeWindow = denote time which thresholds the following question: 
#   for each consecutive # of eventCounts, is the time that elapsed less than or equal to this time window?

countPerTimeWindow <- function(dt, varName, eventCount, timeWindow, newVarName){
        sid <- unique(dt[ , ..varName]) %>% unlist %>% unname
        y <- lapply(sid, 
                    function(x){ 
                            #cat(paste(x, '\n'))
                            time=c(0, diff(dt[get(varName)==x]$Date_Fmt))
                            time_span=roll_sum(time, eventCount)
                            flag=any(time_span<=timeWindow)
                            return(flag)
                    }) %>% unlist
        sid <- data.table(sid)
        names(sid) <- varName
        out <- sid[ , paste(newVarName, 'GTE', eventCount, '_', timeWindow, sep='') := y]
        return(out)
}



#### 3. MED.TXT FUNCTIONS ####
#requires two file inputs: 
# 1) COPD/Asthma Medication File with columns Medication, BrandNames, Disease, Class, TargetMechanism, RxNote
# 2) rxn_ingredient_rpdr_map: File with columns rxn_ingredient, rxn_ingredient_name, c_basecode, c_name

# get rpdr code mappings for the medications associated with disease of interest
ss.rxn.map <- function(ref_meds, map.med.rxn, dx){
        paste(ref_meds[Disease==dx]$Medication, collapse='|') %>%
                c(grep(., map.med.rxn$rxn_ingredient_name), grep(., map.med.rxn$c_name)) %>%
                unique %>%
                as.integer %>%
                .[-which(is.na(.))] %>%
                map.med.rxn[.,] -> map.med.rxn.m
        
        paste(ref_meds[Disease==dx & !is.na(BrandNames)]$BrandNames, collapse='|') %>%
                gsub(', ', '|', .) %>%
                c(grep(., map.med.rxn$rxn_ingredient_name), grep(., map.med.rxn$c_name)) %>%
                unique %>%
                as.integer %>%
                .[-which(is.na(.))] %>%
                map.med.rxn[.,] -> map.med.rxn.bn
        
        map.med.rxn.out <- 
                merge(map.med.rxn.m, map.med.rxn.bn, 
                      by=c('rxn_ingredient', 'rxn_ingredient_name', 'c_basecode', 'c_name', 'Code_Type', 'Code'), all=T)
        
        map.med.rxn.out[ , paste('i', dx, 'med', sep='.') := 1]
        return(map.med.rxn.out)
}


#assign single class of medications to mappings 
assn.rxn.class <- 
        function(ref_meds, # data table with cols: Medication, BrandNames, Disease, Class, Class2, TargetMechanism, RxNote
                 map.med.rxn, # data table with cols: rxn_ingredient, rxn_ingredient_name, c_basecode, c_name
                 rxn.class, # any specific rxn class of interest
                 rxn.class.field='Class' #Class field name: e.g. 'Class', 'Class1', 'Class2'
        ){
                ref_meds[get(rxn.class.field)==rxn.class & !is.na(BrandNames)]$BrandNames %>%
                        paste(., collapse='|') %>%
                        gsub(', ', '|', .) %>% 
                        c(ref_meds[get(rxn.class.field)==rxn.class]$Medication, .) %>%
                        paste(., collapse='|') -> search.string
                
                if(grepl('|$', search.string))
                { search.string <- substr(search.string, 1, nchar(search.string)-1) }
                
                if(nchar(search.string)!=0){ 
                        grepl(search.string, map.med.rxn$rxn_ingredient_name) %>%
                                cbind(., grepl(search.string, map.med.rxn$c_name)) %>%
                                apply(., 1, function(x) ifelse(any(x), 1, 0)) ->
                                i.class
                        out <- map.med.rxn[ , paste('i', rxn.class, sep='.') := i.class]
                } else { out <- map.med.rxn[ , paste('i', rxn.class, sep='.') := 0] }
                
                return(out)
        }


#assign all medication classes to map file
assn.rxn.class.global <- function(rref_meds, mmap.med.rxn, rrxn.class.field='Class'){
        rxn.class.list <- unique(rref_meds[ , get(rrxn.class.field)])
        for(c in rxn.class.list){
                assn.rxn.class(rref_meds, mmap.med.rxn, c, rxn.class.field=rrxn.class.field)
        }
}



#### PFT & PRC DATA FUNCTIONS FOR SPIROMETRY VARIABLES ####



pulCleanUp <- function(dt.pul.id){
        #delete misread lines (e.g. NWH files that don't have good report text delimiters)
        #due to lack of appropriate delimiters in some of the NWH files, the report
        #text is read in as subject lines, in addition to the subject ID lines themselves: remove these.
        if(any(grep('\\s+', dt.pul.id$MRN_Type))) {dt.pul.id <- dt.pul.id[-grep('\\s+', MRN_Type)]}
        
        #format the character date into date format understood by R
        dt.pul.id[ , Date_Fmt := strsplit(dt.pul.id$Report_Date_Time, split=' ') %>% 
                           lapply(., '[[', 1) %>% 
                           unlist %>% 
                           as.Date(., '%m/%d/%Y')
                   ]
        dt.pul.id[ , EMPI := as.integer(EMPI)]
        dt.pul.id[ , MRN := as.integer(EMPI)]
        
        return(dt.pul.id)
}

prcCleanUp <- function(dt.prc){
        if(any( sapply(dt.prc, class) != 'character')){
                ncvarIndex <- which(sapply(dt.prc, class) != 'character')
                paste('Note: The following variable is not of character class: ', 
                      names(dt.prc)[ncvarIndex], '\n', sep='') %>%
                        cat()
        }
        #format date and reorder by EMPI/MRN
        dt.prc[ , Date_Fmt := as.Date(Date, '%m/%d/%Y')]
        dt.prc <- dt.prc[order(EMPI, MRN, Date_Fmt)]
        
        #remove duplicated entries
        dt.prc[ , c(Cs(EMPI, MRN_Type, MRN, Date_Fmt, Procedure_Name, Code_Type, Code)), with=F] %>%
                duplicated %>% 
                which -> t
        dt.prc <- dt.prc[-t]
        
        return(dt.prc)
}



pftSubjIndex <- function(dt.pft.pul, recordLine=1, dlm='\\|', clmns=c(1:9), firstRecordContainsHeader=T){
        #for pft report read in, the header column doesn't have a 'report end' delimiter
        #and is lumped in with the first subject's data. Extract column name info (first line)
        #from the first subject record & save, then remove this line from the first
        #subject record text.
        if(firstRecordContainsHeader){ 
                colNames <- splitLines(dt.pft.pul[[1]][1], dlm=dlm)
                dt.pft.pul[[1]] <- dt.pft.pul[[1]][-1]
        }
        
        lapply(dt.pft.pul, function(x){ splitLines( x[recordLine], dlm=dlm )[clmns] }) %>% 
                do.call('rbind', .) %>% 
                data.table -> out
        
        if(firstRecordContainsHeader) {names(out) <- colNames[clmns]}
        
        return(out)
}


# get the corresponding pul.txt file information per report
getPULIDs <- function(pulFilePath, delim='[|]'){
        pul.id<- processFile(pulFilePath, delim) %>%
                lapply(., function(x) unlist(strsplit(x, split=delim))) %>% 
                lapply(., function(x) x[1:9]) %>% 
                do.call(rbind, .) %>%
                data.table() 
        
        names(pul.id) <- unlist(pul.id[1,])
        pul.id <- pul.id[-1]
        pul.id <- pulCleanUp(pul.id)
        
        return(pul.id)
}

# Function to determine if pdf/rpt/prc code of spirometry is available
pftInfo <- function(dt.pft.id, 
                    dt.prc,
                    pdfDir='/Users/shchu/Documents/Local-HMS-BWH/RPDR/pheno-copd/raw/Pul.zip/',
                    cpt.pft.codes=c("94010", "94060", "94375", "94150", "94200", 
                                    "94200", "94240", "94720", "94070", "94360", 
                                    "94260", "93720", "94621", "94620")
){
        
        pft.fileNames <- list.files(pdfDir)
        ## PFT DATA  
        pft.id.f1 <- dt.pft.id[grep('FILE1', Report_Number)] %>% data.table
        pft.id.f2 <- dt.pft.id[grep('FILE2', Report_Number)] %>% data.table
        pft.id.o <- dt.pft.id[-grep('FILE2|FILE1', Report_Number)] %>% data.table
        
        
        #get columns with file names for Report Text and Report PDFs, where applicable
        pft.id.f1[ , Report_Number_Text := Report_Number]
        pft.id.f2[ , Report_Number_PDF := Report_Number]
        pft.id.o[ , Report_Number_PDF := ifelse(Report_Number %in% pft.fileNames, Report_Number, NA)]
        pft.id.o[ , Report_Number_Text := ifelse(!(Report_Number %in% pft.fileNames), Report_Number, NA)]
        
        #merge 
        pft.id.drv <- merge(pft.id.f1[ , -c('MID'), with=F], pft.id.f2[ , -c('Report_Number', 'MID')],
                            by=c('EMPI', 'MRN_Type', 'MRN', 'Report_Date_Time', 'Date_Fmt','Report_Description', 
                                 'Report_Status', 'Report_Type'), all=T)
        pft.id.drv[ , rpt_flag := 0]
        pft.id.o[ , rpt_flag := 1]
        pft.id.drv <- rbind(pft.id.drv, pft.id.o[ , -c('MID'), with=F]) %>% data.table 
        
        
        ## now incorporate code data
        prc.id.ss <- dt.prc[Code %in% cpt.pft.codes]
        
        prc.pft.ids <- merge(pft.id.drv, prc.id.ss[ , -c('MRN'), with=F], by=Cs(EMPI, MRN_Type, Date_Fmt), all=T)
        prc.pft.ids[ , rpt_prc := ifelse((!is.na(Report_Number_Text) | !is.na(Report_Number_PDF)) & !is.na(Code), 'both', 
                                         ifelse(!is.na(Code), 'no_rpt', 'no_prc'))]
        
        return(prc.pft.ids)
}


getPulReportAttributes <- function(
        x # character string with subject data of outputted data.table that is extracted from getTextReport(<file_path_Pul.txt>)
){
        #read in character string
        strsplit(x[1], split='[|]') %>% 
                unlist %>% 
                {.} -> info_string
        
        #extract report c;ass
        strsplit(info_string[5], split='[[:punct:]]') %>%
                lapply(., '[[', 1) %>%
                unlist %>%
                str_split(., '[[:number:]]') %>%
                lapply(., '[[', 1) %>% 
                unlist%>%
                substr(., start=1, stop=12) %>%
                unlist -> report_class
        
        report_number <- info_string[5]
        
        list(report_class=report_class, info_string=info_string, report_number=report_number)
}


#function to contain strings for spirometry strings
key_pftValueStrings <- function(
        report_class, #report class string, e.g. 'BWHBICSPUL'
        pftMeasure='fev1' #or 'fvc', or 'fvr' (fev1/fvc ratio)
){
        
        #MGPCISPUL 631
        #PHSEPIC 651
        #MGHCOMPASPUL 211
        #BWHBICSPUL 7
        #PHSEPIC 1-6
        #NSMCBREEZE 452, 591
        #NSMCLCR454
        
        #fev1/fvc ratio
        if(pftMeasure=='fvr'){
                pftSearchString <- switch(report_class, 
                                          BWHBICSPUL='^FEV! / FVC|^FEV1/FVC\\s+\\(\\%\\)', #X
                                          MGHCOMPASPPU='^FEV1/FVC %|FEV1/FVC\\s+\\(\\%\\)', # unknown as not represented in this dataset
                                          MGHCOMPASPUL='^FEV! / FVC %|^FEV1/FVC\\s+\\(\\%\\)', #X
                                          MGHPCISPUL="^FEV1/VC|^FEV1/FVC\\s+\\(\\%\\)", #X
                                          NSMCBREEZESU='^FEV1/FVC \\(\\%\\)',#XY
                                          PHSEPIC='^FEV! / FVC \\%', #X
                                          NWHMTOEDEPT='NWHMTOEDEPT',
                                          FHMT='FHMT',
                                          NSMCLCR='NSMCLCR',
                                          DFCIBICSPUL='^FEV1/FVC\\s+\\(%\\)')#XY
        }
        
        
        if(pftMeasure=='fev1'){
                pftSearchString <- switch(report_class, 
                                          BWHBICSPUL='^FEV! L|^FEV1\\s+\\(Lts\\)', #X
                                          MGHCOMPASPPU='^FEV1\\s+\\(L\\)|^FEV1 L|^FEV!', #unknown as not represented in dataset
                                          MGHCOMPASPUL='^FEV!\\s+L|^FEV1\\s+\\(L\\)', # X
                                          MGHPCISPUL="^FEV1\\s+\\(L\\)", #X
                                          NSMCBREEZESU='^FEV1 \\(L\\)', #X
                                          PHSEPIC='^FEV! L', #X
                                          NWHMTOEDEPT='NWHMTOEDEPT',
                                          FHMT='FHMT',
                                          NSMCLCR='NSMCLCR',
                                          DFCIBICSPUL='^FEV1')#XY
        }
        
        if(pftMeasure=='fvc'){
                pftSearchString <- switch(report_class, 
                                          BWHBICSPUL='^FVC\\s+\\(Lts\\)', #X
                                          MGHCOMPASPPU='^FVC\\s+\\(L\\)',  #unknown as not represented in dataset
                                          MGHCOMPASPUL='^FVC L|^FVC\\s+\\(L\\)', #X
                                          MGHPCISPUL="^FVC\\s+\\(L\\)", #X
                                          NSMCBREEZESU='^FVC \\(L\\)',#XY
                                          PHSEPIC='^FVC L', #X
                                          NWHMTOEDEPT='NWHMTOEDEPT',
                                          FHMT='FHMT',
                                          NSMCLCR='NSMCLCR',
                                          DFCIBICSPUL='^FVC') #XY
        }
        
        return(pftSearchString)
}

#### FILE 2 Layouts ####
# BWHBICSPUL Format in Pul.txt [see pul.rpt[[45]] compared to the pdf file]
#  MX | Unit | Predicted Mean | PreBD Actual | PostBD Actual | % Change | PreBD % Pred | PostBD %Pred | Predicted Range 95%
# MGHPCISPUL
#  MX | (Unit) | PreBD Actual | PreBD Predicted | Pre BD %Pred | PreBD CI Range/Eval | PostBD Meas | PostBD %Pred | PostBD %Change
# MGHCOMPASPUL
#  Mx | Unit | PreBD Pred | PreBD Actual | PostBD Actual | PostBD %Change | PreBD %Pred | PostBD %Pred | PreBD CI Range Low | Random | Random | PreBD CI Range High
# PHSEPIC 
#  Mx | Unit | Predicted Mean | PreBD Actual | PostBD Actual | % Change | PreBD %Pred | PostBD %Pred | Predicted Range 95%
# NSMCBREEZE
#  Mx | (Unit) | PreBD Actual | PreBD Predicted | PreBD %Pred | PostBD Actual | PostBD %Pred | PostBD %Change
#####


#function to find correct index in spirometry strings for the corresponding pre/post BD measures
key_pftValueIndices <- function(pft_vals,
                                report_class,
                                report_number){
        if(grepl('FILE2$', report_number)){
                out <- switch(report_class, 
                              BWHBICSPUL=c(4,3,7,5,NA,8),
                              MGHCOMPASPUL=c(3,2,6,4,NA,7),
                              MGHCOMPASPPU=c(3,2,6,4,NA,7),
                              #MGHPCISPUL=c(3,4,5,7,NA,9),
                              NSMCBREEZESU=c(3,4,5,6,NA,8),
                              PHSEPIC=c(3,2,6,4,NA,7))
                if(is.null(out)){
                        cat(paste('Warning: No FILE2 string in Report_Number identified for this Report_Class. Add new index switch to key_pftValueIndices() for this PFT: ', report_number, '.\n', sep=''))}
        }
        if(!grepl('FILE2$', report_number) | grepl('_OBSV$', report_number)){
                out <- switch(report_class, 
                              BWHBICSPUL=c(4, 2, 5, 6, NA, 7),
                              MGHCOMPASPUL=c(2, 3, 4, 7, NA, 8), #XY
                              MGHCOMPASPPU=c(2, 3, 4, 7, NA, 8), #XY
                              MGHPCISPUL=c(2, 3, 4, 7, NA, 8), 
                              NSMCBREEZESU=c(2, 3, 4, 5 , NA, 7),
                              #PHSEPIC=c(3, 2, 6, 4, NA, 6),
                              NSMCLCR=rep(NA, 6),
                              DFCIBICSPUL=c(2, 4, 3, 5, NA, 6))
                if(is.null(out)){
                        cat(paste('Warning: No _OBSV or FILE1 string in Report_Number identified for this Report_Class. Add new index switch to key_pftValueIndices() for this PFT: ', report_number, '.\n', sep=''))}
        }
        if(grepl('PDF|TXT', report_number)){
                out <- switch(report_class, 
                              MGHCOMPASPUL=c(3, 2, 6, 4, NA, 7), #XY
                              NSMCBREEZESU=c(2, 3, 4, 5 , NA, 6)) #XY
                if(is.null(out)){
                        cat(paste('Warning: No PDF|TXT string in Report_Number identified for this Report_Class. Add new index switch to key_pftValueIndices() for this PFT: ', report_number, '.\n', sep=''))}
        }
        return(out)
        
}

# extract the values of interest from all the different pft report types
extractValuesFromReport<- function(pft_vals, pftstring_index){
        if(!all(is.na(pft_vals))){
                i <- pftstring_index
                #if(any(is.na(i)))  pft_vals[] <- NA 
                out <- c(preBD = ifelse(!is.na(pft_vals[i[1]]), pft_vals[i[1]], NA) %>% as.numeric, 
                         preBD_Pred = ifelse(!is.na(pft_vals[i[2]]), pft_vals[i[2]], NA) %>% as.numeric,
                         preBD_PctPred = ifelse(!is.na(pft_vals[i[3]]), pft_vals[i[3]], NA) %>% as.numeric,
                         postBD = ifelse(!is.na(pft_vals[i[4]]), pft_vals[i[4]], NA) %>% as.numeric,
                         postBD_Pred = ifelse(!is.na(pft_vals[i[5]]), pft_vals[i[5]], NA) %>% as.numeric,
                         postBD_PctPred = ifelse(!is.na(pft_vals[i[6]]), pft_vals[i[6]], NA) %>% as.numeric)
                return(out)
        } else { cat('Error: No report values found.\n')}
}




getPFTvals <- function(x, mmeasure='fev1', vrbs=F){ 
        info_string <- getPulReportAttributes(x)$info_string
        report_class <- getPulReportAttributes(x)$report_class
        report_number <- getPulReportAttributes(x)$report_number
        #if(vrbs){print(x[1])}
        
        #search spirometry
        pft_string <- key_pftValueStrings(report_class, pftMeasure=mmeasure)
        #if(vrbs){print(pft_string)}
        if(is.null(pft_string)){cat(paste('Warning: The report class', report_class, 'has a NULL search string. See text of', report_number, 'to identify correct character string for', mmeasure, '\n', sep=' '))}
        
        #get the corresponding line in text file report for our mmeasure of interest
        index <- ifelse(length(pft_string)==0, 
                        NA, 
                        grep(pft_string, x)[1])
        
        #handle spirometry reports from report_class that don't have a format
        if(is.na(index)) { 
                pft_out=c(preBD1 = NA, preBD2=NA, preBD3=NA, postBD1 = NA, postBD2=NA, postBD3=NA)
        } else if(!is.na(index)){
                pft_vals <- str_split(x[index], pft_string) %>% unlist %>% str_split(., '\\s+') %>% .[[2]]
                
                pftstring_index <- key_pftValueIndices(pft_vals, report_class, report_number)
                pft_out <- extractValuesFromReport(pft_vals, pftstring_index)
        }
        
        #if(vrbs){print(pft_out)}
        
        #if(all(is.na(pft_out)) & !grep('FILE2$', report_number)){ print(report_number); print(pftstring_index)}
        #format the date of the report
        fmtdate <- strsplit(info_string[6], split=' ') %>% 
                lapply(., '[[', 1) %>% 
                unlist %>% 
                as.Date(., '%m/%d/%Y')
        
        #format output lines
        cbind(EMPI=info_string[1], EPIC_PMRN=info_string[2], MRN_Type=info_string[3], MRN=info_string[4], 
              Report_Number=info_string[5], Report_Date_Time=info_string[6], Report_Description=info_string[7]) %>% 
                data.table %>%
                cbind(.,  
                      Report_Class=report_class, 
                      Date_Fmt=fmtdate, 
                      preBD1=pft_out[1], 
                      preBD2=pft_out[2],
                      preBD3=pft_out[3],
                      postBD1=pft_out[4],
                      postBD2=pft_out[5],
                      postBD3=pft_out[6]) %>%
                {.} -> out
        
        
        names(out)[(ncol(out)-5):ncol(out)] <- c(paste('preBD', mmeasure, sep='_'),
                                                 paste('preBD_Pred', mmeasure, sep='_'),
                                                 paste('preBD_PctPred', mmeasure, sep='_'),
                                                 paste('postBD', mmeasure, sep='_'),
                                                 paste('postBD_Pred', mmeasure, sep='_'),
                                                 paste('postBD_PctPred', mmeasure, sep='_'))
        #return pft info
        return(out) 
}


dtPFTData <- function(
        pftDat, #pftDat = RPDR pft.txt file format, as a data.table
        measure='fev1',
        vrbs=F
){
        lapply(pftDat, function(x) getPFTvals(x, measure, vrbs)) %>% 
                do.call('rbind', .) -> pft.out
        
        pft.out[ , EMPI := as.integer(EMPI)]
        pft.out[ , MRN := as.integer(MRN)]
        
        return(pft.out)
}



#### PLOT FXNS ####


panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
        usr <- par("usr"); on.exit(par(usr)) 
        par(usr = c(0, 1, 0, 1)) 
        r <- cor(x, y) 
        txt <- format(c(r, 0.123456789), digits=digits)[1] 
        txt <- paste(prefix, txt, sep="") 
        if(missing(cex.cor)) cex <- 0.5/strwidth(txt) 
        
        test <- cor.test(x,y) 
        # borrowed from printCoefmat
        Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " ")) 
        
        #text(0.5, 0.5, txt, cex = cex * abs(r)) 
        text(0.5,0.5,txt,cex=cex)
        text(.8, .8, Signif, cex=cex, col=2) 
}

panel.smooth<- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                         cex = 1, col.smooth = "red", span = 2/3, iter = 3) 
{
        points(x, y, pch = pch, col = col, bg = bg, cex = cex)
        ok <- is.finite(x) & is.finite(y)
        if (any(ok)) 
                lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
                      col = col.smooth)
}

panel.hist<-function(x)
{
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col="cyan")
}



#### performance stats for elnet ####
alg.prfrm <- function(betas, betaCol, dt, th){
        b0 <- betas[1, c(betaCol), with=F] %>% as.numeric
        index.b <- which(!(betas[,betaCol,with=F]==0 | is.na(betas[,betaCol,with=F])))
        logodds <- as.matrix(dt[ , betas$Feature[index.b][-1], with=F]) %*% as.matrix(betas[index.b, betaCol, with=F][-1]) + b0
        phat <- exp(logodds) / (1+exp(logodds))
        dt[ , iAlg := ifelse(phat<th, 0, 1)]
        
        out.counts <- table(truth=dt$iCOPD, test=dt$iAlg)
        out.prfstats <- predictperform(out.counts)
        
        
        list(counts=out.counts, prfstats=out.prfstats)
        
}



predictperform <- function(tbl22, rowTruth=T, case=2){
        if(rowTruth){
                sens=tbl22[2,2]/(tbl22[2,1] + tbl22[2,2])
                ppv=tbl22[2,2]/(tbl22[1,2] + tbl22[2,2])
                npv=tbl22[1,1]/(tbl22[1,1] + tbl22[2,1])
                spec=tbl22[1,1]/(tbl22[1,1] + tbl22[1,2])
        }
        list( sens=sens, spec=spec, ppv=ppv,  npv=npv)
        
}
