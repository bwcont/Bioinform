## DEG on Laurenti et al microarray data
#1.14.2020
#E.C.

setwd("~/Documents/Bioinformatics/LaurentiMicroArray")

### Reading into R ####
HSC1_CB1 <- read.csv("GSM1039578.csv")
HSC1_CB2 <- read.csv("GSM1039579.csv")
HSC1_CB3 <- read.csv("GSM1039580.csv")
HSC1_CB4 <- read.csv("GSM1039581.csv")
HSC1_CB5 <- read.csv("GSM1039582.csv")
MLB_CB1 <- read.csv("GSM1039593.csv")
MLB_CB2 <- read.csv("GSM1039594.csv")
MLB_CB3 <- read.csv("GSM1039595.csv")
MLB_CB5 <- read.csv("GSM1039596.csv")
ETP_Thy1 <- read.csv("GSM1039616.csv")
ETP_Thy2 <- read.csv("GSM1039617.csv")
ETP_Thy3 <- read.csv("GSM1039618.csv")

### Removing first few rows that aren't needed ####
HSC1_CB1 <- HSC1_CB1[43:29374,]
HSC1_CB2 <- HSC1_CB2[43:29374,]
HSC1_CB3 <- HSC1_CB3[43:29374,]
HSC1_CB4 <- HSC1_CB4[43:29374,]
HSC1_CB5 <- HSC1_CB5[43:29374,]
MLB_CB1 <- MLB_CB1[43:29374,]
MLB_CB2 <- MLB_CB2[43:29374,]
MLB_CB3 <- MLB_CB3[43:29374,]
MLB_CB5 <- MLB_CB5[43:29374,]
ETP_Thy1 <- ETP_Thy1[43:29374,]
ETP_Thy2 <- ETP_Thy2[43:29374,]
ETP_Thy3 <- ETP_Thy2[43:29374,]

### Testing my organization methods with one set of cells first####

combined.HSC <- cbind(HSC1_CB1, HSC1_CB2, HSC1_CB3,HSC1_CB4,HSC1_CB5)

### Ensuring all rows are aligned, all genes are ordered the same way across samples####
test.same <- as.character(combined.HSC[,1]) == as.character(combined.HSC[,4]) & 
             as.character(combined.HSC[,1]) == as.character(combined.HSC[,7]) & 
             as.character(combined.HSC[,1]) == as.character(combined.HSC[,10]) & 
             as.character(combined.HSC[,1]) == as.character(combined.HSC[,13])
test.same <- as.data.frame(test.same)
test.same[test.same[,1] == FALSE,]
#All the values are ordered in the same way 

### Renaming####

combined.HSC <- combined.HSC[,c(1,2,3,5,6,8,9,11,12,14,15)] #Remove redundant names
colnames(combined.HSC) <- c("ID_REF", 
                            'HSC1_CB1_VALUE','HSC1_CB1_pval',
                            'HSC1_CB2_VALUE','HSC1_CB2_pval',
                            'HSC1_CB3_VALUE','HSC1_CB3_pval',
                            'HSC1_CB4_VALUE','HSC1_CB4_pval',
                            'HSC1_CB5_VALUE','HSC1_CB5_pval')
combined.HSC <-as.data.frame(combined.HSC)
combined.HSC <- combined.HSC[-1,] #remove top redundant row


#####


