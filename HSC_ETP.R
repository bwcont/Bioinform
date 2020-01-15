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
ETP_Thy3 <- ETP_Thy3[43:29374,]

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
remove(combined.HSC, test.same)


### Go back to the data sets and removing low quality reads####
HSC1_CB1 <- HSC1_CB1[-1,]
colnames(HSC1_CB1) <- c("GeneID","VALUE","QCpval")

HSC1_CB2 <- HSC1_CB2[-1,]
HSC1_CB3 <- HSC1_CB3[-1,]
HSC1_CB4 <- HSC1_CB4[-1,]
HSC1_CB5 <- HSC1_CB5[-1,]
MLB_CB1 <- MLB_CB1[-1,]
MLB_CB2 <- MLB_CB2[-1,]
MLB_CB3 <- MLB_CB3[-1,]
MLB_CB5 <- MLB_CB5[-1,]
ETP_Thy1 <- ETP_Thy1[-1,]
ETP_Thy2 <- ETP_Thy2[-1,]
ETP_Thy3 <- ETP_Thy3[-1,]

colnames(HSC1_CB2) <- c("GeneID","VALUE","QCpval")
colnames(HSC1_CB3) <- c("GeneID","VALUE","QCpval")
colnames(HSC1_CB4) <- c("GeneID","VALUE","QCpval")
colnames(HSC1_CB5) <- c("GeneID","VALUE","QCpval")
colnames(MLB_CB1) <- c("GeneID","VALUE","QCpval")
colnames(MLB_CB2) <- c("GeneID","VALUE","QCpval")
colnames(MLB_CB3) <- c("GeneID","VALUE","QCpval")
colnames(MLB_CB5) <- c("GeneID","VALUE","QCpval")
colnames(ETP_Thy1) <- c("GeneID","VALUE","QCpval")
colnames(ETP_Thy2) <- c("GeneID","VALUE","QCpval")
colnames(ETP_Thy3) <- c("GeneID","VALUE","QCpval")

#ETPs
ETP_Thy1$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(ETP_Thy1$QCpval[i]))> 0.05,
    ETP_Thy1[i,4] <- NA, 
    ETP_Thy1[i,4] <- as.numeric(as.character(ETP_Thy1[i,2]))
    )
}

ETP_Thy2$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(ETP_Thy2$QCpval[i]))> 0.05,
    ETP_Thy2[i,4] <- NA, 
    ETP_Thy2[i,4] <- as.numeric(as.character(ETP_Thy2[i,2]))
  )
}

ETP_Thy3$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(ETP_Thy3$QCpval[i]))> 0.05,
    ETP_Thy3[i,4] <- NA, 
    ETP_Thy3[i,4] <- as.numeric(as.character(ETP_Thy3[i,2]))
  )
}
#MLB
MLB_CB1$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(MLB_CB1$QCpval[i]))> 0.05,
    MLB_CB1[i,4] <- NA, 
    MLB_CB1[i,4] <- as.numeric(as.character(MLB_CB1[i,2]))
  )
}

MLB_CB2$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(MLB_CB2$QCpval[i]))> 0.05,
    MLB_CB2[i,4] <- NA, 
    MLB_CB2[i,4] <- as.numeric(as.character(MLB_CB2[i,2]))
  )
}

MLB_CB3$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(MLB_CB3$QCpval[i]))> 0.05,
    MLB_CB3[i,4] <- NA, 
    MLB_CB3[i,4] <- as.numeric(as.character(MLB_CB3[i,2]))
  )
}

MLB_CB5$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(MLB_CB5$QCpval[i]))> 0.05,
    MLB_CB5[i,4] <- NA, 
    MLB_CB5[i,4] <- as.numeric(as.character(MLB_CB5[i,2]))
  )
}








