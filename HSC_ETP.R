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
#HSC
HSC1_CB1$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(HSC1_CB1$QCpval[i]))> 0.05,
    HSC1_CB1[i,4] <- NA, 
    HSC1_CB1[i,4] <- as.numeric(as.character(HSC1_CB1[i,2]))
  )
}

HSC1_CB2$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(HSC1_CB2$QCpval[i]))> 0.05,
    HSC1_CB2[i,4] <- NA, 
    HSC1_CB2[i,4] <- as.numeric(as.character(HSC1_CB2[i,2]))
  )
}

HSC1_CB3$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(HSC1_CB3$QCpval[i]))> 0.05,
    HSC1_CB3[i,4] <- NA, 
    HSC1_CB3[i,4] <- as.numeric(as.character(HSC1_CB3[i,2]))
  )
}

HSC1_CB4$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(HSC1_CB4$QCpval[i]))> 0.05,
    HSC1_CB4[i,4] <- NA, 
    HSC1_CB4[i,4] <- as.numeric(as.character(HSC1_CB4[i,2]))
  )
}

HSC1_CB5$VALUES <- NA
for (i in 1:29331) {
  ifelse(
    as.numeric(as.character(HSC1_CB5$QCpval[i]))> 0.05,
    HSC1_CB5[i,4] <- NA, 
    HSC1_CB5[i,4] <- as.numeric(as.character(HSC1_CB5[i,2]))
  )
}

#Remove pvals and redundant values
ETP_Thy1 <- ETP_Thy1[,c(-2,-3)]
ETP_Thy2 <- ETP_Thy2[,c(-2,-3)]
ETP_Thy3 <- ETP_Thy3[,c(-2,-3)]
MLB_CB1 <- MLB_CB1[,c(-2,-3)]
MLB_CB2 <- MLB_CB2[,c(-2,-3)]
MLB_CB3 <- MLB_CB3[,c(-2,-3)]
MLB_CB5 <- MLB_CB5[,c(-2,-3)]
HSC1_CB1 <- HSC1_CB1[,c(-2,-3)]
HSC1_CB2 <- HSC1_CB2[,c(-2,-3)]
HSC1_CB3 <- HSC1_CB3[,c(-2,-3)]
HSC1_CB4 <- HSC1_CB4[,c(-2,-3)]
HSC1_CB5 <- HSC1_CB5[,c(-2,-3)]

remove(i)

### Combine into matrix to test MLPs vs ETPs#####
mat.ETP.MLP <- cbind(HSC1_CB1[,1:2], HSC1_CB2[,2], HSC1_CB3[,2], HSC1_CB4[,2], HSC1_CB5[,2],
                 ETP_Thy1[,2], ETP_Thy2[,2], ETP_Thy3[,2],
                 MLB_CB1[,2], MLB_CB2[,2], MLB_CB3[,2],MLB_CB5[,2] )
#~It is at this point i realized i mislabeled MLPs for MLBs. That's ok,
# they are just names, so I'll rename my columns in the matrix.
colnames(mat.ETP.MLP) <- c("GeneID","HSC1","HSC2","HSC3","HSC4","HSC5", 
                           "ETP1", "ETP2","ETP3",
                           "MLP1","MLP2","MLP3","MLP5")
rownames(mat.ETP.MLP) <- 1:length(mat.ETP.MLP$GeneID)

#FC calculation for each gene between EPT and MLP
mat.ETP.MLP$FC.MLP.ETP <- NA

for (z in 1:29331) {
  numF1 <- is.na( mat.ETP.MLP[z, 7:9] ) #ETP NA
  numL1 <- as.numeric(3 - length(numF1[numF1==TRUE])) 
  numF2 <- is.na( mat.ETP.MLP[z, 10:13] ) #MLP NA
  numL2 <- as.numeric(4 - length(numF2[numF2==TRUE]))
  ETP.avg <- ( sum(mat.ETP.MLP$ETP1[z], mat.ETP.MLP$ETP2[z], mat.ETP.MLP$ETP3[z], na.rm = T) / numL1 )
  MLP.avg <- ( sum(mat.ETP.MLP$MLP1[z], mat.ETP.MLP$MLP2[z], mat.ETP.MLP$MLP3[z], mat.ETP.MLP$MLP5[z], na.rm = T) / numL2)
  mat.ETP.MLP$FC.MLP.ETP[z] <- ((ETP.avg - MLP.avg) / MLP.avg)
}
remove(numF1, numF2, z, numL1, numL2, ETP.avg, MLP.avg)

#Wilcox test to calculate for each gene between EPT and MLP
mat.ETP.MLP$pval.ETP.MLP <- NA
for (w in 2:29331) {
  ETP.group <- c(mat.ETP.MLP$ETP1[w], mat.ETP.MLP$ETP2[w], mat.ETP.MLP$ETP3[w]) 
  MLP.group <- c(mat.ETP.MLP$MLP1[w], mat.ETP.MLP$MLP2[w], mat.ETP.MLP$MLP3[w], mat.ETP.MLP$MLP5[w])
  standin <- NA
  ETP.NA <- is.na( ETP.group)
  MLP.NA <- is.na( MLP.group)
  ETP.NA <- as.numeric(length(ETP.NA[ETP.NA == FALSE]))
  MLP.NA <- as.numeric(length(MLP.NA[MLP.NA == FALSE]))
  standin <- ifelse(test = ETP.NA > 0 & MLP.NA > 1, yes = (wilcox.test(ETP.group, MLP.group)[3]), no = NA)
  mat.ETP.MLP$pval.ETP.MLP[w] <- standin
}

mat.ETP.MLP <- as.data.frame(mat.ETP.MLP)

write.table(mat.ETP.MLP, "~/Documents/Bioinformatics/LaurentiMicroArray/matrixFilter.txt", sep="\t")

####Looking for ID2 and E2A (Tcf3) 
#ID2 = NM_002166; E2A = NM_003200 & NM_001136139

####

index <- data.frame()

for (i in 1:33978) {
  probe_seq_name <- xx[i]
  probe_seq_name <- t(as.data.frame(probe_seq_name))
  if (is.null(probe_seq_name) == FALSE) {
    foo <- c(probe_seq_name, row.names(probe_seq_name) )
    index <- rbind(foo, index)
  }
  remove(probe_seq_name, current.probe)
  
} 

