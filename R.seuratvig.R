#Seurat Practice using 10X Sample PBMC data
#E.C. Jan 2021

# loading packages
library(dplyr)
library(Seurat)
library(patchwork)

#Loading sample dataset
pbmb.data <- Read10X(data.dir = "../Documents/scRNA\ Fidanza2020/filtered_gene_bc_matrices/hg19/")

#create objest compatible with Seurat
pbmc <- CreateSeuratObject(counts = pbmb.data, project = "pbmc3k", min.cells = 3, min.features = 200)
#test to see the seurat object was made correctly
pbmc

#brackets add an element (column) to our data structure. This element is the percent of mitochondrial genes found per cell
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#to see the QC metrics that brought about the plots above see:
head(pbmc@meta.data, 10)

#Use FeatureScatter to visualize above metrics in scatter plot:
mt.vs.RNAcount <- FeatureScatter(pbmc, feature2 = "percent.mt", feature1 = "nCount_RNA")
mt.vs.RNAcount

#filter out the cells with less than 200 reads, more than 2500 reads and greater than 5% mitochdnrial genes
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

###Normalization###
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#the resulting notrmalization values are stores in pbmc[["RNA"]]@data. The above parameters are the defaults. This is a verbose example

##Adjusting data of interest##
#Following code subsets the data to only analyze the top 2000 genes with high variability.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#read out the top ten variables genes
top10.var.genes <- head(VariableFeatures(pbmc),10)
top10.var.genes

#quick plot to see variable genes
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot=plot1, points = top10.var.genes, repel = TRUE)

#scale gene expresions to prepare for PCA
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#To take a look at an example of the scaled data (remember, these genes are the most variable):
pbmc[["RNA"]]@scale.data[1:5,1:5]

###Dimensionality Reduction###
pbmc <- RunPCA(pbmc, features = VariableFeatures(object=pbmc))

#Visualizing the PCA
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#Plot to top genes associated with reduction components (in this case only PC1, 2, and 3)
VizDimLoadings(pbmc, dims = 1:3, reduction = "pca")

#PCA plot with PC1 and PC2
DimPlot(pbmc, dims = c(1,2), reduction = "pca")

#Look at heat maps where PC 
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#Loose way of determining how many primary components to include, where does the data level out? In this case it is around 12.
#so our first 10 principle components will be what we use for downstream analysis 
ElbowPlot(pbmc)

#Next we use our PCs (as inputs) and a Seurat algorithm to determine clusters. 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5) #resolution can be changed, it can be higher than 0.5 for data set with more than 3k cells

#optional: see how the cells are labelled by cluster
head(Idents(pbmc), 10) #Note: the # of levels should be the numbner  of clusters you have

#calculate the UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

#Plot the UMAP
DimPlot(pbmc, reduction = "umap")


