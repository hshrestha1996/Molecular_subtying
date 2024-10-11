getwd()
library(WGCNA)
library(ConsensusClusterPlus)
library(MSnSet.utils)
library(ggplot2)

#####load quant data and select the sig proteins #######
DE_result <- read.csv("./input/sig_proteins.csv", row.names = 1)
sig <- DE_result[DE_result$Combined_P.Value<0.01,]

data <- read.csv("./input/Example_proteomicsdata.csv", row.names = 1)
data <- data[,intersect(colnames(data), rownames(sig))]

#####only keep the AD samples for subtypes #######
meta <- read.csv("./input/metadata.csv", row.names = 1)
meta <- meta[meta$Diagnosis == "AD",]
data <- data[intersect(rownames(data), rownames(meta)),]
yan <- data

##-------z-score---------------#
rnames<- rownames(yan)
yan <- apply(t(yan), 1, scale)
rownames(yan) <- rnames

####-----------------Check for outliers and remove them-------------#
sampleTree = hclust(dist(yan), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 120, col = "red") #choose best value for you
#------------keep the cluster below the red line---------------------
clust = cutreeStatic(sampleTree, cutHeight = 120, minSize = 10)
table(clust)
keepSamples = (clust== 1) #Cluster 1 contains the samples we want.
#keepSamples <- clust %in% c(1, 2)
yan = yan[keepSamples, ]
dim(yan)
##--------------------------------------------------------#

###----------run ConsensusClusterPlus--------------------######
yan <- yan[complete.cases(yan), ]
yan <- t(yan)
boxplot(yan)

set.seed(123)
Linkages <- c('average', 'complete', 'ward.D2') #'average', 'complete', 'ward.D2'
distances <- c('spearman', 'euclidean','pearson') #, 'spearman', 'euclidean','pearson'
clusterAlgs <- c('pam','hc', 'km') #'pam','hc', 'km'
for (distance in distances) {
  for (Linkage in Linkages) {
    for (clusterAlg in clusterAlgs) {
      # Create a unique title for each combination
      title <- paste0("AD_Em_sig0.01_Validation_", distance, "_", Linkage, "_", clusterAlg)
      
      # Run the ConsensusClusterPlus function with the current combination
      results <- ConsensusClusterPlus(
        yan,
        maxK = 10,
        reps = 1000,
        pItem = 0.9,
        pFeature = 1,
        title = title,
        distance = distance,
        writeTable = TRUE,
        seed = 123,
        plot = 'png',
        verbose = TRUE,
        clusterAlg = clusterAlg,
        innerLinkage = Linkage,
        finalLinkage = Linkage
      )
      
      icl = calcICL(results,title=title,plot="png", writeTable=TRUE)
      
      data <- yan
      data <- t(data)
      #data <- log2(data)
      data <- data[order(rownames(data)), ]
      
      
      for (i in 2:10) {
        # Create the meta data frame
        meta <- as.data.frame(results[[i]]$consensusClass)
        meta["subtype"] <- paste0("Subtype", meta$`results[[i]]$consensusClass`)
        metadata <- meta
        
        # Clean the metadata
        metadata <- metadata[complete.cases(metadata), ]
        metadata <- metadata[order(rownames(metadata)), ]
        
        # Keep only the common channel
        common_channel <- intersect(rownames(data), rownames(metadata))
        data_filtered <- data[common_channel, ]
        metadata <- metadata[common_channel, ]
        
        # data
        dta <- t(data_filtered)
        dta <- as.matrix(dta)
        
        # Create MSnSet object
        msnset <- MSnSet(exprs = dta)
        annotated_data <- AnnotatedDataFrame(metadata)
        phenoData(msnset) <- annotated_data
        
        # Plot PCA and save the plot
        plot_pca(msnset, phenotype = "subtype", legend_title = "Subtype")
        ggsave(paste0(title,"/km_", i, "subtypes.png"),width = 4,height = 4,dpi = 300)
        #plot_pca(msnset, phenotype = "subtype", legend_title = "Subtype", components = c(1,3))
        #ggsave(paste0(title,"/km_", i, "subtypes_c13.png"),width = 4,height = 4,dpi = 300)
      }
      
    }}
  
}


##--------------------------------------------------------------------------------------#
##--------------------------------------------------------------------------------------#