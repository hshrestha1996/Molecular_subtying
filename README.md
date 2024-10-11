# Molecular_subtying
This pipeline uses proteomics data to identify subtypes.

This R script is designed to analyze proteomics data for identifying distinct subtypes. The script normalizes the data using z-scores to standardize protein expression levels across samples. Outlier detection is performed through hierarchical clustering, allowing for the identification and removal of outlier samples to ensure a more robust analysis.

It then establishes parameters for clustering, including various linkage methods, distance metrics, and clustering algorithms. An iterative process runs through these combinations, utilizing the ConsensusClusterPlus function to identify clusters while calculating the Integrated Cluster Loss (ICL) to assess clustering stability. Finally, PCA plots are generated for each subtype, which are saved as PNG files in designated folders.
