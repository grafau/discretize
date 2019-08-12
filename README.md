# discretize/PAMSil

PAMSil is an R script that clusters (PAM) individuals based on the distance matrix, simulates hybrids, and filters "genetic gradients" using empirical silhouette scores."

## Requirements

You need to have installed Rscript to run PAMSil from your commandline. Additionally, you will need to install following packages in your R instance: "shape", "cluster", "geosphere", "RColorBrewer", and "maps".

## Running

To run PAMSil you need to specify:
 - number of maximum k-populations to be investigated, 
 - path to file containing meta information in four columns (ID, latitude, longitude, country/population_name), without headers
 - path to file with pairwise distance matrix, without headers or row names, should be in the same order as meta file
 - optionally, you can specify your color matrix for each k (columns) and for each population (rows)
 
 
## Outputs

PAMSil will give you multiple output files for each k:
 - "maps" files contain maps with indviduals colored based on clustering/filtering algorithm (filtered inviduals will appear grey)
 - "clust" files contain MDS plots with indviduals colored based on clustering/filtering algorithm, additional pages include piecharts for country of origin, distance matrix heatmap, and clustering barplots.
 - "list2keep" files contain list of indviduals in discrete clusters that can be used in plink to filter out previous files
 - "list2anno" files conatin list of all indviduals with assignation ti discrete clusters

 
 
