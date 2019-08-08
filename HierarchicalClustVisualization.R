####Load Imputed Data####
library("dplyr", lib.loc="~/R/win-library/3.5")
library("tsoutliers", lib.loc="~/R/win-library/3.5")
library("TSclust", lib.loc="~/R/win-library/3.5")
library("factoextra", lib.loc="~/R/win-library/3.5")
library("dendextend", lib.loc="~/R/win-library/3.5")
library("cluster", lib.loc="C:/Program Files/R/R-3.5.1/library")
library("NbClust", lib.loc="~/R/win-library/3.5")
library("dtwclust", lib.loc="~/R/win-library/3.5")
setwd("C:/Users/sorel/Desktop/Paper ICA/Scripts")
MiceImputed=data.frame(read.csv("MiceImputed.csv", header = TRUE, sep = ","))
MiceImputed=MiceImputed[,3:dim(MiceImputed)[2]]
DataFeeders=MiceImputed[,2:dim(MiceImputed)[2]]

#Generate Outlier Free Data Set and store it
#####DATA SET WITHOUT OUTLIER####

#Hierarchical Clust using dendrogram#
#Generate Dendrogram Single Linkage
hcPER.single2=hclust(DissMat2, method="single")
plot(hcPER.single2)
rect.hclust(hcPER.single2, k=3, border=2:12) 
sub_grp.single2=cutree(hcPER.single2, k = 3)
table(sub_grp.single2)
outliers.single2=names(subset(sub_grp.single2, sub_grp.single2!=1))

CleanedDataSet2=DataFeeders%>%dplyr::select(-Alim_10, -Alim_87, -Alim_130, -Alim_158,
                                            -Alim_166, -Alim_173, -Alim_185, -Alim_205, -Alim_224, -Alim_267, -Alim_286,
                                            -Alim_322, -Alim_362, -Alim_366, -Alim_116, -Alim_329, -Alim_329, -Alim_13, -Alim_152,
                                            -Alim_214, -Alim_335, -Alim_346, -Alim_23, -Alim_101, -Alim_333, -Alim_110, -Alim_159, -Alim_167,
                                            -Alim_177)
TSMatrix3=t(as.matrix(CleanedDataSet2))
DissMat3=diss(TSMatrix3, "PER")

hcPER.single3=hclust(DissMat3, method="complete")
plot(hcPER.single3)
rect.hclust(hcPER.single3, k=3, border=2:12) 
sub_grp.single3=cutree(hcPER.single3, k = 3)
table(sub_grp.single3)
names(subset(sub_grp.single3, sub_grp.single3!=1))


#Hierarchical Clust using distance time warping#
CleanedDataSet2=DataFeeders%>%dplyr::select(-Alim_10, -Alim_87, -Alim_130, -Alim_158,
                                            -Alim_166, -Alim_173, -Alim_185, -Alim_205, -Alim_224, -Alim_267, -Alim_286,
                                            -Alim_322, -Alim_362, -Alim_366, -Alim_116, -Alim_329, -Alim_329, -Alim_13, -Alim_152,
                                            -Alim_214, -Alim_335, -Alim_346, -Alim_23, -Alim_101, -Alim_333, -Alim_110, -Alim_159, -Alim_167,
                                            -Alim_177)
#Transform distance matrix into list
TSList=list() 
for (i in 1:dim(CleanedDataSet2)[2]){
  TSList[[i]]=CleanedDataSet2[,i]
} 

#DoParallel routine
require("doParallel")
# Create parallel workers
workers <- makeCluster(6L)
# Preload dtwclust in each worker; not necessary but useful
invisible(clusterEvalQ(workers, library("dtwclust")))
# Register the backend; this step MUST be done
registerDoParallel(workers)

cluster_dtw_h=tsclust(TSList, type = "hierarchical", k = 3,  distance = "dtw_basic", centroid = shape_extraction,
                      control = hierarchical_control(method = "complete"), seed = 390, preproc = NULL, 
                      args = tsclust_args(dist = list(window.size = 5L)))

# Stop parallel workers
stopCluster(workers)
# Go back to sequential computation
registerDoSEQ()

plot(cluster_dtw_h)
rect.hclust(cluster_dtw_h, k=3, border=2:12) 
plot(cluster_dtw_h, type = "sc")
plot(cluster_dtw_h, type = "series", clus = 1)
plot(cluster_dtw_h, type = "series", clus = 2)
plot(cluster_dtw_h, type = "series", clus = 3)
plot(cluster_dtw_h, type = "centroids", clus = 1L)
plot(cluster_dtw_h, type = "centroids", clus = 2L)
plot(cluster_dtw_h, type = "centroids", clus = 3L)


###SELECT MEMBERS OF CLUSTER 1###
clusters_index=cluster_dtw_h@cluster
outlierIndex=c(which(clusters_index==2),which(clusters_index==3)) #Cause I dropped two columns

CleanedDataSet2=DataFeeders%>%dplyr::select(-Alim_10, -Alim_87, -Alim_130, -Alim_158,
                                           -Alim_166, -Alim_173, -Alim_185, -Alim_205, -Alim_224, -Alim_267, -Alim_286,
                                           -Alim_322, -Alim_362, -Alim_366, -Alim13, -Alim145, -Alim324, -Alim201)

#Hierarchical Clust#
#Generate Dendrogram Single Linkage
hcPER.single2=hclust(DissMat2, method="single")
plot(hcPER.single2)
rect.hclust(hcPER.single2, k=3, border=2:12) 
sub_grp.single2=cutree(hcPER.single2, k = 3)
table(sub_grp.single2)
outliers.single2=names(subset(sub_grp.single2, sub_grp.single2!=1))

CleanedDataSet2=DataFeeders%>%dplyr::select(-Alim_10, -Alim_87, -Alim_130, -Alim_158,
                                           -Alim_166, -Alim_173, -Alim_185, -Alim_205, -Alim_224, -Alim_267, -Alim_286,
                                           -Alim_322, -Alim_362, -Alim_366, -Alim_116, -Alim_329, -Alim_329, -Alim_13, -Alim_152,
                                           -Alim_214, -Alim_335, -Alim_346, -Alim_23, -Alim_101, -Alim_333, -Alim_110, -Alim_159, -Alim_167,
                                           -Alim_177)
TSMatrix3=t(as.matrix(CleanedDataSet2))
DissMat3=diss(TSMatrix3, "PER")

hcPER.single3=hclust(DissMat3, method="complete")
plot(hcPER.single3)
rect.hclust(hcPER.single3, k=3, border=2:12) 
sub_grp.single3=cutree(hcPER.single3, k = 3)
table(sub_grp.single3)
names(subset(sub_grp.single3, sub_grp.single3!=1))