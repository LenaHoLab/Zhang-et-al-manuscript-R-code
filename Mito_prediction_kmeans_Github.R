rm(list=ls(all=TRUE));
options(stringsAsFactors=FALSE);

### Load libraries
library(gplots)
library(fgsea)
library(ggplot2)
library(ggfortify)
library(factoextra)
library(dendextend)

### Data locations and data
folder="C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Correlation_analysis/results/NES_Matrices/Heart/"
path_gsea_annot="C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Correlation_analysis/gsea_pathways/genename_annotations_downloaded_MSigDB_v6.2/"
pathways_to_use=c("c2.cp.reactome","c2.cgp","c2.cp.kegg","c5.cc","c5.bp","c5.mf","h.all")
genes_to_use=c("Mitocarta","Negative_control","Candidates")
Colors=c("#FF3333","#33CC33","#3399FF")


### Function to generate the vector of pathway names that contain less than N genes
pathways_N=function(N,path_gsea_annot){
  pathway_N=c()
  for(pathway in pathways_to_use){
    d=c()
    pathway_anno = gmtPathways(paste(path_gsea_annot,pathway,".v6.2.symbols.gmt",sep=""))
    c=names(pathway_anno)
    for (p in pathway_anno) {
      d=rbind(d,as.integer(length(p)))
    }
    pathway_N=rbind(pathway_N,cbind(c,d))
  }
  pathway_N=subset(pathway_N,as.integer(pathway_N[,2])<N)
}

############### MAIN ###############
##### Step 1 : Prepare data #####
D=list()
for (G in genes_to_use) {
  D[[G]]=read.csv(paste(folder,G,"_NES_Matrix.csv",sep=""),header = T,row.names = 1)
}

### Remove any common gene between candidates, negative control and Mitocarta
if(length(genes_to_use==3)){
  D[[1]]=D[[1]][,!colnames(D[[1]])%in%colnames(D[[3]])]
  D[[2]]=D[[2]][,!colnames(D[[2]])%in%colnames(D[[3]])]
}

### Select pathways of interest and create the Dataset names list
all_pathways=rownames(D[[1]])
gsea_N=pathways_N(500,path_gsea_annot = path_gsea_annot)
all_pathways=subset(all_pathways,all_pathways%in%gsea_N[,1])
Pathways_to_keep=subset(all_pathways,all_pathways%in%rownames(D[[1]][rowSums(abs(D[[1]])>3) >100,])) # Criteria for Mito pathway
Dataset=c()
for(d in 1:length(genes_to_use)){
  D[[d]]=D[[d]][rownames(D[[d]])%in%Pathways_to_keep,]
  Dataset=c(Dataset,rep(genes_to_use[d],ncol(D[[d]])))
}

# Merging all matrices
Matrix=as.matrix(merge(D[[1]],D[[2]],by="row.names",all.x=TRUE))
rownames(Matrix)=Matrix[,1]
Matrix=Matrix[,-1]
Matrix=as.matrix(merge(Matrix,D[[3]],by="row.names",all.x=TRUE))
rownames(Matrix)=Matrix[,1]
Matrix=Matrix[,-1]
class(Matrix) <- "numeric"


##### Step 2 : PCA and k-means clusters #####
### Visualize PCA
dat=rbind(Dataset,Matrix)
rownames(dat)[1]="Dataset"
PCA=prcomp(t(Matrix))
autoplot(PCA,data = t(dat),col='Dataset')+scale_color_manual(values = Colors)+theme_classic()

######## Generate clusters and score the candidates based on a training set of Mito-sORFs and non-mito-sORFs
Tr=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/Training_sORFs.csv",header = T)
rownames(Tr)=Tr[,1]
Tr=Tr[Tr[,1]%in%colnames(D[[3]]),]
ltrain=length(Tr[,1])
best_pred_perc=c("Nb of clusters","Best Percentage")
Cand=subset(colnames(D[[3]]),colnames(D[[3]])%in%rownames(Tr))
Cand=cbind(Cand,rep(0,length(Cand)))
rownames(Cand)=Cand[,1]
for (nb.clusters in 2:15) {
  k=kmeans(PCA$x,nb.clusters, nstart = 25)
  
  ### Obtain Mitocarta percentage per cluster
  M=as.data.frame(cbind(k$cluster,as.character(Dataset)))
  cluster_details=data.frame(rep(0,5))
  for (c in 1:nb.clusters) {
    clust=subset(M,M[,1]==c)
    proportion=data.frame(genes_to_use,freq=rowSums(!adist(genes_to_use,clust[,2],partial = T)))
    L=as.integer(proportion$freq[1]+proportion$freq[2])
    percentage=1/L*100*proportion$freq[1:2]
    cluster_details=cbind(cluster_details,c(proportion$freq,percentage))
    N=paste(c)
    colnames(cluster_details)[c+1]=N
    }
  cluster_details=cluster_details[,-1] # Cluster_details contains the number of TF Mito and Candidate genes for each cluster
                                       # It also contains the percentage of Mitocarta and TF genes for each cluster
  cluster_details=cluster_details[,order(cluster_details[4,])] # Sort cluster_details columns by Mitocarta percentage
  
  # Find the best threshold to apply to consider a candidate as predicted to be localized in mitochondria
  Prediction=c("Threshold","Correct") # Contain all the threshold values and their corresponding prediction percentage
  Training=Tr
  for (cl in 1:(nb.clusters-1)) {
    thres=(cluster_details[4,cl]+cluster_details[4,cl+1])/2
    Training=cbind(Training,rep(0,ltrain))
    for (g in rownames(Training)) {
      if (cluster_details[4,as.character(k$cluster[g])]>thres) {
        Training[g,cl+2]=1
      } else Training[g,cl+2]=0
    }
    pred_score=length(which(Training[,cl+2]==Training[,2]))/ltrain
    Prediction=cbind(Prediction,c(thres,pred_score))
  }
  rownames(Prediction)=Prediction[,1]
  # In the last [] of the next 2 lines, different thresholds can lead to the best prediction
  # For lowest false negative put 1
  # For lowest false positive put length(max(Prediction[,2]))
  P=Prediction[2,]
  P=as.numeric(P[-1])
  thres=Prediction[1,1+which(P==max(P))[1]] # Thres with the best prediction and lowest FN
  best_pred_perc=cbind(best_pred_perc,c(nb.clusters,max(P)[1])) # change 1 to length(max(Prediction)) for lowest FP
  
  # Calculate the kmean score pf each candidate
  Cand=cbind(Cand,rep(0,length(Cand[,1])))
  for (ca in rownames(Cand)) {
    if (cluster_details[4,as.character(k$cluster[ca])]>thres) {
      Cand[ca,nb.clusters]=1
    } else Cand[ca,nb.clusters]=0
  }
}
Cand=Cand[,-1]
Cand=Cand[,1:(ncol(Cand)-1)]

Score=as.data.frame(cbind(Tr$Mito,rep(0,ltrain)))
rownames(Score)=rownames(Tr)
for (i in rownames(Score)) {
  cs=as.numeric(Cand[i,])
  Score[i,2]=mean(cs)
}

FP=nrow(subset(Score,Score[,1]==0&Score[,2]>0.75))
FN=nrow(subset(Score,Score[,1]==1&Score[,2]<0.75))



