rm(list=ls(all=TRUE));
options(stringsAsFactors=FALSE);

### Load RNAseq count database and genes to run the gsea data and libraries
RNAseq_Data=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Correlation_analysis/Data/Liver_GSE94660.csv", header=T)# to modify
rownames(RNAseq_Data)=RNAseq_Data$geneID
RNAseq=RNAseq_Data[,-1]
RNAseq=RNAseq[,-1]
Data=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/Mitocarta.csv", header = T) # To modify
Data=Data$geneID
Data=unique(Data)

library(WGCNA)
library(gridExtra)
library(fgsea)

### Correlation function
cor_genomewide_GOI=function(expr_test,gene){
  a=corAndPvalue(t(expr_test),t(expr_test[gene,]),method="spearman")
  a=as.data.frame(cbind(a$cor,a$p))
  colnames(a)=c("correlation","P_value")
  a=a[with(a,order(-correlation)),]
  a=a[-1,] # remove the gene of interest from the list (correlation= 1)
  return(a)
}

MitoMatrix=list()

### gsea pathways selected for the generation of the matrix
pathways_to_use = c("c2.cp.reactome","c2.cgp","c2.cp.kegg","c5.mf","c5.bp","c5.cc","h.all")
pathways_list = list()
path_gsea_annot="C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/gsea_pathways/genename_annotations_downloaded_MSigDB_v6.2/"
for(pathway in pathways_to_use){
  pathway_anno = gmtPathways(paste(path_gsea_annot,pathway,".v6.2.symbols.gmt",sep=""))
  pathways_list[[pathway]] = pathway_anno
  Details=read.csv(paste(path_gsea_annot,pathway,".csv",sep = ""), header = FALSE)
  MitoMatrix[[pathway]]=as.data.frame(Details$V1)
}

P=pathways_list
P_to_run=c(P$c2.cp.reactome,P$c2.cgp,P$c2.cp.kegg,P$c5.mf,P$c5.bp,P$c5.cc,P$h.all)

### Generating the matrix
minsizeconsider=10
maxsizeconsider=500
permutations=10000
genes_of_interest=subset(Data,Data %in% RNAseq_Data$geneID)
path_out="C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Correlation_analysis/results/NES_Matrices/Liver_GSE94660/" # To modify
a=as.data.frame(names(P_to_run))
rownames(a)=a[,1]
i=1
for(gene in genes_of_interest){
  i=i+1
  a=cbind(a,numeric(length(names(P_to_run))))
  names(a)[i]=gene
  cor=cor_genomewide_GOI(expr_test=RNAseq,gene=gene)
  cor2=merge(cor,RNAseq_Data[,1:2],by.x = "row.names", by.y = "geneID")
  cor2=cor2[with(cor2,order(-correlation)),]
  ranked_list=cor2$correlation
  names(ranked_list)=cor2$IDENTIFIER
  Rank_forgsea=ranked_list
  fgseaRes <- fgsea(pathways = P_to_run, stats = Rank_forgsea, minSize=minsizeconsider, maxSize=maxsizeconsider, nperm=permutations)
 # fgseaRes$NES[fgseaRes$pval > 0.05] <- 0   ######  CAREFUL WITH THIS LINE IF YOU WANT TO CONSIDER NOT SIGNIFICANT PATHWAYS  ######
  j=1
  for(path in fgseaRes$pathway){
    a[path,gene]=fgseaRes[j,NES]
    j=j+1
  }
  print(i-1)
}
a=a[,-1]

# Write the matrix in a csv file
write.csv(a,paste(path_out,"Mitocarta_NES_Matrix",".csv",sep="")) # to change