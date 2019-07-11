###
# heart_control_data=readRDS("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Correlation_analysis/Data/1.Controls_vstcounts_data.rds")
# 1.Heart.csv
# 2.DCM.csv
# 9.1.Liver_GSE94660.csv
# 10.SKM.csv
# Homo_sapiens_TF.csv
# Inflammation.csv
# Human.MitoCarta2.0.csv
###############################
###############################

### Clear loaded data
rm(list=ls(all=TRUE));
options(stringsAsFactors=FALSE);

### Load RNAseq count database and genes to run the gsea data and libraries
RNAseq_Data=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Correlation_analysis/Data/10.SKM.csv", header=T)# to change
rownames(RNAseq_Data)=RNAseq_Data$geneID
RNAseq=RNAseq_Data[,-1]
RNAseq=RNAseq[,-1]
Data=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/MSUF.csv", header = T) # To change
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

### gsea pathways
pathways_to_use = c("c2.cp.reactome","c2.cgp","c2.cp.kegg","c5.mf","c5.bp","c5.cc","h.all")
pathways_list = list()
path_gsea_annot="C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Correlation_analysis/gsea_pathways/genename_annotations_downloaded_MSigDB_v6.2/"
for(pathway in pathways_to_use){
  pathway_anno = gmtPathways(paste(path_gsea_annot,pathway,".v6.2.symbols.gmt",sep=""))
  pathways_list[[pathway]] = pathway_anno
  Details=read.csv(paste(path_gsea_annot,pathway,".csv",sep = ""), header = FALSE)
  MitoMatrix[[pathway]]=as.data.frame(Details$V1)
}

P=pathways_list
P_to_run=c(P$c2.cp.reactome,P$c2.cgp,P$c2.cp.kegg,P$c5.mf,P$c5.bp,P$c5.cc,P$h.all)

### Build of matrix
minsizeconsider=10
maxsizeconsider=500
permutations=10000
genes_of_interest=subset(Data,Data %in% RNAseq_Data$geneID) # To change
path_out="C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Correlation_analysis/results/NES_Matrices/SKM/" # To change
a=as.data.frame(names(P_to_run))
#genes_of_interest="ENSG00000178982"
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
 # fgseaRes$NES[fgseaRes$pval > 0.05] <- 0              ##########  CAREFUL WITH THIS LINE   ###########
  j=1
  for(path in fgseaRes$pathway){
    a[path,gene]=fgseaRes[j,NES]
    j=j+1
  }
  print(i-1)
}
a=a[,-1]
write.csv(a,paste(path_out,"RNF8_NES_Matrix",".csv",sep="")) # to change

# csv files to open

Mito=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/Human.MitoCarta2.0.csv", header = T)
TF=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/Homo_sapiens_TF.csv", header = T)
peptides=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/Unique_aa_sequences_V2.csv", header = T)
All_sorfs=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/all_sORFs.csv", header = T)
Random=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/Random_1000.csv", header = T)
Candidates=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/50candidates.csv", header = T)

### Generate a list of random genes which are not in TF or Mitocarta, open csv files

present=Reduce(c,list(Mito$EnsemblGeneID,TF$V1))
rest=subset(heart_control_data$gene_info$Ensembl_ID,!(heart_control_data$gene_info$Ensembl_ID %in% present))

Random=sample(rest,1000)

write.csv(Random,"C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/Random_1000.csv")