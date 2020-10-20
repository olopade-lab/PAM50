exprdata <- read.csv('~/Downloads/BatchQC/salmon_merged_gene_counts.csv') #Change to location of salmon_merged_gene_counts.csv file

for (i in (2:length(colnames(exprdata)))) {
  colnames(exprdata)[i]=substring(colnames(exprdata)[i],2)
}


annotData<-data.frame(
  probe=exprdata$gene_id,
  EntrezGene.ID=rep('na', nrow(exprdata))
  #Gene.Symbol=rep('na', nrow(exprdata))
)
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(
  attributes=c('ensembl_gene_id', "hgnc_symbol","entrezgene_id"),
  mart = mart)
head(genes)

write.csv(genes,'./entrez.csv')

for (i in 1:length(annotData$probe)) {
  index<-which(genes$ensembl_gene_id == annotData$probe[i])
  if (length(index)>0) {
    annotData$EntrezGene.ID[i]<-genes$entrezgene_id[index]
  }
  
  #annotData$Gene.Symbol[i]<-genes$hgnc_symbol[index]
}
row.names(annotData)=annotData$probe

row.names(exprdata) <- exprdata$gene_id
exprdata<-exprdata[2:dim(exprdata)[2]]
exprdata<-t(exprdata)


exprdata <- exprdata[rowSums(exprdata) > 0,]
QuantileExp <-
  apply(exprdata, 2, function(x){quantile(x[x>0], 0.75)});

exprdata <- t(t(exprdata) / QuantileExp);

posRows<-c('18RES10382_1R_T_S_R_RS_T',
           '18RES10385_1R_T_S_R_RS_T',
           '18RES10388_1R_T_S_R_RS_T',
           '18RES10393_1R_T_S_R_RS_T',
           '18RES10394_1R_T_S_R_RS_T',
           '18RES10398_1R_T_S_R_RS_T',
           '18RES10399_1R_T_S_R_RS_T',
           '18RES10414_1R_T_S_R_RS_T',
           '18RES10417_1R_T_S_R_RS_T',
           '18RES10418_1R_T_S_R_RS_T',
           '18RES10428_1R_T_S_R_RS_T',
           '18RES13578_1R_T_S_R_RS_T',
           '18RES13639_1R_T_S_R_RS_T',
           '18RES13652_1R_T_S_R_RS_T',
           '18RES13653_1R_T_S_R_RS_T',
           '18RES13654_1R_T_S_R_RS_T',
           '18RES16043_1R_T_S_R_RS_T',
           '18RES16061_1R_T_S_R_RS_T',
           '18RES16066_1R_T_S_R_RS_T',
           '18RES16067_1R_T_S_R_RS_T',
           '18RES16068_1R_T_S_R_RS_T',
           '18RES16164_1R_T_S_R_RS_T',
           '18RES16179_1R_T_S_R_RS_T',
           '18RES16182_1R_T_S_R_RS_T',
           '18RES04129_1R_T_S_R_RS_T',
           '18RES04372_1R_T_S_R_RS_T',
           '18RES04887_1R_T_S_R_RS_T',
           '18RES04412_1R_T_S_R_RS_T',
           '18RES04367_1R_T_S_R_RS_T',
           '18RES04366_1R_T_S_R_RS_T',
           '18RES04728_1R_T_S_R_RS_T',
           '18RES09403_1R_T_S_R_RS_T',
           '18RES04364_1R_T_S_R_RS_T',
           '18RES04407_1R_T_S_R_RS_T',
           '18RES09402_1R_T_S_R_RS_T')

negRows<-c('18RES10383_1R_T_S_R_RS_T',
           '18RES10389_1R_T_S_R_RS_T',
           '18RES10392_1R_T_S_R_RS_T',
           '18RES10395_1R_T_S_R_RS_T',
           '18RES10401_1R_T_S_R_RS_T',
           '18RES10402_1R_T_S_R_RS_T',
           '18RES10403_1R_T_S_R_RS_T',
           '18RES10405_1R_T_S_R_RS_T',
           '18RES10409_1R_T_S_R_RS_T',
           '18RES10413_1R_T_S_R_RS_T',
           '18RES10426_1R_T_S_R_RS_T',
           '18RES10430_2R_T_S_R_RS_T',
           '18RES13572_1R_T_S_R_RS_T',
           '18RES13580_1R_T_S_R_RS_T',
           '18RES13642_1R_T_S_R_RS_T',
           '18RES13645_1R_T_S_R_RS_T',
           '18RES13651_1R_T_S_R_RS_T',
           '18RES13656_1R_T_S_R_RS_T',
           '18RES13663_1R_T_S_R_RS_T',
           '18RES16069_1R_T_S_R_RS_T',
           '18RES16168_1R_T_S_R_RS_T',
           '18RES16169_1R_T_S_R_RS_T',
           '18RES04413_1R_T_S_R_RS_T',
           '18RES04251_1R_T_S_R_RS_T',
           '18RES04370_1R_T_S_R_RS_T',
           '18RES04260_1R_T_S_R_RS_T',
           '18RES04368_1R_T_S_R_RS_T',
           '18RES04226_1R_T_S_R_RS_T',
           '18RES04403_1R_T_S_R_RS_T',
           '18RES04249_1R_T_S_R_RS_T',
           '18RES04254_1R_T_S_R_RS_T',
           '18RES04139_1R_T_S_R_RS_T',
           '18RES04220_1R_T_S_R_RS_T',
           '18RES04369_1R_T_S_R_RS_T',
           '18RES04371_1R_T_S_R_RS_T',
           '18RES09399_1R_T_S_R_RS_T',
           '18RES04227_1R_T_S_R_RS_T',
           '18RES09394_1R_T_S_R_RS_T',
           '18RES04404_1R_T_S_R_RS_T',
           '18RES04229_1R_T_S_R_RS_T',
           '18RES04885_1R_T_S_R_RS_T',
           '18RES04401_1R_T_S_R_RS_T',
           '18RES04138_1R_T_S_R_RS_T',
           '18RES04410_1R_T_S_R_RS_T',
           '18RES04408_1R_T_S_R_RS_T',
           '18RES04258_1R_T_S_R_RS_T',
           '18RES09404_1R_T_S_R_RS_T',
           '18RES04250_1R_T_S_R_RS_T')

negRows<-sample(negRows, 35, replace = FALSE, prob = NULL)
write.csv(negRows, './negRows.csv')
impRows<-c(negRows, posRows)
exprdata<-exprdata[impRows,]


library(genefu)
data(pam50)

impCols<-c()
for (i in pam50$centroids.map$EntrezGene.ID) {
  a<-which(annotData$EntrezGene.ID == i )
  a<-annotData[a,"probe"]
  impCols<-c(impCols,a)
}

exprdata<-exprdata[,impCols]

rowmed <- apply(df,1,exprdata)
exprdata<-exprdata - rowmed




PAM50Preds<-molecular.subtyping(sbt.model = "pam50", data=exprdata, annot=annotData,do.mapping=TRUE)


Subtypes<-PAM50Preds$subtype
write.csv(Subtypes,'./Subtypes.csv')
probs<-PAM50Preds$subtype.proba
write.csv(probs,'./probs.csv')
BinaryMap<-PAM50Preds$subtype.crisp
write.csv(genes,'./BinaryMap.csv')

Subtypes<-read.csv('./Subtypes.csv')

Basal<-0
Normal<-0
Her2<-0
LumA<-0
LumB<-0

for (k in Subtypes$x) {
  if (k=='Basal') {
    Basal<-Basal+1
  } else if (k=='Normal'){
    Normal<-Normal+1
  } else if (k=='Her2'){
    Her2<-Her2+1
  }
  else if (k=='LumA'){
    LumA<-LumA+1
  }
  else if (k=='LumB'){
    LumB<-LumB+1
  }
}

df <- data.frame(Type=c("Normal", "Basal", "Her2", "LumA", "LumB"),
                 Frequency=c(Normal, Basal, Her2, LumA, LumB))

library(ggplot2)
p<-ggplot(data=df, aes(x=Type, y=Frequency)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal()

ggsave('./pam50freq.pdf',plot=p)
  
