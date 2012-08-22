options(stringsAsFactors=FALSE)

exp_rest <- read.csv(file = "REST_cofactors/expression_restday0_cofactors.csv")

#split xpn into time
day0 <- exp_rest[which(abs(exp_rest[,"day0"])>= 1),"EnsemblID"]
day1 <- exp_rest[which(abs(exp_rest[,"day1"])>= 1),"EnsemblID"]
day2 <- exp_rest[which(abs(exp_rest[,"day2"])>= 1),"EnsemblID"]
day4 <- exp_rest[which(abs(exp_rest[,"day4"])>= 1),"EnsemblID"]

day14 <- unique(c(day1, day2, day4))
day14 <- day14[which(!(day14 %in% day0))]

#get REST data
rest0 <- read.csv(file = "REST_ChIP/results/REST_D0_nearest_peak_to_gene_TSS.csv")
rest4 <- read.csv(file = "REST_ChIP/results/REST_D4_nearest_peak_to_gene_TSS.csv")

## take nearest REST peak
rest0.o <- rest0[order(abs(rest0[,"distancetoFeature"]),decreasing = TRUE),]
rest0.od <- rest0.o[!duplicated(rest0.o[,"Peak"]),]

rest4.o <- rest4[order(abs(rest4[,"distancetoFeature"]),decreasing = TRUE),]
rest4.od <- rest4.o[!duplicated(rest4.o[,"Peak"]),]

##merge in REST day4 data
exp_rest <- exp_rest[,c(2,3,4,5,6,7,8,9,14,16,18,19,20,21,22,23,24,25,26,27,28,29)]
colnames(exp_rest)[8:10] <- paste(colnames(exp_rest)[8:10],"_DAY0",sep="")

rest4 <- rest4.od[,c("EnsemblID","Peak","neg10log10pVal","distancetoFeature")]
colnames(rest4) <- c("EnsemblID","Peak_REST_DAY4","neg10log10pVal_REST_DAY4","distancetoFeature_REST_DAY4"

exp_rest_plus <- merge(exp_rest,rest4, by.x="EnsemblID",by.y="EnsemblID_DAY4", all.x=TRUE)

##merge with DeSeq data for K9/H4ac

ctrl_k9 <- read.csv(file = "Ac_ChIP/results/new_annotation/ctrl_h3k9ac/ctrl_h3k9ac_nearest_peak_to_gene_TSS.csv")
rest_k9 <- read.csv(file = "Ac_ChIP/results/new_annotation/rest_h3k9ac/rest_h3k9ac_nearest_peak_to_gene.csv")
ctrl_h4 <- read.csv(file = "Ac_ChIP/results/new_annotation/ctrl_h4ac/ctrl_h4ac_nearest_peak_to_gene_TSS.csv")
rest_h4 <- read.csv(file = "Ac_ChIP/results/new_annotation/rest_h4ac/rest_h4ac_nearest_peak_to_gene_TSS.csv")

ctrl_k9 <- ctrl_k9[which(ctrl_k9[,"neg10log10pVal"] >=100),]
rest_k9 <- rest_k9[which(rest_k9[,"neg10log10pVal"] >=100),]
ctrl_h4 <- ctrl_h4[which(ctrl_h4[,"neg10log10pVal"] >=100),]
rest_h4 <- rest_h4[which(rest_h4[,"neg10log10pVal"] >=100),]

##chop to size

ctrl_k9 <- ctrl_k9[,c("Peak","Peak_width","neg10log10pVal","EnsemblID","distancetoFeature")]
rest_k9 <- rest_k9[,c("Peak","Peak_width","neg10log10pVal","EnsemblID","distancetoFeature")]
ctrl_h4 <- ctrl_h4[,c("Peak","Peak_width","neg10log10pVal","EnsemblID","distancetoFeature")]
rest_h4 <- rest_h4[,c("Peak","Peak_width","neg10log10pVal","EnsemblID","distancetoFeature")]

colnames(ctrl_k9) <- paste(colnames(ctrl_k9),"_CTRL_K9ac",sep="")
colnames(rest_k9) <- paste(colnames(rest_k9),"_REST_K9ac",sep="")
colnames(ctrl_h4) <- paste(colnames(ctrl_h4),"_CTRL_H4ac",sep="")
colnames(rest_h4) <- paste(colnames(rest_h4),"_REST_H4ac",sep="")

# pull in deseq data

deseq_k9 <- read.csv(file = "Ac_ChIP/results/DeSeq_K9ac/peak_compare.csv")
deseq_h4 <- read.csv(file = "Ac_ChIP/results/DeSeq_H4/peak_compare.csv")

deseq_k9 <- deseq_k9[,c("id","log2FoldChange","padj")]
colnames(deseq_k9) <- c("Peak_K9","log2FoldChange_K9Ac","padj_K9ac")

deseq_h4 <- deseq_h4[,c("id","log2FoldChange","padj")]
colnames(deseq_h4) <- c("Peak_H4","log2FoldChange_H4Ac","padj_H4ac")

## and merge
k9_merge <- merge(ctrl_k9,rest_k9, by.x = "EnsemblID_CTRL_K9ac", by.y = "EnsemblID_REST_K9ac", all.x=T,all.y=T)

h4_merge <- merge(ctrl_h4,rest_h4, by.x = "EnsemblID_CTRL_H4ac", by.y = "EnsemblID_REST_H4ac", all.x=T,all.y=T)

ac_merge <- merge(k9_merge,h4_merge, by.x = "EnsemblID_CTRL_K9ac","EnsemblID_CTRL_H4ac",all.x = T, all.y=T)

##merge with Deseq data

ac_merge_k9deseq <- merge(ac_merge,deseq_k9, by.x = "Peak_REST_K9ac", by.y = "Peak_K9", all.x = T)

ac_merge_h4deseq <- merge(ac_merge_k9deseq,deseq_h4, by.x = "Peak_REST_H4ac", by.y = "Peak_H4", all.x = T)

res_merge <- merge(exp_rest_plus,ac_merge_h4deseq,by.x = "EnsemblID",by.y = "EnsemblID_CTRL_K9ac", all.x = TRUE)

colnames(res_merge)[23:25] <- c("Peak_REST_DAY4","neg10log10pVal_REST_DAY4","distancetoFeature_REST_DAY4")

res_merge <- res_merge[,c(1,2,3,4,5,6,7,8,9,10,23,24,25,11,12,13,14,15,16,17,18,19,20,21,22,28,29,30,31,27,32,33,34,35,36,37,38,26,39,40,41,42,43,44,45)]

write.csv(res_merge, file = "all_xpn_rest_rcor_ac_deseq.csv")


##plot some stuff

#plot change in exp at day 1-4 vs acetylation changes

plot(res_merge[which(res_merge[,"EnsemblID"] %in% day14),"log2FoldChange_K9Ac"],res_merge[which(res_merge[,"EnsemblID"] %in% day14),"day1"])














