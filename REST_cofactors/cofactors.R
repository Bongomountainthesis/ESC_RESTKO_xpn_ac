options(stringsAsFactors=FALSE)

library(IRanges)
library(ChIPpeakAnno)

## import raw data from Yu et al 2011 Genome Res (Hong Bing and Rory's paper)
rest.hb <- read.csv(file = "GSM698696_REST_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)   
rcor1 <- read.csv(file = "GSM698697_RCOR1_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)
rcor2 <- read.csv(file = "GSM698698_RCOR2_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)
rcor3 <- read.csv(file = "GSM698699_RCOR3_ChIP-seqpeaks_mm9.bed", sep = "\t", header = F)
sin3a <- read.csv(file = "GSM698700_SIN3A_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)
sin3b <- read.csv(file = "GSM698701_SIN3B_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)

## add names column for each peak

rest.hb[,5] <- paste(rest.hb[,"V1"],paste(rest.hb[,"V2"],rest.hb[,"V3"],sep="-"),sep=":")
rcor1[,5] <- paste(rcor1[,"V1"],paste(rcor1[,"V2"],rcor1[,"V3"],sep="-"),sep=":")
rcor2[,5] <- paste(rcor2[,"V1"],paste(rcor2[,"V2"],rcor2[,"V3"],sep="-"),sep=":")
rcor3[,5] <- paste(rcor3[,"V1"],paste(rcor3[,"V2"],rcor3[,"V3"],sep="-"),sep=":")
sin3a[,5] <- paste(sin3a[,"V1"],paste(sin3a[,"V2"],sin3a[,"V3"],sep="-"),sep=":")
sin3b[,5] <- paste(sin3b[,"V1"],paste(sin3b[,"V2"],sin3b[,"V3"],sep="-"),sep=":")

## make RD objects from each

rest.hb.rd <- RangedData(ranges=IRanges(
			start=as.numeric(rest.hb$V2),
			end= as.numeric(rest.hb$V3),
			names=as.character(rest.hb$V5)
			),	
			space=as.character(rest.hb$V1)
			)
rcor1.rd <- RangedData(ranges=IRanges(
			start=as.numeric(rcor1$V2),
			end= as.numeric(rcor1$V3),
			names=as.character(rcor1$V5)
			),	
			space=as.character(rcor1$V1)
			)
rcor2.rd <- RangedData(ranges=IRanges(
			start=as.numeric(rcor2$V2),
			end= as.numeric(rcor2$V3),
			names=as.character(rcor2$V5)
			),	
			space=as.character(rcor2$V1)
			)
rcor3.rd <- RangedData(ranges=IRanges(
			start=as.numeric(rcor3$V2),
			end= as.numeric(rcor3$V3),
			names=as.character(rcor3$V5)
			),	
			space=as.character(rcor3$V1)
			)
sin3a.rd <- RangedData(ranges=IRanges(
			start=as.numeric(sin3a$V2),
			end= as.numeric(sin3a$V3),
			names=as.character(sin3a$V5)
			),	
			space=as.character(sin3a$V1)
			)
sin3b.rd <- RangedData(ranges=IRanges(
			start=as.numeric(sin3b$V2),
			end= as.numeric(sin3b$V3),
			names=as.character(sin3b$V5)
			),	
			space=as.character(sin3b$V1)
			)

## import raw data from Yu et al 2011 Genome Res (Hong Bing and Rory's paper)
rest.hb <- read.csv(file = "GSM698696_REST_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)   
rcor1 <- read.csv(file = "GSM698697_RCOR1_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)
rcor2 <- read.csv(file = "GSM698698_RCOR2_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)
rcor3 <- read.csv(file = "GSM698699_RCOR3_ChIP-seqpeaks_mm9.bed", sep = "\t", header = F)
sin3a <- read.csv(file = "GSM698700_SIN3A_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)
sin3b <- read.csv(file = "GSM698701_SIN3B_ChIP-seq_peaks_mm9.bed", sep = "\t", header = F)

## add names column for each peak

rest.hb[,5] <- paste(rest.hb[,"V1"],paste(rest.hb[,"V2"],rest.hb[,"V3"],sep="-"),sep=":")
rcor1[,5] <- paste(rcor1[,"V1"],paste(rcor1[,"V2"],rcor1[,"V3"],sep="-"),sep=":")
rcor2[,5] <- paste(rcor2[,"V1"],paste(rcor2[,"V2"],rcor2[,"V3"],sep="-"),sep=":")
rcor3[,5] <- paste(rcor3[,"V1"],paste(rcor3[,"V2"],rcor3[,"V3"],sep="-"),sep=":")
sin3a[,5] <- paste(sin3a[,"V1"],paste(sin3a[,"V2"],sin3a[,"V3"],sep="-"),sep=":")
sin3b[,5] <- paste(sin3b[,"V1"],paste(sin3b[,"V2"],sin3b[,"V3"],sep="-"),sep=":")

## make RD objects from each

rest.hb.rd <- RangedData(ranges=IRanges(
			start=as.numeric(rest.hb$V2),
			end= as.numeric(rest.hb$V3),
			names=as.character(rest.hb$V5)
			),	
			space=as.character(rest.hb$V1)
			)
rcor1.rd <- RangedData(ranges=IRanges(
			start=as.numeric(rcor1$V2),
			end= as.numeric(rcor1$V3),
			names=as.character(rcor1$V5)
			),	
			space=as.character(rcor1$V1)
			)
rcor2.rd <- RangedData(ranges=IRanges(
			start=as.numeric(rcor2$V2),
			end= as.numeric(rcor2$V3),
			names=as.character(rcor2$V5)
			),	
			space=as.character(rcor2$V1)
			)
rcor3.rd <- RangedData(ranges=IRanges(
			start=as.numeric(rcor3$V2),
			end= as.numeric(rcor3$V3),
			names=as.character(rcor3$V5)
			),	
			space=as.character(rcor3$V1)
			)
sin3a.rd <- RangedData(ranges=IRanges(
			start=as.numeric(sin3a$V2),
			end= as.numeric(sin3a$V3),
			names=as.character(sin3a$V5)
			),	
			space=as.character(sin3a$V1)
			)
sin3b.rd <- RangedData(ranges=IRanges(
			start=as.numeric(sin3b$V2),
			end= as.numeric(sin3b$V3),
			names=as.character(sin3b$V5)
			),	
			space=as.character(sin3b$V1)
			)

## import our REST data
rest <- read.csv(file = "REST_D0_nearest_peak_to_gene_TSS.csv")

## take nearest REST peak
rest.o <- rest[order(abs(rest[,"distancetoFeature"]),decreasing = TRUE),]
rest.od <- rest.o[!duplicated(rest.o[,"Peak"]),]

rest.rd <- RangedData(ranges=IRanges(
			start=as.numeric(rest.od$Peak_start),
			end=as.numeric(rest.od$Peak_end),
			names=as.character(rest.od$Peak)
			),
			space=as.character(rest.od$Chromosome)
			)

rcor1.nearest <- annotatePeakInBatch(rcor1.rd,
                                   AnnotationData=rest.rd,
                                   PeakLocForDistance = "middle",
                                   FeatureLocForDistance = "middle",
                                   output = "both",
				   multiple = TRUE
                                   )

rcor1.nearest <- as.data.frame(rcor1.nearest)
rcor2.nearest <- annotatePeakInBatch(rcor2.rd,
                                   AnnotationData=rest.rd,
                                   PeakLocForDistance = "middle",
                                   FeatureLocForDistance = "middle",
                                   output = "both",
				   multiple = TRUE
                                   )

rcor2.nearest <- as.data.frame(rcor2.nearest)
rcor3.nearest <- annotatePeakInBatch(rcor3.rd,
                                   AnnotationData=rest.rd,
                                   PeakLocForDistance = "middle",
                                   FeatureLocForDistance = "middle",
                                   output = "both",
				   multiple = TRUE
                                   )

rcor3.nearest <- as.data.frame(rcor3.nearest)
sin3a.nearest <- annotatePeakInBatch(sin3a.rd,
                                   AnnotationData=rest.rd,
                                   PeakLocForDistance = "middle",
                                   FeatureLocForDistance = "middle",
                                   output = "both",
				   multiple = TRUE
                                   )

sin3a.nearest <- as.data.frame(sin3a.nearest)
sin3b.nearest <- annotatePeakInBatch(sin3b.rd,
                                   AnnotationData=rest.rd,
                                   PeakLocForDistance = "middle",
                                   FeatureLocForDistance = "middle",
                                   output = "both",
				   multiple = TRUE
                                   )

sin3b.nearest <- as.data.frame(sin3b.nearest)

##then only take cofactor binding events within 100bp of REST site

rcor1.near <- rcor1.nearest[which(abs(rcor1.nearest[,"distancetoFeature"]) <= 100),]
rcor1.near <- rcor1.near[,c("peak","feature","distancetoFeature")]
colnames(rcor1.near) <- c("Peak_RCOR1","feature","distancetoFeature_RCOR1")
rcor2.near <- rcor2.nearest[which(abs(rcor2.nearest[,"distancetoFeature"]) <= 100),]
rcor2.near <- rcor2.near[,c("peak","feature","distancetoFeature")]
colnames(rcor2.near) <- c("Peak_RCOR2","feature","distancetoFeature_RCOR2")
rcor3.near <- rcor3.nearest[which(abs(rcor3.nearest[,"distancetoFeature"]) <= 100),]
rcor3.near <- rcor3.near[,c("peak","feature","distancetoFeature")]
colnames(rcor3.near) <- c("Peak_RCOR3","feature","distancetoFeature_RCOR3")
sin3a.near <- sin3a.nearest[which(abs(sin3a.nearest[,"distancetoFeature"]) <= 100),]
sin3a.near <- sin3a.near[,c("peak","feature","distancetoFeature")]
colnames(sin3a.near) <- c("Peak_SIN3A","feature","distancetoFeature_SIN3A")
sin3b.near <- sin3b.nearest[which(abs(sin3b.nearest[,"distancetoFeature"]) <= 100),]
sin3b.near <- sin3b.near[,c("peak","feature","distancetoFeature")]
colnames(sin3b.near) <- c("Peak_SIN3B","feature","distancetoFeature_SIN3B")

## make dataframe with all REST peaks and then the cofactors they recruit
rest.df <- rest.od[,c("Peak","Chromosome","Peak_start","Peak_end","Peak_width","neg10log10pVal","FoldEnrichment","EnsemblID","distancetoFeature","Symbol")]
colnames(rest.df) <- paste(colnames(rest.df),"_REST",sep="")

rest.cofactors <- merge(rest.df, rcor1.near, by.x="Peak_REST",by.y="feature", all.x = TRUE)
rest.cofactors <- merge(rest.cofactors, rcor2.near, by.x="Peak_REST",by.y="feature", all.x = TRUE)
rest.cofactors <- merge(rest.cofactors, rcor3.near, by.x="Peak_REST",by.y="feature", all.x = TRUE)
rest.cofactors <- merge(rest.cofactors, sin3a.near, by.x="Peak_REST",by.y="feature", all.x = TRUE)
rest.cofactors <- merge(rest.cofactors, sin3b.near, by.x="Peak_REST",by.y="feature", all.x = TRUE)

##count number of cofactors recruited at each REST site
rest.cofactors.count <- cbind(rest.cofactors, seq(1,nrow(rest.cofactors)))
colnames(rest.cofactors.count)[21] <- "Count"

for(i in 1:nrow(rest.cofactors.count)){
rest.cofactors.count[i,"Count"] <- (ifelse(is.na(rest.cofactors[i,"Peak_RCOR1"]), 0, 1) + ifelse(is.na(rest.cofactors[i,"Peak_RCOR2"]), 0, 1) + ifelse(is.na(rest.cofactors[i,"Peak_RCOR3"]), 0, 1) + ifelse(is.na(rest.cofactors[i,"Peak_SIN3A"]), 0, 1) + ifelse(is.na(rest.cofactors[i,"Peak_SIN3B"]), 0, 1))
}

##count occurances of numbers

count <- table(rest.cofactors.count[,"Count"])

##split REST peaks into ALONE and COMPLEX

for(i in 1:nrow(rest.cofactors.count)){
	rest.cofactors.count[i,22] <- ifelse(rest.cofactors.count[i,"Count"] >= 1, "COMPLEX", "ALONE")
	}
colnames(rest.cofactors.count)[22] <- "REST"

##pull in gene expression changes

ftest <- read.csv(file = "../xpn_timecourse/results/f_test.csv")

ftest <- ftest[which(ftest[,"adj.P.Val"] <= 0.05),]
ftest <- ftest[order(abs(ftest[,"day0"]),decreasing = TRUE),]
ftest <- ftest[!duplicated(ftest[,"EnsemblID"]),]
ftest <- ftest[,c("EnsemblID","symbol","day0","day1","day2","day4","adj.P.Val")]

##merge with rest/cofactor sites

exp_rest <- merge(ftest,rest.cofactors.count, by.x = "EnsemblID", by.y = "EnsemblID_REST", all.x = TRUE)

day0 <- exp_rest[which(abs(exp_rest[,"day0"])>= 1),"EnsemblID"]
day1 <- exp_rest[which(abs(exp_rest[,"day1"])>= 1),"EnsemblID"]
day2 <- exp_rest[which(abs(exp_rest[,"day2"])>= 1),"EnsemblID"]
day4 <- exp_rest[which(abs(exp_rest[,"day4"])>= 1),"EnsemblID"]

day14 <- unique(c(day1, day2, day4))
day14 <- day14[which(!(day14 %in% day0))]

##plot REST distance vs gene expression changes
rest.ids <- rest.od[,"EnsemblID"]
exp_rest_restbound <- exp_rest[which(exp_rest[,"EnsemblID"] %in% rest.ids),]

postscript(file = "distancetoFeature_REST_vs_day0_logFC_cofactors.ps",horizontal=FALSE)

plot(exp_rest_restbound[,"distancetoFeature_REST"],exp_rest_restbound[,"day0"],pch=20,xlab="Distance from TSS",ylab="Expression Change after REST KO")

points(exp_rest_restbound[which(exp_rest_restbound[,"REST"] == "COMPLEX"),"distancetoFeature_REST"],exp_rest_restbound[which(exp_rest_restbound[,"REST"] == "COMPLEX"),"day0"],col="red")

dev.off()

## work out how many of the REST genes that change have REST/cofactors bound
day0_rest <- intersect(day0, exp_rest[which(!(is.na(exp_rest[,"Peak_REST"]))),"EnsemblID"])
day0_rest_complex <- exp_rest[which(exp_rest[,"EnsemblID"] %in% day0_rest),]
# 21 out of the 51 genes that change at day 0 are bound by REST, 16 of which have REST complex

day14_rest <- intersect(day14, exp_rest[which(!(is.na(exp_rest[,"Peak_REST"]))),"EnsemblID"])
day14_rest_complex <- exp_rest[which(exp_rest[,"EnsemblID"] %in% day14_rest),]
# 37 out of the 228 genes that change at day 1-4 are bound by REST, 29 of which have REST complexes
##why is there such a difference?? What are the genes that REST binds that do change at day0, could they have an effect at Day1-4?

day0_rest_complex[,"symbol"]
"Ckmt1"   "Car4"    "Crip2"   "Rundc3a" "Chga"    "Scg5"    "Celsr3" 
"Gdap1"   "Stmn3"   "Chrnb2"  "Syp"     "Scg3"    "Galnt9"  "Cplx1"  
"Letmd1"  "Vgf"     "Bex2"    "Tmem145" "Onecut2" "Rab39"   "Ap3b2"  

 plot(exp_rest[which(exp_rest[,"EnsemblID"] %in% rest.ids),"Count"],exp_rest[which(exp_rest[,"EnsemblID"] %in% rest.ids),"FoldEnrichment_REST"])









