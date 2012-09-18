options(stringsAsFactors=FALSE)

library(IRanges)
library(ShortRead)
library(Rsamtools)
library(ChIPpeakAnno)
library(RColorBrewer)
#library(GenomeGraphs)
library(biomaRt)
library(gplots)

##pull in BAM files

d0_rest <- "/mnt/data/REST_CtrlES_D0.clean_sorted_nodups.bam"
d4_rest <- "/mnt/data/REST_CtrlES_D4.clean_sorted_nodups.bam"
rest_k9ac <- "/mnt/data/CME143_GA3R71_export_sorted_nodups.bam"
ctrl_k9ac <- "/mnt/data/CME141_GA3R71_export_sorted_nodups.bam"
rest_h4ac <- "/mnt/data/FloResCre__C18__H4ac_CME117_s_2_export_sorted_nodups.bam"
ctrl_h4ac <- "/mnt/data/FloRes_H4acs_CME118_3_export_sorted_nodups.bam"

d0_rest <- "REST_ChIP/results/alignment/bowtie/REST_CtrlES_D0.clean_sorted_nodups.bam"
d4_rest <- "REST_ChIP/results/alignment/bowtie/REST_CtrlES_D4.clean_sorted_nodups.bam"
rest_k9ac <- "Ac_ChIP/results/alignment/bowtie/CME143_GA3R71_export_sorted_nodups.bam"
ctrl_k9ac <- "Ac_ChIP/results/alignment/bowtie/CME141_GA3R71_export_sorted_nodups.bam"
rest_h4ac <- "Ac_ChIP/results/alignment/bowtie/FloResCre__C18__H4ac_CME117_s_2_export_sorted_nodups.bam"
ctrl_h4ac <- "Ac_ChIP/results/alignment/bowtie/FloRes_H4acs_CME118_3_export_sorted_nodups.bam"

##take REST and K9/K4Ac read files and count overlaps in 200bp bins. Test for correlation between presence of both peaks.

chr_lengths <- as.data.frame(matrix(nrow = 19, ncol = 2))

colnames(chr_lengths) <- c("Chr", "Length")

chr_lengths[1,2] <- as.numeric(197195432)
chr_lengths[2,2] <- as.numeric(181748087)
chr_lengths[3,2] <- as.numeric(159599783)
chr_lengths[4,2] <- as.numeric(155630120)
chr_lengths[5,2] <- as.numeric(152537259)
chr_lengths[6,2] <- as.numeric(152537259)
chr_lengths[7,2] <- as.numeric(152524553)
chr_lengths[8,2] <- as.numeric(131738871)
chr_lengths[9,2] <- as.numeric(124076172)
chr_lengths[10,2] <- as.numeric(129993255)
chr_lengths[11,2] <- as.numeric(121843856)
chr_lengths[12,2] <- as.numeric(121257530)
chr_lengths[13,2] <- as.numeric(120284312)
chr_lengths[14,2] <- as.numeric(125194864)
chr_lengths[15,2] <- as.numeric(103494974)
chr_lengths[16,2] <- as.numeric(98319150)
chr_lengths[17,2] <- as.numeric(95272651)
chr_lengths[18,2] <- as.numeric(90772031)
chr_lengths[19,2] <- as.numeric(61342430)

chr_lengths[,1] <- paste("chr",seq(from = 1, to = 19),sep="")

bin.size <- 500

########## REST DAY 0 - count overlaps

bam_file <- d0_rest

######Chr1

start <- seq(from=0, to = chr_lengths[1,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[1,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[1,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr1_bam <- countBam(bam_file, param=param)

######Chr2

start <- seq(from=0, to = chr_lengths[2,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[2,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[2,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr2_bam <- countBam(bam_file, param=param)

######Chr3

start <- seq(from=0, to = chr_lengths[3,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[3,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[3,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr3_bam <- countBam(bam_file, param=param)

######Chr4

start <- seq(from=0, to = chr_lengths[4,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[4,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[4,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr4_bam <- countBam(bam_file, param=param)

######Chr5

start <- seq(from=0, to = chr_lengths[5,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[5,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[5,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr5_bam <- countBam(bam_file, param=param)

######Chr6

start <- seq(from=0, to = chr_lengths[6,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[6,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[6,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr6_bam <- countBam(bam_file, param=param)

######Chr7

start <- seq(from=0, to = chr_lengths[7,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[7,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[7,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr7_bam <- countBam(bam_file, param=param)

######Chr8

start <- seq(from=0, to = chr_lengths[8,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[8,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[8,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr8_bam <- countBam(bam_file, param=param)

######Chr9

start <- seq(from=0, to = chr_lengths[9,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[9,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[9,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr9_bam <- countBam(bam_file, param=param)

######Chr10

start <- seq(from=0, to = chr_lengths[10,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[10,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[10,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr10_bam <- countBam(bam_file, param=param)

######Chr11

start <- seq(from=0, to = chr_lengths[11,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[11,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[11,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr11_bam <- countBam(bam_file, param=param)

######Chr12

start <- seq(from=0, to = chr_lengths[12,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[12,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[12,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr12_bam <- countBam(bam_file, param=param)

######Chr13

start <- seq(from=0, to = chr_lengths[13,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[13,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[13,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr13_bam <- countBam(bam_file, param=param)

######Chr14

start <- seq(from=0, to = chr_lengths[14,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[14,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[14,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr14_bam <- countBam(bam_file, param=param)

######Chr15

start <- seq(from=0, to = chr_lengths[15,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[15,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[15,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr15_bam <- countBam(bam_file, param=param)

######Chr16

start <- seq(from=0, to = chr_lengths[16,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[16,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[16,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr16_bam <- countBam(bam_file, param=param)

######Chr17

start <- seq(from=0, to = chr_lengths[17,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[17,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[17,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr17_bam <- countBam(bam_file, param=param)

######Chr18

start <- seq(from=0, to = chr_lengths[18,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[18,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[18,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr18_bam <- countBam(bam_file, param=param)

######Chr19

start <- seq(from=0, to = chr_lengths[19,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[19,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[19,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr19_bam <- countBam(bam_file, param=param)

######join all together and save

rest_d0_counts <- c(chr1_bam[,"records"],chr2_bam[,"records"],chr3_bam[,"records"],chr4_bam[,"records"],chr5_bam[,"records"],chr6_bam[,"records"],chr7_bam[,"records"],chr8_bam[,"records"],chr9_bam[,"records"],chr10_bam[,"records"],chr11_bam[,"records"],chr12_bam[,"records"],chr13_bam[,"records"],chr14_bam[,"records"],chr15_bam[,"records"],chr16_bam[,"records"],chr17_bam[,"records"],chr18_bam[,"records"],chr19_bam[,"records"])

save(rest_d0_counts, file = "comparison_results/rest_d0_tag_counts.RData")



########## REST DAY 4 - count overlaps

bam_file <- d4_rest

######Chr1

start <- seq(from=0, to = chr_lengths[1,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[1,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[1,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr1_bam <- countBam(bam_file, param=param)

######Chr2

start <- seq(from=0, to = chr_lengths[2,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[2,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[2,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr2_bam <- countBam(bam_file, param=param)

######Chr3

start <- seq(from=0, to = chr_lengths[3,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[3,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[3,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr3_bam <- countBam(bam_file, param=param)

######Chr4

start <- seq(from=0, to = chr_lengths[4,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[4,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[4,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr4_bam <- countBam(bam_file, param=param)

######Chr5

start <- seq(from=0, to = chr_lengths[5,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[5,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[5,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr5_bam <- countBam(bam_file, param=param)

######Chr6

start <- seq(from=0, to = chr_lengths[6,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[6,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[6,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr6_bam <- countBam(bam_file, param=param)

######Chr7

start <- seq(from=0, to = chr_lengths[7,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[7,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[7,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr7_bam <- countBam(bam_file, param=param)

######Chr8

start <- seq(from=0, to = chr_lengths[8,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[8,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[8,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr8_bam <- countBam(bam_file, param=param)

######Chr9

start <- seq(from=0, to = chr_lengths[9,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[9,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[9,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr9_bam <- countBam(bam_file, param=param)

######Chr10

start <- seq(from=0, to = chr_lengths[10,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[10,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[10,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr10_bam <- countBam(bam_file, param=param)

######Chr11

start <- seq(from=0, to = chr_lengths[11,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[11,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[11,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr11_bam <- countBam(bam_file, param=param)

######Chr12

start <- seq(from=0, to = chr_lengths[12,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[12,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[12,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr12_bam <- countBam(bam_file, param=param)

######Chr13

start <- seq(from=0, to = chr_lengths[13,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[13,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[13,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr13_bam <- countBam(bam_file, param=param)

######Chr14

start <- seq(from=0, to = chr_lengths[14,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[14,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[14,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr14_bam <- countBam(bam_file, param=param)

######Chr15

start <- seq(from=0, to = chr_lengths[15,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[15,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[15,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr15_bam <- countBam(bam_file, param=param)

######Chr16

start <- seq(from=0, to = chr_lengths[16,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[16,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[16,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr16_bam <- countBam(bam_file, param=param)

######Chr17

start <- seq(from=0, to = chr_lengths[17,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[17,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[17,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr17_bam <- countBam(bam_file, param=param)

######Chr18

start <- seq(from=0, to = chr_lengths[18,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[18,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[18,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr18_bam <- countBam(bam_file, param=param)

######Chr19

start <- seq(from=0, to = chr_lengths[19,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[19,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[19,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr19_bam <- countBam(bam_file, param=param)

######join all together and save

rest_d4_counts <- c(chr1_bam[,"records"],chr2_bam[,"records"],chr3_bam[,"records"],chr4_bam[,"records"],chr5_bam[,"records"],chr6_bam[,"records"],chr7_bam[,"records"],chr8_bam[,"records"],chr9_bam[,"records"],chr10_bam[,"records"],chr11_bam[,"records"],chr12_bam[,"records"],chr13_bam[,"records"],chr14_bam[,"records"],chr15_bam[,"records"],chr16_bam[,"records"],chr17_bam[,"records"],chr18_bam[,"records"],chr19_bam[,"records"])

save(rest_d4_counts, file = "/mnt/data/rest_d4_tag_counts.RData")

########## H3K9Ac (Cntl) - count overlaps

bam_file <- ctrl_k9ac

######Chr1

start <- seq(from=0, to = chr_lengths[1,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[1,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[1,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr1_bam <- countBam(bam_file, param=param)

######Chr2

start <- seq(from=0, to = chr_lengths[2,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[2,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[2,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr2_bam <- countBam(bam_file, param=param)

######Chr3

start <- seq(from=0, to = chr_lengths[3,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[3,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[3,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr3_bam <- countBam(bam_file, param=param)

######Chr4

start <- seq(from=0, to = chr_lengths[4,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[4,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[4,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr4_bam <- countBam(bam_file, param=param)

######Chr5

start <- seq(from=0, to = chr_lengths[5,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[5,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[5,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr5_bam <- countBam(bam_file, param=param)

######Chr6

start <- seq(from=0, to = chr_lengths[6,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[6,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[6,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr6_bam <- countBam(bam_file, param=param)

######Chr7

start <- seq(from=0, to = chr_lengths[7,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[7,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[7,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr7_bam <- countBam(bam_file, param=param)

######Chr8

start <- seq(from=0, to = chr_lengths[8,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[8,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[8,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr8_bam <- countBam(bam_file, param=param)

######Chr9

start <- seq(from=0, to = chr_lengths[9,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[9,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[9,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr9_bam <- countBam(bam_file, param=param)

######Chr10

start <- seq(from=0, to = chr_lengths[10,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[10,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[10,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr10_bam <- countBam(bam_file, param=param)

######Chr11

start <- seq(from=0, to = chr_lengths[11,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[11,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[11,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr11_bam <- countBam(bam_file, param=param)

######Chr12

start <- seq(from=0, to = chr_lengths[12,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[12,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[12,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr12_bam <- countBam(bam_file, param=param)

######Chr13

start <- seq(from=0, to = chr_lengths[13,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[13,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[13,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr13_bam <- countBam(bam_file, param=param)

######Chr14

start <- seq(from=0, to = chr_lengths[14,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[14,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[14,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr14_bam <- countBam(bam_file, param=param)

######Chr15

start <- seq(from=0, to = chr_lengths[15,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[15,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[15,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr15_bam <- countBam(bam_file, param=param)

######Chr16

start <- seq(from=0, to = chr_lengths[16,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[16,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[16,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr16_bam <- countBam(bam_file, param=param)

######Chr17

start <- seq(from=0, to = chr_lengths[17,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[17,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[17,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr17_bam <- countBam(bam_file, param=param)

######Chr18

start <- seq(from=0, to = chr_lengths[18,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[18,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[18,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr18_bam <- countBam(bam_file, param=param)

######Chr19

start <- seq(from=0, to = chr_lengths[19,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[19,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[19,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr19_bam <- countBam(bam_file, param=param)

######join all together and save

ctrl_k9ac_counts <- c(chr1_bam[,"records"],chr2_bam[,"records"],chr3_bam[,"records"],chr4_bam[,"records"],chr5_bam[,"records"],chr6_bam[,"records"],chr7_bam[,"records"],chr8_bam[,"records"],chr9_bam[,"records"],chr10_bam[,"records"],chr11_bam[,"records"],chr12_bam[,"records"],chr13_bam[,"records"],chr14_bam[,"records"],chr15_bam[,"records"],chr16_bam[,"records"],chr17_bam[,"records"],chr18_bam[,"records"],chr19_bam[,"records"])

save(ctrl_k9ac_counts, file = "/mnt/data/ctrl_k9ac_tag_counts.RData")


########## ctrl_h4ac - count overlaps

bam_file <- ctrl_h4ac

######Chr1

start <- seq(from=0, to = chr_lengths[1,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[1,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[1,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr1_bam <- countBam(bam_file, param=param)

######Chr2

start <- seq(from=0, to = chr_lengths[2,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[2,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[2,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr2_bam <- countBam(bam_file, param=param)

######Chr3

start <- seq(from=0, to = chr_lengths[3,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[3,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[3,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr3_bam <- countBam(bam_file, param=param)

######Chr4

start <- seq(from=0, to = chr_lengths[4,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[4,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[4,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr4_bam <- countBam(bam_file, param=param)

######Chr5

start <- seq(from=0, to = chr_lengths[5,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[5,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[5,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr5_bam <- countBam(bam_file, param=param)

######Chr6

start <- seq(from=0, to = chr_lengths[6,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[6,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[6,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr6_bam <- countBam(bam_file, param=param)

######Chr7

start <- seq(from=0, to = chr_lengths[7,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[7,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[7,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr7_bam <- countBam(bam_file, param=param)

######Chr8

start <- seq(from=0, to = chr_lengths[8,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[8,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[8,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr8_bam <- countBam(bam_file, param=param)

######Chr9

start <- seq(from=0, to = chr_lengths[9,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[9,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[9,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr9_bam <- countBam(bam_file, param=param)

######Chr10

start <- seq(from=0, to = chr_lengths[10,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[10,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[10,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr10_bam <- countBam(bam_file, param=param)

######Chr11

start <- seq(from=0, to = chr_lengths[11,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[11,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[11,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr11_bam <- countBam(bam_file, param=param)

######Chr12

start <- seq(from=0, to = chr_lengths[12,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[12,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[12,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr12_bam <- countBam(bam_file, param=param)

######Chr13

start <- seq(from=0, to = chr_lengths[13,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[13,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[13,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr13_bam <- countBam(bam_file, param=param)

######Chr14

start <- seq(from=0, to = chr_lengths[14,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[14,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[14,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr14_bam <- countBam(bam_file, param=param)

######Chr15

start <- seq(from=0, to = chr_lengths[15,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[15,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[15,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr15_bam <- countBam(bam_file, param=param)

######Chr16

start <- seq(from=0, to = chr_lengths[16,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[16,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[16,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr16_bam <- countBam(bam_file, param=param)

######Chr17

start <- seq(from=0, to = chr_lengths[17,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[17,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[17,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr17_bam <- countBam(bam_file, param=param)

######Chr18

start <- seq(from=0, to = chr_lengths[18,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[18,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[18,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr18_bam <- countBam(bam_file, param=param)

######Chr19

start <- seq(from=0, to = chr_lengths[19,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[19,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[19,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr19_bam <- countBam(bam_file, param=param)

######join all together and save

ctrl_h4ac_counts <- c(chr1_bam[,"records"],chr2_bam[,"records"],chr3_bam[,"records"],chr4_bam[,"records"],chr5_bam[,"records"],chr6_bam[,"records"],chr7_bam[,"records"],chr8_bam[,"records"],chr9_bam[,"records"],chr10_bam[,"records"],chr11_bam[,"records"],chr12_bam[,"records"],chr13_bam[,"records"],chr14_bam[,"records"],chr15_bam[,"records"],chr16_bam[,"records"],chr17_bam[,"records"],chr18_bam[,"records"],chr19_bam[,"records"])

save(ctrl_h4ac_counts, file = "/mnt/data/ctrl_h4ac_tag_counts.RData")

########## H3K9Ac (rest) - count overlaps

bam_file <- rest_k9ac

######Chr1

start <- seq(from=0, to = chr_lengths[1,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[1,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[1,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr1_bam <- countBam(bam_file, param=param)

######Chr2

start <- seq(from=0, to = chr_lengths[2,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[2,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[2,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr2_bam <- countBam(bam_file, param=param)

######Chr3

start <- seq(from=0, to = chr_lengths[3,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[3,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[3,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr3_bam <- countBam(bam_file, param=param)

######Chr4

start <- seq(from=0, to = chr_lengths[4,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[4,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[4,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr4_bam <- countBam(bam_file, param=param)

######Chr5

start <- seq(from=0, to = chr_lengths[5,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[5,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[5,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr5_bam <- countBam(bam_file, param=param)

######Chr6

start <- seq(from=0, to = chr_lengths[6,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[6,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[6,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr6_bam <- countBam(bam_file, param=param)

######Chr7

start <- seq(from=0, to = chr_lengths[7,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[7,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[7,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr7_bam <- countBam(bam_file, param=param)

######Chr8

start <- seq(from=0, to = chr_lengths[8,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[8,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[8,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr8_bam <- countBam(bam_file, param=param)

######Chr9

start <- seq(from=0, to = chr_lengths[9,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[9,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[9,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr9_bam <- countBam(bam_file, param=param)

######Chr10

start <- seq(from=0, to = chr_lengths[10,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[10,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[10,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr10_bam <- countBam(bam_file, param=param)

######Chr11

start <- seq(from=0, to = chr_lengths[11,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[11,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[11,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr11_bam <- countBam(bam_file, param=param)

######Chr12

start <- seq(from=0, to = chr_lengths[12,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[12,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[12,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr12_bam <- countBam(bam_file, param=param)

######Chr13

start <- seq(from=0, to = chr_lengths[13,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[13,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[13,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr13_bam <- countBam(bam_file, param=param)

######Chr14

start <- seq(from=0, to = chr_lengths[14,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[14,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[14,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr14_bam <- countBam(bam_file, param=param)

######Chr15

start <- seq(from=0, to = chr_lengths[15,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[15,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[15,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr15_bam <- countBam(bam_file, param=param)

######Chr16

start <- seq(from=0, to = chr_lengths[16,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[16,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[16,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr16_bam <- countBam(bam_file, param=param)

######Chr17

start <- seq(from=0, to = chr_lengths[17,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[17,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[17,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr17_bam <- countBam(bam_file, param=param)

######Chr18

start <- seq(from=0, to = chr_lengths[18,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[18,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[18,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr18_bam <- countBam(bam_file, param=param)

######Chr19

start <- seq(from=0, to = chr_lengths[19,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[19,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[19,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr19_bam <- countBam(bam_file, param=param)

######join all together and save

rest_k9ac_counts <- c(chr1_bam[,"records"],chr2_bam[,"records"],chr3_bam[,"records"],chr4_bam[,"records"],chr5_bam[,"records"],chr6_bam[,"records"],chr7_bam[,"records"],chr8_bam[,"records"],chr9_bam[,"records"],chr10_bam[,"records"],chr11_bam[,"records"],chr12_bam[,"records"],chr13_bam[,"records"],chr14_bam[,"records"],chr15_bam[,"records"],chr16_bam[,"records"],chr17_bam[,"records"],chr18_bam[,"records"],chr19_bam[,"records"])

save(rest_k9ac_counts, file = "/mnt/data/rest_k9ac_tag_counts.RData")


########## rest_h4ac - count overlaps

bam_file <- rest_h4ac

######Chr1

start <- seq(from=0, to = chr_lengths[1,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[1,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[1,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr1_bam <- countBam(bam_file, param=param)

######Chr2

start <- seq(from=0, to = chr_lengths[2,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[2,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[2,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr2_bam <- countBam(bam_file, param=param)

######Chr3

start <- seq(from=0, to = chr_lengths[3,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[3,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[3,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr3_bam <- countBam(bam_file, param=param)

######Chr4

start <- seq(from=0, to = chr_lengths[4,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[4,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[4,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr4_bam <- countBam(bam_file, param=param)

######Chr5

start <- seq(from=0, to = chr_lengths[5,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[5,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[5,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr5_bam <- countBam(bam_file, param=param)

######Chr6

start <- seq(from=0, to = chr_lengths[6,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[6,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[6,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr6_bam <- countBam(bam_file, param=param)

######Chr7

start <- seq(from=0, to = chr_lengths[7,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[7,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[7,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr7_bam <- countBam(bam_file, param=param)

######Chr8

start <- seq(from=0, to = chr_lengths[8,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[8,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[8,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr8_bam <- countBam(bam_file, param=param)

######Chr9

start <- seq(from=0, to = chr_lengths[9,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[9,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[9,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr9_bam <- countBam(bam_file, param=param)

######Chr10

start <- seq(from=0, to = chr_lengths[10,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[10,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[10,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr10_bam <- countBam(bam_file, param=param)

######Chr11

start <- seq(from=0, to = chr_lengths[11,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[11,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[11,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr11_bam <- countBam(bam_file, param=param)

######Chr12

start <- seq(from=0, to = chr_lengths[12,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[12,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[12,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr12_bam <- countBam(bam_file, param=param)

######Chr13

start <- seq(from=0, to = chr_lengths[13,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[13,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[13,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr13_bam <- countBam(bam_file, param=param)

######Chr14

start <- seq(from=0, to = chr_lengths[14,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[14,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[14,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr14_bam <- countBam(bam_file, param=param)

######Chr15

start <- seq(from=0, to = chr_lengths[15,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[15,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[15,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr15_bam <- countBam(bam_file, param=param)

######Chr16

start <- seq(from=0, to = chr_lengths[16,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[16,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[16,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr16_bam <- countBam(bam_file, param=param)

######Chr17

start <- seq(from=0, to = chr_lengths[17,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[17,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[17,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr17_bam <- countBam(bam_file, param=param)

######Chr18

start <- seq(from=0, to = chr_lengths[18,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[18,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[18,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr18_bam <- countBam(bam_file, param=param)

######Chr19

start <- seq(from=0, to = chr_lengths[19,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[19,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[19,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr19_bam <- countBam(bam_file, param=param)

######join all together and save

rest_h4ac_counts <- c(chr1_bam[,"records"],chr2_bam[,"records"],chr3_bam[,"records"],chr4_bam[,"records"],chr5_bam[,"records"],chr6_bam[,"records"],chr7_bam[,"records"],chr8_bam[,"records"],chr9_bam[,"records"],chr10_bam[,"records"],chr11_bam[,"records"],chr12_bam[,"records"],chr13_bam[,"records"],chr14_bam[,"records"],chr15_bam[,"records"],chr16_bam[,"records"],chr17_bam[,"records"],chr18_bam[,"records"],chr19_bam[,"records"])

save(rest_h4ac_counts, file = "/mnt/data/rest_h4ac_tag_counts.RData")



###############################################

######would be good to train this on H3K4me3 modification to see if see any correlation

#take from MLA NS K4me3

bam_file <- "/space/MLA2_MLA2dNeuron_ChIP/results/alignment/bowtie/MLA_NS_H3K4me3_CMN054_s_2_export_sorted_nodups.bam"

######Chr1

start <- seq(from=0, to = chr_lengths[1,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[1,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[1,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr1_bam <- countBam(bam_file, param=param)

######Chr2

start <- seq(from=0, to = chr_lengths[2,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[2,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[2,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr2_bam <- countBam(bam_file, param=param)

######Chr3

start <- seq(from=0, to = chr_lengths[3,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[3,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[3,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr3_bam <- countBam(bam_file, param=param)

######Chr4

start <- seq(from=0, to = chr_lengths[4,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[4,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[4,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr4_bam <- countBam(bam_file, param=param)

######Chr5

start <- seq(from=0, to = chr_lengths[5,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[5,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[5,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr5_bam <- countBam(bam_file, param=param)

######Chr6

start <- seq(from=0, to = chr_lengths[6,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[6,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[6,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr6_bam <- countBam(bam_file, param=param)

######Chr7

start <- seq(from=0, to = chr_lengths[7,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[7,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[7,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr7_bam <- countBam(bam_file, param=param)

######Chr8

start <- seq(from=0, to = chr_lengths[8,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[8,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[8,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr8_bam <- countBam(bam_file, param=param)

######Chr9

start <- seq(from=0, to = chr_lengths[9,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[9,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[9,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr9_bam <- countBam(bam_file, param=param)

######Chr10

start <- seq(from=0, to = chr_lengths[10,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[10,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[10,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr10_bam <- countBam(bam_file, param=param)

######Chr11

start <- seq(from=0, to = chr_lengths[11,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[11,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[11,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr11_bam <- countBam(bam_file, param=param)

######Chr12

start <- seq(from=0, to = chr_lengths[12,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[12,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[12,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr12_bam <- countBam(bam_file, param=param)

######Chr13

start <- seq(from=0, to = chr_lengths[13,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[13,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[13,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr13_bam <- countBam(bam_file, param=param)

######Chr14

start <- seq(from=0, to = chr_lengths[14,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[14,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[14,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr14_bam <- countBam(bam_file, param=param)

######Chr15

start <- seq(from=0, to = chr_lengths[15,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[15,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[15,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr15_bam <- countBam(bam_file, param=param)

######Chr16

start <- seq(from=0, to = chr_lengths[16,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[16,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[16,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr16_bam <- countBam(bam_file, param=param)

######Chr17

start <- seq(from=0, to = chr_lengths[17,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[17,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[17,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr17_bam <- countBam(bam_file, param=param)

######Chr18

start <- seq(from=0, to = chr_lengths[18,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[18,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[18,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr18_bam <- countBam(bam_file, param=param)

######Chr19

start <- seq(from=0, to = chr_lengths[19,2] - bin.size, by=bin.size)
end   <- seq(from=0 + bin.size, to = chr_lengths[19,2], by=bin.size)

chr_bins <- RangedData(ranges = IRanges(
                   start= start,
                   end = end
                   ),
		   space = chr_lengths[19,1]
                 )

## get counts from BAM file
param <- ScanBamParam(which=chr_bins)

chr19_bam <- countBam(bam_file, param=param)

######join all together and save

mla_k4me3_counts <- c(chr1_bam[,"records"],chr2_bam[,"records"],chr3_bam[,"records"],chr4_bam[,"records"],chr5_bam[,"records"],chr6_bam[,"records"],chr7_bam[,"records"],chr8_bam[,"records"],chr9_bam[,"records"],chr10_bam[,"records"],chr11_bam[,"records"],chr12_bam[,"records"],chr13_bam[,"records"],chr14_bam[,"records"],chr15_bam[,"records"],chr16_bam[,"records"],chr17_bam[,"records"],chr18_bam[,"records"],chr19_bam[,"records"])

save(mla_k4me3_counts, file = "/mnt/data/mla_k4me3_tag_counts.RData")


#########################

#reload files in


rest_d0_counts <- get(load("comparison_results/rest_d0_tag_counts.RData"))
rest_d4_counts <- get(load("comparison_results/rest_d4_tag_counts.RData"))

ctrl_k9ac_counts <- get(load("comparison_results/ctrl_k9ac_tag_counts.RData"))
rest_k9ac_counts <- get(load("comparison_results/rest_k9ac_tag_counts.RData"))

ctrl_h4ac_counts <- get(load("comparison_results/ctrl_h4ac_tag_counts.RData"))
rest_h4ac_counts <- get(load("comparison_results/rest_h4ac_tag_counts.RData"))

mla_k4me3_counts <- get(load("comparison_results/mla_k4me3_tag_counts.RData"))

###need to take into account peaks that aren't really peaks - filter on minimum number of reads in peaks?

rest_d0_peaks <- read.csv(file = "REST_ChIP/results/REST_D0_nearest_peak_to_gene_TSS.csv")
rest_d4_peaks <- read.csv(file = "REST_ChIP/results/REST_D4_nearest_peak_to_gene_TSS.csv")

cntl_k9ac_peaks <- read.csv(file = "Ac_ChIP/results/new_annotation/ctrl_h3k9ac/ctrl_h3k4ac_nearest_or_overlapping_peak_to_gene.csv")
rest_k9ac_peaks <- read.csv(file = "Ac_ChIP/results/new_annotation/rest_h3k9ac/rest_h3k9ac_nearest_or_overlapping_peak_to_gene.csv")

cntl_h4ac_peaks <- read.csv(file = "Ac_ChIP/results/new_annotation/ctrl_h4ac/ctrl_h4ac_nearest_or_overlapping_peak_to_gene.csv")
rest_h4ac_peaks <- read.csv(file = "Ac_ChIP/results/new_annotation/rest_h4ac/rest_h4ac_nearest_or_overlapping_peak_to_gene.csv")


##normalise to read counts

rest_d0_depth <- as.numeric(18488744)
rest_d4_depth <- as.numeric(17510073)

ctrl_k9ac_depth <- as.numeric(13960865)
rest_k9ac_depth <- as.numeric(19301764)

ctrl_h4ac_depth <- as.numeric(19561815)
rest_h4ac_depth <- as.numeric(17039692)

mla_k4me3_depth <- as.numeric(18769907)

total_depth <- rest_d0_depth + rest_d4_depth + ctrl_k9ac_depth + rest_k9ac_depth + ctrl_h4ac_depth + rest_h4ac_depth + mla_k4me3_depth

rest_d0_depth <- rest_d0_depth/total_depth
rest_d4_depth <- rest_d4_depth/total_depth

ctrl_k9ac_depth <- ctrl_k9ac_depth/total_depth
rest_k9ac_depth <- rest_k9ac_depth/total_depth

ctrl_h4ac_depth <- ctrl_h4ac_depth/total_depth
rest_h4ac_depth <- rest_h4ac_depth/total_depth

mla_k4me3_depth <- mla_k4me3_depth/total_depth

## normalise

rest_d0_counts <- rest_d0_counts / rest_d0_depth
rest_d4_counts <- rest_d4_counts / rest_d4_depth

ctrl_k9ac_counts <- ctrl_k9ac_counts / ctrl_k9ac_depth
rest_k9ac_counts <- rest_k9ac_counts / rest_k9ac_depth

ctrl_h4ac_counts <- ctrl_h4ac_counts / ctrl_h4ac_depth
rest_h4ac_counts <- rest_h4ac_counts / rest_h4ac_depth

mla_k4me3_counts <- mla_k4me3_counts / mla_k4me3_depth


##plot together

#plot(rest_d0_counts, rest_d4_counts)

#plot(rest_d0_counts, ctrl_k9ac_counts)
#plot(rest_d0_counts, ctrl_h4ac_counts)

##vs mla k4

#plot(rest_d0_counts, mla_k4me3_counts)

#################

##make a plot clustered together to show overlap of REST at D0 and K9ac with/without REST

#make matrix to cluster

peak_matrix <- cbind(rest_d0_counts, ctrl_k9ac_counts, rest_k9ac_counts)

#remove lines that are 0 everywhere
peak_matrix_sum <- matrix(nrow = nrow(peak_matrix), ncol = 1)

for(i in 1:nrow(peak_matrix)){
	peak_matrix_sum[i] <- sum(peak_matrix[i,])
	}
	
peak_matrix_counts <- peak_matrix[which(peak_matrix_sum >3),]

#scale peaks
#peak_matrix_scale <- scale(peak_matrix_counts)


##save and export to cluster 3 as R can't process this many lines...
name <- seq(1,nrow(peak_matrix_counts))

peak_matrix_name <- cbind(name,peak_matrix_counts)

write.table(peak_matrix_name, file = "comparison_results/tag_count_matrix.txt", sep = "\t",row.names = F, col.names = F)


##command for Cluster 3 (to be run from command line)
# cluster -f comparison_results/tag_count_matrix.txt -ng -ca a -g 2

## make differential peak




heatmap.2(peak_matrix_counts)

plot(peak_matrix_counts[seq(1,nrow(peak_matrix_counts)),1])








