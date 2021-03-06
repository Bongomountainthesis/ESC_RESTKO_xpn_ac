options(stringsAsFactors=FALSE)

rest_h3k9ac <- read.csv("rest_h3k9ac_nearest_peak_to_gene_TSS.csv")
rest_h4ac <- read.csv("rest_h4ac_nearest_peak_to_gene_TSS.csv")
ctrl_h3k9ac <- read.csv("ctrl_h3k9ac_nearest_peak_to_gene_TSS.csv")
ctrl_h4ac <- read.csv("ctrl_h4ac_nearest_peak_to_gene_TSS.csv")

###apply cut offs from neg10log10pval > 100

ctrl_h3k9ac <- ctrl_h3k9ac[which(ctrl_h3k9ac[,"neg10log10pVal"]>=100),]
rest_h3k9ac <- rest_h3k9ac[which(rest_h3k9ac[,"neg10log10pVal"]>=100),]
ctrl_h4ac <- ctrl_h4ac[which(ctrl_h4ac[,"neg10log10pVal"]>=100),]
rest_h4ac <- rest_h4ac[which(rest_h4ac[,"neg10log10pVal"]>=100),]

##get IDs
r_k9 <- rest_h3k9ac[,"EnsemblID"]
c_k9 <- ctrl_h3k9ac[,"EnsemblID"]

r_h4 <- rest_h4ac[,"EnsemblID"]
c_h4 <- ctrl_h4ac[,"EnsemblID"]

#shared?

shared_k9 <- intersect(r_k9,c_k9)
shared_h4 <- intersect(r_h4,c_h4)

shared_r_k9 <- rest_h3k9ac[which(rest_h3k9ac[,"EnsemblID"] %in% c_k9),]
shared_c_k9 <- ctrl_h3k9ac[which(rest_h3k9ac[,"EnsemblID"] %in% r_k9),]

shared_r_h4 <- rest_h4ac[which(rest_h4ac[,"EnsemblID"] %in% c_h4),]
shared_c_h4 <- ctrl_h4ac[which(rest_h4ac[,"EnsemblID"] %in% r_h4),]

#unique?
unique_r_k9 <- rest_h3k9ac[which(!(rest_h3k9ac[,"EnsemblID"] %in% c_k9)),]
unique_c_k9 <- ctrl_h3k9ac[which(!(ctrl_h3k9ac[,"EnsemblID"] %in% r_k9)),]

unique_r_h4 <- rest_h4ac[which(!(rest_h4ac[,"EnsemblID"] %in% c_h4)),]
unique_c_h4 <- ctrl_h4ac[which(!(ctrl_h4ac[,"EnsemblID"] %in% r_h4)),]

shared.df <- merge(ctrl_h3k9ac,rest_h3k9ac, by.x = "EnsemblID", by.y = "EnsemblID", suffixes = c("_ctrl","restko"))
write.csv(shared.df,file = "H3K9ac_peaks_shared_between_ctrl_and_restko_ESCs.csv")

#how many of these have REST binding sites??

rest <- read.csv(file="../../REST_ChIP/results/REST_D0_nearest_peak_to_gene_TSS.csv")
rest <- rest[which(rest[,"neg10log10pVal"] >=100),]

rest.ids <- rest[,"EnsemblID"]

rest_shared_k9 <- intersect(rest.ids,shared_k9)
rest_shared_h4 <- intersect(rest.ids,shared_h4)

rest_unique_r_k9 <- intersect(rest.ids,unique_r_k9[,"EnsemblID"])
rest_unique_c_k9 <- intersect(rest.ids,unique_c_k9[,"EnsemblID"])

rest_unique_r_h4 <- intersect(rest.ids,unique_r_h4[,"EnsemblID"])
rest_unique_c_h4 <- intersect(rest.ids,unique_c_h4[,"EnsemblID"])





unique_binding_restko <- merge(unique_restko, rest_sites, by.x = "EnsemblID", by.y = "EnsemblID", all.x = TRUE, suffixes = c("_H3K9ac","_REST"))
unique_binding_cnt <- merge(unique_cnt, rest_sites, by.x = "EnsemblID", by.y = "EnsemblID", all.x = TRUE, suffixes = c("_H3K9ac","_REST"))

tidy_restko <- unique_binding_restko[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,35,36,37,38)]
tidy_cnt <- unique_binding_cnt[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,35,36,37,38)]

write.csv(tidy_restko, file = "unique_RESTKO_H3K9ac_and_nearest_REST_site.csv")
write.csv(tidy_cnt, file = "unique_control_H3K9ac_and_nearest_REST_site.csv")

length(which(!is.na(unique_binding_restko[,"EnsemblID"])))
length(which(!is.na(unique_binding_cnt[,"EnsemblID"])))

#969 K9Ac peaks have a REST binding site in rest ko

unique_binding_REST <- unique_binding_restko[which(!is.na(unique_binding_restko[,"Peak_REST"])),]
unique_binding_cnt <- unique_binding_cnt[which(!is.na(unique_binding_cnt[,"Peak_REST"])),]

dim(unique_binding_REST)
#441
dim(unique_binding_cnt)
#23

foo_rest <- unique_binding_REST[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,35,36,37,38)]
foo_cnt <- unique_binding_cnt[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,35,36,37,38)]

write.csv(foo_rest, file = "unique_RESTKO_H3K9ac_with_nearest_REST_site.csv")
write.csv(foo_cnt, file = "unique_control_H3K9ac_with_nearest_REST_site.csv")


###########same again for H4ac

# rest k4ac = 14883
# ctrl k4ac = 14187

r_h4 <- rest_h4ac[,"EnsemblID"]
c_h4 <- ctrl_h4ac[,"EnsemblID"]

shared_h4 <- intersect(r_h4, c_h4)

shared_h4.df <- merge(ctrl_h4ac,rest_h4ac, by.x = "EnsemblID", by.y = "EnsemblID", suffixes = c("_ctrl","restko"))
write.csv(shared_h4.df,file = "H4ac_peaks_shared_between_ctrl_and_restko_ESCs.csv")


# shared_h4 = 11134

unique_restko <- rest_h4ac[which(!(rest_h4ac[,"EnsemblID"] %in% c_h4)),]

# new k4ac after rest ko = 3749

#how many of these have REST binding sites??

rest_sites <- read.csv(file="../results/REST_binding_sites.csv")

unique_binding <- merge(unique_restko, rest_sites, by.x = "EnsemblID", by.y = "EnsemblID", all.x = TRUE, suffixes = c("_H4ac","_REST"))

fooba <- unique_binding[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,35,36,37,38)]

write.csv(fooba, file = "unique_RESTKO_H4ac_and_nearest_REST_site.csv")

length(which(!is.na(unique_binding[,"EnsemblID"])))

#497 new K9Ac peaks have a REST binding site

unique_binding_REST <- unique_binding[which(!is.na(unique_binding[,"Peak_REST"])),]

foo <- unique_binding_REST[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,35,36,37,38)]

write.csv(foo, file = "unique_RESTKO_H4ac_with_nearest_REST_site.csv")

##plot width versus size w/wo REST
#remove syt4 as its a weird peak
ctrl_h3k9ac <- ctrl_h3k9ac[which(!(ctrl_h3k9ac[,"Symbol"] == "Syt4")),]
rest_h3k9ac <- rest_h3k9ac[which(!(rest_h3k9ac[,"Symbol"] == "Syt4")),]

par(mfrow=c(2,2))
plot(ctrl_h3k9ac[,"Peak_width"],ctrl_h3k9ac[,"neg10log10pVal"],pch=20,xlim = c(0,40000),ylim=c(0,2750),xlab="Peak width",ylab="neg10log10pVal")
plot(rest_h3k9ac[,"Peak_width"],rest_h3k9ac[,"neg10log10pVal"],pch=20,col = "red",xlim = c(0,40000),ylim=c(0,2750),xlab="Peak width",ylab="neg10log10pVal")

plot(ctrl_h4ac[,"Peak_width"],ctrl_h4ac[,"neg10log10pVal"],pch=20,xlim = c(0,80000),ylim=c(0,1750),xlab="Peak width",ylab="neg10log10pVal")
plot(rest_h4ac[,"Peak_width"],rest_h4ac[,"neg10log10pVal"],pch=20,col = "red",xlim = c(0,80000),ylim=c(0,1750),xlab="Peak width",ylab="neg10log10pVal")

## overlay points for those bound by REST and/or cofactors

rest <- read.csv(file="../../REST_ChIP/results/REST_D0_nearest_peak_to_gene_TSS.csv")
rest <- rest[which(rest[,"neg10log10pVal"] >=100),]

rest.ids <- rest[,"EnsemblID"]

par(mfrow=c(2,2))
plot(ctrl_h3k9ac[,"Peak_width"],ctrl_h3k9ac[,"neg10log10pVal"],pch=20,xlim = c(0,40000),ylim=c(0,2750),xlab="Peak width",ylab="neg10log10pVal")
points(ctrl_h3k9ac[which(ctrl_h3k9ac[,"EnsemblID"] %in% rest.ids),"Peak_width"],ctrl_h3k9ac[which(ctrl_h3k9ac[,"EnsemblID"] %in% rest.ids),"neg10log10pVal"],col="blue")

plot(rest_h3k9ac[,"Peak_width"],rest_h3k9ac[,"neg10log10pVal"],pch=20,col = "red",xlim = c(0,40000),ylim=c(0,2750),xlab="Peak width",ylab="neg10log10pVal")
points(rest_h3k9ac[which(rest_h3k9ac[,"EnsemblID"] %in% rest.ids),"Peak_width"],rest_h3k9ac[which(rest_h3k9ac[,"EnsemblID"] %in% rest.ids),"neg10log10pVal"],col="blue")

plot(ctrl_h4ac[,"Peak_width"],ctrl_h4ac[,"neg10log10pVal"],pch=20,xlim = c(0,80000),ylim=c(0,1750),xlab="Peak width",ylab="neg10log10pVal")
points(ctrl_h4ac[which(ctrl_h4ac[,"EnsemblID"] %in% rest.ids),"Peak_width"],ctrl_h4ac[which(ctrl_h4ac[,"EnsemblID"] %in% rest.ids),"neg10log10pVal"],col="blue")

plot(rest_h4ac[,"Peak_width"],rest_h4ac[,"neg10log10pVal"],pch=20,col = "red",xlim = c(0,80000),ylim=c(0,1750),xlab="Peak width",ylab="neg10log10pVal")
points(rest_h4ac[which(rest_h4ac[,"EnsemblID"] %in% rest.ids),"Peak_width"],rest_h4ac[which(rest_h4ac[,"EnsemblID"] %in% rest.ids),"neg10log10pVal"],col="blue")













