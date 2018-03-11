##########################################################################################
########### ANNOTATE DIFFERENTIALLY BOUND (DB) REGIONS and GENERATE FIGS ################
##########################################################################################
qlogin -l h_vmem=10G -l h_rt=05:00:00 #or run locally
R

################# INSTALL/LOAD REQUIRED PACKAGES
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db"); biocLite("ReactomePA"); biocLite("DOSE"); biocLite("ChIPseeker")
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R") 
zzz<-c("pheatmap","grid","gplots","ggplot2","export","devtools","DESeq2","pasilla","Biobase","EBSeq","dplyr","data.table",
       "genefilter","FactoMineR","VennDiagram","DOSE","ReactomePA","org.Hs.eg.db","clusterProfiler","pathview",
       "ChIPseeker","TxDb.Hsapiens.UCSC.hg19.knownGene","GO.db")
lapply(zzz, require, character.only = TRUE)
library("DOSE")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

################# READ IN BED FILE(S) & ASSIGN PROMOTER REGIONS
file1 = "/Users/jchap12/Desktop/temp/ATACseq/DBA_ATAC_DB-peaks-UP.bed" 
file2 = "/Users/jchap12/Desktop/temp/ATACseq/DBA_ATAC_DB-peaks-DOWN.bed"
fileslist = list(up=file1, down=file2)
peak1 <- readPeakFile(fileslist[[1]])
peak2 <- readPeakFile(fileslist[[2]])
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

################# HEATMAP OF BINDING TO TSS REGIONS FOR INDIVIDUAL PEAKSETS
# tagMatrix1 <- getTagMatrix(peak1, windows=promoter)
# png('UP.DBs.heatmap.png', width = 2000, height = 4000, res = 300, pointsize = 14)
# tagHeatmap(tagMatrix1, xlim=c(-3000, 3000), color="red")
# dev.off()

# tagMatrix2 <- getTagMatrix(peak2, windows=promoter)
# png('DOWN.DBs.heatmap.png', width = 2000, height = 4000, res = 300, pointsize = 14)
# tagHeatmap(tagMatrix2, xlim=c(-3000, 3000), color="red")
# dev.off()

################# AVERAGE PROFILE OF BINDING TO TSS REGIONS FOR INDIVIDUAL PEAKSETS
# png('UP.DBs.avgProfile-TSS.png', width = 2000, height = 4000, res = 300, pointsize = 14)
# plotAvgProf(tagMatrix1, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
# dev.off()

# png('DOWN.DBs.avgProfile-TSS.png', width = 2000, height = 4000, res = 300, pointsize = 14)
# plotAvgProf(tagMatrix2, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
# dev.off()

######################### PEAK ANNOTATION FOR INDIVIDUAL PEAKSETS
peakAnno1 <- annotatePeak(fileslist[[1]], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakGRanges1 <- as.GRanges(peakAnno1)
peakDataFrame1 <- as.data.frame(peakAnno1)
write.table(peakDataFrame1, "UP.DBs.peakAnno.txt", sep="\t", row.names=FALSE)

peakAnno2 <- annotatePeak(fileslist[[2]], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakGRanges2 <- as.GRanges(peakAnno2)
peakDataFrame2 <- as.data.frame(peakAnno2)
write.table(peakDataFrame2, "DOWN.DBs.peakAnno.txt", sep="\t", row.names=FALSE)

########################### VISUALIZE GENOMIC ANNOTATION
# png('UP.DBs.Annotation-pie.png', width = 4000, height = 4000, res = 300, pointsize = 14)
# plotAnnoPie(peakAnno1) ## generates genome annotation pie chart
# dev.off()
# png('DOWN.DBs.Annotation-pie.png', width = 4000, height = 4000, res = 300, pointsize = 14)
# plotAnnoPie(peakAnno2) ## generates genome annotation pie chart
# dev.off()

# png('UP.DBs.Annotation-bar.png', width = 4000, height = 1500, res = 300, pointsize = 14)
# plotAnnoBar(peakAnno1) ## generates genome annotation bar chart
# dev.off()
# png('DOWN.DBs.Annotation-bar.png', width = 4000, height = 1500, res = 300, pointsize = 14)
# plotAnnoBar(peakAnno2) ## generates genome annotation bar chart
# dev.off()

png('UP.DBs.Annotation-vennPie.png', width = 4000, height = 4000, res = 300, pointsize = 14)
vennpie(peakAnno1) ## generates genome annotation vennpie chart
dev.off()
png('DOWN.DBs.Annotation-vennPie.png', width = 4000, height = 4000, res = 300, pointsize = 14)
vennpie(peakAnno2) ## generates genome annotation vennpie chart
dev.off()

################### VISUALIZE DISTRIBUTION OF BINDING LOCI RELATIVE TO TSS
# png('UP.DBs.distToTSS-bar.png', width = 4000, height = 1500, res = 300, pointsize = 14)
# plotDistToTSS(peakAnno1, title="Distribution of transcription factor-binding loci\nrelative to TSS")
# dev.off()
# png('DOWN.DBs.distToTSS-bar.png', width = 4000, height = 1500, res = 300, pointsize = 14)
# plotDistToTSS(peakAnno2, title="Distribution of transcription factor-binding loci\nrelative to TSS")
# dev.off()

##########################################################################################
########################### COMPARISON OF MULTIPLE PEAKSETS ##############################
tagMatrixList <- lapply(fileslist, getTagMatrix, windows=promoter)

################################### AVERAGE PROFILES
png('DBs.avg-profiles.png', width = 2000, height = 2000, res = 300, pointsize = 14)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()
################################### PEAK HEATMAPS
png('DBs.peak-heatmaps.png', width = 3000, height = 3000, res = 300, pointsize = 14)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

################################### PEAK HEATMAPS and AVERAGE PROFILES  -- long range
long_range <- getPromoters(TxDb=txdb, upstream=25000, downstream=25000)
tagMatrixList_long <- lapply(fileslist, getTagMatrix, windows=long_range)

png('DBs.peak-longrange-heatmaps.png', width = 3000, height = 3000, res = 300, pointsize = 14)
tagHeatmap(tagMatrixList_long, xlim=c(-25000, 25000), color=NULL)
dev.off()

png('DBs.longrange-avg-profiles.png', width = 2000, height = 2000, res = 300, pointsize = 14)
plotAvgProf(tagMatrixList_long, xlim=c(-25000, 25000))
dev.off()

##########################################################################################
########################### PEAK ANNOTATION COMPARISION
peakAnnoList <- lapply(fileslist, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

png('DBs.annotate-compare-bar.png', width = 3000, height = 3000, res = 300, pointsize = 14)
plotAnnoBar(peakAnnoList)
dev.off()

png('DBs.annotateTSS-compare-bar.png', width = 3000, height = 3000, res = 300, pointsize = 14)
plotDistToTSS(peakAnnoList)
dev.off()

##########################################################################################
############################ FUNCTIONAL PROFILE COMPARISONS ##############################
### TAKE ALL PEAKS FOR THIS ANALYSIS
# genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# names(genes) = sub("_", "\n", names(genes))

### TAKE A SUBSET OF ONLY PROMOTER PROXIMAL PEAKS
proximal_genes_UP = as.data.frame(peakAnnoList$up)
proximal_genes_UP2 <- proximal_genes_UP[proximal_genes_UP$distanceToTSS<=3000,]
proximal_genes_UP3 <- proximal_genes_UP2[proximal_genes_UP2$distanceToTSS>=-3000,]
# write.table(proximal_genes_UP3, "proximal_genes_UP.txt", sep="\t", row.names=FALSE)

proximal_genes_DOWN = as.data.frame(peakAnnoList$down)
proximal_genes_DOWN2 <- proximal_genes_DOWN[proximal_genes_DOWN$distanceToTSS<=3000,]
proximal_genes_DOWN3 <- proximal_genes_DOWN2[proximal_genes_DOWN2$distanceToTSS>=-3000,]
# write.table(proximal_genes_DOWN3, "proximal_genes_DOWN.txt", sep="\t", row.names=FALSE)

prox_genes_UP_vector <- as.vector(proximal_genes_UP3$geneId)
prox_genes_DOWN_vector <- as.vector(proximal_genes_DOWN3$geneId)

genes <- list(up=c(prox_genes_UP_vector), down=c(prox_genes_DOWN_vector))
names(genes) = sub("_", "\n", names(genes))

########################################## GO MF
compMF <- compareCluster(geneCluster=genes, fun="enrichGO", ont="MF",pvalueCutoff=0.05, pAdjustMethod="BH")
png('DBs.MF.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compMF, showCategory=10, title="Molecular Function Enrichment")
dev.off()

compMF_drop1 <- clusterProfiler::dropGO(compMF, level = 1, term = NULL)
png('DBs.MF_drop1.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compMF_drop1, showCategory=10, title="Molecular Function Enrichment")
dev.off()

compMF_drop2 <- clusterProfiler::dropGO(compMF_drop1, level = 2, term = NULL)
png('DBs.MF_drop2.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compMF_drop2, showCategory=10, title="Molecular Function Enrichment")
dev.off()

compMF_drop3 <- clusterProfiler::dropGO(compMF_drop2, level = 3, term = NULL)
png('DBs.MF_drop3.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compMF_drop3, showCategory=10, title="Molecular Function Enrichment")
dev.off()

compMF_drop4 <- clusterProfiler::dropGO(compMF_drop3, level = 4, term = NULL)
png('DBs.MF_drop4.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compMF_drop4, showCategory=10, title="Molecular Function Enrichment")
dev.off()

compMF_drop5 <- clusterProfiler::dropGO(compMF_drop4, level = 5, term = NULL)
png('DBs.MF_drop5.compCluster.png', width=3500, height=2000, res=300, pointsize=14)
plot(compMF_drop5, showCategory=10, title="Molecular Function Enrichment")
dev.off()

#### WRITE OUT SUMMARY TO A TAB-DELIMITED TEXT FILE
compMF_drop5_summary <- (summary(compMF_drop5))
write.table(compMF_drop5_summary, "DBs.MF_drop5.compCluster.txt", sep="\t")

########################################## GO BP
compBP <- compareCluster(geneCluster=genes, fun="enrichGO", ont="BP",pvalueCutoff=0.05, pAdjustMethod="BH")
png('DBs.BP.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compBP, showCategory=10, title="Biological Process Enrichment")
dev.off()

compBP_drop1 <- clusterProfiler::dropGO(compBP, level = 1, term = NULL)
png('DBs.BP_drop1.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compBP_drop1, showCategory=10, title="Biological Process Enrichment")
dev.off()

compBP_drop2 <- clusterProfiler::dropGO(compBP_drop1, level = 2, term = NULL)
png('DBs.BP_drop2.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compBP_drop2, showCategory=10, title="Biological Process Enrichment")
dev.off()

compBP_drop3 <- clusterProfiler::dropGO(compBP_drop2, level = 3, term = NULL)
png('DBs.BP_drop3.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compBP_drop3, showCategory=10, title="Biological Process Enrichment")
dev.off()

compBP_drop4 <- clusterProfiler::dropGO(compBP_drop3, level = 4, term = NULL)
png('DBs.BP_drop4.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(compBP_drop4, showCategory=10, title="Biological Process Enrichment")
dev.off()

compBP_drop5 <- clusterProfiler::dropGO(compBP_drop4, level = 5, term = NULL)
png('DBs.BP_drop5.compCluster.png', width=2500, height=2000, res=300, pointsize=14)
plot(compBP_drop5, showCategory=10, title="Biological Process Enrichment")
dev.off()

compBP_drop5_summary <- (summary(compBP_drop5))
write.table(compBP_drop5_summary, "DBs.BP_drop5.compCluster.txt", sep="\t")

########################################## DISEASE ONTOLOGY
cc_eDO <- clusterProfiler::compareCluster(geneCluster=genes, fun="enrichDO", pvalueCutoff= 0.1)
png('DBs.DO.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(cc_eDO, showCategory=10)
dev.off()

cc_eDO_summary <- (summary(cc_eDO))
write.table(cc_eDO_summary, "DBs.DO.compCluster.txt", sep="\t")

########################################## REACTOME
cc_reactome <- clusterProfiler::compareCluster(geneCluster=genes, fun= "ReactomePA::enrichPathway",
												organism= "human", pvalueCutoff= 0.1)
png('DBs.Reactome.compCluster.png', width=1500, height=1500, res=300, pointsize=14)
plot(cc_reactome, showCategory=10)
dev.off()

cc_reactome_summary <- (summary(cc_reactome))
write.table(cc_reactome_summary, "DBs.reactome.compCluster.txt", sep="\t", row.names=FALSE)

########################################## KEGG
cc_KEGG <- clusterProfiler::compareCluster(geneCluster=genes, fun= "enrichKEGG", 
											use_internal_data= FALSE, organism= "human",
											pvalueCutoff= 0.1)
png('DBs.KEGG.compCluster.png', width=2000, height=2000, res=300, pointsize=14)
plot(cc_KEGG, showCategory=10)
dev.off()

cc_KEGG_summary <- (summary(cc_KEGG))
write.table(cc_KEGG_summary, "DBs.KEGG.compCluster.txt", sep="\t", row.names=FALSE)

##########################################################################################