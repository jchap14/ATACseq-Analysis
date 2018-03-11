########### DETERMINE DIFFERENTIALLY BOUND (DB) REGIONS and GENERATE FIGS ################
qlogin -l h_vmem=10G -l h_rt=05:00:00
# change back to working directory
module load r
R

##### install Diffbind
# source("http://bioconductor.org/biocLite.R")
# biocLite("DiffBind")
# # browseVignettes("DiffBind")

##### set experiment title
Title <- "siKLF_ETS"

################### READ IN PEAKSETS with *.csv containing metadata #####
library(DiffBind)
samples   <- read.csv(paste(Title, ".metadata.csv", sep=''))
DiffBind_Obj <- dba(sampleSheet=paste(Title, ".metadata.csv", sep=''))

##### look at some peak correlations prior to signal consideration

## generate a correlation heatmap (initial clustering of peaks)
png(paste(Title, ".occupancy.correlation.heatmap.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
plot(DiffBind_Obj)
dev.off()

## generate plot to determine how many peaks overlap between samples
olap.rate <- dba.overlap(DiffBind_Obj,mode=DBA_OLAP_RATE)
png(paste(Title, ".Peak-overlap-rate.plot.png", sep=''), type="cairo",
 width = 2000, height = 2000, res = 300, pointsize = 14)
plot(olap.rate,type='b',ylab='# peaks',xlab='Overlap at least this many peaksets')
dev.off()

## generate Venn-diagrams of overlapping peaks within $conditions, $replicates, etc.
png(paste(Title, ".overlapping-peaks-STATIC.Venn.png", sep=''), type="cairo",
 width = 2000, height = 2000, res = 300, pointsize = 14)
dba.plotVenn(DiffBind_Obj, DiffBind_Obj$masks$Static)
dev.off()

png(paste(Title, ".overlapping-peaks-FLOW.Venn.png", sep=''), type="cairo",
 width = 2000, height = 2000, res = 300, pointsize = 14)
dba.plotVenn(DiffBind_Obj, DiffBind_Obj$masks$Flow)
dev.off()

##########################################################################################
################# COUNT READS UNDER MERGED (or CONSENSUS) PEAKSET ########################
DiffBind_Obj <- dba.count(DiffBind_Obj, minOverlap=1,bParallel=TRUE)

## generate a correlation heatmap (affinity clustering takes reads into account)
png(paste(Title, ".affinity.correlation.heatmap.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
plot(DiffBind_Obj)
dev.off()

##########################################################################################
################# DIFFERENTIAL BINDING ANALYSIS (NO BLOCKING) ############################
## ESTABLISH A CONTRAST (what to compare for DB analysis)
DiffBind_Obj <- dba.contrast(DiffBind_Obj, categories=DBA_CONDITION)

## PERFORM DB ANALYSIS
DiffBind_Obj <- dba.analyze(DiffBind_Obj)

## GENERATE A CORRELATION HEATMAP (affinity clustering using only DB sites)
png(paste(Title, ".diffBound.correlation.heatmap.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
plot(DiffBind_Obj, contrast=1)
dev.off()

## RETRIEVE DB SITES
DiffBind_Obj.DB <- dba.report(DiffBind_Obj)
# write.table(DiffBind_Obj.DB, "DiffBind_Obj.DB_sites.txt", sep="\t") ## doesnt work yet

## GENERATE MA PLOT: RED = DB SITE
png(paste(Title, ".diffBound.MAplot.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
dba.plotMA(DiffBind_Obj,bXY=FALSE)
dev.off()

## GENERATE PCA PLOT USING AFFINITY DATA AT ALL SITES
png(paste(Title, ".allsites.PCA.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
dba.plotPCA(DiffBind_Obj,attributes=DBA_CONDITION,label=DBA_ID)
dev.off()

## GENERATE PCA PLOT USING AFFINITY DATA AT DB SITES
png(paste(Title, ".diffBound.PCA.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300)
dba.plotPCA(DiffBind_Obj,attributes=DBA_CONDITION,contrast=1,th=.1,label=DBA_ID)
dev.off()

## GENERATE BOXPLOT: shows distribution of reads over all conditions
png(paste(Title, ".read_distribution.boxplot.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
dba.plotBox(DiffBind_Obj)
dev.off()

## GENERATE HEATMAP OF DB SITES
png(paste(Title, ".diffBound.heatmap.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
dba.plotHeatmap(DiffBind_Obj, contrast=1, correlations=FALSE,scale="row")
dev.off()

##########################################################################################
####################### DB ANALYSIS USING A BLOCKING FACTOR ##############################
## ESTABLISH A CONTRAST (what to compare for DB analysis, specify blocking factor)
DiffBind_Obj <- dba.contrast(DiffBind_Obj,categories=DBA_CONDITION, block=DBA_TISSUE)

## PERFORM DB ANALYSIS
DiffBind_Obj <- dba.analyze(DiffBind_Obj)

## GENERATE MA PLOT: RED = DB SITE
png(paste(Title, ".diffBound.blocked.MAplot.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
dba.plotMA(DiffBind_Obj,method=DBA_EDGER_BLOCK)
dev.off()

## GENERATE A CORRELATION HEATMAP (affinity clustering using only DB sites)
png(paste(Title, ".diffBound.blocked.correlation.heatmap.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
dba.plotHeatmap(DiffBind_Obj,contrast=1,method=DBA_EDGER_BLOCK,attributes=c(DBA_TISSUE,DBA_CONDITION,DBA_REPLICATE))
dev.off()

## GENERATE A HEATMAP OF DB SITES
png(paste(Title, ".diffBound.blocked.heatmap.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
dba.plotHeatmap(DiffBind_Obj,contrast=1,method=DBA_EDGER_BLOCK, correlations=FALSE,scale="row")
dev.off()

## GENERATE PCA PLOT USING AFFINITY DATA AT DB SITES
png(paste(Title, ".diffBound.blocked.PCA.png", sep=''), type="cairo",
 width = 4000, height = 2000, res = 300, pointsize = 14)
dba.plotPCA(DiffBind_Obj,contrast=1,method=DBA_EDGER_BLOCK,attributes=DBA_CONDITION,label=DBA_ID)
dev.off()

## Export DB SITES for local loading
DiffBind_Obj.DB <- dba.report(DiffBind_Obj,method=DBA_EDGER_BLOCK,contrast=1,th=.1,
bNormalized=TRUE, file="ATAC_DB-peaks")

##########################################################################################
########################### DB ANALYSIS METHODS COMPARISON ###############################

## compare the sites identified using edgeR and DESeq2 (go back and do DESeq2 also)
treat.block <- dba.report(DiffBind_Obj,method=DBA_ALL_METHODS_BLOCK,bDB=TRUE,bAll=TRUE)
png(paste(Title, ".diffBound.DB-method.Venn.png", sep=''), type="cairo",
 width = 4000, height = 4000, res = 300, pointsize = 14)
dba.plotVenn(treat.block,1:2,label1="edgeR", label2="edgeR Blocked")
dev.off()

##########################################################################################
q(save = "no")

## convert .csv peak report to .bedgraph (really .bed) for IGV loading
cat DBA_ATAC_DB-peaks.csv | tr ',' '\t' | cut -f1,2,3,7 | tr -d '"' | tail -n+2> a-temp.bed
cat a-temp.bed | sort -k1,1 -k2,2n > c.temp.txt
cat > b.temp.txt << EOF
track name="Flow_vs_Static_Peaks" color=255,0,0 altColor=0,0,255
EOF
cat b.temp.txt c.temp.txt > DBA_ATAC_DB-peaks.bedgraph; rm a-temp.bed b.temp.txt c.temp.txt

track name="Flow_vs_Static_Peaks" color=255,0,0 altColor=0,0,255
###################### GO TO CHIPSEEKER PACKAGE TO ANNOTATE PEAKS ########################
