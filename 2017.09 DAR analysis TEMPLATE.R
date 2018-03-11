##########################################################################################
########### DETERMINE DIFFERENTIALLY ACCESSIBLE (DA) REGIONS and GENERATE FIGS ###########
##########################################################################################

############################# Step 1: Always do this on the cluster ####
screen -S R_session
qlogin -l h_vmem=10G -l h_rt=24:00:00
module add r
cd /srv/gsfs0/projects/snyder/chappell/JR/SLO_newATAC/siKLF_ETS_analysis #or CWD
R
library("DiffBind")

##### READ IN PEAKSETS with *.csv containing metadata
## set experiment title
Title <- "siKLF_ETS"
metadata <- read.csv(paste(Title, ".metadata.csv", sep='')) #metadata has to match exact format as example

# have to set minOverlap to <2 (default) if wanting to test a lesser overlap
# this may cause a memory failure if qlogin not working well
treatment <- dba(sampleSheet=metadata, minOverlap=3)

## generate a correlation heatmap (initial clustering of peaks)
pdf(paste(Title,".peaks.correlation.heatmap.pdf", sep=''), width= 8, height= 8, pointsize= 14)
plot(treatment, main= "occupancy.correlation.heatmap")
dev.off()

##### generate plot to determine how many peaks overlap between samples
olap.rate <- dba.overlap(treatment,mode=DBA_OLAP_RATE)
pdf(paste(Title,".Peak-overlap-rate.plot.pdf", sep=''), width= 5, height= 5, pointsize= 14)
plot(olap.rate,type='b',ylab='# peaks',xlab='Overlap at least this many peaksets')
dev.off()

# ##### get a dataframe of the full peakset used for DA testing
# full_overlapped_peaks<- dba.peakset(treatment, peaks= NULL, bRetrieve=T,
#                                     minOverlap=1, DataType= DBA_DATA_FRAME)
## export for annotation
# write.table(full_overlapped_peaks, "SLO.full_peakset", sep="\t", row.names= F)

# export pre-counts DBA object for local manipulation. This may cause issues locally.
# dba.save(treatment, file='SLO_ATAC', dir='.', pre='dba_', ext='RData', bMinimize=FALSE)

##### COUNT READS UNDER MERGED (or CONSENSUS) PEAKSET
treatment <- dba.count(treatment, minOverlap=3, bParallel=T)
## get the full matrix of counts & export them for local manipulation
counts.matrix <- dba.peakset(treatment, bRetrieve=T, DataType=DBA_DATA_FRAME)
write.table(counts.matrix, paste(Title,".counts.matrix",sep=''), sep="\t", row.names= F)

## generate a correlation heatmap (affinity clustering takes reads into account)
pdf(paste(Title,".affinity.correlation.heatmap.pdf", sep=''), width= 8, height= 8, pointsize = 14)
plot(treatment)
dev.off()

############### DIFFERENTIAL accessibility analysis locally ####
## set experiment title
Title <- "siKLF_ETS"

## read in metadata
metadata <- read.delim(paste(Title, ".metadata.csv", sep=''), quote="\"'", sep = ",")

##### Read in the count matrix
y <- read.delim(paste(Title, ".counts.matrix", sep=''), quote="\"'")
## ID columns containing count info
CountCols <- as.character(metadata$SampleID)
## create matrix from count containing columns
y.matrix <- as.matrix(y[,CountCols])
class(y.matrix) <- "integer"

# ##### construct a DESeq2 DataSet w/ count matrix for ATAC without blocking batch
require("DESeq2")
dds <- DESeqDataSetFromMatrix(countData= y.matrix, colData= metadata,
                              design= ~Condition) # ~batch + Condition (if desired)
dds <- dds[ rowSums(counts(dds)) > 1, ] ## pre-filter rows that have only 0 or 1 read
dds <- DESeq(dds) #run differential test
##### perform rlog transform of counts matrix to adjust for PCA/heatmap
rld <- rlog(dds)
log2.rlog.counts <- assay(rld); x <- log2.rlog.counts

## If desired: remove batch effect, merge corrected rLog counts w/ gene names for later incorporation
# x <- removeBatchEffect(x, batch= as.character(metadata$Treatment)) #batch by date of prep
# x.symbols <- cbind(y[,c(1:3)], x)

##### PCA: PC2 vs PC1 ####
## assign a numeric matrix to mtrx
mtrx <- x

## Calculations
require("genefilter")
rv <- rowVars(mtrx)
## select # of genes to consider for PCA (change max 2 min to specify top 500)
select <- order(rv, decreasing= T)[seq_len(min(10000, length(rv)))]
pca <- prcomp(t(mtrx[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], metadata)

## specify plot aesthetics
require("ggplot2")
color     <- "Condition"
label     <- d$SampleID
shape     <- "Treatment"
mainTitle <- "PCA"
textSize  <- element_text(size= 14)

## Plot
a <- ggplot(d, aes_string(x= "PC1", y= "PC2", color=color, shape=shape, label=label)) +
  geom_point(size= 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() +
  geom_text(aes(label=label),hjust=0, vjust=0, size= 5) + ggtitle(mainTitle) +
  theme(plot.title= element_text(size= 14, face= "bold"), axis.text= textSize,
        legend.text= textSize, legend.title= textSize, axis.title= textSize, strip.text= textSize)
plot(a)

## Export to powerpoint
require("export")
graph2ppt(file=paste(Title, ".PCA.pptx",sep=''), width=10, height=7, append=T)

##### Transcriptome PCA: PC3 vs PC2 ####
## Plot
a <- ggplot(d, aes_string(x= "PC3", y= "PC4", color=color, shape=shape, label=label)) +
  geom_point(size= 3) + xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() +
  geom_text(aes(label=label),hjust=0, vjust=0, size= 5) + ggtitle(mainTitle) +
  theme(plot.title= element_text(size= 14, face= "bold"), axis.text= textSize,
        legend.text= textSize, legend.title= textSize, axis.title=textSize, strip.text=textSize)
plot(a)

## Export to powerpoint
graph2ppt(file=paste(Title, ".PCA.pptx",sep=''), width=10, height=7, append=T)

##### Heatmap of top variable genes ####
x <- x[select, ]         #set x to the top variable genes from above
## cluster rows by pearson, complete
hr <- hclust(as.dist(1-cor(t(x), method='pearson')), method='complete')
hc <- hclust(as.dist(1-cor(x, method='pearson')), method='complete')
## Heatmap2 w/ color bar. h & k modify cutree (k overrides h). Specify margins here as well.
mycl <- cutree(hr, h=max(hr$height)/3, k = 5);
mycol <- sample(rainbow(256)); mycol <- mycol[as.vector(mycl)]
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)
## generate heatmap with rows clustered, but not columns
require(gplots)
png('10ktopVarPeaks.heatmap_r.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=NA, dendrogram= c("row"), col=my_palette, scale="row",
          density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## generate a heatmap with rows & columns clustered
png('10ktopVarPeaks.heatmap_rc.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram= c("both"), col=my_palette,
          scale="row", density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## plot & export column dendrograms
require(dendextend)
hc.dend <- hc %>% as.dendrogram #convert to dendrogram
plot(hc.dend, main = "Column Dendrogram")
graph2ppt(file="10ktopVarPeaks.dendrograms.pptx", width=10, height=7, append=T)
## plot & export row dendrograms
hr.dend <- hr %>% as.dendrogram #convert to dendrogram
hr.dend.cut <- cut(hr.dend, h= 1.9) #cut at desired height
plot(hr.dend.cut$upper, horiz = T, main = "Row Dendrogram")
graph2ppt(file="10ktopVarPeaks.dendrograms.pptx", width=10, height=7, append=T)

#################################################################################################
##### DA Testing: LS_vs_ST ####
res      <- results(dds, contrast=c("Condition","LS","ST"))
summary(res)
resMF_df <- as.data.frame(res)
## merge significant results w/ gene names & rlog counts
resMF_counts   <- cbind(resMF_df, x.symbols)
res.all        <- resMF_counts[order(resMF_counts$padj),] #order the results by the smallest adj pval
DARs.LS_vs_ST.all_res <- res.all

##### CHIPSEEKER PACKAGE TO ANNOTATE PEAKS
FDR <- 0.01 #set desired FDR for remaining analysis here
DARs.LS_vs_ST <- subset(res.all, padj < FDR) 

##### convert to gRanges objects
DARs.LS_vs_ST.gRanges <- makeGRangesFromDataFrame(DARs.LS_vs_ST,keep.extra.columns=T)
promoter              <- getPromoters(TxDb=txdb, upstream=2500, downstream=500)

##### ANNOTATE DARs, subset BEDs for Motif search, & export
DARs.LS_vs_ST.anno    <- annotatePeak(DARs.LS_vs_ST.gRanges, tssRegion=c(-2500, 500), TxDb=txdb, annoDb="org.Hs.eg.db")
DARs.LS_vs_ST.anno.df <- as.data.frame(DARs.LS_vs_ST.anno)
## calculate average read counts for LS & ST
DARs.LS_vs_ST.anno.df$LS_avg  <- rowMeans(DARs.LS_vs_ST.anno.df[,c(12:14)])
DARs.LS_vs_ST.anno.df$ST_avg  <- rowMeans(DARs.LS_vs_ST.anno.df[,c(17:20)])
## reorder columns for easy interpretation
DARs.LS_vs_ST.anno.df <- DARs.LS_vs_ST.anno.df[,c(1:4,7,11,31,29,33:34,12:13,15:16,19,14,17,18,20,26,22:28,30,32)]
## export annotated DARs
write.table(DARs.LS_vs_ST.anno.df, "DARs.LS_vs_ST.FDR1.txt", sep="\t", row.names=F)

##### Subset DARs by LS- or ST- enrichment & make BEDs for Motif search
## LS-enriched DARs
DARs.LS_vs_ST.LS_enriched     <- subset(DARs.LS_vs_ST.anno.df, log2FoldChange > 0)
DARs.LS_vs_ST.LS_enriched.bed <- DARs.LS_vs_ST.LS_enriched[,c(1:3)]
## add a name_distToTSS column to bed for ID after motif finding
DARs.LS_vs_ST.LS_enriched.bed$namedist <- paste(DARs.LS_vs_ST.LS_enriched$SYMBOL,
                                                DARs.LS_vs_ST.LS_enriched$distanceToTSS, sep='_')
write.table(DARs.LS_vs_ST.LS_enriched.bed, "DARs.LS_vs_ST.LS_enriched.bed", sep="\t", row.names=F, col.names=F, quote=F)

## ST-enriched DARs
DARs.LS_vs_ST.ST_enriched     <- subset(DARs.LS_vs_ST.anno.df, log2FoldChange < -0)
DARs.LS_vs_ST.ST_enriched.bed <-  DARs.LS_vs_ST.ST_enriched[,c(1:3)]
## add a name_distToTSS column to bed for ID after motif finding
DARs.LS_vs_ST.ST_enriched.bed$namedist <- paste(DARs.LS_vs_ST.ST_enriched$SYMBOL,
                                                DARs.LS_vs_ST.ST_enriched$distanceToTSS, sep='_')
write.table(DARs.LS_vs_ST.ST_enriched.bed, "DARs.LS_vs_ST.ST_enriched.bed", sep="\t", row.names=F, col.names=F, quote=F)

##### Subset DARs by Proximal or Distal (>3kb) & make BEDs for Motif search
## subset TSS_proximal (within -/+ 3kb) DARs
DARs.LS_vs_ST.LS_enriched.TSSprox <- subset(DARs.LS_vs_ST.LS_enriched, abs(distanceToTSS) < 3000)
DARs.LS_vs_ST.LS_enriched.TSSdist <- subset(DARs.LS_vs_ST.LS_enriched, abs(distanceToTSS) > 3000)
DARs.LS_vs_ST.ST_enriched.TSSdist <- subset(DARs.LS_vs_ST.ST_enriched, abs(distanceToTSS) > 3000)
DARs.LS_vs_ST.ST_enriched.TSSprox <- subset(DARs.LS_vs_ST.ST_enriched, abs(distanceToTSS) < 3000)

##### VISUALIZE GENOMIC ANNOTATION
## create LS- & ST- annotated objects
DARs.LS_vs_ST.LS.gRanges <- makeGRangesFromDataFrame(DARs.LS_vs_ST.LS_enriched.bed)
DARs.LS_vs_ST.ST.gRanges <- makeGRangesFromDataFrame(DARs.LS_vs_ST.ST_enriched.bed)
DARs.LS_vs_ST.LS.anno    <- annotatePeak(DARs.LS_vs_ST.LS.gRanges, tssRegion=c(-2500, 500),
                                         TxDb=txdb, annoDb="org.Hs.eg.db")
DARs.LS_vs_ST.ST.anno    <- annotatePeak(DARs.LS_vs_ST.ST.gRanges, tssRegion=c(-2500, 500),
                                         TxDb=txdb, annoDb="org.Hs.eg.db")

## make a VennPie graphic
vennpie(DARs.LS_vs_ST.LS.anno) ## generates genome annotation vennpie chart
graph2ppt(file="DARs.LS_vs_ST.LS_genome_anno.ppt", width=7, height=7, append=T); dev.off()
vennpie(DARs.LS_vs_ST.ST.anno) ## generates genome annotation vennpie chart
graph2ppt(file="DARs.LS_vs_ST.ST_genome_anno.ppt", width=7, height=7, append=T); dev.off()

##### COMPARISON OF MULTIPLE PEAKSETS
fileslist     <- list(DARs.LS_vs_ST.LS=DARs.LS_vs_ST.LS.gRanges,
                      DARs.LS_vs_ST.ST=DARs.LS_vs_ST.ST.gRanges)
promoter      <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(fileslist, getTagMatrix, windows=promoter)

##### AVERAGE PROFILES
pdf('DARs.LS_vs_ST.avg_binding.pdf', width= 10, height= 6, pointsize= 14)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()

##### PEAK HEATMAPS
png('DARs.LS_vs_ST.binding_heatmaps.png', width= 6, height= 10, units= "in", res= 300, pointsize= 14)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

##### PEAK ANNOTATION COMPARISION
peakAnnoList <- lapply(fileslist, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=F)
#
plotAnnoBar(peakAnnoList)
graph2ppt(file="DARs.LS_vs_ST.annotate_compare_bar.ppt", width=10, height=6, append=T)
dev.off()
#
plotDistToTSS(peakAnnoList)
graph2ppt(file="DARs.LS_vs_ST.annotate_compare_bar.ppt", width=10, height=6, append=T)
dev.off()

##### DAR Heatmap
x <- as.matrix(DARs.LS_vs_ST.anno.df[,c(11:19)])
## cluster rows by pearson, complete
hr <- hclust(as.dist(1-cor(t(x), method='pearson')), method='complete')
hc <- hclust(as.dist(1-cor(x, method='pearson')), method='complete')
## Heatmap2 w/ color bar. h & k modify cutree (k overrides h). Specify margins here as well.
mycl       <- cutree(hr, h=max(hr$height)/3, k = 5);
mycol      <- sample(rainbow(256)); mycol <- mycol[as.vector(mycl)]
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)
## generate heatmap with rows clustered, but not columns
png('DARs.LS_vs_ST.heatmap_r.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=NA, dendrogram= c("row"), col=my_palette, scale="row",
          density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## generate a heatmap with rows & columns clustered
png('DARs.LS_vs_ST.heatmap_rc.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram= c("both"), col=my_palette,
          scale="row", density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## plot & export column dendrograms
hc.dend <- hc %>% as.dendrogram #convert to dendrogram
plot(hc.dend, main = "Column Dendrogram")
graph2ppt(file="DARs.LS_vs_ST.dendrograms.pptx", width=10, height=7, append=T)
## plot & export row dendrograms
hr.dend <- hr %>% as.dendrogram #convert to dendrogram
hr.dend.cut <- cut(hr.dend, h= 1.9) #cut at desired height
plot(hr.dend.cut$upper, horiz = T, main = "Row Dendrogram")
graph2ppt(file="DARs.LS_vs_ST.dendrograms.pptx", width=10, height=7, append=T)

##### DA Testing: LS_vs_OS ####
res      <- results(dds, contrast=c("Condition","LS","OS"))
summary(res)
resMF_df <- as.data.frame(res)
## merge significant results w/ gene names & rlog counts
resMF_counts   <- cbind(resMF_df, x.symbols)
res.all        <- resMF_counts[order(resMF_counts$padj),] #order the results by the smallest adj pval
DARs.LS_vs_OS.all_res <- res.all

##### CHIPSEEKER PACKAGE TO ANNOTATE PEAKS
FDR <- 0.01 #set desired FDR for remaining analysis here
DARs.LS_vs_OS <- subset(res.all, padj < FDR) 

##### convert to gRanges objects
DARs.LS_vs_OS.gRanges <- makeGRangesFromDataFrame(DARs.LS_vs_OS,keep.extra.columns=T)
promoter              <- getPromoters(TxDb=txdb, upstream=2500, downstream=500)

##### ANNOTATE DARs, subset BEDs for Motif search, & export
DARs.LS_vs_OS.anno    <- annotatePeak(DARs.LS_vs_OS.gRanges, tssRegion=c(-2500, 500), TxDb=txdb, annoDb="org.Hs.eg.db")
DARs.LS_vs_OS.anno.df <- as.data.frame(DARs.LS_vs_OS.anno)
# calculate average read counts for LS & OS
DARs.LS_vs_OS.anno.df$LS_avg  <- rowMeans(DARs.LS_vs_OS.anno.df[,c(12:14)])
DARs.LS_vs_OS.anno.df$OS_avg  <- rowMeans(DARs.LS_vs_OS.anno.df[,c(15:16)])
#reorder columns for easy interpretation
DARs.LS_vs_OS.anno.df <- DARs.LS_vs_OS.anno.df[,c(1:4,7,11,31,29,33:34,12:13,15:16,19,14,17,18,20,26,22:28,30,32)]
#export annotated DARs
write.table(DARs.LS_vs_OS.anno.df, "DARs.LS_vs_OS.FDR1.txt", sep="\t", row.names=F)

##### Subset DARs by LS- or OS- enrichment & make BEDs for Motif search
## LS-enriched DARs
DARs.LS_vs_OS.LS_enriched     <- subset(DARs.LS_vs_OS.anno.df, log2FoldChange > 0)
DARs.LS_vs_OS.LS_enriched.bed <- DARs.LS_vs_OS.LS_enriched[,c(1:3)]
## add a name_distToTSS column to bed for ID after motif finding
DARs.LS_vs_OS.LS_enriched.bed$namedist <- paste(DARs.LS_vs_OS.LS_enriched$SYMBOL,
                                                DARs.LS_vs_OS.LS_enriched$distanceToTSS, sep='_')
write.table(DARs.LS_vs_OS.LS_enriched.bed, "DARs.LS_vs_OS.LS_enriched.bed", sep="\t", row.names=F, col.names=F, quote=F)

## OS-enriched DARs
DARs.LS_vs_OS.OS_enriched     <- subset(DARs.LS_vs_OS.anno.df, log2FoldChange < -0)
DARs.LS_vs_OS.OS_enriched.bed <-  DARs.LS_vs_OS.OS_enriched[,c(1:3)]
## add a name_distToTSS column to bed for ID after motif finding
DARs.LS_vs_OS.OS_enriched.bed$namedist <- paste(DARs.LS_vs_OS.OS_enriched$SYMBOL,
                                                DARs.LS_vs_OS.OS_enriched$distanceToTSS, sep='_')
write.table(DARs.LS_vs_OS.OS_enriched.bed, "DARs.LS_vs_OS.OS_enriched.bed", sep="\t", row.names=F, col.names=F, quote=F)

##### Subset DARs by Proximal or Distal (>3kb) & make BEDs for Motif search
## subset TSS_proximal (within -/+ 3kb) DARs
DARs.LS_vs_OS.LS_enriched.TSSprox <- subset(DARs.LS_vs_OS.LS_enriched, abs(distanceToTSS) < 3000)
DARs.LS_vs_OS.LS_enriched.TSSdist <- subset(DARs.LS_vs_OS.LS_enriched, abs(distanceToTSS) > 3000)
DARs.LS_vs_OS.OS_enriched.TSSdist <- subset(DARs.LS_vs_OS.OS_enriched, abs(distanceToTSS) > 3000)
DARs.LS_vs_OS.OS_enriched.TSSprox <- subset(DARs.LS_vs_OS.OS_enriched, abs(distanceToTSS) < 3000)

##### VISUALIZE GENOMIC ANNOTATION
## create LS- & OS- annotated objects
DARs.LS_vs_OS.LS.gRanges <- makeGRangesFromDataFrame(DARs.LS_vs_OS.LS_enriched.bed)
DARs.LS_vs_OS.OS.gRanges <- makeGRangesFromDataFrame(DARs.LS_vs_OS.OS_enriched.bed)
DARs.LS_vs_OS.LS.anno    <- annotatePeak(DARs.LS_vs_OS.LS.gRanges, tssRegion=c(-2500, 500),
                                         TxDb=txdb, annoDb="org.Hs.eg.db")
DARs.LS_vs_OS.OS.anno    <- annotatePeak(DARs.LS_vs_OS.OS.gRanges, tssRegion=c(-2500, 500),
                                         TxDb=txdb, annoDb="org.Hs.eg.db")

## make a VennPie graphic
vennpie(DARs.LS_vs_OS.LS.anno) ## generates genome annotation vennpie chart
graph2ppt(file="DARs.LS_vs_OS.LS_genome_anno.ppt", width=7, height=7, append=T); dev.off()
vennpie(DARs.LS_vs_OS.OS.anno) ## generates genome annotation vennpie chart
graph2ppt(file="DARs.LS_vs_OS.OS_genome_anno.ppt", width=7, height=7, append=T); dev.off()

##### COMPARISON OF MULTIPLE PEAKSETS
fileslist     <- list(DARs.LS_vs_OS.LS=DARs.LS_vs_OS.LS.gRanges,
                      DARs.LS_vs_OS.OS=DARs.LS_vs_OS.OS.gRanges)
promoter      <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(fileslist, getTagMatrix, windows=promoter)

##### AVERAGE PROFILES
pdf('DARs.LS_vs_OS.avg_binding.pdf', width= 10, height= 6, pointsize= 14)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()

##### PEAK HEATMAPS
png('DARs.LS_vs_OS.binding_heatmaps.png', width= 6, height= 10, units= "in", res= 300, pointsize= 14)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

##### PEAK ANNOTATION COMPARISION
peakAnnoList <- lapply(fileslist, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=F)
#
plotAnnoBar(peakAnnoList)
graph2ppt(file="DARs.LS_vs_OS.annotate_compare_bar.ppt", width=10, height=6, append=T)
dev.off()
#
plotDistToTSS(peakAnnoList)
graph2ppt(file="DARs.LS_vs_OS.annotate_compare_bar.ppt", width=10, height=6, append=T)
dev.off()

##### DAR Heatmap
x <- as.matrix(DARs.LS_vs_OS.anno.df[,c(11:19)])
## cluster rows by pearson, complete
hr <- hclust(as.dist(1-cor(t(x), method='pearson')), method='complete')
hc <- hclust(as.dist(1-cor(x, method='pearson')), method='complete')
## Heatmap2 w/ color bar. h & k modify cutree (k overrides h). Specify margins here as well.
mycl       <- cutree(hr, h=max(hr$height)/3, k = 5);
mycol      <- sample(rainbow(256)); mycol <- mycol[as.vector(mycl)]
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)
## generate heatmap with rows clustered, but not columns
png('DARs.LS_vs_OS.heatmap_r.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=NA, dendrogram= c("row"), col=my_palette, scale="row",
          density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## generate a heatmap with rows & columns clustered
png('DARs.LS_vs_OS.heatmap_rc.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram= c("both"), col=my_palette,
          scale="row", density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## plot & export column dendrograms
hc.dend <- hc %>% as.dendrogram #convert to dendrogram
plot(hc.dend, main = "Column Dendrogram")
graph2ppt(file="DARs.LS_vs_OS.dendrograms.pptx", width=10, height=7, append=T)
## plot & export row dendrograms
hr.dend <- hr %>% as.dendrogram #convert to dendrogram
hr.dend.cut <- cut(hr.dend, h= 1.9) #cut at desired height
plot(hr.dend.cut$upper, horiz = T, main = "Row Dendrogram")
graph2ppt(file="DARs.LS_vs_OS.dendrograms.pptx", width=10, height=7, append=T)


##### DA Testing: ST_vs_OS # finds no DARs with FDR=10% ####
res      <- results(dds, contrast=c("Condition","ST","OS"))
summary(res) # finds no DARs with FDR=10%, thus stop here for now

##### Comparisons of the DARs found in LS_vs_ST & LS_vs_OS: DARs.LS_vs_dF.intersection ####

## Read them in again to avoid confusion
DARs.LS_vs_ST.df <- as.data.table(read.delim('DARs.LS_vs_ST.FDR1.txt', quote="\"'", sep = "\t"))
DARs.LS_vs_OS.df <- as.data.table(read.delim('DARs.LS_vs_OS.FDR1.txt', quote="\"'", sep = "\t"))
## Set DAR coords as key
setkey(DARs.LS_vs_ST.df, seqnames, start, end)
setkey(DARs.LS_vs_OS.df, seqnames, start, end)
## intersect DARs
DARs.LS_vs_dF.intersection <- foverlaps(DARs.LS_vs_ST.df, DARs.LS_vs_OS.df, type="any", nomatch=0)
write.table(DARs.LS_vs_dF.intersection, "DARs.LS_vs_df.intersection.txt", sep="\t", row.names=F, col.names=T, quote=F)

##### draw a venn diagram to illustrate overlap
draw.pairwise.venn(area1=4178 , area2=1283 , cross.area= 1187)
graph2ppt(file="DARs.LS_vs_dF.intersection.venn.pptx", width=3, height=3, append=T)
dev.off()

##### scatterplot fold changes & linear regress
ggplot(DARs.LS_vs_dF.intersection, aes(x=DARs.LS_vs_dF.intersection$log2FoldChange,
  y=DARs.LS_vs_dF.intersection$i.log2FoldChange)) + geom_point(shape=1) +
  geom_smooth(method=lm) + xlab("LS vs OS log2FC") + ylab("LS vs ST log2FC") + 
  ggtitle("DA of LS vs dF intersection") +  theme(plot.title= element_text(size= 14,
  face= "bold"), axis.text= element_text(size= 14), legend.text= element_text(size= 14),
  legend.title= element_text(size= 14), axis.title= element_text(size= 14),
  strip.text= element_text(size= 14))
graph2ppt(file="LS_vs_dF.intersection.scatterplot.ppt", width=3, height=3.15, append=T)
## fit linear model, get r2 (r is Pearson's correlation coefficient)
summary(lm(DARs.LS_vs_dF.intersection$i.log2FoldChange ~ DARs.LS_vs_dF.intersection$log2FoldChange))

##### make a heatmap
x <- as.matrix(DARs.LS_vs_dF.intersection[,c(11:19)])
## cluster rows by pearson, complete
hr <- hclust(as.dist(1-cor(t(x), method='pearson')), method='complete')
hc <- hclust(as.dist(1-cor(x, method='pearson')), method='complete')
## Heatmap2 w/ color bar. h & k modify cutree (k overrides h). Specify margins here as well.
mycl       <- cutree(hr, h=max(hr$height)/3, k = 5);
mycol      <- sample(rainbow(256)); mycol <- mycol[as.vector(mycl)]
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)
## generate a heatmap with rows & columns clustered
png('DARs.LS_vs_dF.intersection.heatmap_rc.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram= c("both"), col=my_palette,
          scale="row", density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## plot & export column dendrograms
hc.dend <- hc %>% as.dendrogram #convert to dendrogram
plot(hc.dend, main = "Column Dendrogram")
graph2ppt(file="DARs.LS_vs_dF.intersection.dendrograms.pptx", width=10, height=7, append=T)
## plot & export row dendrograms
hr.dend <- hr %>% as.dendrogram #convert to dendrogram
hr.dend.cut <- cut(hr.dend, h= 1.9) #cut at desired height
plot(hr.dend.cut$upper, horiz = T, main = "Row Dendrogram")
graph2ppt(file="DARs.LS_vs_dF.intersection.dendrograms.pptx", width=10, height=7, append=T)

##### Make an condition averaged heatmap from DARs.LS_vs_dF.intersection
## make a heatmap from LS_vs_dF
x <- as.matrix(DARs.LS_vs_dF.intersection[,c(9,10,38)])
## cluster rows by pearson, complete
hr <- hclust(as.dist(1-cor(t(x), method='pearson')), method='complete')
hc <- hclust(as.dist(1-cor(x, method='pearson')), method='complete')
## Heatmap2 w/ color bar. h & k modify cutree (k overrides h). Specify margins here as well.
mycl       <- cutree(hr, h=max(hr$height)/3, k = 5);
mycol      <- sample(rainbow(256)); mycol <- mycol[as.vector(mycl)]
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)
## generate a heatmap with rows & columns clustered
png('DARs.LS_vs_dF.intersection.avg.heatmap_rc.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram= c("both"), col=my_palette,
          scale="row", density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## plot & export column dendrograms
hc.dend <- hc %>% as.dendrogram #convert to dendrogram
plot(hc.dend, main = "Column Dendrogram")
graph2ppt(file="DARs.LS_vs_dF.intersection.avg.dendrograms.pptx", width=10, height=7, append=T)
## plot & export row dendrograms
hr.dend <- hr %>% as.dendrogram #convert to dendrogram
hr.dend.cut <- cut(hr.dend, h= 1.9) #cut at desired height
plot(hr.dend.cut$upper, horiz = T, main = "Row Dendrogram")
graph2ppt(file="DARs.LS_vs_dF.intersection.avg.dendrograms.pptx", width=10, height=7, append=T)

##### Split DARs.LS_vs_dF.intersection by LS or dF enriched for MOTIF finding ####
##### LS-enriched
DARs.LS_vs_dF.intersection.LS_enriched     <- subset(DARs.LS_vs_dF.intersection, log2FoldChange > 0)
DARs.LS_vs_dF.intersection.LS_enriched.bed <-  DARs.LS_vs_dF.intersection.LS_enriched[,c(1:3,7,8)]
## add a name_distToTSS column to bed for ID after motif finding
DARs.LS_vs_dF.intersection.LS_enriched.bed$namedist <- paste(DARs.LS_vs_dF.intersection.LS_enriched$SYMBOL,
                                                DARs.LS_vs_dF.intersection.LS_enriched$distanceToTSS, sep='_')
write.table(DARs.LS_vs_dF.intersection.LS_enriched.bed[,c(1:3,6)], "DARs.LS_vs_dF.intersection.LS_enriched.bed", sep="\t", row.names=F, col.names=F, quote=F)

##### dF-enriched
DARs.LS_vs_dF.intersection.dF_enriched     <- subset(DARs.LS_vs_dF.intersection, log2FoldChange < 0)
DARs.LS_vs_dF.intersection.dF_enriched.bed <-  DARs.LS_vs_dF.intersection.dF_enriched[,c(1:3,7,8)]
## add a name_distToTSS column to bed for ID after motif finding
DARs.LS_vs_dF.intersection.dF_enriched.bed$namedist <- paste(DARs.LS_vs_dF.intersection.dF_enriched$SYMBOL,
                                                             DARs.LS_vs_dF.intersection.dF_enriched$distanceToTSS, sep='_')
write.table(DARs.LS_vs_dF.intersection.dF_enriched.bed[,c(1:3,6)], "DARs.LS_vs_dF.intersection.dF_enriched.bed", sep="\t", row.names=F, col.names=F, quote=F)

##### Comparisons of the DARs found in either LS_vs_ST or LS_vs_OS: DARs.LS_vs_dF.union ####
DARs.LS_vs_ST.df_sub <- DARs.LS_vs_ST.df[,c(1:3,11:19)]
DARs.LS_vs_OS.df_sub <- DARs.LS_vs_OS.df[,c(1:3,11:19)]
DARs.LS_vs_dF.union  <- unique(rbind(DARs.LS_vs_ST.df_sub,DARs.LS_vs_OS.df_sub))
## calculate average values for each condition
DARs.LS_vs_dF.union$LS_avg <- rowMeans(DARs.LS_vs_dF.union[,c(4,5,9)])
DARs.LS_vs_dF.union$OS_avg <- rowMeans(DARs.LS_vs_dF.union[,c(6,7)])
DARs.LS_vs_dF.union$ST_avg <- rowMeans(DARs.LS_vs_dF.union[,c(8,10:12)])
## make a heatmap
x <- as.matrix(DARs.LS_vs_dF.union[,c(4:12)])
## cluster rows by pearson, complete
hr <- hclust(as.dist(1-cor(t(x), method='pearson')), method='complete')
hc <- hclust(as.dist(1-cor(x, method='pearson')), method='complete')
## Heatmap2 w/ color bar. h & k modify cutree (k overrides h). Specify margins here as well.
mycl       <- cutree(hr, h=max(hr$height)/3, k = 5);
mycol      <- sample(rainbow(256)); mycol <- mycol[as.vector(mycl)]
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)
## generate a heatmap with rows & columns clustered
png('DARs.LS_vs_dF.union.heatmap_rc.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram= c("both"), col=my_palette,
          scale="row", density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## plot & export column dendrograms
hc.dend <- hc %>% as.dendrogram #convert to dendrogram
plot(hc.dend, main = "Column Dendrogram")
graph2ppt(file="DARs.LS_vs_dF.union.dendrograms.pptx", width=10, height=7, append=T)
## plot & export row dendrograms
hr.dend <- hr %>% as.dendrogram #convert to dendrogram
hr.dend.cut <- cut(hr.dend, h= 1.9) #cut at desired height
plot(hr.dend.cut$upper, horiz = T, main = "Row Dendrogram")
graph2ppt(file="DARs.LS_vs_dF.union.dendrograms.pptx", width=10, height=7, append=T)

##### Make an condition averaged heatmap from DARs.LS_vs_dF.union
## make a heatmap from LS_vs_dF
x <- as.matrix(DARs.LS_vs_dF.union[,c(13:15)])
## cluster rows by pearson, complete
hr <- hclust(as.dist(1-cor(t(x), method='pearson')), method='complete')
hc <- hclust(as.dist(1-cor(x, method='pearson')), method='complete')
## Heatmap2 w/ color bar. h & k modify cutree (k overrides h). Specify margins here as well.
mycl       <- cutree(hr, h=max(hr$height)/3, k = 5);
mycol      <- sample(rainbow(256)); mycol <- mycol[as.vector(mycl)]
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)
## generate a heatmap with rows & columns clustered
png('DARs.LS_vs_dF.union.avg.heatmap_rc.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(x,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram= c("both"), col=my_palette,
          scale="row", density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## plot & export column dendrograms
hc.dend <- hc %>% as.dendrogram #convert to dendrogram
plot(hc.dend, main = "Column Dendrogram")
graph2ppt(file="DARs.LS_vs_dF.union.avg.dendrograms.pptx", width=10, height=7, append=T)
## plot & export row dendrograms
hr.dend <- hr %>% as.dendrogram #convert to dendrogram
hr.dend.cut <- cut(hr.dend, h= 1.9) #cut at desired height
plot(hr.dend.cut$upper, horiz = T, main = "Row Dendrogram")
graph2ppt(file="DARs.LS_vs_dF.union.avg.dendrograms.pptx", width=10, height=7, append=T)

##### Generate the union of DARs to test some correlations ####
## import DA test results & fix names
LS_v_ST.fullDAtest <- DARs.LS_vs_ST.all_res[,c("CHR","START","END","log2FoldChange","padj")]
colnames(LS_v_ST.fullDAtest)[4:5] <- c("log2FC.LS_vs_ST","padj.LS_vs_ST")
LS_v_ST.fullDAtest$coords <- paste(LS_v_ST.fullDAtest$CHR, LS_v_ST.fullDAtest$START, LS_v_ST.fullDAtest$END)

## import DA test results, fix names, subset to necessary columns
LS_v_OS.fullDAtest <- DARs.LS_vs_OS.all_res[,c("CHR","START","END","log2FoldChange","padj")]
colnames(LS_v_OS.fullDAtest)[4:5] <- c("log2FC.LS_vs_OS","padj.LS_vs_OS")
LS_v_OS.fullDAtest$coords <- paste(LS_v_OS.fullDAtest$CHR, LS_v_OS.fullDAtest$START, LS_v_OS.fullDAtest$END)

## merge the 2
union2 <- merge(LS_v_ST.fullDAtest, LS_v_OS.fullDAtest, by.x= "coords", by.y= "coords")
## subset to DARs from either LS vs dF comparison
DARs.union.df <- subset(union2, padj.LS_vs_ST < 0.1 | padj.LS_vs_OS < 0.1)

##### make correlation plots from the union of DEGs ####
##### Test/plot correlation between DEGs: LS vs ST & LS vs OS ####
Title  <- "Fold Change of Union DARs" #set title for exported filenames
df     <- DARs.union.df          #set dataframe of interest to df
dfy    <- df$log2FC.LS_vs_ST #set y feature
y_axis <- "LS_vs_ST (log2FC)"    #set y axis
dfx    <- df$log2FC.LS_vs_OS #set x feature
x_axis <- "log2FC (log2FC)"  #set x axis
Size   <- element_text(size= 14) #set text size for plot
## scatterplot with best fit line
ggplot(df, aes(x=dfx, y=dfy)) + geom_point(shape=1) + geom_smooth(method=lm) + xlab(x_axis)+
  ylab(y_axis) + ggtitle(Title) + theme(plot.title= element_text(size= 14, face= "bold"),
                                        axis.text= Size, legend.text= Size, legend.title= Size,
                                        axis.title= Size, strip.text= Size)
## fit linear model, get r2 (r is Pearson's correlation coefficient)
summary(lm(dfy ~ dfx))
graph2ppt(file=paste(Title,".correlations.pptx",sep=''), width=3, height=3.15, append=T)

##### Install/load required libraries ####
# source("http://bioconductor.org/biocLite.R")
# biocLite(pkgs = c("DESeq2","data.table","pasilla","DESeq","limma","ReportingTools",
# "GenomicRanges"))
# install.packages(pkgs= c("rJava","ReporteRs","ReporteRsjars","ggplot2","rtable","xtable",
#                          "VennDiagram","taRifx","devtools","dplyr","dendextend"))
# devtools::install_github('tomwenseleers/export',local=F)
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
zzz<-c("pheatmap","grid","gplots","ggplot2","export","devtools","DESeq2","pasilla","Biobase",
       "EBSeq","dplyr","data.table", "genefilter","FactoMineR","VennDiagram","DOSE","ReactomePA",
       "org.Hs.eg.db","clusterProfiler","pathview","DiffBind","dendextend","limma","ReportingTools",
       "TxDb.Hsapiens.UCSC.hg19.knownGene","GO.db","ChIPseeker","GenomicRanges")
lapply(zzz, require, character.only= T)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
