source("R_scripts/00_importation.R")

ngl.abund.clusters.cast <- readRDS("R_scripts/ngl.abund.clusters.cast.RDS")
mgCST.hclust <- readRDS("R_scripts/mgCST.hclust.RDS")
mgCST.dist <- readRDS("R_scripts/mgCST.dist.RDS")

relabund<-ngl.abund.clusters.cast/rowSums(ngl.abund.clusters.cast)

mgCSTs.sort.df <- read_csv("Data/mgCST_sort_color.csv")
mgCSTs.samples.df <-read_csv("Data/mgCST_samples_color.csv")
subspecies.with.colors <- read_csv("Data/subspecies_with_colors.csv")
subspecies.with.colors$Subspecies <- gsub("\\.", "_", subspecies.with.colors$Subspecies)

# Define the conditions
parameters <- read.csv("Data/mgCSTs_parameters_streamlit.csv")
deepsplit <- parameters$deepsplit
mincluster <- parameters$minClusterSize

# Use subset to filter the data
mgCSTs.samples <- subset(mgCSTs.samples.df, deepSplit == deepsplit & minClusterSize == mincluster)
mgCSTs.sort <- subset(mgCSTs.sort.df, deepSplit == deepsplit & minClusterSize == mincluster)

## order relabund table by most abundant taxon plus mgss's
test<-names(relabund[order(colSums(relabund), decreasing = TRUE)])
test.s<- unique(substr(test,1,nchar(test)-2))
taxon.order<-c()
for (i in test.s)
{
  l<-names(relabund[grep(pattern = i, names(relabund))])
  taxon.order<-append(taxon.order, l, after = length(taxon.order))
}
taxon.order<-unique(taxon.order)
relabund.s<-relabund[taxon.order]
l<-as.vector(table(sort(as.numeric(mgCST.hclust$order))))
colfunc <- colorRampPalette(c("khaki", "limegreen", "darkslategray1", "mediumblue", "magenta", "red"))
names(relabund.s)<-gsub(x=names(relabund.s), pattern = "_", replacement = " ")

dtc.df <- as.data.frame(mgCSTs.samples[c('dtc', 'sampleID')])


mgCST <- mgCSTs.sort[c('mgCST', "color_mgCST")]


color_mgcst <- character()
for (element in mgCSTs.samples$mgCST){
  color_mgcst <- append(color_mgcst, as.character(mgCST[mgCST$mgCST == element, "color_mgCST"]))
}
color_mgcst
mgCSTs.samples$color_mgCST <- color_mgcst

z<-rownames(relabund)
z.1<-z[order.dendrogram(as.dendrogram(mgCST.hclust))]
rownames(mgCSTs.samples)<-mgCSTs.samples$sampleID
mgCSTs.samples<-mgCSTs.samples[z.1,]
mgCSTs.samples<-mgCSTs.samples[order(mgCSTs.samples$mgCST),]
relabund.s<-relabund[mgCSTs.samples$sampleID, taxon.order]
#relabund.s<-relabund.s[order(mgCSTs.samples$mgCST),]
# l<-as.vector(table(mgCSTs.samples$mgCST))
# #l<-prepend(0, l)
# n<-cumsum(l)


colsidecolor <- mgCSTs.samples$color_mgCST # to display the color of samples by mgCST
# colsidecolor <- mgCSTs.samples$color_domTaxa # to display color of each sample's domTaxa

pdf("Medias/mgCST_heatmap.pdf", width=10, height=10)
gplots::heatmap.2(t(as.matrix(relabund.s[,1:200])),
                  Colv = FALSE,
                  Rowv = FALSE,
                  col=colfunc(300),
                  labRow = NULL,
                  labCol = FALSE,
                  keysize= 1.0,
                  densadj=0,
                  density.info='none',
                  key = FALSE,
                  key.ylab=NA,
                  key.title=NA,
                  key.ytickfun=NA,
                  key.xlab="Relative Abundance",
                  trace="none",
                  cexRow = 0.5, cexCol = 0.1, adjRow = c(1, NA), offsetRow = -56.5,
                  main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s))),
                  ColSideColors = colsidecolor,
                  title(main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s))), line = -2))

legend("top", ncol = 5,
       # legend = mgCST$mgCST,
       legend = mgCSTs.sort$domTaxa,
       col = mgCSTs.sort$color_mgCST,
       lty= 1, lwd = 4, cex=.5, bty="n", title = "mgCSTs")

dev.off()

print("Rscript has run")