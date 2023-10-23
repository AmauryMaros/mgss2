setwd("~/Desktop/mgss2")

source("00_importation.R")

ngl.abund.clusters.cast <- readRDS("RDS_files/ngl.abund.clusters.cast.RDS")
mgCST.hclust <- readRDS("RDS_files/mgCST.hclust.RDS")
mgCST.dist <- readRDS("RDS_files/mgCST.dist.RDS")

relabund <- ngl.abund.clusters.cast/rowSums(ngl.abund.clusters.cast)

mgCSTs.sort.df <- read_csv("/Users/amaros/Desktop/mgss2/streamlit_app/mgCST_sort_color.csv")
mgCSTs.samples.df <- read_csv("/Users/amaros/Desktop/mgss2/streamlit_app/Data/mgCSTs.samples.df.csv")
subspecies.with.colors <- read_csv("/Users/amaros/Desktop/mgss2/streamlit_app/Data/subspecies_with_colors.csv")
subspecies.with.colors$Subspecies <- gsub("\\.", "_", subspecies.with.colors$Subspecies)


# Define the conditions
parameters <- read.csv("/Users/amaros/Desktop/mgss2/streamlit_app/R_scripts/mgCSTs_parameters_streamlit.csv")
deepsplit <- parameters$deepsplit
mincluster <- parameters$minClusterSize

# Use subset to filter the data
mgCSTs.samples <- subset(mgCSTs.samples.df, deepSplit == deepsplit & minClusterSize == mincluster)
mgCSTs.sort <- subset(mgCSTs.sort.df, deepSplit == deepsplit & minClusterSize == mincluster)



# Create an empty vector for 'color'
color <- character(0)
# Loop through elements in data2$domTaxa
for (element in mgCSTs.samples$mgCST) {

  # Subset the subspecies_with_colors data frame based on the modified 'element'
  subset_data <- mgCSTs.sort[mgCSTs.sort$mgCST == element,]
  
  # Check if 'subset_data' is not empty
  if (nrow(subset_data) > 0) {
    # Extract the 'Color' from the first row of 'subset_data'
    color <- c(color, subset_data$color[1])
  } else {
    # If 'subset_data' is empty, use the default color
    color <- c(color, "#8c8c8c")
  }
}

# Add the 'color' vector as a new column in data2
mgCSTs.samples$color <- color


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


mgCST <- mgCSTs.sort[c('mgCST', "color")]


# mgcsts.cols <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8,"Set2"), RColorBrewer::brewer.pal(12,"Set3"), RColorBrewer::brewer.pal(12, "Paired"))
mgcsts.cols <- mgCST$color
# mgcsts.cols <- mgCSTs.samples$color


# colored banner under the dendrogram
colsidecolors2 <- mgcsts.cols[as.factor(dtc.df[match(rownames(relabund.s), dtc.df$sampleID), "dtc"])]

pdf("/Users/amaros/Desktop/mgss2/streamlit_app/Medias/new_mgCST_heatmap_30.pdf", width=10, height=10)
gplots::heatmap.2(t(as.matrix(relabund.s[,1:200])),
                  Colv = FALSE,#as.dendrogram(mgCST.hclust),
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
                  dendrogram = "column",
                  main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s))),
                  ColSideColors = colsidecolors2,
                  title(main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s))), line = -2))


legend("topleft", ncol = 2,
       legend = mgCST$mgCST,
       col = mgCST$color,
       # legend = sort(unique(as.factor(dtc.df[match(rownames(relabund.s), dtc.df$sampleID), "dtc"]))),
       # col = sort(unique(mgCST$color[as.factor(dtc.df[match(rownames(relabund.s), dtc.df$sampleID), "dtc"])])),
       # legend = sort(unique(mgCSTs.sort$mgCST)),
       # col = mgCSTs.sort$color,
       lty= 1, lwd = 4, cex=.5, bty="n", title = "mgCSTs")

dev.off()

print("09_heatmap.R has run")
