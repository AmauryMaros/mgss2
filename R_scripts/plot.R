setwd("~/Desktop/mgss2")

source("00_importation.R")

ngl.abund.clusters.cast <- readRDS("RDS_files/ngl.abund.clusters.cast.RDS")
mgCST.hclust <- readRDS("RDS_files/mgCST.hclust.RDS")
mgCST.dist <- readRDS("RDS_files/mgCST.dist.RDS")

relabund<-ngl.abund.clusters.cast/rowSums(ngl.abund.clusters.cast)

rownames(mgCST.dist)<-rownames(relabund)
colnames(mgCST.dist)<-rownames(relabund)

dtc<-dynamicTreeCut::cutreeDynamic(mgCST.hclust, distM = mgCST.dist, method="hybrid", minClusterSize = 10, deepSplit = 3)
dtc.df<-as.data.frame(dtc)
dtc.df$sampleID<-rownames(dtc.df)

mgCSTs.samples<-cbind(dtc.df, domTaxa=colnames(relabund)[max.col(relabund, ties.method="first")], relabund=apply(relabund, 1, max))
mgCSTs<-as.data.frame(mgCSTs.samples %>% group_by(dtc) %>% dplyr::summarise(meanRelabund=mean(relabund)))

as.data.frame(mgCSTs.samples %>% count(dtc, domTaxa) %>% group_by(dtc) %>% slice(which.max(n)))

mgCSTs$domTaxa<-as.data.frame(mgCSTs.samples %>% count(dtc, domTaxa) %>% group_by(dtc) %>% slice(which.max(n)))[,2]
# for.sort<-c("Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_iners", "Lactobacillus_jensenii",
#             "BVAB1", "Gardnerella_vaginalis", "Gardnerella_swidsinki", "Gardnerella_piotii", "Gardnerella_leopoldi",
#             "Bifidobacterium_breve", "UBA629_sp005465875")
# 
# # for.sort <- c("Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_iners", "Lactobacillus_jensenii", "BVAB1", "Gardnerella_vaginalis", "Bifidobacterium_breve")
# v<-vector()
# for (i in for.sort){
#   v<-append(v, sort(as.vector(mgCSTs[grepl(pattern = i, x = mgCSTs$domTaxa), "domTaxa"])))
# }
# v
# mgCSTs.sort<-mgCSTs[pmatch(v, mgCSTs$domTaxa), ]

# mgCSTs.sort$mgCST<-c(1:46)

mgCSTs <- mgCSTs[order(mgCSTs$domTaxa),]

mgCST <- as.data.frame(
  rbind(
    c("1", "#FE0308"), c("2", "#F54C5E"), c("3", "#F07084"), c("4", "#EC94A5"),
    c("5", "#F0BCCC"), c("6", "#F6D3DA"), c("7", "#86C61A"), c("8", "#B4DB29"),
    c("9", "#DBEA77"), c("10", "#FF7200"), c("11", "#F68A11"), c("12", "#F8A40E"),
    c("13", "#F3BC11"), c("14", "#f7d15a"), c("15", "#FAE50D"), c("16", "#F3F46E"),
    c("17", "#448A73"), c("18", "#89BEAB"), c("19", "#BCD6CD"), c("20", "#221886"),
    c("21", "#3E3792"), c("22", "#5D579E"), c("23", "#7C76AC"), c("24", "#9A98BF"),
    c("25", "#C9C8D8"), c("26", "#98C999"), c("27", "#989898"), c("28", "#123456"),
    c("29", "#ABCDEF"), c("30", "#FF5733"), c("31", "#00FF00"), c("32", "#FFFF00"), 
    c("33", "#FF00FF"), c("34", "#00FFFF"), c("35", "#AABBCC"), c("36", "#112233"),
    c("37", "#445566"), c("38", "#778899"), c("39", "#99AABB"), c("40", "#CCDDEE"),
    c("41", "#556677"), c("42", "#8899AA"), c("43", "#BBCCEE"), c("44", "#FFEEDD"),
    c("45", "#765432"), c("46", "#DCBA98"),c("47", "#AAB2CC"), c("", "white"), c("NA", "black"))
)
names(mgCST) = c('dtc','color')
# write_csv(mgCST, "mgCST_color.csv")
# 
mgCST<-as.data.frame(rbind(c("1", "#FE0308"),c("2", "#F54C5E"), c("3", "#F07084"), c("4", "#EC94A5"),
                           c("5", "#F0BCCC"),c("6", "#F6D3DA"),c("7", "#86C61A"), c("8", "#B4DB29"),
                           c("9", "#DBEA77"), c("10", "#FF7200"), c("11", "#F68A11"),c("12", "#F8A40E"),
                           c("13", "#F3BC11"),c("14", "#f7d15a"), c("15", "#FAE50D"),c("16", "#F3F46E"),
                           c("17", "#448A73"),c("18", "#89BEAB"), c("19", "#BCD6CD") ,c("20", "#221886")
                           ,c("21", "#3E3792"),c("22", "#5D579E"),c("23", "#7C76AC"),c("24", "#9A98BF"),
                           c("25", "#C9C8D8"),c("26", "#98C999"), c("27", "#989898"), c("", "white"), c("NA", "black")))

names(mgCST)<-c("mgCST", "color")
mgCSTs$mgCST<-c(1:31)
mgCSTs$color<-mgCST[match(mgCSTs$dtc, mgCST$mgCST), "color"]
mgCSTs.samples<-merge(mgCSTs.samples, mgCSTs[,c("dtc", "mgCST", "color")], all.x=TRUE)

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

### Color by Project
# 
sources<-read.csv("SOURCE_DATA/VIRGO2_projects.csv", header=T)
# proj.cols<-c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(5,"Set2"))
# colsidecolors = proj.cols[as.factor(sources[match(rownames(relabund.s), sources$sampleID), "Project"])]
# 
# pdf("FIGURES/mgCST_heatmap_test_by_project.pdf", width=10, height=10)
# gplots::heatmap.2(t(as.matrix(relabund.s[,1:100:200])), Colv = as.dendrogram(mgCST.hclust), Rowv = FALSE, col=colfunc(100), labRow = NULL, labCol = FALSE, keysize= 1.0, densadj=0, density.info='none', key = FALSE, key.ylab=NA, key.title=NA, key.ytickfun=NA, key.xlab="Relative Abundance", trace="none", cexRow = 0.5, cexCol = 0.1, adjRow = c(1, NA), offsetRow = -56.5, dendrogram = "column", main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s))), ColSideColors = colsidecolors, title(main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s))), line = -2))
# legend("topleft", ncol = 2, legend = unique(as.factor(sources[match(rownames(relabund.s), sources$sampleID), "Project"])), col = unique(proj.cols[as.factor(sources[match(rownames(relabund.s), sources$sampleID), "Project"])]), lty= 1, lwd = 4, cex=.5, bty="n", title = "Source Project")
# dev.off()
# 
# relabund.s.s<-relabund.s[, !grepl("MultiGenera", names(relabund.s))]
# pdf("FIGURES/mgCST_heatmap_test_by_project_noMultigenera.pdf", width=10, height=10)
# gplots::heatmap.2(t(as.matrix(relabund.s.s[,1:100:200])), Colv = as.dendrogram(mgCST.hclust), Rowv = FALSE, col=colfunc(100), labRow = NULL, labCol = FALSE, keysize= 1.0, densadj=0, density.info='none', key = FALSE, key.ylab=NA, key.title=NA, key.ytickfun=NA, key.xlab="Relative Abundance", trace="none", cexRow = 0.5, cexCol = 0.1, adjRow = c(1, NA), offsetRow = -56.5, dendrogram = "column", main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s.s))), ColSideColors = colsidecolors, title(main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s.s))), line = -2))
# legend("topleft", ncol = 2, legend = unique(as.factor(sources[match(rownames(relabund.s.s), sources$sampleID), "Project"])), col = unique(proj.cols[as.factor(sources[match(rownames(relabund.s.s), sources$sampleID), "Project"])]), lty= 1, lwd = 4, cex=.5, bty="n", title = "Source Project")
# dev.off()

### Color by mgCSTs

dtc.df$sampleID <- sources$sampleID
mgcsts.cols <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8,"Set2"), RColorBrewer::brewer.pal(12,"Set3"), RColorBrewer::brewer.pal(12, "Paired"))
colsidecolors2 = mgcsts.cols[as.factor(dtc.df[match(rownames(relabund.s), dtc.df$sampleID), "dtc"])]

pdf("/Users/amaros/Desktop/mgss2/streamlit_app/Medias/mgCST_heatmap_by_mgcst.pdf", width=10, height=10)
gplots::heatmap.2(t(as.matrix(relabund.s[,1:80])), Colv = as.dendrogram(mgCST.hclust), Rowv = FALSE, col=colfunc(100), labRow = NULL, labCol = FALSE, keysize= 1.0, densadj=0, density.info='none', key = FALSE, key.ylab=NA, key.title=NA, key.ytickfun=NA, key.xlab="Relative Abundance", trace="none", cexRow = 0.5, cexCol = 0.1, adjRow = c(1, NA), offsetRow = -56.5, dendrogram = "column", main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s))), ColSideColors = colsidecolors2, title(main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s))), line = -2))
legend("topleft", ncol = 2, legend = sort(unique(as.factor(dtc.df[match(rownames(relabund.s), dtc.df$sampleID), "dtc"]))), col = unique(mgcsts.cols[as.factor(dtc.df[match(rownames(relabund.s), dtc.df$sampleID), "dtc"])]), lty= 1, lwd = 4, cex=.5, bty="n", title = "mgCSTs")
dev.off()

relabund.s.s<-relabund.s[, !grepl("MultiGenera", names(relabund.s))]
pdf("/Users/amaros/Desktop/mgss2/streamlit_app/Medias/mgCST_heatmap_by_mgcst_noMultigenera.pdf", width=10, height=10)
gplots::heatmap.2(t(as.matrix(relabund.s.s[,1:80])), Colv = as.dendrogram(mgCST.hclust), Rowv = FALSE, col=colfunc(100), labRow = NULL, labCol = FALSE, keysize= 1.0, densadj=0, density.info='none', key = FALSE, key.ylab=NA, key.title=NA, key.ytickfun=NA, key.xlab="Relative Abundance", trace="none", cexRow = 0.5, cexCol = 0.1, adjRow = c(1, NA), offsetRow = -56.5, dendrogram = "column", main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s.s))), ColSideColors = colsidecolors2, title(main = paste("mgCSTs JSD, Ward linkage HC, nSamples=", paste(nrow(relabund.s.s))), line = -2))
legend("topleft", ncol = 2, legend = sort(unique(as.factor(dtc.df[match(rownames(relabund.s.s), dtc.df$sampleID), "dtc"]))), col = unique(mgcsts.cols[as.factor(dtc.df[match(rownames(relabund.s.s), dtc.df$sampleID), "dtc"])]), lty= 1, lwd = 4, cex=.5, bty="n", title = "mgCSTs")
dev.off()

print("09_heatmap.R has run")