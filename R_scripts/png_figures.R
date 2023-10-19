setwd("~/Desktop/mgss2")

# mgSs.stats <- readRDS("RDS_files/mgSs.stats.RDS")

# for (i in names(mgSs.stats)) {
#   ggplot(mgSs.stats[[i]], aes(x=No_Genes, y=log10(species_coverage), color=sample_cluster))+geom_point()+scale_color_brewer(palette="Set3") + 
#     theme_classic() + 
#     xlab("No. Genes") + 
#     ylab("log10(Coverage)") + 
#     theme(legend.position = "none") + 
#     ggtitle(i)
#   ggsave(paste("streamlit_app/Medias/mgSs_coverage/", i, "_subspecies_coverage_by_NoGenes.png", sep = ""), width=4, height=3)
#   print(i)

#   ggplot(mgSs.stats[[i]], aes(x=as.factor(sample_cluster), y=log10(species_coverage),fill=as.factor(sample_cluster), group=as.factor(sample_cluster),color=as.factor(sample_cluster))) + 
#     geom_violin()+scale_fill_brewer(palette="Set3") + 
#     scale_color_brewer(palette="Set3") + 
#     theme_classic() + 
#     theme(legend.position = "none") + 
#     ggtitle(i) + 
#     xlab("mgSs") + 
#     ylab("log10(Coverage)") +
#     geom_jitter(size=0.05, color="black", width=0.2)
  
#   ggsave(paste("streamlit_app/Medias/mgSs_coverage/", i, "_subspecies_coverage_boxplot.png", sep = ""), width=4, height=3)
  
#   print(i)
# }

install.packages("feather")
library(feather)
samples.clusters <- readRDS("samples.clusters.RDS")
sapply(names(samples.clusters),
       function(i) write_feather(samples.clusters[[i]], paste0("samples_clusters_df/",i,"_samples_clusters.feather")))