# script for paper reproducibility  ####

getwd()
setwd("D:/Box Sync/02 - RNA seq work/RNA-seq-CABBI/2020_end/")

options(stringsAsFactors = F)
library(edgeR)
library(WGCNA)
library(affycoretools)
library(limma)
library(gplots)
library(rgl)
library(Glimma)
library(rtracklayer)
library(ggplot2)
library(PCAtools)

rm(list = ls())

# Input read counts ####
load("input/counts_file.Rdata")
write.table(d$counts, file = "results/Read_counts_raw.txt", sep = "\t")

d$counts[apply(d$counts,1,sum)==0,] %>%nrow()
# 45 - genes with zero expression in ALL  conditions 
d <- calcNormFactors(d) 
logCPM <- cpm(d, log=T, prior.count = 3)

# filtering out low reads  #### 
filter_index <- rowSums(logCPM>log2(1))>=3
d.filt <- d[filter_index,,keep.lib.size=F]
d.filt <- calcNormFactors(d.filt)
logCPM.filt <- cpm(d.filt, log = T, prior.count = 3)
write.table(logCPM.filt, file = "results/CPM_normalized_log2.txt", sep = "\t")

# PCA #### 
x11()
plotPCA(logCPM.filt, pch = 16, col = d.filt$samples$group, 
        groupnames = levels(d.filt$samples$group), main = "RT - PCA")

test_new <- pca(logCPM.filt)
postscript("results/PCA_2020_12_20.eps", width = 6, height = 6)
biplot(test_new, xlim = c(-100,100), ylim = c(-75,75), 
       lab = d.filt$samples$group, 
       pointSize = 4, gridlines.major = FALSE, gridlines.minor = FALSE)
dev.off()

jpeg("results/PCA_2020_12_20.jpeg", width = 6, height = 6, units = "in", res = 300)
biplot(test_new, xlim = c(-100,100), ylim = c(-75,75), 
       lab = d.filt$samples$group, 
       pointSize = 4, gridlines.major = FALSE, gridlines.minor = FALSE)
dev.off() 
# one-way anova #### 
group <- factor(d$samples$group, levels=c("YP","Glc", "Xyl", "Ac", "So" ))
design <- model.matrix(~group)
colnames(design)[-1] <- paste0(levels(group)[-1], "_vs_", levels(group)[1])

fit <- eBayes(lmFit(logCPM.filt, design = design), trend = T)
fit$genes <- d.filt$genes

cont.matrix <- makeContrasts(Glc_vs_YP, Xyl_vs_YP, Ac_vs_YP, So_vs_YP,
                             Xyl_vs_Glc = Xyl_vs_YP - Glc_vs_YP,
                             Ac_vs_Glc = Ac_vs_YP - Glc_vs_YP,
                             So_vs_Glc = So_vs_YP - Glc_vs_YP,
                             levels = design)
fit2 <- contrasts.fit(fit, contrasts = cont.matrix)
fit2 <- eBayes(fit2, trend = T)
summary(decideTests(fit2))
# Glc_vs_YP Xyl_vs_YP Ac_vs_YP So_vs_YP Xyl_vs_Glc Ac_vs_Glc So_vs_Glc
# Down        2630      1983     2192     1073       1283      2734      2419
# NotSig      2920      4199     4096     5833       5021      2821      3250
# Up          2522      1890     1784     1166       1768      2517      2403

Sig_index <- topTable(fit, n=Inf, sort.by = "none")$adj.P.Val < 0.05
table(Sig_index)
# Sig_index
# FALSE  TRUE 
# 964  7108 

# Get pairwise ####

source("input/summarizeFit.R")

out.conts <- summarizeFit(fit2, addAnova = TRUE)
names(out.conts)
colnames(out.conts)[27:29] <- c("Fstat.ANOVA","rawP.ANOVA","FDR.ANOVA")

# Heatmaps ####

logCPM.heat <- t(scale(t(logCPM.filt[Sig_index,])))
col.pan <- colorpanel(100, "blue", "white", "red")

x.cluster.h1 <- hclust(dist(logCPM.heat))

x11(w = 12, h = 6)
plot(x.cluster.h1, labels = FALSE)
abline(h = 7 )

rowcols.h1 <- cutree(x.cluster.h1, h = 7)

setEPS()
postscript("results/Heatmap_Rhoto_GeneBlocks.eps", width = 6, height = 10)
heatmap.2(logCPM.heat, col=col.pan, Rowv=as.dendrogram(x.cluster.h1),Colv=FALSE, scale="row",
          trace="none", dendrogram="none", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3,6,9,12), sepcolor = "black",
          margins = c(7,2),density.info = "none", RowSideColors = labels2colors(rowcols.h1),
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat), " genes" ))
dev.off()

jpeg("results/Heatmap_Rhoto_GeneBlocks.jpeg", width = 6, height = 10, units = "in", res = 600)
heatmap.2(logCPM.heat, col=col.pan, Rowv=as.dendrogram(x.cluster.h1),Colv=FALSE, scale="row",
          trace="none", dendrogram="none", cexRow=1, cexCol=0.9,labRow ="",
          keysize = 0.8,lwid = c(1.4,4), lhei = c(0.8,4), colsep = c(3,6,9,12), sepcolor = "black",
          margins = c(7,2),density.info = "none", RowSideColors = labels2colors(rowcols.h1),
          main = paste0("1-way ANOVA\nFDR < 0.05","\n",nrow(logCPM.heat), " genes" ))
dev.off()

# Combine data #### 
all.equal(out.conts$gene_id, rownames(logCPM.filt))
out.conts <- cbind(out.conts, logCPM.filt)
names(out.conts)
#  [1] "gene_id"         "proteinId"       "EC"              "KEGG"           
#  [5] "Amean"           "FC.Glc_vs_YP"    "rawP.Glc_vs_YP"  "FDR.Glc_vs_YP"  
#  [9] "FC.Xyl_vs_YP"    "rawP.Xyl_vs_YP"  "FDR.Xyl_vs_YP"   "FC.Ac_vs_YP"    
# [13] "rawP.Ac_vs_YP"   "FDR.Ac_vs_YP"    "FC.So_vs_YP"     "rawP.So_vs_YP"  
# [17] "FDR.So_vs_YP"    "FC.Xyl_vs_Glc"   "rawP.Xyl_vs_Glc" "FDR.Xyl_vs_Glc" 
# [21] "FC.Ac_vs_Glc"    "rawP.Ac_vs_Glc"  "FDR.Ac_vs_Glc"   "FC.So_vs_Glc"   
# [25] "rawP.So_vs_Glc"  "FDR.So_vs_Glc"   "Fstat.ANOVA"     "rawP.ANOVA"     
# [29] "FDR.ANOVA"       "heat_color"      "heat_order"      "Rhoto.Ac.1"     
# [33] "Rhoto.Ac.2"      "Rhoto.Ac.3"      "Rhoto.Glc.1"     "Rhoto.Glc.2"    
# [37] "Rhoto.Glc.3"     "Rhoto.So.1"      "Rhoto.So.2"      "Rhoto.So.3"     
# [41] "Rhoto.Xyl.1"     "Rhoto.Xyl.2"     "Rhoto.Xyl.3"     "Rhoto.YP.1"     
# [45] "Rhoto.YP.2"      "Rhoto.YP.3"
