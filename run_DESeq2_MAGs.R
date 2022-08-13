library("DESeq2")
library("Biobase")
library("ggplot2")
library("data.table")
library("pheatmap")
library("tidyverse")
library("tibble")
library("vegan")
library("gplots")
library("compositions")
library("ape")
library("cowplot")
library("plotly")
library("purrr")
library("reshape2")
library("MASS")
library("ggpubr")
library("stringr")
library(rstatix)
library(broom)
library("RColorBrewer")

working_directory = "~/Dropbox/PostDoc/4th_wedge/MAGs_disrat0.33/"
setwd(working_directory)

my_palette <- colorRampPalette(c('#ffffff','#8856a7'))(n = 299)

# Load MAG taxonomy file
tax <- read.csv("MAG_tax.csv", header=T, stringsAsFactors = F)
tax$Bin <- gsub(".fa", "", tax$Bin)

# Load Rsubread output for DESeq
infile = "Rsubread_output.csv"
mtx_counts = read.csv(infile, header = TRUE, row.names = 1, stringsAsFactors = F)
colnames(mtx_counts) <- gsub(".unf.bam", "", colnames(mtx_counts))
colnames(mtx_counts) <- gsub("12C16O", "", colnames(mtx_counts))
mtx_counts$bin <- gsub("_k[0-9].*", "", rownames(mtx_counts))
# Keep only significantly enriched MAGs (CI.L>0)
enr_MAGs <- c()
enr_MAGs_100 <- read.csv("AFE_MAGs_CIL_tax_100.csv", stringsAsFactors = F, header=T)
enr_MAGs <- unique(enr_MAGs_100$bin)
enr_MAGs_50 <- read.csv("AFE_MAGs_CIL_tax_50.csv", stringsAsFactors = F, header=T)
enr_MAGs <- c(enr_MAGs, unique(enr_MAGs_50$bin))
enr_MAGs_tax <- tax[tax$Bin %in% unique(enr_MAGs),]
enr_MAGs_tax$phylum <- gsub(".*p__", "", gsub(";c.*", "", enr_MAGs_tax$classification))
enr_MAGs_tax$class <- gsub(".*c__", "", gsub(";o.*", "", enr_MAGs_tax$classification))
write.csv(enr_MAGs_tax, "enr_MAGs_all_tax.csv")[,-1]
enr_MAGs_tax$Bin <- gsub("\\.", "_", gsub("-", "_", enr_MAGs_tax$Bin))

# Pie chart of MAG taxonomy
ggcols <- read.csv("~/Dropbox/PostDoc/4th_wedge/MAGs_disrat0.33/class_ggcolors.csv", stringsAsFactors = F)
enr_MAGs_tax <- enr_MAGs_tax[order(enr_MAGs_tax$phylum),]
p_freq <- as.data.frame(matrix(ncol=2, nrow=length(unique(enr_MAGs_tax$phylum))))
colnames(p_freq) <- c("Phylum","count")
i=1
for (p in unique(enr_MAGs_tax$phylum)) {
  p_freq$Phylum[i] <- p
  p_freq$count[i] <- nrow(tax[enr_MAGs_tax$phylum==p,])
  i=i+1
}
c_freq <- as.data.frame(matrix(ncol=2, nrow=length(unique(enr_MAGs_tax$class))))
colnames(c_freq) <- c("Class","count")
i=1
for (c in unique(enr_MAGs_tax$class)) {
  c_freq$Class[i] <- c
  c_freq$count[i] <- nrow(enr_MAGs_tax[enr_MAGs_tax$class==c,])
  i=i+1
}
c_freq <- c_freq %>% mutate(rand=1)
c_freq <- merge(c_freq, unique(enr_MAGs_tax[,3:4]), by.x="Class",by.y="class")
c_freq <- merge(c_freq, ggcols, by="Class", all.x=T)
c_freq <- c_freq[order(c_freq$phylum),]
library(plotrix)
pie(p_freq$count, labels=p_freq$Phylum, col=rainbow(nrow(p_freq)), main="Enriched MAG collection by phylum")
ggplot(c_freq, aes(fill=Class, y=count, x=rand)) + 
  geom_col(position="stack", color="black")

# Analyze only genes from isotopically enriched MAGs
mtx_counts_enr <- mtx_counts[mtx_counts$bin %in% enr_MAGs_tax$Bin,]
mtx_counts <- mtx_counts_enr[,-ncol(mtx_counts_enr)]

# DRAM annotations
# CAZy distilled functions (use flag --distillate_gene_names when running DRAM distill)
ann_cazy <- read.csv("DRAM/distill_MAGs_enr_all_metabolism_CAZy.csv",header=T,stringsAsFactors = F)

# kofam annotations
ann <- read.delim("DRAM/annotations.tsv", sep="\t", header=T, stringsAsFactors = F, strip.white=TRUE)
ann <- ann[ann$kegg_id != "" | ann$peptidase_id != "" | ann$cazy_hits != "",]
# DRAM adds the bin name again - remove it
ann$X <- gsub("^12.*_12", "12", ann$X)
ann$X <- gsub("^H.*_H", "H", ann$X)

mtx_matrix = as.matrix(sapply(mtx_counts, as.integer)) 
rownames(mtx_matrix) = rownames(mtx_counts)
mtx_matrix[is.na(mtx_matrix)] = 0
mtx_matrix = mtx_matrix[,order(colnames(mtx_matrix))]

designfile = "mtx_design_ecofun_genomes4.txt" 
# Without time point 4
#designfile = "From_Erin/mtx_design_ecofun_genomes3.txt" 
#mtx_matrix <- mtx_matrix[,-which(colnames(mtx_matrix) %like% "H4_")]

mtx_design = read.delim(designfile, header = TRUE)
mtx_design <- mtx_design[!(rownames(mtx_design) %in% c("H3_Rhizo_39", "H1_RhizoLitter_2", "H2_RhizoLitter_9")),] # Remove these because these columns are mostly 0s, was messing up the geometric mean 

# Running DESeq2: comparing by condition (combination of time+prec)
dds2 = DESeqDataSetFromMatrix(countData = mtx_matrix, colData = mtx_design, design = ~ condition)
dds2 = dds2[ rowSums(counts(dds2)) > 5, ] # Filter out rows with < 5 counts
design(dds2) = ~ 0+condition
dds2 = DESeq(dds2)
res2 = results(dds2)
resultsNames(dds2)
saveRDS(dds2, "dds_traits_4thwedge_enr.rds")
dds2_counts <- counts(dds2, normalized=TRUE)
write.csv(dds2_counts, "dds2_counts_enr.csv")
rownames(dds2_counts) <- rownames(res2)
dds2_ann <- merge(dds2_counts, ann, by.x="row.names", by.y="X")
write.csv(dds2_ann, "dds2_counts_ann_enr.csv")

# Beta diversity matrix
dds2_bray_dm = vegdist(t(dds2_counts), method="bray")
dds2_dist = as.dist(dds2_bray_dm)

# PCoA
pcoa_bray <- pcoa(dds2_dist)
jac_variances <- data.frame(pcoa_bray$values$Relative_eig) %>% 
  rownames_to_column(var = "PCaxis") %>% 
  data.frame
pcoa_bray_df <- data.frame(pcoa_bray$vectors) %>% 
  rownames_to_column(var = "sample") %>% 
  data.frame

pcoa_bray_df$sample <- gsub("_unf_init.bam", "", gsub("X_", "", pcoa_bray_df$sample))
pcoa_bray_df$grp <- gsub("_P[0-9].*$", "", pcoa_bray_df$sample)
pcoa_bray_df$plotnum <- gsub(".*P","P", pcoa_bray_df$sample)
pcoa_bray_df$time <- gsub("h.*", "", pcoa_bray_df$grp)
pcoa_bray_df$time <- factor(pcoa_bray_df$time, levels=c(0,24,48,72,168))
eigenvalues<-round(jac_variances[,2], digits = 4)*100

ggplot(pcoa_bray_df) +
  geom_point(aes(x = Axis.1, y = Axis.2, color=time, shape=plotnum), size = 4) +
  ylab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 1 ',eigenvalues[1],'%')) + #Extract x axis value from variance
  ggtitle('Log-Ratio PCoA Ordination') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_manual(values=rev(brewer.pal(name="Paired", n=6))) +
  scale_shape_manual(values = c(15, 0, 1, 2, 16, 17))

# Compare between consecutive time points
res = res2
resultsNames(dds2)
res2 = results(dds2, contrast = c("condition", "24h_100", "0h_100")) # group, numerator, denominator
comparison = "24h_100_vs_0h_100"
res_sig = which(res2$padj < 0.05)
length(res_sig)
baseMean_df = data.frame(res2$baseMean)
betas_df = data.frame(res2$log2FoldChange)
stat_df = data.frame(res2$stat)
padj_df = data.frame(res2$padj)
setnames(baseMean_df, "res2.baseMean", paste0("baseMean_", comparison))
setnames(betas_df, "res2.log2FoldChange", paste0("betas_", comparison))
setnames(stat_df, "res2.stat", paste0("stat_", comparison))
setnames(padj_df, "res2.padj", paste0("padj_", comparison))

res2 = results(dds2, contrast = c("condition", "48h_100", "24h_100")) # group, numerator, denominator
comparison = "48h_100_vs_24h_100"
res_sig = which(res2$padj < 0.05)
length(res_sig)
baseMean_df = data.frame(res2$baseMean)
betas_df = data.frame(betas_df, res2$log2FoldChange)
stat_df = data.frame(stat_df, res2$stat)
padj_df = data.frame(padj_df, res2$padj)
setnames(baseMean_df, "res2.baseMean", paste0("baseMean_", comparison))
setnames(betas_df, "res2.log2FoldChange", paste0("betas_", comparison))
setnames(stat_df, "res2.stat", paste0("stat_", comparison))
setnames(padj_df, "res2.padj", paste0("padj_", comparison))

res2 = results(dds2, contrast = c("condition", "72h_100", "48h_100")) # group, numerator, denominator
comparison = "72h_100_vs_48h_100"
res_sig = which(res2$padj < 0.05)
length(res_sig)
baseMean_df = data.frame(res2$baseMean)
betas_df = data.frame(betas_df, res2$log2FoldChange)
stat_df = data.frame(stat_df, res2$stat)
padj_df = data.frame(padj_df, res2$padj)
setnames(baseMean_df, "res2.baseMean", paste0("baseMean_", comparison))
setnames(betas_df, "res2.log2FoldChange", paste0("betas_", comparison))
setnames(stat_df, "res2.stat", paste0("stat_", comparison))
setnames(padj_df, "res2.padj", paste0("padj_", comparison))

res2 = results(dds2, contrast = c("condition", "168h_100", "72h_100"))
comparison = "168h_100_vs_72h_100"
res_sig = which(res2$padj < 0.05)
length(res_sig)
baseMean_df = data.frame(res2$baseMean)
betas_df = data.frame(betas_df, res2$log2FoldChange)
stat_df = data.frame(stat_df, res2$stat)
padj_df = data.frame(padj_df, res2$padj)
setnames(baseMean_df, "res2.baseMean", paste0("baseMean_", comparison))
setnames(betas_df, "res2.log2FoldChange", paste0("betas_", comparison))
setnames(stat_df, "res2.stat", paste0("stat_", comparison))
setnames(padj_df, "res2.padj", paste0("padj_", comparison))

res2 = results(dds2, contrast = c("condition", "48h_50", "0h_50"))
comparison = "48h_50_vs_0h_50"
res_sig = which(res2$padj < 0.05)
length(res_sig)
baseMean_df = data.frame(res2$baseMean)
betas_df = data.frame(betas_df, res2$log2FoldChange)
stat_df = data.frame(stat_df, res2$stat)
padj_df = data.frame(padj_df, res2$padj)
setnames(baseMean_df, "res2.baseMean", paste0("baseMean_", comparison))
setnames(betas_df, "res2.log2FoldChange", paste0("betas_", comparison))
setnames(stat_df, "res2.stat",paste0("stat_", comparison))
setnames(padj_df, "res2.padj", paste0("padj_", comparison))

res2 = results(dds2, contrast = c("condition", "168h_50", "48h_50"))
comparison = "168h_50_vs_48h_50"
res_sig = which(res2$padj < 0.05)
length(res_sig)
baseMean_df = data.frame(res2$baseMean)
betas_df = data.frame(betas_df, res2$log2FoldChange)
stat_df = data.frame(stat_df, res2$stat)
padj_df = data.frame(padj_df, res2$padj)
setnames(baseMean_df, "res2.baseMean", paste0("baseMean_", comparison))
setnames(betas_df, "res2.log2FoldChange", paste0("betas_", comparison))
setnames(stat_df, "res2.stat",paste0("stat_", comparison))
setnames(padj_df, "res2.padj", paste0("padj_", comparison))

res2 = results(dds2, contrast = c("condition", "0h_50", "0h_100"))
comparison = "0h_50_vs_0h_100"
res_sig = which(res2$padj < 0.05)
length(res_sig)
baseMean_df = data.frame(res2$baseMean)
betas_df = data.frame(betas_df, (res2$log2FoldChange))
stat_df = data.frame(stat_df, res2$stat)
padj_df = data.frame(padj_df, res2$padj)
setnames(baseMean_df, "res2.baseMean", "baseMean_0h_50_vs_0h_100")
setnames(betas_df, colnames(betas_df)[ncol(betas_df)], "betas_0h_50_vs_0h_100")
setnames(stat_df, "res2.stat","stat_0h_50_vs_0h_100")
setnames(padj_df, "res2.padj", "padj_0h_50_vs_0h_100")

res2 = results(dds2, contrast = c("condition", "48h_100", "48h_50"))
comparison = "48h_100_vs_48h_50"
res_sig = which(res2$padj < 0.05)
length(res_sig)
baseMean_df = data.frame(res2$baseMean)
betas_df = data.frame(betas_df, res2$log2FoldChange)
stat_df = data.frame(stat_df, res2$stat)
padj_df = data.frame(padj_df, res2$padj)
setnames(baseMean_df, "res2.baseMean", paste0("baseMean_", comparison))
setnames(betas_df, "res2.log2FoldChange", paste0("betas_", comparison))
setnames(stat_df, "res2.stat",paste0("stat_", comparison))
setnames(padj_df, "res2.padj", paste0("padj_", comparison))

res2 = results(dds2, contrast = c("condition", "168h_100", "168h_50"))
comparison = "168h_100_vs_168h_50"
res_sig = which(res2$padj < 0.05)
length(res_sig)
baseMean_df = data.frame(res2$baseMean)
betas_df = data.frame(betas_df, res2$log2FoldChange)
stat_df = data.frame(stat_df, res2$stat)
padj_df = data.frame(padj_df, res2$padj)
setnames(baseMean_df, "res2.baseMean", paste0("baseMean_", comparison))
setnames(betas_df, "res2.log2FoldChange", paste0("betas_", comparison))
setnames(stat_df, "res2.stat",paste0("stat_", comparison))
setnames(padj_df, "res2.padj", paste0("padj_", comparison))

res_all_data_bulkcomp = data.frame(betas_df, stat_df, padj_df)
rownames(res_all_data_bulkcomp) = rownames(res)
colnames(res_all_data_bulkcomp)[7] <- "betas_0h_50_vs_0h_100"
res_all_data_bulkcomp$Bin <- gsub("_k[0-9].*", "", rownames(res_all_data_bulkcomp))
saveRDS(res_all_data_bulkcomp, file = paste(working_directory,"/res_enrMAGs.rds", sep = ""))

# To reload the saved RDS:
res_all_data_bulkcomp <- readRDS("res_enrMAGs.rds")

# Only keep p values where the comparison by Wald test was significantly different
padj_matrix = as.matrix(res_all_data_bulkcomp[,19:27]) # This is essentially padj_df, but only saved res_all_data, so working from that
padj_matrix[is.na(padj_matrix)] = 1 # Change NA to 1 so will get filtered out in next step
padj_matrix[which(padj_matrix > 0.05)] = NA # Convert all p values > 0.05 to NA
res_all_data_bulkcomp_sig = data.frame(res_all_data_bulkcomp[,c(1:9,28)], as.data.frame(padj_matrix)) # Remake res_all_data with only significant p values, all non-significant p values are NA
res_all_data_bulkcomp_sig$X <- rownames(res_all_data_bulkcomp_sig)

saveRDS(res_all_data_bulkcomp_sig, file = paste(working_directory,"/res_enrMAGs_sig.rds", sep = ""))

tax$Bin <- gsub("\\.", "_", gsub("-","_", tax$Bin))
res_all_data_bulkcomp_sig <- merge(res_all_data_bulkcomp_sig, tax, by = "Bin")
write.csv(res_all_data_bulkcomp_sig, file = paste0(working_directory,"results_enrMAGs_sig.csv"))

# To load the saved RDS:
res_all_data_bulkcomp_sig <- readRDS("res_enrMAGs_sig.rds")
res_all_data_bulkcomp_sig$X <- rownames(res_all_data_bulkcomp_sig)

# subset to rows that have at least one significant p value (<0.05)
res_all_data_bulkcomp_sig_only <- subset(res_all_data_bulkcomp_sig, !(is.na(padj_24h_100_vs_0h_100) & (is.na(padj_48h_100_vs_24h_100)) & (is.na(padj_72h_100_vs_48h_100)) & (is.na(padj_168h_100_vs_72h_100)) & (is.na(padj_168h_100_vs_72h_100)) & (is.na(padj_48h_50_vs_0h_50)) & (is.na(padj_168h_50_vs_48h_50)) & (is.na(padj_0h_50_vs_0h_100)) & (is.na(padj_48h_100_vs_48h_50)) & (is.na(padj_168h_100_vs_168h_50))))
nrow(res_all_data_bulkcomp_sig_only)
write.csv(res_all_data_bulkcomp_sig_only, "DESeq2_sigonly_MAGs_enr.csv")
res_all_data_bulkcomp_sig_only_ann <- merge(res_all_data_bulkcomp_sig_only, ann, by="X", all.x=T)
write.csv(res_all_data_bulkcomp_sig_only_ann, "res_all_data_bulkcomp_sig_only_MAGs_enr.csv")

# To load the saved dataframes of genes with significant differential abundance
res_all_data_bulkcomp_sig_only <- read.csv("DESeq2_sigonly_MAGs_enr.csv",stringsAsFactors = F, strip.white=TRUE)
res_all_data_bulkcomp_sig_only_ann <- read.csv("res_all_data_bulkcomp_sig_only_MAGs_enr.csv",stringsAsFactors = F, strip.white=TRUE)[,-1]

# Create Venn diagrams of significant genes from time points / treatments
library("ggvenn")
x <- list(
  as.character(res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_24h_100_vs_0h_100),1]),
  as.character(res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_48h_100_vs_24h_100),1]),
  as.character(res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_72h_100_vs_48h_100),1]),
  as.character(res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_168h_100_vs_72h_100),1])
)
names(x) <- c("24h", "48h", "72h", "168h")
ggvenn(x)
x <- list(
  res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_48h_50_vs_0h_50),1],
  res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_168h_50_vs_48h_50),1]
)
names(x) <- c("48h", "168h")
ggvenn(x)
x <- list(
  res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_0h_50_vs_0h_100),1],
  res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_48h_100_vs_48h_50),1],
  res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_168h_100_vs_1),1]
)
names(x) <- c("0h", "48h", "168h")
ggvenn(x)

# subset to rows that have at least one significant p value (<0.05) in the 100% treatment
res_all_data_bulkcomp_sig_only_100 <- subset(res_all_data_bulkcomp_sig_only, !(is.na(padj_24h_100_vs_0h_100) & (is.na(padj_48h_100_vs_24h_100)) & (is.na(padj_72h_100_vs_48h_100)) & (is.na(padj_168h_100_vs_72h_100)) & (is.na(padj_168h_100_vs_72h_100))))
nrow(res_all_data_bulkcomp_sig_only_100)
effect_size_100 <- c(
nrow(subset(res_all_data_bulkcomp_sig_only_100, !(is.na(padj_24h_100_vs_0h_100)))),
nrow(subset(res_all_data_bulkcomp_sig_only_100, !(is.na(padj_48h_100_vs_24h_100)))),
nrow(subset(res_all_data_bulkcomp_sig_only_100, !(is.na(padj_72h_100_vs_48h_100)))),
nrow(subset(res_all_data_bulkcomp_sig_only_100, !(is.na(padj_168h_100_vs_72h_100)))))
names(effect_size_100) <- c("24_vs_0", "48_vs_24", "72_vs_48", "168_vs_72")

# subset to rows that have at least one significant p value (<0.05) in the 50% treatment
res_all_data_bulkcomp_sig_only_50 <- subset(res_all_data_bulkcomp_sig_only, !((is.na(padj_48h_50_vs_0h_50)) & (is.na(padj_168h_50_vs_48h_50))))
nrow(res_all_data_bulkcomp_sig_only_50)
effect_size_50 <- c(
  nrow(subset(res_all_data_bulkcomp_sig_only_100, !(is.na(padj_48h_50_vs_0h_50)))),
  nrow(subset(res_all_data_bulkcomp_sig_only_100, !(is.na(padj_168h_50_vs_48h_50)))))
names(effect_size_50) <- c("48_vs_0", "168_vs_48")

# subset to rows that have at least one significant p value (<0.05) when comparing treatments
res_all_data_bulkcomp_sig_only_treat <- subset(res_all_data_bulkcomp_sig_only, !((is.na(padj_0h_50_vs_0h_100)) & (is.na(padj_48h_100_vs_48h_50)) & (is.na(padj_168h_100_vs_168h_50))))
nrow(res_all_data_bulkcomp_sig_only_treat)
effect_size_treat <- c(
nrow(subset(res_all_data_bulkcomp_sig_only_treat, !((is.na(padj_0h_50_vs_0h_100))))),
nrow(subset(res_all_data_bulkcomp_sig_only_treat, !((is.na(padj_48h_100_vs_48h_50))))),
nrow(subset(res_all_data_bulkcomp_sig_only_treat, !((is.na(padj_168h_100_vs_168h_50))))))
names(effect_size_treat) <- c("0", "48", "168")

# Figure 3A
effect_size_all <- c(effect_size_100, effect_size_50, effect_size_treat)
barplot(log10(1+effect_size_all), ylab="Log 10 number of MAGs with p-value < 0.05", las=2)

# Plot the log2-fold change compared to t0 
betas <- res_all_data_bulkcomp_sig
betas <- betas %>% filter_at(vars(starts_with("padj")), any_vars(.<0.05))
#betas <- merge(betas, ann, by.x="row.names",by.y="X")
pivot <- betas[,c(1:10, 12:ncol(betas))]
colnames(pivot)[2:10] <- gsub("betas_", "up_", colnames(pivot)[2:10])
colnames(pivot)[11:19] <- gsub("padj_", "down_", colnames(pivot)[11:19])
for (c in 2:10) {
  pivot[is.na(pivot[,c+9]),c] <- 0
  pivot[which(pivot[,c] > 0), c] <- 1
  pivot[which(pivot[,c] > 0), c+9] <- 0
  pivot[which(pivot[,c] < 0), c+9] <- 1
  pivot[which(pivot[,c] < 0), c] <- 0
}
pivot[is.na(pivot)] <- 0
colSums(pivot[2:18])
write.csv(pivot, "DESeq_pivot_enrMAGs.csv")

# Volcano plots of differential abundance
p24 <- ggplot(data=betas, aes(x=betas_24h_100_vs_0h_100, y=-log10(padj_24h_100_vs_0h_100))) + 
  geom_point() + theme_classic() + ggtitle("24h") + 
  ylab("-log10 adjujsted p-value") + xlab("log2 fold change") + 
  geom_vline(xintercept=c(-1, 1), col="red")
p48 <- ggplot(data=betas, aes(x=betas_48h_100_vs_24h_100, y=-log10(padj_48h_100_vs_24h_100))) + 
  geom_point() + theme_classic() + ggtitle("48h") + 
  ylab("-log10 adjujsted p-value") + xlab("log2 fold change") + 
  geom_vline(xintercept=c(-1, 1), col="red")
p72 <- ggplot(data=betas, aes(x=betas_72h_100_vs_48h_100, y=-log10(padj_72h_100_vs_48h_100))) + 
  geom_point() + theme_classic() + ggtitle("72h") + 
  ylab("-log10 adjujsted p-value") + xlab("log2 fold change") + 
  geom_vline(xintercept=c(-1, 1), col="red")
p168 <- ggplot(data=betas, aes(x=betas_168h_100_vs_72h_100, y=-log10(padj_168h_100_vs_72h_100))) + 
  geom_point() + theme_classic() + ggtitle("168h") + 
  ylab("-log10 adjujsted p-value") + xlab("log2 fold change") + 
  geom_vline(xintercept=c(-1, 1), col="red")
ggarrange(p24, p48, p72, p168, nrow=2, ncol=2)

p48 <- ggplot(data=betas, aes(x=betas_48h_50_vs_0h_50, y=-log10(padj_48h_50_vs_0h_50))) + 
  geom_point() + theme_classic() + ggtitle("48h 50%") + 
  ylab("-log10 adjujsted p-value") + xlab("log2 fold change") + 
  geom_vline(xintercept=c(-1, 1), col="red")
p168 <- ggplot(data=betas, aes(x=betas_168h_50_vs_48h_50, y=-log10(padj_168h_50_vs_48h_50))) + 
  geom_point() + theme_classic() + ggtitle("168h 50%") + 
  ylab("-log10 adjujsted p-value") + xlab("log2 fold change") + 
  geom_vline(xintercept=c(-1, 1), col="red")
ggarrange(p48, p168, nrow=1, ncol=2)

p0 <- ggplot(data=betas, aes(x=betas_0h_50_vs_0h_100, y=-log10(padj_0h_50_vs_0h_100))) + 
  geom_point() + theme_classic() + ggtitle("0h legacy") + 
  ylab("-log10 adjujsted p-value") + xlab("log2 fold change") + 
  geom_vline(xintercept=c(-1, 1), col="red")
p48 <- ggplot(data=betas, aes(x=betas_48h_100_vs_48h_50, y=-log10(padj_48h_100_vs_48h_50))) + 
  geom_point() + theme_classic() + ggtitle("48h legacy") + 
  ylab("-log10 adjujsted p-value") + xlab("log2 fold change") + 
  geom_vline(xintercept=c(-1, 1), col="red")
p168 <- ggplot(data=betas, aes(x=betas_168h_100_vs_168h_50, y=-log10(padj_168h_100_vs_168h_50))) + 
  geom_point() + theme_classic() + ggtitle("168h legacy") + 
  ylab("-log10 adjujsted p-value") + xlab("log2 fold change") + 
  geom_vline(xintercept=c(-1, 1), col="red")
ggarrange(p0, p48, p168, nrow=1, ncol=3)

# Violin plots of significant betas - sup. fig. 2
plot.new()
par(mfrow=c(2, 4))
min_betas <- min(betas$betas_0h_50_vs_0h_100)
max_betas <- max(betas$betas_48h_100_vs_48h_50)
violin_plot(betas$betas_24h_100_vs_0h_100[!is.na(betas$padj_24h_100_vs_0h_100)], main="24h 100%", ylim=c(min_betas, max_betas))
violin_plot(betas$betas_48h_100_vs_24h_100[!is.na(betas$padj_48h_100_vs_24h_100)], main="48h 100%", ylim=c(min_betas, max_betas))
violin_plot(betas$betas_72h_100_vs_48h_100[!is.na(betas$padj_72h_100_vs_48h_100)], main="72h 100%", ylim=c(min_betas, max_betas))
violin_plot(betas$betas_168h_100_vs_72h_100[!is.na(betas$padj_168h_100_vs_72h_100)], main="168h 100%", ylim=c(min_betas, max_betas))
violin_plot(betas$betas_48h_50_vs_0h_50[!is.na(betas$padj_48h_50_vs_0h_50)], main="48h 50%", ylim=c(min_betas, max_betas))
violin_plot(betas$betas_168h_50_vs_48h_50[!is.na(betas$padj_168h_50_vs_48h_50)], main="168h 50%", ylim=c(min_betas, max_betas))
violin_plot(betas$betas_0h_50_vs_0h_100[!is.na(betas$padj_0h_50_vs_0h_100)], main="0h legacy", ylim=c(min_betas, max_betas))
violin_plot(betas$betas_48h_100_vs_48h_50[!is.na(betas$padj_48h_100_vs_48h_50)], main="48h legacy", ylim=c(min_betas, max_betas))

# CAZy genes and pathways
ann_cazy_long <- melt(ann_cazy, id.vars=1:5)
ann_cazy_long <- ann_cazy_long %>% filter(value!="")
# DRAM adds the bin name to the beginning of the fasta header. Remove that.
ann_cazy_long$value <- sub("^.*_12", "12", ann_cazy_long$value)
ann_cazy_long$value <- sub("^.*_H", "H", ann_cazy_long$value)
colnames(ann_cazy_long)[6] <- "Bin"
colnames(ann_cazy_long)[7] <- "X"
ann_cazy_long$Bin <- gsub("-", "_", gsub("\\.", "_",ann_cazy_long$Bin))
betas_CAZy <- merge(betas, ann_cazy_long)
colnames(betas_CAZy)[2] <- "ORF"
betas_grouped_CAZy <- betas_CAZy %>% group_by(subheader) %>% summarize(across(starts_with("betas"), ~ mean(.x, na.rm = TRUE)))

nHalf <- 50
Min <- min(betas_grouped_CAZy[,8:10])
Max <- max(betas_grouped_CAZy[,8:10])
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white","lightgoldenrod1","red"), space="Lab")(nHalf)
# Make blue-yellow-red palette for heatmap
rampcols <- c(rc1, rc2)

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)
heatmap.2(as.matrix(betas_grouped_CAZy[,8:10]), cexCol=0.7, cexRow = 0.4, density.info="none", trace="none", labRow = betas_grouped_CAZy$subheader, Rowv=TRUE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)

# Change columns based on subset (100% 2:5, 50% 6:7, legacy 8:10)
Min <- min(betas_grouped_CAZy[,2:10])
Max <- max(betas_grouped_CAZy[,2:10])
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white","lightgoldenrod1","red"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

# 100% only
heatmap.2(as.matrix(betas_grouped_CAZy[,2:5]), cexCol=0.7, cexRow = 0.4, density.info="none", trace="none", labRow = betas_grouped_CAZy$subheader, Rowv=TRUE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)
# keep only values >1.2 or < 1.2
betas_grouped_CAZy_1.2 <- betas_grouped_CAZy %>% filter_at(vars(2:5), any_vars(abs(.)>=1.2))
heatmap.2(as.matrix(betas_grouped_CAZy_1.2[,2:5]), cexCol=0.7, cexRow = 0.2, density.info="none", trace="none", labRow = betas_grouped_CAZy_1.2$subheader, Rowv=TRUE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)

# 50% only
heatmap.2(as.matrix(betas_grouped_CAZy[,6:7]), cexCol=0.7, cexRow = 0.4, density.info="none", trace="none", labRow = betas_grouped_CAZy$subheader, Rowv=TRUE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)
betas_grouped_CAZy_1.2 <- betas_grouped_CAZy %>% filter_at(vars(6:7), any_vars(abs(.)>=1.2))
heatmap.2(as.matrix(betas_grouped_CAZy_1.2[,6:7]), cexCol=0.7, cexRow = 0.2, density.info="none", trace="none", labRow = betas_grouped_CAZy_1.2$subheader, Rowv=TRUE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)

# Precipitation treatment
heatmap.2(as.matrix(betas_grouped_CAZy[,8:10]), cexCol=0.7, cexRow = 0.4, density.info="none", trace="none", labRow = betas_grouped_CAZy$subheader, Rowv=TRUE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)
betas_grouped_CAZy_1.2 <- betas_grouped_CAZy %>% filter_at(vars(8:10), any_vars(abs(.)>=1.2))
heatmap.2(as.matrix(betas_grouped_CAZy_1.2[,8:10]), cexCol=0.7, cexRow = 0.1, density.info="none", trace="none", labRow = betas_grouped_CAZy_1.2$subheader, Rowv=TRUE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)

# Nitrogen cycle genes
betas <- merge(betas, ann, by="X")
nitrogen_ko <- c("K10944", "K10945", "K10946", "K10535","K00370", "K00371", "K00374", "K02567", "K02568", "K00368", "K15864", "K04561", "K02305", "K15877", "K00376", "K00362", "K00363", "K03385", "K15876", "K00367", "K10534", "K00372", "K00360", "K17877", "K00366", "K01915", "K00260", "K00264", "K00265", "K00266", "K15371", "K01183", "K03791", "K17525", "K00990", "K04751", "K07708", "K07712", "K03320", "K02586", "K02588", "K02591", "K00531", "K05916", "K01725", "K20932", "K20933", "K20934", "K20935", "K04752", "K07004", "K07273", "K01428", "K01429", "K01430")
nitrogen_gene <- c("amoA", "amoB", "amoC", "hao", "narG", "narH", "narI", "napA", "napB", "nirK", "nirS", "norB", "norC", "cyp55", "nosZ", "nirB", "nirD", "nrfA", "nrfH", "narB", "NR", "nasA", "nasB", "nit-6", "nirA", "glnA", "gudB", "glt1", "gltB", "gltD", "gdh2", "chit1", "putative_chitinase","CHID1", "glnD", "glnB", "glnL", "glnG", "amt", "nifD", "nifK", "nifH", "anfG", "hmp", "cynS", "hzsC", "hzsB", "hzsA", "hdh", "glnK","Xds", "lys", "ureA", "ureB", "ureC")
ko_to_gene_N <- read.csv("ko_to_gene_N.csv", header=T, stringsAsFactors = F)[,1:4]
betas_N <- betas[betas$kegg_id %in% ko_to_gene_N$KO,]
betas_grouped_N <- betas_N %>% group_by(kegg_id) %>% summarise(across(starts_with("betas"), ~ mean(.x, na.rm = TRUE)))
betas_grouped_N <- merge(betas_grouped_N, ko_to_gene_N, by.x="kegg_id", by.y="KO")
betas_cols <- which(colnames(betas_grouped_N) %like% "betas")
# Order by pathway
betas_grouped_N <- betas_grouped_N[order(betas_grouped_N$pathway),]
colmin <- min(betas_grouped_N[,2:10]) # -5
colmax <- max(betas_grouped_N[,2:10]) # 3.4
nHalf <- 50

Min <- colmin
Max <- colmax
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white","lightgoldenrod1","red"), space="Lab")(nHalf)
# Make heatmap palette
rampcols <- c(rc1, rc2)

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)
heatmap.2(as.matrix(betas_grouped_N[,betas_cols]), cexRow = 0.4, density.info="none", trace="none", labRow = betas_grouped_N$gene, Rowv=TRUE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)
heatmap.2(as.matrix(betas_grouped_N[,2:5]), cexRow = 0.4, density.info="none", trace="none", labRow = betas_grouped_N$gene, Rowv=FALSE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)
heatmap.2(as.matrix(betas_grouped_N[,6:7]), cexRow = 0.4, density.info="none", trace="none", labRow = betas_grouped_N$gene, Rowv=FALSE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)
heatmap.2(as.matrix(betas_grouped_N[,8:10]), cexRow = 0.4, density.info="none", trace="none", labRow = betas_grouped_N$gene, Rowv=FALSE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)
heatmap.2(as.matrix(betas_grouped_N[,9:10]), cexRow = 0.4, density.info="none", trace="none", labRow = betas_grouped_N$gene, Rowv=FALSE, Colv=FALSE, symbreaks = FALSE, col=rampcols, breaks = rampbreaks)

# Subset annotated significant genes by time point
res_all_data_bulkcomp_sig_only_ann_24 <- res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_24h_100_vs_0h_100),]
res_all_data_bulkcomp_sig_only_ann_24 <- res_all_data_bulkcomp_sig_only_ann_24[order(res_all_data_bulkcomp_sig_only_ann_24$betas_24h_100_vs_0h_100, decreasing=T),]
write.csv(res_all_data_bulkcomp_sig_only_ann_24, "res_all_data_bulkcomp_sig_only_ann_24_100.csv")
res_all_data_bulkcomp_sig_only_ann_48_100 <- res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_48h_100_vs_24h_100),]
res_all_data_bulkcomp_sig_only_ann_48_100 <- res_all_data_bulkcomp_sig_only_ann_48_100[order(res_all_data_bulkcomp_sig_only_ann_48_100$betas_48h_100_vs_24h_100, decreasing=T),]
write.csv(res_all_data_bulkcomp_sig_only_ann_48_100, "res_all_data_bulkcomp_sig_only_ann_48_100.csv")
res_all_data_bulkcomp_sig_only_ann_72 <- res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_72h_100_vs_48h_100),]
res_all_data_bulkcomp_sig_only_ann_72 <- res_all_data_bulkcomp_sig_only_ann_72[order(res_all_data_bulkcomp_sig_only_ann_24$betas_72h_100_vs_48h_100, decreasing=T),]
write.csv(res_all_data_bulkcomp_sig_only_ann_72, "res_all_data_bulkcomp_sig_only_ann_72_100.csv")
res_all_data_bulkcomp_sig_only_ann_168_100 <- res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_168h_100_vs_72h_100),]
res_all_data_bulkcomp_sig_only_ann_168_100 <- res_all_data_bulkcomp_sig_only_ann_168_100[order(res_all_data_bulkcomp_sig_only_ann_168_100$betas_168h_100_vs_72h_100, decreasing=T),]
write.csv(res_all_data_bulkcomp_sig_only_ann_168_100, "res_all_data_bulkcomp_sig_only_ann_168_100.csv")
res_all_data_bulkcomp_sig_only_ann_48_50 <- res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_48h_50_vs_0h_50),]
res_all_data_bulkcomp_sig_only_ann_48_50 <- res_all_data_bulkcomp_sig_only_ann_48_50[order(res_all_data_bulkcomp_sig_only_ann_48_50$betas_48h_100_vs_48h_50, decreasing=T),]
write.csv(res_all_data_bulkcomp_sig_only_ann_48_50, "res_all_data_bulkcomp_sig_only_ann_48_50.csv")
res_all_data_bulkcomp_sig_only_ann_168_50 <- res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_168h_50_vs_48h_50),]
res_all_data_bulkcomp_sig_only_ann_168_50 <- res_all_data_bulkcomp_sig_only_ann_168_50[order(res_all_data_bulkcomp_sig_only_ann_168_50$betas_168h_50_vs_48h_50, decreasing=T),]
write.csv(res_all_data_bulkcomp_sig_only_ann_168_50, "res_all_data_bulkcomp_sig_only_ann_168_50.csv")
res_all_data_bulkcomp_sig_only_ann_0_legacy <- res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_0h_100_vs_0h_50),]
res_all_data_bulkcomp_sig_only_ann_0_legacy <- res_all_data_bulkcomp_sig_only_ann_0_legacy[order(res_all_data_bulkcomp_sig_only_ann_0_legacy$padj_0h_100_vs_0h_50, decreasing=T),]
write.csv(res_all_data_bulkcomp_sig_only_ann_0_legacy, "res_all_data_bulkcomp_sig_only_ann_0_legacy.csv")
res_all_data_bulkcomp_sig_only_ann_48_legacy <- res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_48h_100_vs_48h_50),]
res_all_data_bulkcomp_sig_only_ann_48_legacy <- res_all_data_bulkcomp_sig_only_ann_48_legacy[order(res_all_data_bulkcomp_sig_only_ann_48_legacy$padj_48h_100_vs_48h_50, decreasing=T),]
write.csv(res_all_data_bulkcomp_sig_only_ann_48_legacy, "res_all_data_bulkcomp_sig_only_ann_48_legacy.csv")
res_all_data_bulkcomp_sig_only_ann_168_legacy <- res_all_data_bulkcomp_sig_only_ann[!is.na(res_all_data_bulkcomp_sig_only_ann$padj_168h_100_vs_168h_50),]
res_all_data_bulkcomp_sig_only_ann_168_legacy <- res_all_data_bulkcomp_sig_only_ann_168_legacy[order(res_all_data_bulkcomp_sig_only_ann_168_legacy$padj_168h_100_vs_168h_50, decreasing=T),]
write.csv(res_all_data_bulkcomp_sig_only_ann_168_legacy, "res_all_data_bulkcomp_sig_only_ann_168_legacy.csv")

# Get MAGs with lower confidence interval limit > 0 either in 100% or 50% treatment for figure 3
enr_MAGs <- c()
enr_MAGs_100 <- read.csv("AFE_MAGs_CIL_tax_100.csv", stringsAsFactors = F, header=T)
enr_MAGs <- unique(enr_MAGs_100$bin)
enr_MAGs_50 <- read.csv("AFE_MAGs_CIL_tax_50.csv", stringsAsFactors = F, header=T)
enr_MAGs <- c(enr_MAGs, unique(enr_MAGs_50$bin))
afe100 <- enr_MAGs_100[,c(which(colnames(enr_MAGs_100) %like% "AFE"), which(colnames(enr_MAGs_100) == "bin"))]
afe100 <- unique(afe100)
afe50 <- enr_MAGs_50[,c(which(colnames(enr_MAGs_50) %like% "AFE"), which(colnames(enr_MAGs_50) == "bin"))]
afe50 <- unique(afe50)
afe <- merge(afe100, afe50, by="bin", all=T)
afe[is.na(afe)] <- 0
colnames(afe) <- c("bin", "AFE_24h", "AFE_48h", "AFE_72h", "AFE_168h", "AFE_48h_50", "AFE_168h_50")
nHalf = 50
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    

# Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white","lightgoldenrod1","red"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)
afe <- merge(afe, tax, by.x="bin", by.y="Bin")
afe$Phylum <- gsub(".*;p__", "", gsub(";c.*", "", afe$classification))
afe$Class <- gsub(".*;c__", "", gsub(";o.*", "", afe$classification))
afe <- merge(afe, ggcols, by="Class")
afe_cols <- which(colnames(afe) %like% "AFE")
afe <- afe[order(afe$Phylum, afe$Color,decreasing = T),]
# Draw AFE heatmap
heatmap.2(as.matrix(afe[,afe_cols]), col=rampcols, cexCol = 1, cexRow=0.2, density.info="none", trace="none", labRow = afe$bin, Rowv = FALSE, Colv=FALSE, RowSideColors = afe$Color)
# Create RPKM heatmap for figure 3
rpkm <- read.delim("coverm_RPKM.txt", header=T, stringsAsFactors = F, sep="\t")
rpkm <- rpkm[rowSums(rpkm[,2:25])>0,]
rpkm_melt <- melt(rpkm)
rpkm_melt$variable <- gsub("_unf_init.RPKM", "", rpkm_melt$variable)
rpkm_melt$grp <- gsub("_P.*", "", rpkm_melt$variable)
rpkm_melt <- rpkm_melt %>% group_by(Genome, grp) %>% summarize_at("value", mean)
rpkm <- dcast(rpkm_melt, Genome~grp)
rpkm_enr <- rpkm[rpkm$Genome %in% afe$bin,]
rpkm_enr <- rpkm_enr[match(afe$bin, rpkm_enr$Genome),]
ord <- c("X12C16O_0h_100", "X12C16O_24h_100", "X12C16O_48h_100", "X12C16O_72h_100", "X12C16O_168h_100", "X12C16O_0h_50", "X12C16O_48h_50", "X12C16O_168h_50")
heatmap.2(as.matrix(rpkm_enr[,ord]), col=rc2, cexCol = 1, cexRow=0.2, density.info="none", trace="none", labRow = rpkm_enr$Genome, Rowv = FALSE, Colv=FALSE)

leg <- unique(afe[,c("Class", "Color")])
plot.new()
legend("topleft", legend=leg$Class, fill=leg$Color, xpd=T, inset=c(0,-0.2))

# Functions in enriched MAGs - 79 with AFE confidence interval lower limit>0
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
library(RColorBrewer)
enr_ann <- read.delim("DRAM/distill_MAGs_enr_all.tsv",sep="\t", header=T, stringsAsFactors = F)
enr_ann_afe <- merge(enr_ann, afe, by.x="genome", by.y="bin")
colnames(enr_ann_afe)
enr_ann_afe[enr_ann_afe == "True"] <- 1
enr_ann_afe[enr_ann_afe == "False"] <- 0
cazy_cols <- which(colnames(enr_ann_afe) %like% "CAZ")
n_cols <- which(colnames(enr_ann_afe) %like% "Nitrogen")
s_cols <- which(colnames(enr_ann_afe) %like% "Sulfur")
enr_ann_afe[enr_ann_afe$AFE_24h>0,ncol(enr_ann_afe)+1] <- "#b3cde3"
colnames(enr_ann_afe)[ncol(enr_ann_afe)] <- "AFE_24h_100_col"
enr_ann_afe[enr_ann_afe$AFE_48h_100>0,ncol(enr_ann_afe)+1] <- "#8c96c6"
colnames(enr_ann_afe)[ncol(enr_ann_afe)] <- "AFE_48h_100_col"
enr_ann_afe[enr_ann_afe$AFE_72h>0,ncol(enr_ann_afe)+1] <- "#8856a7"
colnames(enr_ann_afe)[ncol(enr_ann_afe)] <- "AFE_72h_100_col"
enr_ann_afe[enr_ann_afe$AFE_168h_100>0,ncol(enr_ann_afe)+1] <- "#810f7c"
colnames(enr_ann_afe)[ncol(enr_ann_afe)] <- "AFE_168h_100_col"
enr_ann_afe$AFE_24h_100_col[is.na(enr_ann_afe$AFE_24h_100_col)] <- "lightgrey"
enr_ann_afe$AFE_48h_100_col[is.na(enr_ann_afe$AFE_48h_100_col)] <- "lightgrey"
enr_ann_afe$AFE_72h_100_col[is.na(enr_ann_afe$AFE_72h_100_col)] <- "lightgrey"
enr_ann_afe$AFE_168h_100_col[is.na(enr_ann_afe$AFE_168h_100_col)] <- "lightgrey"
enr_ann_afe[enr_ann_afe$AFE_48h_50 > 0,ncol(enr_ann_afe)+1] <- "#8c96c6"
colnames(enr_ann_afe)[ncol(enr_ann_afe)] <- "AFE_48h_50_col"
enr_ann_afe[enr_ann_afe$AFE_168h_50 > 0,ncol(enr_ann_afe)+1] <- "#810f7c"
colnames(enr_ann_afe)[ncol(enr_ann_afe)] <- "AFE_168h_50_col"
enr_ann_afe$AFE_48h_50_col[is.na(enr_ann_afe$AFE_48h_50_col)] <- "lightgrey"
enr_ann_afe$AFE_168h_50_col[is.na(enr_ann_afe$AFE_168h_50_col)] <- "lightgrey"
enr_ann_afe[,34:99] <- apply(enr_ann_afe[,34:99], 1, as.numeric)
pal <- brewer.pal(8, "Blues")

enr_ann_afe <- enr_ann_afe[order(enr_ann_afe$AFE_24h_100, enr_ann_afe$AFE_48h_100, enr_ann_afe$AFE_72h_100, enr_ann_afe$AFE_168h_100),]
heatmap.3(as.matrix(enr_ann_afe[,cazy_cols]), scale="none", dendrogram="none", margins=c(6,12),
          Rowv=FALSE, Colv=FALSE, RowSideColors=t(enr_ann_afe[,(ncol(enr_ann_afe)-5):ncol(enr_ann_afe)]), symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none", main="CAZy", labCol=colnames(enr_ann_afe)[cazy_cols], 
          labRow=enr_ann_afe$genome, cexRow=0.4, col=pal, RowSideColorsSize=2)

# N cycling genes by HMMs and peptidases identified by hotpep
enr_ann_hmm <- ann[ann$fasta %in% enr_ann$genome,]
enr_ann_hmm_N <- enr_ann_hmm %>% filter(kegg_id %in% nitrogen_ko & peptidase_id == "")
temp <- enr_ann_afe[,c(1, (ncol(enr_ann_afe)-12):ncol(enr_ann_afe))]
enr_ann_hmm_N <- merge(enr_ann_hmm_N, temp, by.y="genome", by.x="fasta")
enr_ann_hmm_pep <- enr_ann_hmm %>% filter(peptidase_id != "")
enr_ann_hmm_pep <- merge(enr_ann_hmm_pep, temp, by.y="genome", by.x="fasta")
enr_ann_hmm_N <- enr_ann_hmm_N[, -c(2:8, 11:23)]
enr_ann_hmm_N <- unique(enr_ann_hmm_N)
enr_ann_hmm_N$count24_100 <- 1
enr_ann_hmm_N$count48_100 <- 1
enr_ann_hmm_N$count72_100 <- 1
enr_ann_hmm_N$count168_100 <- 1
enr_ann_hmm_N$count48_50 <- 1
enr_ann_hmm_N$count168_50 <- 1
enr_ann_hmm_N$count24_100[enr_ann_hmm_N$AFE_24h==0] <- 0
enr_ann_hmm_N$count48_100[enr_ann_hmm_N$AFE_48h==0] <- 0
enr_ann_hmm_N$count72_100[enr_ann_hmm_N$AFE_72h==0] <- 0
enr_ann_hmm_N$count168_100[enr_ann_hmm_N$AFE_168h==0] <- 0
enr_ann_hmm_N$count48_50[enr_ann_hmm_N$AFE_48h_50==0] <- 0
enr_ann_hmm_N$count168_50[enr_ann_hmm_N$AFE_168h_50==0] <- 0

temp <- enr_ann_hmm_N[enr_ann_hmm_N$AFE_24h_100_col != "lightgrey",] %>% group_by(kegg_id) %>% tally(count24_100, name="h24")
N_df <- data.frame(matrix(nrow=nrow(temp), ncol=2))
N_df <- temp
temp <- enr_ann_hmm_N[enr_ann_hmm_N$AFE_72h_100_col != "lightgrey",] %>% group_by(kegg_id) %>% tally(count72_100, name="h72")
N_df <- merge(N_df, temp)
temp <- enr_ann_hmm_N[enr_ann_hmm_N$AFE_48h_100_col != "lightgrey",] %>% group_by(kegg_id) %>% tally(count48_100, name="h48")
N_df <- merge(N_df, temp)
temp <- enr_ann_hmm_N[enr_ann_hmm_N$AFE_168h_100_col != "lightgrey",] %>% group_by(kegg_id) %>% tally(count168_100, name="h168")
N_df <- merge(N_df, temp)
temp <- enr_ann_hmm_N[enr_ann_hmm_N$AFE_48h_50_col != "lightgrey",] %>% group_by(kegg_id) %>% tally(count48_50, name="h48_50")
N_df <- merge(N_df, temp)
temp <- enr_ann_hmm_N[enr_ann_hmm_N$AFE_168h_50_col != "lightgrey",] %>% group_by(kegg_id) %>% tally(count168_50, name="h168_50")
N_df <- merge(N_df, temp)
N_df <- N_df[order(-N_df$h24),]

par(mfrow=c(3,2))
barplot(N_df$h24, las=2, main="N pathways in MAGs enriched at 24h", ylim=c(0,80))
barplot(N_df$h48, las=2, main="N pathways in MAGs enriched at 72h", ylim=c(0,80))
barplot(N_df$h72, las=2, main="N pathways in MAGs enriched at 48h", ylim=c(0,80))
barplot(N_df$h168, las=2, main="N pathways in MAGs enriched at 168h", ylim=c(0,80))
barplot(N_df$h48_50, las=2, main="N pathways in MAGs enriched at 48h 50%", ylim=c(0,80))
barplot(N_df$h168_50, las=2, main="N pathways in MAGs enriched at 168h 50%", ylim=c(0,80))

N_df <- data.frame(N_df)
N_df_melt <- melt(N_df)
kruskal.test(value~variable, data=N_df_melt) # p=0.72

# Peptidases per genome per time point
enr_ann_hmm_pep$count24_100 <- 1
enr_ann_hmm_pep$count48_100 <- 1
enr_ann_hmm_pep$count72_100 <- 1
enr_ann_hmm_pep$count168_100 <- 1
enr_ann_hmm_pep$count48_50 <- 1
enr_ann_hmm_pep$count168_50 <- 1
enr_ann_hmm_pep$count24_100[enr_ann_hmm_N$AFE_24h==0] <- 0
enr_ann_hmm_pep$count48_100[enr_ann_hmm_N$AFE_48h==0] <- 0
enr_ann_hmm_pep$count72_100[enr_ann_hmm_N$AFE_72h==0] <- 0
enr_ann_hmm_pep$count168_100[enr_ann_hmm_N$AFE_168h==0] <- 0
enr_ann_hmm_pep$count48_50[enr_ann_hmm_N$AFE_48h_50==0] <- 0
enr_ann_hmm_pep$count168_50[enr_ann_hmm_N$AFE_168h_50==0] <- 0

pep_tally <- melt(enr_ann_hmm_pep[, c(1,35:40)])
pep_tally <- pep_tally %>% group_by(fasta, variable) %>% tally(value)
pep_tally <- data.frame(pep_tally)

# figure 7 D
ggboxplot(data=pep_tally, x="variable", y="n") + ylab("Peptidase genes per MAGs") + xlab("")
pep_tally$variable <- as.character(pep_tally$variable)
aovpep <- aov(n~variable, data=pep_tally)
summary(aovpep) # p=0.0005, F=4.5
tukeypep <- tukey_hsd(aovpep) # only significant p<0.05 differences are 168_50~48_100, 48_50~48_100, 72_100~48_50
# More peptidases per genome in the 50% samples

# CAZy abundance changes between time points
cazy_df <- c()
plot.new()
par(mfrow=c(3,2))
temp <- colSums(enr_ann_afe[enr_ann_afe$AFE_24h_100_col != "lightgrey",cazy_cols])
cazy_df$h24 <- temp

ord <- order(-temp)
temp <- colSums(enr_ann_afe[enr_ann_afe$AFE_72h_100_col != "lightgrey",cazy_cols])
cazy_df$h72 <- temp

temp <- colSums(enr_ann_afe[enr_ann_afe$AFE_48h_100_col != "lightgrey",cazy_cols])
cazy_df$h48 <- temp

temp <- colSums(enr_ann_afe[enr_ann_afe$AFE_168h_100_col != "lightgrey",cazy_cols])
cazy_df$h168 <- temp

temp <- colSums(enr_ann_afe[enr_ann_afe$AFE_48h_50_col != "lightgrey",cazy_cols])
cazy_df$h48_50 <- temp

temp <- colSums(enr_ann_afe[enr_ann_afe$AFE_168h_50_col != "lightgrey",cazy_cols])
cazy_df$h168_50 <- temp

cazy_df <- data.frame(cazy_df)
cazy_df$gene <- gsub("CAZy  ", "", gsub("\\.", " ", rownames(cazy_df)))
cazy_df_melt <- melt(cazy_df)
kruskal.test(value~variable, data=cazy_df_melt) # there is at least one statistically significant group
k24_48 <- kruskal.test(value~variable, data=cazy_df_melt[cazy_df_melt$variable %in% c("h24", "h48"),]) #p=0.02
k48_72 <- kruskal.test(value~variable, data=cazy_df_melt[cazy_df_melt$variable %in% c("h48", "h72"),]) #p=8.9e-05
k72_168 <- kruskal.test(value~variable, data=cazy_df_melt[cazy_df_melt$variable %in% c("h72", "h168"),]) #p=0.04
k48_168 <- kruskal.test(value~variable, data=cazy_df_melt[cazy_df_melt$variable %in% c("h48_50", "h168_50"),]) #p=0.32
k48 <- kruskal.test(value~variable, data=cazy_df_melt[cazy_df_melt$variable %in% c("h48", "h48_50"),]) #p=8.9e-05
k168 <- kruskal.test(value~variable, data=cazy_df_melt[cazy_df_melt$variable %in% c("h168", "h168_50"),]) #p=0.25
# Extracting p-values of all comparisons
p_vals <- c(k24_48$p.value, k48_72$p.value, k72_168$p.value, k48_168$p.value, k48$p.value, k168$p.value)
# Correcting for multiple comparisons with Benjamini-Hochberg
p_vals_BH <- p.adjust(p_vals, "BH", length(p_vals)) # 0.002, 0.00009, 0.04, 0.32, 0.00009, 0.25

# CAZy abundance barplots with manual colors - figure 6
enr_ann_afe <- enr_ann_afe %>%
  mutate(AFE_24h_100_col = ifelse(AFE_24h>0, "#b3cde3", "lightgrey")) %>%
  mutate(AFE_48h_100_col = ifelse(AFE_48h>0, "#8c96c6", "lightgrey")) %>%
  mutate(AFE_72h_100_col = ifelse(AFE_72h>0, "#8856a7", "lightgrey")) %>%
  mutate(AFE_168h_100_col = ifelse(AFE_168h>0, "#810f7c", "lightgrey")) %>%
  mutate(AFE_48h_50_col = ifelse(AFE_48h_50>0, "#8c96c6", "lightgrey")) %>%
  mutate(AFE_168h_50_col = ifelse(AFE_168h_50>0, "#810f7c", "lightgrey"))

temp1 <- enr_ann_afe %>%
  filter_at(vars(starts_with("AFE")), any_vars(. > 0)) %>%
  select(matches('CAZy|col$')) %>%
  pivot_longer(!ends_with("col"), names_to = "pathway", values_to = "count") %>%
  mutate(pathway = gsub("CAZy\\.\\.", "", pathway)) %>%
  mutate(pathway = gsub("\\.*$", "", pathway)) %>%
  mutate(pathway_short = case_when(
    pathway == "Alpha.galactans" ~ "Ag",
    pathway == "Alpha.mannan" ~ "Am",
    pathway == "Amorphous.Cellulose" ~ "Ac",
    pathway == "Arabinan" ~ "Ar",
    pathway == "Arabinose.cleavage" ~ "Ae",
    pathway == "Beta.galactan..pectic.galactan" ~ "Bg",
    pathway == "Beta.mannan" ~ "Bm",
    pathway == "Chitin" ~ "Ch",
    pathway == "Crystalline.Cellulose" ~ "Cc",
    pathway == "Fucose.Cleavage" ~ "Fu",
    pathway == "Mixed.Linkage.glucans" ~ "Ml",
    pathway == "Mucin" ~ "Mu",
    pathway == "Pectin" ~ "Pe",
    pathway == "Polyphenolics" ~ "Po",
    pathway == "Rhamnose.cleavage" ~ "Rh",
    pathway == "Starch" ~ "St",
    pathway == "Sulf.Polysachharides" ~ "Sp",
    pathway == "Xylans" ~ "Xs",
    pathway == "Xyloglucan" ~ "Xn"
  )) %>%
  mutate(pathway_short = factor(pathway_short, levels=c("Ch", "Am", "Ae", "Rh", "Fu", "Mu", "Po", "Ar", "Pe", "Xs", 
                                                        "Xn", "St", "Bm", "Cc", "Ac", "Bg", "Ag", "Ml", "Sp"))) %>%
  mutate(grp = factor(case_when(
    pathway_short %in% c("Am", "Ch") ~ "Fungi",
    pathway_short %in% c("Ae", "Rh", "Fu") ~ "Mono",
    pathway_short == "Mu" ~ "Mucin",
    pathway_short %in% c("Ag", "Ml", "Sp", "Bg") ~ "Poly",
    TRUE ~ "Plant"
  ), levels = c("Fungi", "Mono", "Mucin", "Plant", "Poly"))) %>%
  mutate(col = case_when(
    pathway_short == "Am" ~ "#FEE391",
    pathway_short == "Ch" ~ "#FEC44F",
    pathway_short == "Ae" ~ "#2171B5",
    pathway_short == "Rh" ~ "#6BAED6",
    pathway_short == "Fu" ~ "#BDD7E7",
    pathway_short == "Mu" ~ "#FF00FF",
    pathway_short == "Po" ~ "#00441B",
    pathway_short == "Ar" ~ "#006D2C",
    pathway_short == "Pe" ~ "#238B45",
    pathway_short == "Xs" ~ "#41AB5D",
    pathway_short == "Xn" ~ "#74C476",
    pathway_short == "St" ~ "#A1D99B",
    pathway_short == "Bm" ~ "#C7E9C0",
    pathway_short == "Cc" ~ "#E5F5E0",
    pathway_short == "Ac" ~ "#F7FCF5",
    pathway_short == "Bg" ~ "#A50F15",
    pathway_short == "Ag" ~ "#DE2D26",
    pathway_short == "Ml" ~ "#FB6A4A",
    pathway_short == "Sp" ~ "#FCAE91"
  ))


temp <- temp1 %>%
  filter(AFE_24h_100_col != "lightgrey") %>%
  group_by(pathway_short, grp, col) %>%
  summarise_at("count", sum) %>%
  mutate(count = count*100/nrow(enr_ann_afe[enr_ann_afe$AFE_24h_100_col != "lightgrey",])) %>%
  arrange(grp, pathway_short)
p24 <- ggplot(temp, aes(x=pathway_short, y=count, fill=col)) + 
  geom_bar(stat = "identity", color="black", width=0.8, position = position_dodge(width=0.2)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_identity() +
  ylim(c(0,41)) + xlab(NULL) + ylab(NULL) + ggtitle("24h")

temp <- temp1 %>%
  filter(AFE_48h_100_col != "lightgrey") %>%
  group_by(pathway_short, grp, col) %>%
  summarise_at("count", sum) %>%
  mutate(count = count*100/nrow(enr_ann_afe[enr_ann_afe$AFE_48h_100_col != "lightgrey",])) %>%
  arrange(grp, pathway_short)
p48 <- ggplot(temp, aes(x=pathway_short, y=count, fill=col)) + 
  geom_bar(stat = "identity", color="black", width=0.8, position = position_dodge(width=0.2)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_identity() +
  ylim(c(0,41)) + xlab(NULL) + ylab(NULL) + ggtitle("48h")

temp <- temp1 %>%
  filter(AFE_72h_100_col != "lightgrey") %>%
  group_by(pathway_short, grp, col) %>%
  summarise_at("count", sum) %>%
  mutate(count = count*100/nrow(enr_ann_afe[enr_ann_afe$AFE_72h_100_col != "lightgrey",])) %>%
  arrange(grp, pathway_short) 
p72 <- ggplot(temp, aes(x=pathway_short, y=count, fill=col)) + 
  geom_bar(stat = "identity", color="black", width=0.8, position = position_dodge(width=0.2)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_identity() +
  ylim(c(0,41)) + xlab(NULL) + ylab(NULL) + ggtitle("72h")

temp <- temp1 %>%
  filter(AFE_168h_100_col != "lightgrey") %>%
  group_by(pathway_short, grp, col) %>%
  summarise_at("count", sum) %>%
  mutate(count = count*100/nrow(enr_ann_afe[enr_ann_afe$AFE_168h_100_col != "lightgrey",])) %>%
  arrange(grp, pathway_short)
p168 <- ggplot(temp, aes(x=pathway_short, y=count, fill=col)) + 
  geom_bar(stat = "identity", color="black", width=0.8, position = position_dodge(width=0.2)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_identity() +
  ylim(c(0,41)) + xlab(NULL) + ylab(NULL) + ggtitle("168h")

temp <- temp1 %>%
  filter(AFE_48h_50_col != "lightgrey") %>%
  group_by(pathway_short, grp, col) %>%
  summarise_at("count", sum) %>%
  mutate(count = count*100/nrow(enr_ann_afe[enr_ann_afe$AFE_48h_50_col != "lightgrey",])) %>%
  arrange(grp, pathway_short) 
p48_50 <- ggplot(temp, aes(x=pathway_short, y=count, fill=col)) + 
  geom_bar(stat = "identity", color="black", width=0.8, position = position_dodge(width=0.2)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_identity() +
  ylim(c(0,41)) + xlab(NULL) + ylab(NULL) + ggtitle("48h 50%")

temp <- temp1 %>%
  filter(AFE_168h_50_col != "lightgrey") %>%
  group_by(pathway_short, grp, col) %>%
  summarise_at("count", sum) %>%
  mutate(count = count*100/nrow(enr_ann_afe[enr_ann_afe$AFE_168h_50_col != "lightgrey",])) %>%
  arrange(grp, pathway_short) 
p168_50 <- ggplot(temp, aes(x=pathway_short, y=count, fill=col)) + 
  geom_bar(stat = "identity", color="black", width=0.8, position = position_dodge(width=0.2)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_identity() +
  ylim(c(0,41)) + xlab(NULL) + ylab(NULL) + ggtitle("168h 50%")

ggarrange(p24, p72, p48, p168, p48_50, p168_50, ncol=2, nrow=3)
ggsave(filename="fig6_CAZy_manual_ord.pdf", width=7.5, height=6)

plot.new()
par(mfrow=c(1,1))
legend("topright", legend=temp$pathway_short[order(temp$grp, temp$pathway_short)], fill=temp$col[order(temp$grp, temp$pathway_short)], xpd=TRUE)

# CAZy counts and differential abundance in enriched, unenriched but detected and undetected MAGs - supplementary figure 4
# The DRAM output has names of CAZy genes that I can pull from betas rows
enr_ann <- read.delim("DRAM/distill_MAGs_enr_all.tsv", sep="\t", header=T, stringsAsFactors = F)
unenr_ann <- read.delim("DRAM/detected_MAGs_product.tsv", sep="\t", header=T, stringsAsFactors = F)
all_ann <- read.delim("DRAM/all_MAGs_product.tsv", sep="\t", header=T, stringsAsFactors = F)
cazy_cols <- which(colnames(enr_ann) %like% "CAZ")
enr_ann[enr_ann == "True"] <- 1
enr_ann[enr_ann == "False"] <- 0
unenr_ann[unenr_ann == "True"] <- 1
unenr_ann[unenr_ann == "False"] <- 0
all_ann[all_ann == "True"] <- 1
all_ann[all_ann == "False"] <- 0
enr_ann[,cazy_cols] <- apply(enr_ann[,cazy_cols], 1, as.numeric)
unenr_ann[,cazy_cols] <- apply(unenr_ann[,cazy_cols], 1, as.numeric)
all_ann[,cazy_cols] <- apply(all_ann[,cazy_cols], 1, as.numeric)

cazy_df <- data.frame(rbind(colSums(enr_ann[,cazy_cols]), 
                            colSums(unenr_ann[,cazy_cols]), 
                            colSums(all_ann[,cazy_cols])))
rownames(cazy_df) <- c("enr", "unenr", "all")
colnames(cazy_df) <- gsub("CAZy  ", "", gsub("\\.", " ", colnames(cazy_df)))
wilcox.test(t(cazy_df[1,]), t(cazy_df[2,]), paired=TRUE) 
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(cazy_df[c(1,3),order(-cazy_df[3,])]), las=2, beside=T, xlab="Substrate targetted", ylab="Number of MAGs with function")

# Supplementary figure 3 - CAZy per MAG

# Compare completeness between growers and detected
mag_info <- read.csv("genomeInformation_derep.csv", stringsAsFactors = F)
mag_info$MAG <- gsub("\\.fa", "", mag_info$MAG)
enr_mags <- read.csv("enr_MAGs_all_tax.csv")[,-1] %>% select("Bin")
unenr_mags <- unenr_ann$genome
enr_comp <- mag_info$completeness[mag_info$MAG %in% enr_mags$Bin]
unenr_comp <- mag_info$completeness[mag_info$MAG %in% unenr_mags]
t.test(enr_comp, unenr_comp)

# Compare completeness between fast and slow growers
fast_grow <- scan("../Traits_paper/fast_growers.txt", sep="\n", what = "charater")
slow_grow <- enr_mags$Bin[!(enr_mags$Bin %in% fast_grow)]
fast_comp <- mag_info$completeness[mag_info$MAG %in% fast_grow]
slow_comp <- mag_info$completeness[mag_info$MAG %in% slow_grow]
t.test(fast_comp, slow_comp)

all_cazy <- read.csv("DRAM/all_MAGs_metabolism_summary_cazy.csv") %>%
  filter(header=="CAZY") %>%
  select(!c("gene_description", "module", "header", "subheader", "fasta")) %>%
  column_to_rownames("gene_id") %>%
  t() %>% data.frame() %>% rownames_to_column("genome")
all_cazy$genome <- gsub("^X", "", all_cazy$genome)
all_cazy$genome <- gsub("_H.*", "", all_cazy$genome)
fast_grow <- gsub("-", "\\.", fast_grow)
slow_grow <- gsub("-", "\\.", slow_grow)

fast_grow_cazy <- all_cazy %>%
  filter(genome %in% fast_grow) %>%
  column_to_rownames("genome") 
fast_grow_cazy[fast_grow_cazy != ""] <- 1
fast_grow_cazy[fast_grow_cazy == ""] <- 0 
fast_grow_cazy <- fast_grow_cazy %>% mutate_all(as.numeric)

slow_grow_cazy <- all_cazy %>%
  filter(genome %in% slow_grow) %>%
  column_to_rownames("genome") 
slow_grow_cazy[slow_grow_cazy != ""] <- 1
slow_grow_cazy[slow_grow_cazy == ""] <- 0 
slow_grow_cazy <- slow_grow_cazy %>% mutate_all(as.numeric)

fgc_sum <- rowSums(fast_grow_cazy)
fsc_sum <- rowSums(slow_grow_cazy)
t.test(fgc_sum, fsc_sum, alternative="less") #p=0.02, t=-2
mean(fgc_sum) #33
mean(fsc_sum) #41

# CAZy differential abundance from enriched MAGs
betas_CAZy_enr <- betas_CAZy[betas_CAZy$Bin %in% enr_MAGs$enr_MAGs_,] #817 genes
# Get rid of central carbon cycle genes
unique(betas_CAZy_enr$subheader)
betas_CAZy_enr <- betas_CAZy_enr[!(betas_CAZy_enr$subheader %in% c("glycolysis", "TCA", "", "pentose pathway", "galacturonic acid degradation")),]
plant <- c("Xylan", "Mixed-Linkage glucans", "Beta-mannan", "Xyloglucan", "Polyphenolics", "Arabinan", "Pectin", "Amorphous Cellulose", "Starch", "Alpha-mannan")
my_palette <- colorRampPalette(c("blue","lightblue", "white","lightgoldenrod1","red"))(n = 100)
for (i in 1:length(plant)) {
  temp <- betas_CAZy_enr[betas_CAZy_enr$subheader %like% plant[i],]
  for (j in 3:11) {
    temp[is.na(temp[,j+10]), j] <- 0
  }
  print(plant[i])
  print(temp$variable[temp$betas_48h_100_vs_24h_100>0])
  if (nrow(temp) > 1) {
    pdf(paste0("../Traits_paper/CAZy_enr_betas/sigonly_", plant[i], ".pdf"))
    heatmap.2(as.matrix(unique(temp[, 3:11])), trace="none", density.info="none", Colv=FALSE, col=my_palette ,labRow = temp$Bin, cexRow = 0.5, main=plant[i])
    dev.off()
  }
}
temp <- betas_CAZy_enr[betas_CAZy_enr$subheader %like% "Chitin",]
for (j in 3:11) {
  temp[is.na(temp[,j+10]), j] <- 0
}
temp <- betas_CAZy_enr[betas_CAZy_enr$subheader %like% "Mucin",]
for (j in 3:11) {
  temp[is.na(temp[,j+10]), j] <- 0
}
mono <- c("Arabinose", "Fucose", "Rhamnose", "Chitin", "Mucin")
for (i in 1:length(mono)) {
  temp <- betas_CAZy_enr[betas_CAZy_enr$subheader %like% mono[i],]
  for (j in 3:11) {
    temp[is.na(temp[,j+10]), j] <- 0
  }
  print(mono[i])
  print(temp$variable[temp$betas_48h_100_vs_24h_100>0])
  if (nrow(temp) > 1) {
    pdf(paste0("../Traits_paper/CAZy_enr_betas/sigonly_", mono[i], ".pdf"))
    heatmap.2(as.matrix(unique(temp[, 3:11])), trace="none", density.info="none", Colv=FALSE, col=my_palette ,labRow = temp$Bin, cexRow = 0.5, main=mono[i])
    dev.off()
  }
}
poly <- c("Sulf", "Alpha-galactans", "Beta-galactans")
for (i in 1:length(poly)) {
  temp <- betas_CAZy_enr[betas_CAZy_enr$subheader %like% poly[i],]
  for (j in 3:11) {
    temp[is.na(temp[,j+10]), j] <- 0
  }
  print(poly[i])
  print(temp$variable[temp$betas_48h_100_vs_24h_100>0])
  if (nrow(temp) > 1) {
    pdf(paste0("../Traits_paper/CAZy_enr_betas/sigonly_", mono[i], ".pdf"))
    heatmap.2(as.matrix(unique(temp[, 3:11])), trace="none", density.info="none", Colv=FALSE, col=my_palette ,labRow = temp$Bin, cexRow = 0.5, main=poly[i])
    dev.off()
  }
}

# N cycling gene counts and differential abundance in enriched, unenriched but detected and undetected MAGs
# The DRAM output has names of CAZy genes that I can pull from betas rows 435...
enr_ann <- read.delim("DRAM/enr_MAGs_ann.tsv", sep="\t", header=T, stringsAsFactors = F)
#unenr_ann <- read.delim("DRAM/unenr_MAGs_product.tsv", sep="\t", header=T, stringsAsFactors = F)
all_ann <- read.delim("DRAM/annotations.tsv", sep="\t", header=T, stringsAsFactors = F)
pep_cols <- which(colnames(enr_ann) %like% "pept")
kegg_cols <- which(colnames(enr_ann) %like% "kegg")
enr_ann_N <- enr_ann[enr_ann$kegg_id %in% nitrogen_ko | enr_ann$cazy_hits %like% "chitinase",]
enr_ann_N <- merge(enr_ann_N, ko_to_gene_N, by.x="kegg_id", by.y="KO", all.x=T)
enr_ann_N$gene[enr_ann_N$kegg_id == ""] <- "Chitinase"
enr_ann_pep <- enr_ann[enr_ann$peptidase_id!="" & !(enr_ann$kegg_id %in% nitrogen_ko),]
all_ann_N <- all_ann[all_ann$kegg_id %in% nitrogen_ko | all_ann$cazy_hits %like% "chitinase",]
all_ann_N <- merge(all_ann_N, ko_to_gene_N, by.x="kegg_id", by.y="KO", all.x=T)
all_ann_N$gene[all_ann_N$kegg_id == ""] <- "Chitinase"
all_ann_pep <- all_ann[all_ann$peptidase_id!="" & !(all_ann$kegg_id %in% nitrogen_ko),]
inorN <- as.character(unique(all_ann_N$gene))
inorN <- inorN[which(!(inorN %in% c("lys", "ureA", "ureB", "ureC", "Xds", "chit1", "putative_chitinase")))]
N_df <- data.frame(matrix(ncol=length(inorN), nrow=2))
rownames(N_df) <- c("enr", "all")
colnames(N_df) <- inorN
# collapse multiple gene copies in the same genomes
inorN_enr <- enr_ann_N[enr_ann_N$gene %in% inorN, c(1,3,10,23:25)]
inorN_enr <- unique(inorN_enr)
inorN_all <- all_ann_N[all_ann_N$gene %in% inorN, c(1,3,10,23:25)]
inorN_all <- unique(inorN_all)
for (i in 1:ncol(N_df)) {
  N_df[1, i] <- nrow(inorN_enr[inorN_enr$gene == colnames(N_df[i]),])/1.14 #calculate % of enriched MAGs that have the gene at least once
  N_df[2, i] <- nrow(inorN_all[inorN_all$gene == colnames(N_df[i]),])/5.03 #calculate % of MAGs that have the gene at least once
}
t.test(t(N_df[1,]), t(N_df[2,]), paired=TRUE) 
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(N_df[,order(-N_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("Growing organisms", "Whole community"), fill=c("black", "lightgrey"))

orN <- c("lys", "ureA", "ureB", "ureC", "Xds", "chit1", "putative_chitinase", "Chitinase")
N_df <- data.frame(matrix(ncol=length(orN), nrow=2))
rownames(N_df) <- c("enr", "all")
colnames(N_df) <- orN
# collapse multiple gene copies in the same genomes
orN_enr <- enr_ann_N[enr_ann_N$gene %in% orN, c(1,3,10,23:25)]
orN_enr <- unique(orN_enr)
orN_all <- all_ann_N[all_ann_N$gene %in% orN, c(1,3,10,23:25)]
orN_all <- unique(orN_all)
for (i in 1:ncol(N_df)) {
  N_df[1, i] <- nrow(orN_enr[orN_enr$gene == colnames(N_df[i]),])/1.14 #calculate % of enriched MAGs that have the gene at least once
  N_df[2, i] <- nrow(orN_all[orN_all$gene == colnames(N_df[i]),])/5.03 #calculate % of MAGs that have the gene at least once
}
t.test(t(N_df[1,]), t(N_df[2,]), paired=TRUE) 
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(N_df[,order(-N_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("Growing organisms", "Whole community"), fill=c("black", "lightgrey"))

# Number of extracellular peptidase genes per genome - because almost all genomes have a peptidase
length(unique(enr_ann_pep$fasta)) #446/503 and 113/114
enr_ann_pep$peptidase_grp <- substring(enr_ann_pep$peptidase_family, 1, 1)
all_ann_pep$peptidase_grp <- substring(all_ann_pep$peptidase_family, 1, 1)
pep_per_mag_enr <- enr_ann_pep %>% group_by(fasta, peptidase_grp) %>% count()
pep_per_mag_all <- all_ann_pep %>% group_by(fasta, peptidase_grp) %>% count()
pep_per_mag_enr_grouped <- pep_per_mag_enr %>% group_by(peptidase_grp) %>% summarize_at("freq", sum)
pep_per_mag_all_grouped <- pep_per_mag_all %>% group_by(peptidase_grp) %>% summarize_at("freq", sum)
pep_df <- merge(pep_per_mag_enr_grouped, pep_per_mag_all_grouped, by="peptidase_grp", all=T)
pep_df[is.na(pep_df)] <- 0
pep_df <- pep_df[-which(pep_df$peptidase_grp %in% c("p", "I")),]
colnames(pep_df)[2:3] <- c("n_enr", "n_all")
pep_df$n_enr <- pep_df$n_enr/114
pep_df$n_all <- pep_df$n_all/503
t.test(pep_df$n_enr, pep_df$n_all, paired=TRUE) 
plot.new()
par(mfrow=c(1,1))
pep_df <- pep_df[order(-pep_df$n_all),]
rownames(pep_df) <- pep_df$peptidase_grp
barplot(as.matrix(t(pep_df[,-1])), las=2, beside=T, xlab="Peptidase Group", ylab="Mean number of peptidates per MAGs")
legend("topright", legend=c("Growing organisms", "Whole community"), fill=c("black", "lightgrey"))

# Replicate Alex G's figure 3 (also Nunan 2020)
# requires # of MAGs with specific pathway based on DRAM out of total MAGs growing, 
# mean AFE of these MAGs per time point
# Actually this is pretty much the comparison I already did except I used growing/not growing as a binary operator instead of AFE
dram_misc <- read.csv("DRAM/all_MAGs_metabolism_summary_misc.csv", header=T, stringsAsFactors = F)
dram_pep <- read.csv("DRAM/all_MAGs_metabolism_summary_pep.csv", header=T, stringsAsFactors = F)
dram_cazy <- read.csv("DRAM/all_MAGs_metabolism_summary_cazy.csv", header=T, stringsAsFactors = F)
dram_trans <- read.csv("DRAM/all_MAGs_metabolism_summary_trans.csv", header=T, stringsAsFactors = F)
dram_res <- data.frame(rbind(dram_misc, dram_pep, dram_cazy, dram_trans))
dram_res <- dram_res[, -7] # removing column "fasta"
afe[afe < 0]<- 0
dram_res[dram_res == ""] <- 0
temp <- dram_res[, -(1:5)]
temp[temp != 0] <- 1
temp <- data.frame(t(apply(temp, 1, as.numeric)))
dram_res[, 6:ncol(dram_res)] <- temp
enr_t_org <- as.character(afe$bin)
enr_t <- gsub("-", "\\.", enr_t_org)
enr_t <- gsub("^1", "X1", enr_t)
dram_res <- dram_res[,c(1:5, which(colnames(dram_res) %in% enr_t))] # Keep only traits of MAGs that grow at some time point
dram_res <- dram_res[which(rowSums(dram_res[, 6:ncol(dram_res)]) > 0),] # keep only traits represented in MAGs
# I'm now left with 1523 traits in 97 MAGs

traits_afe <- data.frame(matrix(ncol=17, nrow=nrow(dram_res)))
traits_afe[, 1:5] <- dram_res[, 1:5]
colnames(traits_afe)[1:5] <- colnames(dram_res)[1:5]
colnames(traits_afe)[6:ncol(traits_afe)] <- c(gsub("AFE", "fMAGs", colnames(afe)[-(1:2)]), colnames(afe)[-(1:2)])

ks <- c()
# Get fraction of growing MAGs that have a trait
for (c in 3:(ncol(afe))) {
  enr_t_org <- afe$bin[afe[,c] > 0] # which MAGs are enriched at a time point
  enr_t <- gsub("-", "\\.", enr_t_org) # MAG names are slightly different in the dram results
  enr_t <- gsub("^1", "X1", enr_t)
  traits_t <- dram_res[, which(colnames(dram_res) %in% enr_t)] # get traits for the growing MAGs at the time point
  traits_afe[, c+3] <- rowSums(traits_t)/length(enr_t) # get their fraction out of total MAGs growing at that time point
  # Get the mean AFE of MAGs enriched at this time point and also possessing a trait per trait
  mags_w_trait <- apply(traits_t, 1, function(x) {which(x == 1)})
  l <- lapply(mags_w_trait, function (x) {which(!is.na(match(enr_t, names(x))))})
  # Run a Kolmogorov-Smirnoff test comparing the afe distriubtions of growing MAGs with a trait and all growing MAGs per time point
  ks[[c-1]] <- sapply(l, function (x) {if (length(x) > 3) ks.test(afe[afe$bin %in% enr_t_org[x], c], afe[,c])$p.value})
  traits_afe[, c+9] <- data.frame(sapply(l, function(x) {mean(afe[afe$bin %in% enr_t_org[x], c])}))
}
traits_afe[traits_afe == "NaN"] <- NA

# Find significant traits per time point
ks_sig <- data.frame(matrix(nrow=length(ks[[1]]), ncol=(length(ks)-1)))
for (i in 1:(length(ks)-1)) {
  for (j in 1:length(ks[[i+1]])) {
    if (!(is_null(ks[[i+1]][[j]]))) {
      ks_sig[j, i] <- ks[[i+1]][[j]]
    }
    else ks_sig[j, i] <- NA
  }
}
# Adjusting p-values for multiple comparisons with False Discovery Rate
ks_sig <- data.frame(apply(ks_sig, 2, function (x) {p.adjust(x, method="fdr", length(x))}))
ks_sig <- data.frame(dram_res[, 1:5], ks_sig)
colnames(ks_sig)[6:11] <- paste0("pval_", colnames(afe)[3:8])
ks_sig$uniq_id <- 1:nrow(ks_sig)
# Throwing out traits with all NAs - left with 1148 out of 1496
ks_sig1 <- ks_sig[!is.na(ks_sig$pval_AFE_24h) | !is.na(ks_sig$pval_AFE_48h) | !is.na(ks_sig$pval_AFE_72h) | !is.na(ks_sig$pval_AFE_168h) | !is.na(ks_sig$pval_AFE_48h_50) | !is.na(ks_sig$pval_AFE_168h_50), ]
# Throwing out rows where all p values are > 0.05 - left with 1007
ks_sig_only <- ks_sig1[ks_sig1$pval_AFE_24h <= 0.05 | ks_sig1$pval_AFE_48h <= 0.05 | ks_sig1$pval_AFE_72h <= 0.05 | ks_sig1$pval_AFE_168h <= 0.05 | ks_sig1$pval_AFE_48h_50 <= 0.05 | ks_sig1$pval_AFE_168h_50 <= 0.05, ]
ks_sig_only <- merge(ks_sig_only, traits_afe[, -(1:5)], by.x="uniq_id", by.y="row.names")
write.csv(ks_sig_only, "../Traits_paper/Traits_by_AFE_AlexStyle.csv")

nrow(ks_sig_only[ks_sig_only$header=="CAZY",]) #107 CAZymes with significantly different afe distribution
temp <- ks_sig_only[ks_sig_only$header=="CAZY",]
temp <- temp[order(temp$subheader),]
temp[is.na(temp)] <- 1
# Changing p-values to be binary: 0 is under 0.05, 1 is over 0.05
temp$pval_AFE_24h[temp$pval_AFE_24h < 0.05] <- 0
temp$pval_AFE_24h[temp$pval_AFE_24h >= 0.05] <- 1
temp$pval_AFE_48h[temp$pval_AFE_48h < 0.05] <- 0
temp$pval_AFE_48h[temp$pval_AFE_48h >= 0.05] <- 1
temp$pval_AFE_72h[temp$pval_AFE_72h < 0.05] <- 0
temp$pval_AFE_72h[temp$pval_AFE_72h >= 0.05] <- 1
temp$pval_AFE_168h[temp$pval_AFE_168h < 0.05] <- 0
temp$pval_AFE_168h[temp$pval_AFE_168h >= 0.05] <- 1
temp$pval_AFE_48h_50[temp$pval_AFE_48h_50 < 0.05] <- 0
temp$pval_AFE_48h_50[temp$pval_AFE_48h_50 >= 0.05] <- 1
temp$pval_AFE_168h_50[temp$pval_AFE_168h_50 < 0.05] <- 0
temp$pval_AFE_168h_50[temp$pval_AFE_168h_50 >= 0.05] <- 1

# Bubble plot: each row is a trait, columns are samples, size by relative abundance, color by AFE, shape by p-value?
ggplot(data=temp, aes(x="24h", y=gene_id)) + geom_point(aes(size=fMAGs_24h, color=AFE_24h)) + 
  scale_color_gradient(low="white", high="purple") + theme_classic()
ggplot(data=temp, aes(x="24h", y=gene_id)) + geom_point(aes(size=fMAGs_48h, color=AFE_48h)) + 
  scale_color_gradient(low="white", high="purple") + theme_classic()

# Trehalose synthesis gene counts and differential abundance in enriched, unenriched but detected and undetected MAGs
trehalose <- c("K00697", "K16055", "K23377", "K02777", "K02819", "K01087", "K13057", "K01236", "K05343")
enr_ann <- read.delim("DRAM/enr_MAGs_ann.tsv", sep="\t", header=T, stringsAsFactors = F)
all_ann <- read.delim("DRAM/annotations.tsv", sep="\t", header=T, stringsAsFactors = F)
kegg_cols <- which(colnames(enr_ann) %like% "kegg")
enr_ann_T <- enr_ann[enr_ann$kegg_id %in% trehalose,]
all_ann_T <- all_ann[all_ann$kegg_id %in% trehalose,]
T_df <- data.frame(matrix(ncol=length(trehalose), nrow=2))
rownames(T_df) <- c("enr", "all")
colnames(T_df) <- trehalose
# collapse multiple gene copies in the same genomes
T_enr <- enr_ann_T[enr_ann_T$kegg_id %in% trehalose, c(2,9,10)]
T_enr <- unique(T_enr)
T_enr_50 <- T_enr[T_enr$fasta %in% afe50$bin,]
T_enr_100 <- T_enr[T_enr$fasta %in% afe100$bin,]
T_all <- all_ann_T[all_ann_T$kegg_id %in% trehalose, c(2,9,10)]
T_all <- unique(T_all)
for (i in 1:ncol(T_df)) {
  T_df[1, i] <- nrow(T_enr_50[T_enr_50$kegg_id == colnames(T_df[i]),])/(nrow(T_enr_50)/100) #calculate % of enriched MAGs that have the gene at least once
  T_df[2, i] <- nrow(T_enr_100[T_enr_100$kegg_id == colnames(T_df[i]),])/(nrow(T_enr_100)/100) #calculate % of MAGs that have the gene at least once
}
colnames(T_df)
t.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE, alternative = "greater") 
wilcox.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE, alternative = "greater")
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(T_df[,order(-T_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("Growing at 50%", "Growing at 100%"), fill=c("black", "lightgrey"))

# EPS synthesis gene counts and differential abundance in enriched, unenriched but detected and undetected MAGs
eps <- c("PF02563", "PF01451", "PF02706", "PF13807", "PF13614", "PF04932",  "PF01270", "PF06055", 
           "PF01943", "PF14667", "PF13440", "PF01554", "PF01061", "PF12698", "PF13414", "PF14524", 
           "PF05159", "PF02348", "PF01380", "PF00571", "PF13641", "PF07238", "PF13437", "PF08238", 
           "PF13372", "PF05048", "PF03062", "PF05426", "PF03170", "PF13432", "PF14559")
enr_ann$pfam <- gsub(".*\\[", "", gsub("].*", "", gsub("\\..*","", enr_ann$pfam_hits)))
all_ann$pfam <- gsub(".*\\[", "", gsub("].*", "", gsub("\\..*","", all_ann$pfam_hits)))
enr_ann_T <- enr_ann[enr_ann$pfam %in% eps,]
all_ann_T <- all_ann[all_ann$pfam %in% eps,]
T_df <- data.frame(matrix(ncol=length(eps), nrow=2))
rownames(T_df) <- c("enr", "all")
colnames(T_df) <- eps
# collapse multiple gene copies in the same genomes
T_enr <- enr_ann_T[enr_ann_T$pfam %in% eps, c(2,23)]
T_enr <- unique(T_enr)
T_all <- all_ann_T[all_ann_T$pfam %in% eps, c(2,23)]
T_all <- unique(T_all)
for (i in 1:ncol(T_df)) {
  T_df[1, i] <- nrow(T_enr[T_enr$pfam == colnames(T_df[i]),])/1.13 #calculate % of enriched MAGs that have the gene at least once
  T_df[2, i] <- nrow(T_all[T_all$pfam == colnames(T_df[i]),])/5.03 #calculate % of MAGs that have the gene at least once
}
t.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE, alternative="greater") 
wilcox.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE, alternative="greater")
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(T_df[,order(-T_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("Growing organisms", "Whole community"), fill=c("black", "lightgrey"))
T_enr_50 <- T_enr[T_enr$fasta %in% afe50$bin,]
T_enr_100 <- T_enr[T_enr$fasta %in% afe100$bin,]
for (i in 1:ncol(T_df)) {
  T_df[1, i] <- nrow(T_enr_50[T_enr_50$pfam == colnames(T_df[i]),])/(nrow(T_enr_50)/100) #calculate % of enriched MAGs that have the gene at least once
  T_df[2, i] <- nrow(T_enr_100[T_enr_100$pfam == colnames(T_df[i]),])/(nrow(T_enr_100)/100) #calculate % of MAGs that have the gene at least once
}
t.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE) 
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(T_df[,order(-T_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("50%", "100%"), fill=c("black", "lightgrey"))
# Compare time points
T_df <- data.frame(matrix(ncol=length(eps), nrow=6))
rownames(T_df) <- c("24h", "48h", "72h", "168h", "48h_50", "168h_50")
colnames(T_df) <- eps
# Get the gene content for MAGs that are actually growing at each time point
T_24_100 <- T_enr[T_enr$fasta %in% afe100$bin[afe100$AFE_24h>0],]
T_48_100 <- T_enr[T_enr$fasta %in% afe100$bin[afe100$AFE_48h>0],]
T_72_100 <- T_enr[T_enr$fasta %in% afe100$bin[afe100$AFE_72h>0],]
T_168_100 <- T_enr[T_enr$fasta %in% afe100$bin[afe100$AFE_168h>0],]
T_48_50 <- T_enr[T_enr$fasta %in% afe50$bin[afe100$AFE_48h>0],]
T_168_50 <- T_enr[T_enr$fasta %in% afe50$bin[afe100$AFE_168h>0],]
for (i in 1:ncol(T_df)) {
  T_df[1, i] <- nrow(T_24_100[T_24_100$pfam == colnames(T_df[i]),])/(length(unique(T_24_100$fasta))/100) #calculate % of enriched MAGs that have the gene at least once at 24h 100%
  T_df[2, i] <- nrow(T_48_100[T_48_100$pfam == colnames(T_df[i]),])/(length(unique(T_48_100$fasta))/100) #calculate % of enriched MAGs that have the gene at least once at 48h 100%
  T_df[3, i] <- nrow(T_72_100[T_72_100$pfam == colnames(T_df[i]),])/(length(unique(T_72_100$fasta))/100) #calculate % of enriched MAGs that have the gene at least once at 72h 100%
  T_df[4, i] <- nrow(T_168_100[T_168_100$pfam == colnames(T_df[i]),])/(length(unique(T_168_100$fasta))/100) #calculate % of enriched MAGs that have the gene at least once at 168h 100%
  T_df[5, i] <- nrow(T_48_50[T_48_50$pfam == colnames(T_df[i]),])/(length(unique(T_48_50$fasta))/100) #calculate % of enriched MAGs that have the gene at least once at 48h 50%
  T_df[6, i] <- nrow(T_168_50[T_168_50$pfam == colnames(T_df[i]),])/(length(unique(T_168_50$fasta))/100) #calculate % of enriched MAGs that have the gene at least once at 168h 50%
}
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(T_df[,order(-T_df[1,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function", col=brewer.pal(6, "Paired"))
legend("topright", legend=c("Growing at 24h 100%", "Growing at 48h 100%", "Growing at 72h 100%", 
                            "Growing at 168h 100%", "Growing at 48h 50%", "Growing at 168h 50%"), 
       fill=brewer.pal(6, "Paired"))
eps_melt <- melt(t(T_df))
colnames(eps_melt) <- c("PFAM", "Time", "pct_MAGs")
eps_melt2 <- eps_melt[eps_melt$PFAM %in% c("PF13614", "PF01061", "PF13641", "PF00571", "PF01380"),]
ggboxplot(data=eps_melt2, x="Time", y="pct_MAGs", fill="Time")
eps_aov <- anova_test(pct_MAGs~Time, data=eps_melt)
summary(eps_aov) # Not significant
coefficients(eps_aov)
eps_aov <- aov(pct_MAGs~Time, data=eps_melt)
summary(eps_aov) # Significant, not sure why...is it because now it's doing paired comparisons?
# Maybe try in MAGs that are detected at 0h vs. whole community?
mags_t0 <- rpkm[rpkm$X12C16O_0h_100_P08>0 | rpkm$X12C16O_0h_100_P14>0 | rpkm$X12C16O_0h_100_P16>0, "Genome"]
t0_ann_T <- all_ann[all_ann$fasta %in% mags_t0 & all_ann$pfam %in% eps,]
all_ann_T <- all_ann[all_ann$pfam %in% eps,]
T_df <- data.frame(matrix(ncol=length(eps), nrow=2))
rownames(T_df) <- c("t0", "all")
colnames(T_df) <- eps
# collapse multiple gene copies in the same genomes
T_enr <- t0_ann_T[t0_ann_T$pfam %in% eps, c(2,23)]
T_enr <- unique(T_enr)
T_all <- all_ann_T[all_ann_T$pfam %in% eps, c(2,23)]
T_all <- unique(T_all)
for (i in 1:ncol(T_df)) {
  T_df[1, i] <- nrow(T_enr[T_enr$pfam == colnames(T_df[i]),])/(nrow(T_enr)/100) #calculate % of enriched MAGs that have the gene at least once
  T_df[2, i] <- nrow(T_all[T_all$pfam == colnames(T_df[i]),])/(nrow(T_all)/100) #calculate % of MAGs that have the gene at least once
}
t.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE) 
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(T_df[,order(-T_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("Organisms detected at 0h", "Whole community"), fill=c("black", "lightgrey"))
# Following Kegg EPS biosynthesis pathway (two component system map)
eps <- c("K19667", "K01991", "K25307", "K16692", "K01791", "K02472", "K16712", "K16713")
eps_gene <- c("XpsR", "EpsA", "EpsP", "EpsB", "EpsC", "EpsD", "EpsE", "EpsF")
enr_ann_T <- enr_ann[enr_ann$kegg_id %in% eps,]
all_ann_T <- all_ann[all_ann$kegg_id %in% eps,]
T_df <- data.frame(matrix(ncol=length(eps), nrow=2))
rownames(T_df) <- c("enr", "all")
colnames(T_df) <- eps_gene
# collapse multiple gene copies in the same genomes
T_enr <- enr_ann_T[enr_ann_T$kegg_id %in% eps, c(2,9,10)]
T_enr <- unique(T_enr)
T_all <- all_ann_T[all_ann_T$kegg_id %in% eps, c(2,9,10)]
T_all <- unique(T_all)
for (i in 1:length(eps)) {
  T_df[1, i] <- nrow(T_enr[T_enr$kegg_id == eps[i],])/1.14 #calculate % of enriched MAGs that have the gene at least once
  T_df[2, i] <- nrow(T_all[T_all$kegg_id == eps[i],])/5.03 #calculate % of MAGs that have the gene at least once
}
t.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE, alternative="greater")
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(T_df[,order(-T_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("Growing organisms", "Whole community"), fill=c("black", "lightgrey"))
T_enr_50 <- T_enr[T_enr$fasta %in% afe50$bin,]
T_enr_100 <- T_enr[T_enr$fasta %in% afe100$bin,]
for (i in 1:ncol(T_df)) {
  T_df[1, i] <- nrow(T_enr_50[T_enr_50$kegg_id == eps[i],])/(nrow(T_enr_50)/100) #calculate % of enriched MAGs that have the gene at least once
  T_df[2, i] <- nrow(T_enr_100[T_enr_100$kegg_id == eps[i],])/(nrow(T_enr_100)/100) #calculate % of MAGs that have the gene at least once
}
t.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE) 
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(T_df[,order(-T_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("50%", "100%"), fill=c("black", "lightgrey"))
# Maybe try in MAGs that are detected at 0h vs. whole community?
mags_t0 <- rpkm[rpkm$X12C16O_0h_100_P08>0 | rpkm$X12C16O_0h_100_P14>0 | rpkm$X12C16O_0h_100_P16>0, "Genome"]
t0_ann_T <- all_ann[all_ann$fasta %in% mags_t0 & all_ann$kegg_id %in% eps,]
all_ann_T <- all_ann[all_ann$kegg_id %in% eps,]
T_df <- data.frame(matrix(ncol=length(eps), nrow=2))
rownames(T_df) <- c("t0", "all")
colnames(T_df) <- eps_gene
# collapse multiple gene copies in the same genomes
T_enr <- t0_ann_T[t0_ann_T$kegg_id %in% eps, c(2,9,10)]
T_enr <- unique(T_enr)
T_all <- all_ann_T[all_ann_T$kegg_id %in% eps, c(2,9,10)]
T_all <- unique(T_all)
for (i in 1:ncol(T_df)) {
  T_df[1, i] <- nrow(T_enr[T_enr$kegg_id == eps[i],])/(nrow(T_enr)/100) #calculate % of enriched MAGs that have the gene at least once
  T_df[2, i] <- nrow(T_all[T_all$kegg_id == eps[i],])/(nrow(T_all)/100) #calculate % of MAGs that have the gene at least once
}
t.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE) 
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(T_df[,order(-T_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("Organisms detected at 0h", "Whole community"), fill=c("black", "lightgrey"))

# Trehalose synthesis gene counts and differential abundance in enriched, unenriched but detected and undetected MAGs
trehalose <- c("K00697", "K16055", "K23377", "K02777", "K02819", "K01087", "K13057", "K01236", "K05343")
enr_ann <- read.delim("DRAM/enr_MAGs_ann.tsv", sep="\t", header=T, stringsAsFactors = F)
all_ann <- read.delim("DRAM/annotations.tsv", sep="\t", header=T, stringsAsFactors = F)
kegg_cols <- which(colnames(enr_ann) %like% "kegg")
enr_ann_T <- enr_ann[enr_ann$kegg_id %in% trehalose,]
all_ann_T <- all_ann[all_ann$kegg_id %in% trehalose,]
T_df <- data.frame(matrix(ncol=length(trehalose), nrow=2))
rownames(T_df) <- c("enr", "all")
colnames(T_df) <- trehalose
# collapse multiple gene copies in the same genomes
T_enr <- enr_ann_T[enr_ann_T$kegg_id %in% trehalose, c(2,9,10)]
T_enr <- unique(T_enr)
T_enr_50 <- T_enr[T_enr$fasta %in% afe50$bin,]
T_enr_100 <- T_enr[T_enr$fasta %in% afe100$bin,]
T_all <- all_ann_T[all_ann_T$kegg_id %in% trehalose, c(2,9,10)]
T_all <- unique(T_all)
# Compare 50% to 100%
for (i in 1:ncol(T_df)) {
  T_df[1, i] <- nrow(T_enr_50[T_enr_50$kegg_id == colnames(T_df[i]),])/(nrow(T_enr_50)/100) #calculate % of enriched MAGs that have the gene at least once
  T_df[2, i] <- nrow(T_enr_100[T_enr_100$kegg_id == colnames(T_df[i]),])/(nrow(T_enr_100)/100) #calculate % of MAGs that have the gene at least once
}
colnames(T_df)
t.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE, alternative = "greater") 
wilcox.test(t(T_df[1,]), t(T_df[2,]), paired=TRUE, alternative = "greater")
plot.new()
par(mfrow=c(1,1))
barplot(as.matrix(T_df[,order(-T_df[2,])]), las=2, beside=T, xlab="Gene name", ylab="% of MAGs with function")
legend("topright", legend=c("Growing at 50%", "Growing at 100%"), fill=c("black", "lightgrey"))

# Diversity measures for 50% vs. 100%
enr_ann_50 <- enr_ann[enr_ann$fasta %in% afe50$bin,]
enr_ann_100 <- enr_ann[enr_ann$fasta %in% afe100$bin,]
length(unique(enr_ann_50$pfam))
length(unique(enr_ann_100$pfam))

# Detected at t0
mags_t0_50 <- rpkm[rpkm$X12C16O_0h_50_P09>0 | rpkm$X12C16O_0h_50_P10>0 | rpkm$X12C16O_0h_50_P13>0, "Genome"]
mags_t0_100 <- rpkm[rpkm$X12C16O_0h_100_P08>0 | rpkm$X12C16O_0h_100_P14>0 | rpkm$X12C16O_0h_100_P16>0, "Genome"]
all_ann_50_0 <- all_ann[all_ann$fasta %in% mags_t0_50,]
all_ann_100_0 <- all_ann[all_ann$fasta %in% mags_t0_100,]
length(unique(all_ann_50_0$pfam))
length(unique(all_ann_100_0$pfam))
pfam_0_50 <- all_ann_50_0 %>% count(pfam)
pfam_0_100 <- all_ann_100_0 %>% count(pfam)
H_0_50 <- diversity(pfam_0_50$n, "invsimpson")
H_0_100 <- diversity(pfam_0_100$n, "invsimpson")
# Detected at t48
mags_t48_50 <- rpkm[rpkm$X12C16O_48h_50_P09>0 | rpkm$X12C16O_48h_50_P10>0 | rpkm$X12C16O_48h_50_P13>0, "Genome"]
mags_t48_100 <- rpkm[rpkm$X12C16O_48h_100_P08>0 | rpkm$X12C16O_48h_100_P14>0 | rpkm$X12C16O_48h_100_P16>0, "Genome"]
all_ann_50_48 <- all_ann[all_ann$fasta %in% mags_t48_50,]
all_ann_100_48 <- all_ann[all_ann$fasta %in% mags_t48_100,]
length(unique(all_ann_50_48$pfam))
length(unique(all_ann_100_48$pfam))
pfam_48_50 <- all_ann_50_48 %>% count(pfam)
pfam_48_100 <- all_ann_100_48 %>% count(pfam)
H_48_50 <- diversity(pfam_48_50$n, "invsimpson")
H_48_100 <- diversity(pfam_48_100$n, "invsimpson")
# Detected at t168
mags_t168_50 <- rpkm[rpkm$X12C16O_168h_50_P09>0 | rpkm$X12C16O_168h_50_P10>0 | rpkm$X12C16O_168h_50_P13>0, "Genome"]
mags_t168_100 <- rpkm[rpkm$X12C16O_168h_100_P08>0 | rpkm$X12C16O_168h_100_P14>0 | rpkm$X12C16O_168h_100_P16>0, "Genome"]
all_ann_50_168 <- all_ann[all_ann$fasta %in% mags_t168_50,]
all_ann_100_168 <- all_ann[all_ann$fasta %in% mags_t168_100,]
length(unique(all_ann_50_168$pfam))
length(unique(all_ann_100_168$pfam))
pfam_168_50 <- all_ann_50_168 %>% count(pfam)
pfam_168_100 <- all_ann_100_168 %>% count(pfam)
H_168_50 <- diversity(pfam_168_50$n, "invsimpson")
H_168_100 <- diversity(pfam_168_100$n, "invsimpson")

# What characterizes an organism that grows at a certain time point?
# find the overlap in KOs and plot on a kegg map
enr_ann <- read.delim("DRAM/enr_MAGs_ann.tsv", sep="\t", header=T, stringsAsFactors = F)

enr_ann_24 <- enr_ann[enr_ann$fasta %in% afe$bin[afe$AFE_24h>0], c("fasta", "kegg_id")]
enr_ann_24_df <- dcast(enr_ann_24[enr_ann_24$kegg_id!="",], kegg_id~fasta)
temp <- enr_ann_24_df[,-1]
temp[temp>1] <- 1
enr_ann_24_df[,-1] <- temp
kos_24 <- enr_ann_24_df$kegg_id[which(rowSums(enr_ann_24_df[,-1])>(ncol(enr_ann_24_df)/2))]
kos_24 <- data.frame(paste0("cons24_", seq(from=1, to=length(kos_24))), kos_24, stringsAsFactors = F)
colnames(kos_24)[1] <- "id"
write.table(kos_24, file="../Traits_paper/kos_enr_24.txt", sep="\t", quote=F, row.names=F, col.names=F)

enr_ann_48 <- enr_ann[enr_ann$fasta %in% afe$bin[afe$AFE_48h>0], c("fasta", "kegg_id")]
enr_ann_48_df <- dcast(enr_ann_48[enr_ann_48$kegg_id!="",], kegg_id~fasta)
temp <- enr_ann_48_df[,-1]
temp[temp>1] <- 1
enr_ann_48_df[,-1] <- temp
kos_48 <- enr_ann_48_df$kegg_id[which(rowSums(enr_ann_48_df[,-1])>(ncol(enr_ann_48_df)/2))]
kos_48 <- data.frame(paste0("cons48_", seq(from=1, to=length(kos_48))), kos_48, stringsAsFactors = F)
colnames(kos_48)[1] <- "id"
write.table(kos_48, file="../Traits_paper/kos_enr_48.txt", sep="\t", quote=F, row.names=F, col.names=F)

enr_ann_72 <- enr_ann[enr_ann$fasta %in% afe$bin[afe$AFE_72h>0], c("fasta", "kegg_id")]
enr_ann_72_df <- dcast(enr_ann_72[enr_ann_72$kegg_id!="",], kegg_id~fasta)
temp <- enr_ann_72_df[,-1]
temp[temp>1] <- 1
enr_ann_72_df[,-1] <- temp
kos_72 <- enr_ann_72_df$kegg_id[which(rowSums(enr_ann_72_df[,-1])>(ncol(enr_ann_72_df)/2))]
kos_72 <- data.frame(paste0("cons72_", seq(from=1, to=length(kos_72))), kos_72, stringsAsFactors = F)
colnames(kos_72)[1] <- "id"
write.table(kos_72, file="../Traits_paper/kos_enr_72.txt", sep="\t", quote=F, row.names=F, col.names=F)

enr_ann_168 <- enr_ann[enr_ann$fasta %in% afe$bin[afe$AFE_168h>0], c("fasta", "kegg_id")]
enr_ann_168_df <- dcast(enr_ann_168[enr_ann_168$kegg_id!="",], kegg_id~fasta)
temp <- enr_ann_168_df[,-1]
temp[temp>1] <- 1
enr_ann_168_df[,-1] <- temp
kos_168 <- enr_ann_168_df$kegg_id[which(rowSums(enr_ann_168_df[,-1])>(ncol(enr_ann_168_df)/2))]
kos_168 <- data.frame(paste0("cons168_", seq(from=1, to=length(kos_168))), kos_168, stringsAsFactors = F)
colnames(kos_168)[1] <- "id"
write.table(kos_168, file="../Traits_paper/kos_enr_168.txt", sep="\t", quote=F, row.names=F, col.names=F)

enr_ann_48_50 <- enr_ann[enr_ann$fasta %in% afe$bin[afe$AFE_48h_50>0], c("fasta", "kegg_id")]
enr_ann_48_50_df <- dcast(enr_ann_48_50[enr_ann_48_50$kegg_id!="",], kegg_id~fasta)
temp <- enr_ann_48_50_df[,-1]
temp[temp>1] <- 1
enr_ann_48_50_df[,-1] <- temp
kos_48_50 <- enr_ann_48_50_df$kegg_id[which(rowSums(enr_ann_48_50_df[,-1])>(ncol(enr_ann_48_50_df)/2))]
kos_48_50 <- data.frame(paste0("cons4850_", seq(from=1, to=length(kos_48_50))), kos_48_50, stringsAsFactors = F)
colnames(kos_48_50)[1] <- "id"
write.table(kos_48_50, file="../Traits_paper/kos_enr_48_50.txt", sep="\t", quote=F, row.names=F, col.names=F)

enr_ann_168_50 <- enr_ann[enr_ann$fasta %in% afe$bin[afe$AFE_168h_50>0], c("fasta", "kegg_id")]
enr_ann_168_50_df <- dcast(enr_ann_168_50[enr_ann_168_50$kegg_id!="",], kegg_id~fasta)
temp <- enr_ann_168_50_df[,-1]
temp[temp>1] <- 1
enr_ann_168_50_df[,-1] <- temp
kos_168_50 <- enr_ann_168_50_df$kegg_id[which(rowSums(enr_ann_168_50_df[,-1])>(ncol(enr_ann_168_50_df)/2))]
kos_168_50 <- data.frame(paste0("cons16850_", seq(from=1, to=length(kos_168_50))), kos_168_50, stringsAsFactors = F)
colnames(kos_168_50)[1] <- "id"
write.table(kos_168_50, file="../Traits_paper/kos_enr_168_50.txt", sep="\t", quote=F, row.names=F, col.names=F)

