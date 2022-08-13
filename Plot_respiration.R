setwd("~/Dropbox/PostDoc/4th_wedge/Traits_paper/")
library(tidyverse)

resp13 <- read.csv("resp13C_forR.csv")

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(resp13, varname="d13CO2", 
                    groupnames=c("Time", "Group"))
df2$Group <- factor(df2$Group)

p<- ggplot(df2, aes(x=Time, y=d13CO2, color=Group, group=Group)) + 
  geom_line(stat="identity", color="black", 
           position=position_dodge(width=3)) +
  geom_errorbar(aes(ymin=d13CO2-sd, ymax=d13CO2+sd), width=.2,
                position=position_dodge(2)) +
  ylab("delta 13CO2 (per mil)") + xlab("Time (h)") +
  theme_classic()
ggsave(p, filename="resp13C.pdf")

resp <- read.csv("resp_total_forR.csv")
df2 <- data_summary(resp, varname="adj_CO2", 
                    groupnames=c("Time", "Group"))
df2$Group <- factor(df2$Group, levels=c("100%", "50%"))

p<- ggplot(df2, aes(x=Time, y=adj_CO2, color=Group, group=Group)) + 
  geom_line(stat="identity", color="black", 
            position=position_dodge(width=3)) +
  geom_errorbar(aes(ymin=adj_CO2-sd, ymax=adj_CO2+sd), width=.2,
                position=position_dodge(2)) +
  ylab("CO2 respired (ppm)") + xlab("Time (h)") +
  theme_classic() + scale_y_continuous(labels = scales::scientific)
p
ggsave(p, filename="resp_total.pdf")


# Plot cumulative growth and mortality
mort <- read.csv("Cumulative_mortality_summary_data.csv", header=T)
mort$prec <- factor(mort$prec)
ggplot(mort, aes(x=time, y=mort, group=prec, color=prec)) +
  geom_point() + 
  geom_smooth(method = "lm", fill=NA) + 
  theme_classic() + xlab("Time (h)") +
  ylab("Cumulative mortality (16S rRNA gene copies)")
  
summary(lm(mort$mort[mort$prec=="100"]~mort$time[mort$prec=="100"]))

# Supplemetary table 1 - coverage, AFE and MAG classification
setwd("~/Dropbox/PostDoc/4th_wedge/MAGs_disrat0.33/")
rpkm <- read.delim("coverm_rpkm_wFirmicutes.tsv", sep="\t", header=T, stringsAsFactors = F)
afe <- read.csv("../Traits_paper/AFE_all.csv")[,-1]

colnames(rpkm)[1]
colnames(afe)[1]

colnames(rpkm)[1] <- "MAG"
colnames(afe)[1] <- "MAG"

library(tidyverse)

sup_tbl <- left_join(rpkm, afe)

colnames(sup_tbl) <- gsub("X12C16O_", "t", gsub("\\.bam\\.sort\\.", "_", colnames(sup_tbl)))
colnames(sup_tbl) <- gsub("18O", "AFE", colnames(sup_tbl))
colnames(sup_tbl) <- gsub("_unf", "", colnames(sup_tbl))

sup_tbl$Growing <- FALSE
afe_cols <- which(colnames(sup_tbl) %like% "AFE")
sup_tbl$Growing[rowSums(sup_tbl[,afe_cols], na.rm = T)>0] <- TRUE

sup_tbl$Detected <- FALSE
rpkm_cols <- which(colnames(sup_tbl) %like% "RPKM")
sup_tbl$Detected[rowSums(sup_tbl[,rpkm_cols], na.rm = T)>0] <- TRUE

sup_tbl <- sup_tbl %>% select(-classification)

tax <- read.csv("../MAGs_disrat0.33/MAG_tax.csv", header=T, stringsAsFactors = F)
colnames(tax)[1] <- "MAG"

sup_tbl <- sup_tbl %>%
  left_join(tax)

drep <- read.csv("../MAGs_disrat0.33/genomeInformation_derep.csv", header=T, stringsAsFactors = F) %>%
  mutate(MAG = gsub("\\.fa", "", MAG))

sup_tbl <- sup_tbl %>%
  left_join(drep)

write.csv(sup_tbl, "../Traits_paper/sup_tbl_S1_MAGinfo.csv")

# How much of the community did growing and detected organisms account for?

counts <- read.delim("4W_unf_coverm_count_wFirmicutes.tsv", sep="\t", header=T, stringsAsFactors = F)
rps <- read.csv("../Unfractionated_numOfReads.csv", header=T, stringsAsFactors = F) %>%
  mutate(Sample = gsub("_unf", "", Sample)) %>%
  mutate(Sample = paste0("X", Sample))

colnames(counts)
colnames(rps)
rps$Sample

growing <- sup_tbl$MAG[sup_tbl$Growing]
detected <- sup_tbl$MAG[sup_tbl$Detected]

counts_growing <- colSums(counts[counts$Genome %in% growing,-1])
counts_detected <- colSums(counts[counts$Genome %in% detected,-1]) - counts_growing

pct_growing <- counts_growing*100/rps$PostQC.fwd.reads
pct_detected <- counts_detected*100/rps$PostQC.fwd.reads
mean(pct_growing)
sd(pct_growing)
mean(pct_detected)
sd(pct_detected)
