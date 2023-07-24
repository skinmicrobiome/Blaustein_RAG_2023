###########################
##### RAG METAGENOMES #####
###########################

# ryan.blaustein@nih.gov; rblauste@umd.edu #
# last edited: 07/16/2023 #


##############################
### Read Mapping: Overview ###
##############################

###
### load modules
###

library(ggplot2)
library(gtools)
library(reshape2)
library(gridExtra)
library(vegan)
library(lme4)
library(lmerTest)
library(gplots)
library(ape)
library(ggtext)
library(phangorn)

###
### load data
###

## metadata
rag_hv_meta = read.csv("data/metadata/rag_hv_metadata.csv", h=T, row.names=1)
rag_hv_meta = rag_hv_meta[grep("DNA", rag_hv_meta$Metagenome_type),]
#rag_hv_meta = rag_hv_meta[grep("RNA", rag_hv_meta$Metagenome_type),]
head(rag_hv_meta)

## taxa file
taxa_map = read.csv("data/MAGs_read_mapping/read_map_taxa_list.csv", row.names=1, h=T)
taxa_map = taxa_map[order(taxa_map$MAG),]
head(taxa_map)

## read mapping: unique
map_unique = read.table("data/MAGs_read_mapping/read_map_unique.txt", h=T, row.names=1, comment.char = "")
map_unique = map_unique[order(rownames(map_unique)),]
# check order
which(rownames(map_unique) != paste0(taxa_map$MAG, ".fa"))
## correct for read mapping coverage (>30% prokaryotes; >75% viruses)
map_unique_cov = map_unique[,grep("Coverage", colnames(map_unique))]
map_unique_counts = map_unique[,grep("Counts", colnames(map_unique))]
# prokaryotes
map_unique_counts_PE = as.matrix(map_unique_counts[grep("Euk|Arch|Bac", taxa_map$Domain),])
dim(map_unique_counts_PE) # should be 5421 rows
map_unique_counts_PE[which(as.matrix(map_unique_cov[grep("Euk|Arch|Bac", taxa_map$Domain),])<0.30)] = 0
#virus
map_unique_counts_V = as.matrix(map_unique_counts[grep("Vir", taxa_map$Domain),])
dim(map_unique_counts_V) # should be 9378 rows
map_unique_counts_V[which(as.matrix(map_unique_cov[grep("Vir", taxa_map$Domain),])<0.75)] = 0
map_unique = as.data.frame(rbind(map_unique_counts_PE, map_unique_counts_V))
map_unique = map_unique[order(rownames(map_unique)),]
head(map_unique)
# check order
which(rownames(map_unique) != paste0(taxa_map$MAG, ".fa"))

## read mapping stats
rm_stats = read.csv("data/MAGs_read_mapping/read_map_stats.csv", row.names=1, h=T)
head(rm_stats)


#######################
### SAMPLE OVERVIEW ###
#######################

## set metagenome table
samples_dna = data.frame(subject = rag_hv_meta$Subject_ID,
                         group = rag_hv_meta$Subject_pop,
                         site = rag_hv_meta$Body_site_ID,
                         age = rag_hv_meta$Subject_age,
                         pres = rep(1, nrow(rag_hv_meta)))
samples_dna = samples_dna[-c(grep("Ct", samples_dna$site)),]
samples_dna$site = factor(samples_dna$site,
                          levels = c("Hp", "Vf",
                                     "Ac", "Ic", "Pc", "Ph",
                                     "Mb", "Ra", "N", "St"))

# factor for subject order
samples_dna$subject = factor(samples_dna$subject,
                             levels = c("RAG_pt01", "RAG_pt02", "RAG_pt03", "RAG_pt04", "RAG_pt05", "RAG_pt06", "RAG_pt07", "RAG_pt08",
                                        "HC01", "HC02", "HC03", "HC04", "HV01", "HV02", "HV03", "HV04", "HV05", "HV06", "HV07", "HV08", "HV09", "HV10", "HV11", "HV12"))

levels(factor(samples_dna$subject))
samples_dna$group_age = paste(substr(samples_dna$group, start=1, stop=2), 
                              samples_dna$age, sep = "_")
samples_dna$group_age = factor(samples_dna$group_age,
                               levels = c("ra_child", "ra_adult","hv_child", "hv_adult"))

# plot (Figure S1-A)
ggplot(samples_dna,
       aes(x=subject, y=site, fill=group_age)) +
  geom_point(pch=21, #fill = "red",
             size=5) +
  xlab("Subject ID") +
  ylab("Body Site") +
  scale_fill_manual(values=c("pink", "red", "skyblue", "blue"),
                    labels=c("RAG (child)", "RAG (adult)", "HC (child)", "HV (adult)")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle=90, hjust=1, vjust=0.5
        ),
        axis.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  theme(legend.position = "none")


###
### set RNA-metagenome table
###

## metadata
rag_hv_meta_RNA = read.csv("data/metadata/rag_hv_metadata.csv", h=T, row.names=1)
rag_hv_meta_RNA = rag_hv_meta_RNA[grep("RNA", rag_hv_meta_RNA$Metagenome_type),]
head(rag_hv_meta_RNA)

## subset
samples_rna = data.frame(subject = rag_hv_meta_RNA$Subject_ID,
                         group = rag_hv_meta_RNA$Subject_pop,
                         site = rag_hv_meta_RNA$Body_site_ID,
                         age = rag_hv_meta_RNA$Subject_age,
                         pres = rep(1, nrow(rag_hv_meta_RNA)))
samples_rna = samples_rna[-c(grep("Ct", samples_rna$site)),]
samples_rna$site = factor(samples_rna$site,
                          levels = c("N", "St"))

# factor for subject order
levels(factor(samples_rna$subject))
samples_rna$group_age = paste(substr(samples_rna$group, start=1, stop=2), 
                              samples_rna$age, sep = "_")
samples_rna$group_age = factor(samples_rna$group_age,
                               levels = c("ra_child", "ra_adult"))

# plot (Figure S1-B)
ggplot(samples_rna,
       aes(x=subject, y=site, fill=group_age)) +
  geom_point(pch=21, #fill = "red",
             size=5) +
  xlab("Subject ID") +
  ylab("Body Site") +
  scale_fill_manual(values=c("pink", "red"),
                    labels=c("RAG (child)", "RAG (adult)")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle=60, hjust=1
        ),
        axis.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  theme(legend.position = "none")


##########################################
### Compare BWA read mapping to Kraken ###
##########################################

## summary stats (skin/nares)
# bwa-read mapping
100*tapply(rm_stats$rm_properly_paired[grep("rag_skin", rm_stats$Subject_group)], rm_stats$Subject_group[grep("rag_skin", rm_stats$Subject_group)], mean)
100*tapply(rm_stats$rm_properly_paired[grep("rag_skin", rm_stats$Subject_group)], rm_stats$Subject_group[grep("rag_skin", rm_stats$Subject_group)], sd)
# kraken
100*tapply(rm_stats$kraken_root[grep("rag_skin", rm_stats$Subject_group)], rm_stats$Subject_group[grep("rag_skin", rm_stats$Subject_group)], mean)
100*tapply(rm_stats$kraken_root[grep("rag_skin", rm_stats$Subject_group)], rm_stats$Subject_group[grep("rag_skin", rm_stats$Subject_group)], sd)
# comparison
wilcox.test(rm_stats$rm_properly_paired[grep("rag_skin", rm_stats$Subject_group)], rm_stats$kraken_root[grep("rag_skin", rm_stats$Subject_group)])

## summary stats (stool)
# bwa-read mapping
100*tapply(rm_stats$rm_properly_paired[grep("rag_stool", rm_stats$Subject_group)], rm_stats$Subject_group[grep("rag_stool", rm_stats$Subject_group)], mean)
100*tapply(rm_stats$rm_properly_paired[grep("rag_stool", rm_stats$Subject_group)], rm_stats$Subject_group[grep("rag_stool", rm_stats$Subject_group)], sd)
# kraken
100*tapply(rm_stats$kraken_root[grep("rag_stool", rm_stats$Subject_group)], rm_stats$Subject_group[grep("rag_stool", rm_stats$Subject_group)], mean)
100*tapply(rm_stats$kraken_root[grep("rag_stool", rm_stats$Subject_group)], rm_stats$Subject_group[grep("rag_stool", rm_stats$Subject_group)], sd)
# comparison
wilcox.test(rm_stats$rm_properly_paired[grep("rag_stool", rm_stats$Subject_group)], rm_stats$kraken_root[grep("rag_stool", rm_stats$Subject_group)])


###
### plot RAG-deficient patient data 
###

## RAG vs. HV

# set data
rm_stats_plot = melt(rm_stats[-c(grep("control", rm_stats$Subject_group)),])

# boxplots (Figure S3)
ggplot(rm_stats_plot[grep("rag_s", rm_stats_plot$Subject_group),],
       aes(x=Subject_group, y=value,
           fill=factor(variable, levels = c("kraken_root", "rm_properly_paired")))) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0, outlier.colour = "white", alpha=0.6, size=1.3) +
  geom_jitter(position = position_jitterdodge(dodge.width=0.8, jitter.width=0.4),
              aes(fill=factor(variable, levels = c("kraken_root", "rm_properly_paired"))), 
              size=1.8, alpha=0.6, 
              pch=21, 
              color="black") +
  theme_bw() +
  scale_fill_manual(values = c("white", "black"),
                    labels = c("Kraken", "BWA-MEM")) +
  scale_x_discrete(labels = c("RAG-Skin", "RAG-Stool")) +
  ylim(0,1) +
  ylab("Reads Classified") +
  theme(axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust=1),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank())


#####################################
### Read Mapping: Skin Microbiota ###
#####################################

## relative abundance by MAG
map_unique_counts = map_unique
map_unique_ra = as.data.frame(matrix(nrow=nrow(map_unique_counts), ncol=ncol(map_unique_counts)))
for (i in 1:159) {
  map_unique_ra[,i] = map_unique_counts[,i]/sum(map_unique_counts[,i])
}
colnames(map_unique_ra) = colnames(map_unique_counts)
rownames(map_unique_ra) = rownames(map_unique_counts)
head(map_unique_ra)
apply(map_unique_ra, 2, sum)

## check averages
# skin/nares
map_unique_ra$avg_skin = apply(map_unique_ra[, grep("hv_skin|rag_skin", rag_hv_meta$Subject_group)], 1, mean)
map_unique_ra$avg_hv_skin = apply(map_unique_ra[, grep("hv_skin", rag_hv_meta$Subject_group)], 1, mean)
map_unique_ra$avg_hv_skin_child = apply(map_unique_ra[, intersect(grep("hv_skin", rag_hv_meta$Subject_group), grep("child", rag_hv_meta$Subject_age))], 1, mean)
map_unique_ra$avg_hv_skin_adult = apply(map_unique_ra[, intersect(grep("hv_skin", rag_hv_meta$Subject_group), grep("adult", rag_hv_meta$Subject_age))], 1, mean)
map_unique_ra$avg_rag_skin = apply(map_unique_ra[, grep("rag_skin", rag_hv_meta$Subject_group)], 1, mean)
map_unique_ra$avg_rag_skin_child = apply(map_unique_ra[, intersect(grep("rag_skin", rag_hv_meta$Subject_group), grep("child", rag_hv_meta$Subject_age))], 1, mean)
map_unique_ra$avg_rag_skin_adult = apply(map_unique_ra[, intersect(grep("rag_skin", rag_hv_meta$Subject_group), grep("adult", rag_hv_meta$Subject_age))], 1, mean)
# stool
map_unique_ra$avg_stool = apply(map_unique_ra[, grep("stool", rag_hv_meta$Subject_group)], 1, mean)
map_unique_ra$avg_hv_stool = apply(map_unique_ra[, grep("hv_stool", rag_hv_meta$Subject_group)], 1, mean)
map_unique_ra$avg_rag_stool = apply(map_unique_ra[, grep("rag_stool", rag_hv_meta$Subject_group)], 1, mean)
# all
map_unique_ra$avg_all = apply(map_unique_ra, 1, mean)

## add taxa
map_unique_ra = data.frame(map_unique_ra, taxa_map)
head(map_unique_ra[order(map_unique_ra$avg_all, decreasing = T), c(160:170)], 20)
head(map_unique_ra[order(map_unique_ra$avg_all, decreasing = T), c(160:172)], 20)

##
## subset microbial species (and higher level taxonomic groups for MAGs with non-defined species)
##

# set data frame
map_ra_species = as.data.frame(matrix(nrow=length(tapply(map_unique_ra[,8], 
                                                         paste(map_unique_ra$Domain, map_unique_ra$Phylum, map_unique_ra$Class, map_unique_ra$Order, 
                                                               map_unique_ra$Family, map_unique_ra$Genus, map_unique_ra$Species), sum)),
                                      ncol=170))
# run for every row
for (i in 1:170) {
  map_ra_species[,i] = tapply(map_unique_ra[,i], 
                              paste(map_unique_ra$Domain, map_unique_ra$Phylum, map_unique_ra$Class, map_unique_ra$Order, 
                                    map_unique_ra$Family, map_unique_ra$Genus, map_unique_ra$Species), sum)
}
colnames(map_ra_species) = colnames(map_unique_ra[,1:170])
rownames(map_ra_species) = names(tapply(map_unique_ra[,8], 
                                        paste(map_unique_ra$Domain, map_unique_ra$Phylum, map_unique_ra$Class, map_unique_ra$Order, 
                                              map_unique_ra$Family, map_unique_ra$Genus, map_unique_ra$Species), sum))
head(map_ra_species[order(map_ra_species$avg_skin, decreasing = T), c(167:170)], 20)

# re-order
map_ra_species = map_ra_species[order(map_ra_species$avg_skin, decreasing = T),]
map_ra_species[1:12, 163:166]
dim(map_ra_species)

## Table S5
#write.table(map_ra_species, "supp_table_5.txt", sep = "\t", quote = FALSE)


###
### Beta diversity
###

## subset skin taxa
skin_taxa = map_ra_species[,grep("hv_skin|rag_skin", rag_hv_meta$Subject_group)]
skin_taxa$avg = apply(skin_taxa, 1, mean)
skin_taxa = skin_taxa[order(skin_taxa$avg, decreasing=T),]
skin_taxa[1:30,78:79]
# meta
skin_meta = rag_hv_meta[grep("hv_skin|rag_skin", rag_hv_meta$Subject_group),]

## HV skin microbiota: associations with subject, skin site, skin niche, and age group
# beta-diversity
beta <- vegdist(t(skin_taxa[,grep("hv_", skin_meta$Subject_group)]))
# PERMANOVAs (Table S4)
adonis2(beta ~ Subject_ID, skin_meta[grep("hv_", skin_meta$Subject_group),], permutations=999)
adonis2(beta ~ Body_site_ID, skin_meta[grep("hv_", skin_meta$Subject_group),], permutations=999)
adonis2(beta ~ Subject_age, skin_meta[grep("hv_", skin_meta$Subject_group),], permutations=999)

## RAG skin microbiota: associations with subject, skin site, skin niche, and age group
# beta-diversity
beta <- vegdist(t(skin_taxa[,grep("rag", skin_meta$Subject_group)]))
# PERMANOVAs (Table S4)
adonis2(beta ~ Subject_ID, skin_meta[grep("rag", skin_meta$Subject_group),], permutations=999)
adonis2(beta ~ Body_site_ID, skin_meta[grep("rag", skin_meta$Subject_group),], permutations=999)
adonis2(beta ~ Body_site_niche, skin_meta[grep("rag", skin_meta$Subject_group),], permutations=999)
adonis2(beta ~ Subject_age, skin_meta[grep("rag", skin_meta$Subject_group),], permutations=999)

## RAG clinical phenotype
beta <- vegdist(t(skin_taxa[,intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID),
                                       grep("rag", skin_meta$Subject_group))]))
pcoa <- cmdscale(beta, k=4, eig=T)
ord <- as.data.frame(pcoa$points)
names(ord) <- c("pcoa1", "pcoa2", "pcoa3", "pcoa4")
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# metadata
ord$group_id = skin_meta[intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID),
                                   grep("rag", skin_meta$Subject_group)),]$Subject_group
ord$Site_ID = skin_meta[intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID),
                                  grep("rag", skin_meta$Subject_group)),]$Body_site_ID
ord$Site_ID = factor(ord$Site_ID, levels=c("Vf", "Ac", "Ra", "N"))
ord$patient_ID = skin_meta[intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID),
                                     grep("rag", skin_meta$Subject_group)),]$Subject_ID
ord$patient_age = skin_meta[intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID),
                                      grep("rag", skin_meta$Subject_group)),]$Subject_age
ord$age_site = paste(ord$patient_age, ord$Site_ID, sep = "_")
ord$age_site = factor(ord$age_site, 
                      levels=c("adult_Vf", "adult_Ac", "adult_Ra", "adult_N",
                               "child_Vf", "child_Ac", "child_Ra", "child_N"))
ord$RAG_type = "SCID"
ord$RAG_type[grep("01", ord$patient_ID)] = "LS"
ord$RAG_type[grep("02", ord$patient_ID)] = "CID-G/AI"
ord$RAG_type[grep("03", ord$patient_ID)] = "SCID"
ord$RAG_type[grep("04", ord$patient_ID)] = "CID"
ord$RAG_type[grep("05", ord$patient_ID)] = "CID-G/AI"
ord$RAG_type[grep("06", ord$patient_ID)] = "CID-G/AI"
ord$RAG_type[grep("07", ord$patient_ID)] = "CID-G/AI"
ord$RAG_type[grep("08", ord$patient_ID)] = "CID"

## plot pcoa (Figure S4)
ggplot(data=ord, 
       aes(x=pcoa1, y=-pcoa2, fill=RAG_type, col=RAG_type,
           shape=age_site)) +
  geom_vline(xintercept=0, lty=2, color="darkgray") +
  geom_hline(yintercept=0, lty=2, color="darkgray") +
  geom_point(size=4, alpha=0.9, stroke=1.3) +
  theme_bw() +
  scale_fill_manual(values=c("tan", "brown", "gray", "darkgreen")) +
  scale_color_manual(values=c("tan", "brown", "gray", "darkgreen")) +
  scale_shape_manual(values = c(25, 22, 24, 23, 6, 0, 2, 5)) +
  xlab("PCoA1 (26.6%)") +
  ylab("PCoA2 (16.2%)") +
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=13),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13))

#PERMANOVAs
adonis2(beta ~ RAG_type, ord, permutations=999)
adonis2(beta ~ Site_ID, ord, permutations=999)
adonis2(beta ~ patient_age, ord, permutations=999)
adonis2(beta ~ patient_ID, ord, permutations=999)


###
### RAG vs. HV skin microbiota: site match
###

beta <- vegdist(t(skin_taxa[,grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID)]))
pcoa <- cmdscale(beta, k=4, eig=T)
ord <- as.data.frame(pcoa$points)
names(ord) <- c("pcoa1", "pcoa2", "pcoa3", "pcoa4")
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# metadata
ord$group_id = skin_meta[grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID),]$Subject_group
ord$Site_ID = skin_meta[grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID),]$Body_site_ID
ord$Site_ID = factor(ord$Site_ID, levels=c("Vf", "Ac", "Ra", "N"))
ord$patient_age = skin_meta[grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID),]$Subject_age
ord$age_site = paste(ord$patient_age, ord$Site_ID, sep = "_")
ord$age_site = factor(ord$age_site, 
                      levels=c("adult_Vf", "adult_Ac", "adult_Ra", "adult_N",
                               "child_Vf", "child_Ac", "child_Ra", "child_N"))
## plot pcoa (Figure 1-C)
ggplot(data=ord, 
       aes(x=pcoa1, y=pcoa2, fill=group_id, col=group_id, 
           shape=age_site)) +
  geom_vline(xintercept=0, lty=2, color="darkgray") +
  geom_hline(yintercept=0, lty=2, color="darkgray") +
  geom_point(size=3.5, alpha=0.6, stroke=1) +
  theme_bw() +
  scale_fill_manual(values=c("blue", "red")) +
  scale_color_manual(values=c("blue", "red")) +
  scale_shape_manual(values = c(25, 22, 24, 23, 6, 0, 2, 5)) +
  xlab("PCoA1 (19.0%)") +
  ylab("PCoA2 (15.7%)") +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))

## PERMANOVAs (Table S4)

#RAG-adult
beta <- vegdist(t(skin_taxa[,intersect(intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID), grep("rag", skin_meta$Subject_group)),
                                       grep("adult", skin_meta$Subject_age))]))
adonis2(beta ~ Subject_ID, skin_meta[intersect(intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID), grep("rag", skin_meta$Subject_group)),
                                              grep("adult", skin_meta$Subject_age)),], permutations=999)
adonis2(beta ~ Body_site_ID, skin_meta[intersect(intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID), grep("rag", skin_meta$Subject_group)),
                                                 grep("adult", skin_meta$Subject_age)),], permutations=999)
# RAG-child
beta <- vegdist(t(skin_taxa[,intersect(intersect(grep("Ac|Vf|Ra", skin_meta$Body_site_ID), grep("rag", skin_meta$Subject_group)),
                                       grep("child", skin_meta$Subject_age))]))
adonis2(beta ~ Subject_ID, skin_meta[intersect(intersect(grep("Ac|Vf|Ra", skin_meta$Body_site_ID), grep("rag", skin_meta$Subject_group)),
                                              grep("child", skin_meta$Subject_age)),], permutations=999)
adonis2(beta ~ Body_site_ID, skin_meta[intersect(intersect(grep("Ac|Vf|Ra", skin_meta$Body_site_ID), grep("rag", skin_meta$Subject_group)),
                                                 grep("child", skin_meta$Subject_age)),], permutations=999)
# HV
beta <- vegdist(t(skin_taxa[,intersect(intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID), grep("hv_", skin_meta$Subject_group)),
                                       grep("adult", skin_meta$Subject_age))]))
adonis2(beta ~ Subject_ID, skin_meta[intersect(intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID), grep("hv_", skin_meta$Subject_group)),
                                              grep("adult", skin_meta$Subject_age)),], permutations=999)
adonis2(beta ~ Body_site_ID, skin_meta[intersect(intersect(grep("Ac|Vf|Ra|N", skin_meta$Body_site_ID), grep("hv_", skin_meta$Subject_group)),
                                                 grep("adult", skin_meta$Subject_age)),], permutations=999)

# HC
beta <- vegdist(t(skin_taxa[,intersect(intersect(grep("Ac|Vf|Ra", skin_meta$Body_site_ID), grep("hv_", skin_meta$Subject_group)),
                                       grep("child", skin_meta$Subject_age))]))
adonis2(beta ~ Subject_ID, skin_meta[intersect(intersect(grep("Ac|Vf|Ra", skin_meta$Body_site_ID), grep("hv_", skin_meta$Subject_group)),
                                              grep("child", skin_meta$Subject_age)),], permutations=999)
adonis2(beta ~ Body_site_ID, skin_meta[intersect(intersect(grep("Ac|Vf|Ra", skin_meta$Body_site_ID), grep("hv_", skin_meta$Subject_group)),
                                                 grep("child", skin_meta$Subject_age)),], permutations=999)


###
### Taxon enrichments
###

## SKIN sites (Ac, Ra, Vf)
## LME model -- RAG vs. HV (run adult and child comparisons separately)
test = data.frame(t(map_ra_species[,1:159]), rag_hv_meta)
p_vals = c(1:nrow(map_ra_species))
for (i in 1:nrow(map_ra_species)) {
  p_vals[i] = summary(lmer(test[intersect(grep("hv_|rag_", test$Subject_group), intersect(grep("Ac|Vf|Ra", test$Body_site_ID), grep("adult", test$Subject_age))),i] ~ Subject_group + (1|Body_site_ID),
                           test[intersect(grep("hv_|rag_", test$Subject_group), intersect(grep("Ac|Vf|Ra", test$Body_site_ID), grep("adult", test$Subject_age))),]))$coefficients[2,5]
  #p_vals[i] = summary(lmer(test[intersect(grep("Ac|Vf|Ra", test$Body_site_ID), grep("child", test$Subject_age)),i] ~ Subject_group + (1|Body_site_ID),
  #                         test[intersect(grep("Ac|Vf|Ra", test$Body_site_ID), grep("child", test$Subject_age)),]))$coefficients[2,5]
}
names(p_vals) = rownames(map_ra_species)
# result
length(which(p_vals < 0.05))
names(p_vals)[which(p_vals < 0.05)]
p_vals[which(p_vals < 0.05)][1:10]
# result - bh correction
length(which(p.adjust(p_vals, method = "BH") < 0.1))
head(names(p_vals)[which(p.adjust(p_vals, method = "BH") < 0.1)], 10)
head(p_vals[which(p.adjust(p_vals, method = "BH") < 0.1)], 10)
p.adjust(p_vals, method = "BH")[which(p.adjust(p_vals, method = "BH") < 0.1)][1:10]
p.adjust(p_vals, method = "BH")[which(p_vals < 0.1)][1:40]

# papillomaviridae
mean(as.numeric(map_ra_species[2,intersect(grep("rag_", test$Subject_group), intersect(grep("Ac|Vf|Ra", test$Body_site_ID), grep("adult", test$Subject_age)))]))
sd(as.numeric(map_ra_species[2,intersect(grep("rag_", test$Subject_group), intersect(grep("Ac|Vf|Ra", test$Body_site_ID), grep("adult", test$Subject_age)))]))/sqrt(18)

## Nares (N)
## MW test & t-test for significance
p_vals = as.data.frame(matrix(nrow=nrow(map_ra_species), ncol=2))
for (i in 1:nrow(map_ra_species)) {
  p_vals[i, 1] = t.test(as.numeric(map_ra_species[i, intersect(intersect(grep("hv_skin", rag_hv_meta$Subject_group), grep("adult", rag_hv_meta$Subject_age)), grep("N", rag_hv_meta$Body_site_ID))]),
                        as.numeric(map_ra_species[i, intersect(intersect(grep("rag_skin", rag_hv_meta$Subject_group), grep("adult", rag_hv_meta$Subject_age)), grep("N", rag_hv_meta$Body_site_ID))]))$p.val
  p_vals[i, 2] = wilcox.test(as.numeric(map_ra_species[i, intersect(intersect(grep("hv_skin", rag_hv_meta$Subject_group), grep("adult", rag_hv_meta$Subject_age)), grep("N", rag_hv_meta$Body_site_ID))]),
                             as.numeric(map_ra_species[i, intersect(intersect(grep("rag_skin", rag_hv_meta$Subject_group), grep("adult", rag_hv_meta$Subject_age)), grep("N", rag_hv_meta$Body_site_ID))]))$p.val
}
colnames(p_vals) = c("t_test", "mw_test")
rownames(p_vals) = rownames(map_ra_species)
# result
length(which(p_vals[,2] < 0.05))
p_vals[which(p_vals[,2] < 0.05),]
# result - bh correction
length(which(p.adjust(p_vals[,2], method = "BH") < 0.1))
rownames(p_vals)[which(p.adjust(p_vals[,2], method = "BH") < 0.1)]


###
### Select taxa for plot
###

apply(map_ra_species, 2, sum)
map_ra_species[c(1:15), 160:163]
rownames(map_ra_species)[1:15]

## subset
map_taxa_select = rbind(apply(map_ra_species[grep("d__Bacteria", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("d__Eukaryota", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("d__Virus", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Actinobacteriota", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Bacteroidota", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Firmicutes", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Proteobacteria", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("g__Corynebacterium", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("g__Streptococcus", rownames(map_ra_species)),], 2, sum),
                        map_ra_species[c(1:3,6:10,12),]) # top taxa, with Corynebacterium and Streptococcus at genus-level
rownames(map_taxa_select) = c("d__Bacteria", "d__Eukaryota", "d__Virus",
                              "d__Bacteria p__Actinobacteriota", "d__Bacteria p__Bacteroidota", 
                              "d__Bacteria p__Firmicutes", "d__Bacteria p__Proteobacteria", 
                              "d__Bacteria p__Actinobacteriota c__Actinomycetia o__Mycobacteriales f__Mycobacteriaceae g__Corynebacterium",
                              "d__Bacteria p__Firmicutes c__Bacilli o__Lactobacillales f__Streptococcaceae g__Streptococcus",
                              rownames(map_ra_species[c(1:3,6:10,12),]))
head(map_taxa_select)
map_taxa_select
rownames(map_taxa_select)

# fungi
map_taxa_select[grep("^d__Eukaryota$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Eukaryota$", rownames(map_taxa_select)),] - map_taxa_select[grep("Malass", rownames(map_taxa_select)),]
# virus
map_taxa_select[grep("^d__Virus$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Virus$", rownames(map_taxa_select)),] - map_taxa_select[grep("Papill", rownames(map_taxa_select)),]
# phyla
map_taxa_select[grep("^d__Bacteria p__Actinobacteriota$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria p__Actinobacteriota$", rownames(map_taxa_select)),] - apply(map_taxa_select[grep("Cuti|Cory|Micr", rownames(map_taxa_select)),], 2, sum)
map_taxa_select[grep("^d__Bacteria p__Firmicutes$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria p__Firmicutes$", rownames(map_taxa_select)),] - apply(map_taxa_select[grep("Staph|Strep", rownames(map_taxa_select)),], 2, sum)
map_taxa_select[grep("^d__Bacteria p__Proteobacteria$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria p__Proteobacteria$", rownames(map_taxa_select)),] - apply(map_taxa_select[grep("QFNR01|Escher|Haem", rownames(map_taxa_select)),], 2, sum)
# bacteria
map_taxa_select[grep("^d__Bacteria$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria$", rownames(map_taxa_select)),] - apply(map_taxa_select[-c(grep("^d__Bacteria$|d__Eukaryota|d__Virus|d__Archaea", rownames(map_taxa_select))),], 2, sum)

# check sums
apply(map_taxa_select, 2, sum)

## RA by subject
mpa_ra = map_taxa_select
mpa_ra = mpa_ra[order(rownames(mpa_ra), decreasing=T),]
rownames(mpa_ra)
dim(mpa_ra)

###
### By subject: averages for matching skin sites
###

mpa_bar = data.frame(melt(as.data.frame(rbind(tapply(as.numeric(mpa_ra[1, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[2, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[3, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[4, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[5, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[6, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[7, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[8, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[9, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[10, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[11, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[12, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[13, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[14, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[15, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[16, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[17, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean),
                                              tapply(as.numeric(mpa_ra[18, intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))]), rag_hv_meta$Subject_ID[intersect(grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk", rag_hv_meta$Subject_group))], mean)
                                              ))),
                     taxon = rep(c("     Papillomaviridae", "Virus/Phage", 
                                   "     Malassezia restricta", "Fungi", 
                                   "     Haemophilus influenzae",  "     Escherichia coli",  "     QFRN01/Candidatus pellibacterium", "  Proteobacteria", 
                                   "     Staphylococcus hominis", "     Staphylococcus epidermidis", "     Streptococcus sp.", "  Firmicutes", 
                                   "  Bacteroidota",
                                   "     Cutibacterium acnes", "     Corynebacterium sp.","     Micrococcus luteus", "  Actinobacteriota", 
                                   "Bacteria"), 18))

# add group type and age
mpa_bar$group_id = rep("RAG (child)", nrow(mpa_bar))
mpa_bar$group_id[grep("RAG_pt05|RAG_pt06|RAG_pt07|RAG_pt08", mpa_bar$variable)] = "RAG (adult)"
mpa_bar$group_id[grep("HC", mpa_bar$variable)] = "HC (child)"
mpa_bar$group_id[grep("HV", mpa_bar$variable)] = "HV (adult)"
mpa_bar$group_id = factor(mpa_bar$group_id,
                          levels = c("RAG (adult)", "HV (adult)", "RAG (child)", "HC (child)"))

# taxa factor
mpa_bar$taxon = factor(mpa_bar$taxon,
                       levels = rev(c("     Papillomaviridae", "Virus/Phage", "     Malassezia restricta", "Fungi", 
                                      "     Haemophilus influenzae",  "     Escherichia coli",  "     QFRN01/Candidatus pellibacterium", "  Proteobacteria", 
                                      "     Streptococcus sp.", "     Staphylococcus hominis", "     Staphylococcus epidermidis", "  Firmicutes", 
                                      "  Bacteroidota",
                                      "     Micrococcus luteus", "     Corynebacterium sp.","     Cutibacterium acnes", "  Actinobacteriota", 
                                      "Bacteria")))

# rename factors for axis
mpa_bar$variable = as.character(mpa_bar$variable)

# plot (Figure 1-A)
ggplot(mpa_bar,
       aes(x=variable, y=value, fill=taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Subject ID") +
  facet_grid(~group_id, scales="free", space="free"
  ) +
  scale_fill_manual(values=rev(c("#CCCCFF", "#9933FF", "#FFCCFF", "#FF33FF", 
                                 "#CCFFCC", "#66FF66",  "#339933", "#006600", # greens
                                 "#CCFFFF", "#66CCFF", "#3399FF", "blue", # blues 
                                 "#FF3333", #reds
                                 "#FFCC99", "#CC6633", "#993300", "#663300", 
                                 "lightgray")),
                    labels = rev(c("*Papillomaviridae*", "Virus/Phage", "*Malassezia restricta*", "Fungi", 
                                   "*Haemophilus influenzae*",  "*Escherichia coli*",  "*Candidatus pellibacterium*/QFRN01", "Proteobacteria", 
                                   "*Streptococcus* sp.", "*Staphylococcus hominis*", "*Staphylococcus epidermidis*", "Firmicutes", 
                                   "Bacteroidota",
                                   "*Micrococcus luteus*", "*Corynebacterium* sp.","*Cutibacterium acnes*", "Actinobacteriota", 
                                   "Bacteria"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(size=11, angle=90, vjust=0.5),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 13)) +
  guides(fill=guide_legend(ncol=1)) #+
#theme(legend.position = "none")

###
### By subject: averages for nares
###

mpa_bar = data.frame(melt(as.data.frame(rbind(tapply(as.numeric(mpa_ra[1, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[2, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[3, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[4, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[5, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[6, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[7, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[8, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[9, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[10, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[11, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[12, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[13, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[14, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[15, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[16, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[17, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean),
                                              tapply(as.numeric(mpa_ra[18, grep("N", rag_hv_meta$Body_site_ID)]), rag_hv_meta$Subject_ID[grep("N", rag_hv_meta$Body_site_ID)], mean)))),
                     taxon = rep(c("     Papillomaviridae", "Virus/Phage", "     Malassezia restricta", "Fungi", 
                                   "     Haemophilus influenzae",  "     Escherichia coli",  "     QFRN01/Candidatus pellibacterium", "  Proteobacteria", 
                                   "     Staphylococcus hominis", "     Staphylococcus epidermidis", "     Streptococcus sp.", "  Firmicutes", 
                                   "  Bacteroidota",
                                   "     Cutibacterium acnes", "     Corynebacterium sp.","     Micrococcus luteus", "  Actinobacteriota", 
                                   "Bacteria"), 11))

# add group type and age
mpa_bar$group_id = rep("RAG (child)", nrow(mpa_bar))
mpa_bar$group_id[grep("RAG_pt05|RAG_pt06|RAG_pt07|RAG_pt08", mpa_bar$variable)] = "RAG (adult)"
mpa_bar$group_id[grep("HC", mpa_bar$variable)] = "HC (child)"
mpa_bar$group_id[grep("HV", mpa_bar$variable)] = "HV (adult)"
mpa_bar$group_id = factor(mpa_bar$group_id,
                          levels = c("RAG (adult)", "HV (adult)", "RAG (child)", "HC (child)"))

# taxa factor
mpa_bar$taxon = factor(mpa_bar$taxon,
                       levels = rev(c("     Papillomaviridae", "Virus/Phage", "     Malassezia restricta", "Fungi", 
                                      "     Haemophilus influenzae",  "     Escherichia coli",  "     QFRN01/Candidatus pellibacterium", "  Proteobacteria", 
                                      "     Streptococcus sp.", "     Staphylococcus hominis", "     Staphylococcus epidermidis", "  Firmicutes", 
                                      "  Bacteroidota",
                                      "     Micrococcus luteus", "     Corynebacterium sp.","     Cutibacterium acnes", "  Actinobacteriota", 
                                      "Bacteria")))

# rename factors for axis
mpa_bar$variable = as.character(mpa_bar$variable)

# add blank for Nares
mpa_bar = data.frame(rbind(mpa_bar[1:18,], mpa_bar))
mpa_bar[1:18,]$variable = "HV06"
mpa_bar[1:18,]$value = 0

# plot (Figure 1-B)
ggplot(mpa_bar,
       aes(x=variable, y=value, fill=taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Subject ID") +
  facet_grid(~group_id, scales="free", space="free"
  ) +
  scale_fill_manual(values=rev(c("#CCCCFF", "#9933FF", "#FFCCFF", "#FF33FF", 
                                 "#CCFFCC", "#66FF66",  "#339933", "#006600", # greens
                                 "#CCFFFF", "#66CCFF", "#3399FF", "blue", # blues 
                                 "#FF3333", #reds
                                 "#FFCC99", "#CC6633", "#993300", "#663300", 
                                 "lightgray")),
                    labels = rev(c("*Papillomaviridae*", "Virus/Phage", "*Malassezia restricta*", "Fungi", 
                                   "*Haemophilus influenzae*",  "*Escherichia coli*",  "*Candidatus pellibacterium*/QFRN01", "Proteobacteria", 
                                   "*Streptococcus* sp.", "*Staphylococcus hominis*", "*Staphylococcus epidermidis*", "Firmicutes", 
                                   "Bacteroidota",
                                   "*Micrococcus luteus*", "*Corynebacterium* sp.","*Cutibacterium acnes*", "Actinobacteriota", 
                                   "Bacteria"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(size=11, angle=90, vjust=0.5),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 13)) +
  guides(fill=guide_legend(ncol=1)) #+
  #theme(legend.position = "none")


#########################
### RAG Skin: Novelty ###
#########################

## set rag skin novelty data frame

# call RAG-patient-derived MAGs
rsn = map_unique_ra[grep("d__Bacteria", map_unique_ra$Domain),]
rsn = rsn[order(rsn$avg_all, decreasing = T), ]
##rsn = rsn[-c(grep("SMGC|hv_skin|hv_stool", rownames(rsn))),]
#rsn = rsn[c(grep("rag", rownames(rsn))),]

# subset to matching sites with HV/HC
heat_plot = as.matrix(rsn[, intersect(grep("N|Ac|Vf|Ra|Ct", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk|con", rag_hv_meta$Subject_group))])
head(heat_plot)
dim(heat_plot)
length(tapply(heat_plot[1,], rag_hv_meta$Subject_ID[intersect(grep("N|Ac|Vf|Ra|Ct", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk|con", rag_hv_meta$Subject_group))], mean))

# calculate average
heat_plot2 = matrix(nrow = nrow(heat_plot), ncol = 20)
for (i in 1:dim(heat_plot2)[1]) {
  heat_plot2[i,] = tapply(heat_plot[i,], rag_hv_meta$Subject_ID[intersect(grep("N|Ac|Vf|Ra|Ct", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk|con", rag_hv_meta$Subject_group))], mean)
}
colnames(heat_plot2) = names(tapply(heat_plot[1,], rag_hv_meta$Subject_ID[intersect(grep("N|Ac|Vf|Ra|Ct", rag_hv_meta$Body_site_ID), grep("hv_sk|rag_sk|con", rag_hv_meta$Subject_group))], mean))
rownames(heat_plot2) = rownames(heat_plot)

# phylum text colors
levels(factor(rsn$Phylum))
rsn$phy_col = rep("purple", length(rsn$Phylum))
rsn$phy_col[grep("Actin", rsn$Phylum)] = "#663300"
rsn$phy_col[grep("Bactero", rsn$Phylum)] = "#FF3333"
rsn$phy_col[grep("Verr", rsn$Phylum)] = "#FFCC33"
rsn$phy_col[grep("Firm", rsn$Phylum)] = "#3399FF"
rsn$phy_col[grep("Prot", rsn$Phylum)] = "#006600"

# adjust row removal
rownames(rsn) == rownames(heat_plot2)
colnames(heat_plot2)

# remove species that are >0.1% abundant in HV/HC and controls
rsn = rsn[which(apply(heat_plot2[,c(2:11)], 1, mean)<0.001),]
heat_plot2 = heat_plot2[which(apply(heat_plot2[,c(2:11)], 1, mean)<0.001),]

rsn = rsn[which(apply(heat_plot2[,c(1, 12)], 1, mean)<0.001),]
heat_plot2 = heat_plot2[which(apply(heat_plot2[,c(1, 12)], 1, mean)<0.001),]

# remove species that are <0.1% abundant in RAG-skin samples
rsn = rsn[which(apply(heat_plot2[,c(13:20)], 1, mean)>0.001),]
heat_plot2 = heat_plot2[which(apply(heat_plot2[,c(13:20)], 1, mean)>0.001),]

# re-order plot
colnames(heat_plot2)
heat_plot2 = heat_plot2[,13:20]

heat_plot2 = heat_plot2[order(paste(rsn$Genus, rsn$Species, sep = "|")),]
rsn = rsn[order(paste(rsn$Genus, rsn$Species, sep = "|")),]

heat_plot2 = heat_plot2[order(rsn$Phylum),]
rsn = rsn[order(rsn$Phylum),]

# set color for MAG source
rsn$source_col = rep("blue", length(rsn$phy_col))
rsn$source_col[grep("rag", rownames(rsn))] = "red"

# set taxa
taxa = substr(rsn$Species, start=4, stop=nchar(rsn$Species))
taxa[which(taxa == "")] = paste(substr(rsn$Genus, start=4, stop=nchar(rsn$Genus))[which(taxa == "")], "sp.")

# set legend
max(heat_plot2)
quantile(heat_plot2)
sort(heat_plot2, decreasing = T)
mypalette = colorRampPalette(c("#FFCCFF", "#9900CC", "#660099", "#330066", "#000033", "black"))(n=370)
mypalette[1] = "white"
mycolors = seq(0, 0.37, 0.001)

## plot heatmap (Figure 3-A)
heatmap.2(heat_plot2,
          Colv = NA,
          Rowv = NA,
          dendrogram = "none",
          col=mypalette,
          breaks=mycolors,
          cexRow = 1,
          cexCol = 1.1,
          scale = "none",
          margins = c(5, 20),
          density.info = c("none"),
          labRow=as.expression(lapply(taxa, function(a) bquote(italic(.(a))))),
          colRow = rsn$phy_col,
          adjCol = 0.8,
          key.title = NA,
          tracecol = NA,
          linecol= "gray")

# legend
dev.off()
plot(1:50, col="white")
legend(1,49, cex = 1, legend = c("Actinobacteriota", "Bacteroidota", "Cyanobacteria", "Firmicutes", "Proteobacteria"),
       c("#663300", "#FF3333", "purple", "#3399FF", "#006600")) 

# plot key for heat colors
dev.off()
hist(c(0:350), breaks = 350, col = mypalette, 
     border = mypalette, ylim = c(0,0.5),
     xlab = NULL, ylab = NULL, 
     axes = FALSE,
     main = NA)
axis(side = 1, labels = c("0", "7", "14", "21", "28", "35"),
     at = c(0, 70, 140, 210, 280, 350), las = 1, cex = 2)


### Check ABX controls

subset = map_unique_ra[grep(paste(rownames(heat_plot2), collapse="|"), rownames(map_unique_ra)), grep("hvskin_abx", rag_hv_meta$Subject_group)]
apply(subset, 1, max)[order(apply(subset, 1, max))]

sub_avg = as.data.frame(matrix(nrow=nrow(subset), ncol=8))
for (i in 1:nrow(subset)) {
  sub_avg[i,] = tapply(as.numeric(subset[i,]), paste(rag_hv_meta$Subject_ID, rag_hv_meta$Longitudinal_ABX_or_HSCT, sep="_")[grep("hvskin_abx", rag_hv_meta$Subject_group)], mean)

}
colnames(sub_avg) = names(tapply(as.numeric(subset[1,]), paste(rag_hv_meta$Subject_ID, rag_hv_meta$Longitudinal_ABX_or_HSCT, sep="_")[grep("hvskin_abx", rag_hv_meta$Subject_group)], mean))
rownames(sub_avg) = rownames(subset)

sub_avg[which(apply(sub_avg, 1, max) > 0.001),]
sub_avg[which(apply(sub_avg[,c(1,2,4,5,7,8)], 1, max) > 0.001),]
sub_avg[which(apply(sub_avg, 1, mean) > 0.001),]
sub_avg[which(apply(sub_avg[,c(1,2,4,5,7,8)], 1, mean) > 0.001),]

sub_avg[which(apply(sub_avg[,c(1,4,7)], 1, mean) > 0.001),]
sub_avg[which(apply(sub_avg[,c(2,5,8)], 1, mean) > 0.001),]

map_unique_ra[grep(paste(rownames(sub_avg[which(apply(sub_avg[,c(1,4,7)], 1, mean) > 0.001),]), collapse = "|"), rownames(map_unique_ra)),171:182]
map_unique_ra[grep(paste(rownames(sub_avg[which(apply(sub_avg[,c(1,2,4,5,7,8)], 1, mean) > 0.001),]), collapse = "|"), rownames(map_unique_ra)),171:182]


######################################
### Read Mapping: Stool Microbiota ###
######################################

###
### Select taxa for plot
###

map_ra_species$avg_stool_B = apply(map_ra_species[,intersect(which(is.na(rag_hv_meta$HSCT_post_month) == TRUE),
                                                             grep("St", rag_hv_meta$Body_site_ID))], 1, mean)

map_ra_species = map_ra_species[order(map_ra_species$avg_stool_B, decreasing = T),]
map_ra_species[c(1:15), c(141:145)]
map_ra_species[c(1:6, 11, 13), c(141:145)]
grep("Escheric", rownames(map_ra_species))

## subset
map_taxa_select = rbind(apply(map_ra_species[grep("d__Archaea", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("d__Bacteria", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("d__Eukaryota", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("d__Virus", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Actinobacteriota", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Bacteroidota", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Firmicutes", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Proteobacteria", rownames(map_ra_species)),], 2, sum),
                        map_ra_species[c(1:6, 11, 13),])
#map_ra_species[c(1:9),])

rownames(map_taxa_select) = c("d__Archaea", 
                              "d__Bacteria", 
                              "d__Eukaryota", 
                              "d__Virus",
                              "d__Bacteria p__Actinobacteriota", "d__Bacteria p__Bacteroidota", 
                              "d__Bacteria p__Firmicutes", "d__Bacteria p__Proteobacteria", 
                              rownames(map_ra_species[c(1:6, 11, 13),]))
head(map_taxa_select)
map_taxa_select
rownames(map_taxa_select)

# virus
map_taxa_select[grep("^d__Virus$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Virus$", rownames(map_taxa_select)),] - map_taxa_select[grep("Sipho", rownames(map_taxa_select)),]
# phyla
map_taxa_select[grep("^d__Bacteria p__Bacteroidota$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria p__Bacteroidota$", rownames(map_taxa_select)),] - apply(map_taxa_select[grep("f__Bacteroidaceae|f__Tannerellaceae|f__Rikenell", rownames(map_taxa_select)),], 2, sum)
map_taxa_select[grep("^d__Bacteria p__Firmicutes$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria p__Firmicutes$", rownames(map_taxa_select)),] - apply(map_taxa_select[grep("f__Ruminococcaceae|f__Lachnospiraceae|Enterococcus", rownames(map_taxa_select)),], 2, sum)
map_taxa_select[grep("^d__Bacteria p__Proteobacteria$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria p__Proteobacteria$", rownames(map_taxa_select)),] - apply(map_taxa_select[grep("Escher", rownames(map_taxa_select)),], 2, sum)
# bacteria
map_taxa_select[grep("^d__Bacteria$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria$", rownames(map_taxa_select)),] - apply(map_taxa_select[-c(grep("^d__Bacteria$|d__Eukaryota|d__Virus|d__Archaea", rownames(map_taxa_select))),], 2, sum)

# check sums
apply(map_taxa_select, 2, sum)

## RA by subject
mpa_ra = map_taxa_select
mpa_ra = mpa_ra[order(rownames(mpa_ra), decreasing=T),]
rownames(mpa_ra)
# combine bac, euks and archaea
mpa_ra[15,] = mpa_ra[15,] + mpa_ra[16,] + mpa_ra[3,]
mpa_ra = mpa_ra[-c(16,3),]
apply(mpa_ra, 2, sum)
dim(mpa_ra)
rownames(mpa_ra)

mpa_bar = mpa_ra[,intersect(grep("St", rag_hv_meta$Body_site_ID), 
                            which(is.na(rag_hv_meta$HSCT_post_month) == TRUE))]
mpa_bar$taxa = c("  Siphoviridae", "Virus/Phage", 
                 "Proteobacteria", 
                 "  Faecalibacterium prausnitzii_G", "  Blautia_A wexlerae", "Firmicutes", 
                 "  Parabacteroides distonasis", "  Alistipes onderdonkii", "  Phocaeicola dorei", 
                 "  Bacteroides xylanisolvens", "  Bacteroides uniformis", "Bacteroidota",
                 "Actinobacteriota", 
                 "Other")
colnames(mpa_bar) = rag_hv_meta$Subject_ID[intersect(grep("St", rag_hv_meta$Body_site_ID), 
                                                     which(is.na(rag_hv_meta$HSCT_post_month) == TRUE))]
mpa_bar = data.frame(melt(mpa_bar))
colnames(mpa_bar)[1] = "taxon"
head(mpa_bar, 30)

# add group type
mpa_bar$group_id = rep("HV", nrow(mpa_bar))
mpa_bar$group_id[grep("R", mpa_bar$variable)] = "RAG"

# taxa factor
mpa_bar$taxon = factor(mpa_bar$taxon,
                       levels = rev(c("  Siphoviridae", "Virus/Phage", 
                                      "Proteobacteria", 
                                      "  Faecalibacterium prausnitzii_G", "  Blautia_A wexlerae", "Firmicutes", 
                                      "  Alistipes onderdonkii", "  Parabacteroides distonasis", "  Phocaeicola dorei",
                                      "  Bacteroides xylanisolvens", "  Bacteroides uniformis", "Bacteroidota",
                                      "Actinobacteriota", 
                                      "Other")))

mpa_bar$group_id = factor(mpa_bar$group_id,
                          levels = c("RAG", "HV"))

# rename factors for axis
mpa_bar$variable = as.character(mpa_bar$variable)

# plot (Figure 2-B)
head(mpa_bar)
ggplot(mpa_bar,
       aes(x=variable, y=value, fill=taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Subject ID") +
  facet_grid(~group_id, scales="free", space="free") +
  scale_fill_manual(values=rev(c("#CCCCFF", "#9933FF",
                                 "#006600", # greens
                                 "#CCFFFF", "#66CCFF", "blue", # blues 
                                 "#FFFFCC", "#FFCCCC", "#FF9999", "#FF3333", "#CC0033", "darkred", #reds
                                 #"#663300", #ET
                                 "#996633",
                                 "lightgray")),
                    labels=rev(c("*Siphoviridae*", "Virus/Phage", 
                                 "Proteobacteria", 
                                 "*Faecalibacterium prausnitzii_G*", "*Blautia_A wexlerae*", "Firmicutes", 
                                 "*Alistipes onderdonkii*", "*Parabacteroides distonasis*", "*Phocaeicola dorei*",
                                 "*Bacteroides xylanisolvens*", "*Bacteroides uniformis*", "Bacteroidota",
                                 "Actinobacteriota", 
                                 "Other"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90, vjust=0.5),
        axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 13)) +
  guides(fill=guide_legend(ncol=1)) +
  theme(legend.position = "none")


###
### Time Series
###

###
### Select taxa for plot (HSCT)
###

map_ra_species$avg_stool_H = apply(map_ra_species[,intersect(which(is.na(rag_hv_meta$HSCT_post_month) == FALSE),
                                                             grep("St", rag_hv_meta$Body_site_ID))], 1, mean)

map_ra_species = map_ra_species[order(map_ra_species$avg_stool_H, decreasing = T),]
map_ra_species[c(1:9), c(141:145)]
grep("Escheric", rownames(map_ra_species))

## subset
map_taxa_select = rbind(apply(map_ra_species[grep("d__Archaea", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("d__Bacteria", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("d__Eukaryota", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("d__Virus", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Actinobacteriota", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Bacteroidota", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Firmicutes", rownames(map_ra_species)),], 2, sum),
                        apply(map_ra_species[grep("p__Proteobacteria", rownames(map_ra_species)),], 2, sum),
                        map_ra_species[c(1:9),])

rownames(map_taxa_select) = c("d__Archaea", 
                              "d__Bacteria", 
                              "d__Eukaryota", 
                              "d__Virus",
                              "d__Bacteria p__Actinobacteriota", "d__Bacteria p__Bacteroidota", 
                              "d__Bacteria p__Firmicutes", "d__Bacteria p__Proteobacteria", 
                              rownames(map_ra_species[c(1:9),]))
head(map_taxa_select)
map_taxa_select
rownames(map_taxa_select)

# virus
map_taxa_select[grep("^d__Virus$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Virus$", rownames(map_taxa_select)),] - map_taxa_select[grep("Sipho", rownames(map_taxa_select)),]
# phyla
map_taxa_select[grep("^d__Bacteria p__Bacteroidota$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria p__Bacteroidota$", rownames(map_taxa_select)),] - apply(map_taxa_select[grep("f__Bacteroidaceae|f__Tannerellaceae|f__Rikenell", rownames(map_taxa_select)),], 2, sum)
map_taxa_select[grep("^d__Bacteria p__Firmicutes$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria p__Firmicutes$", rownames(map_taxa_select)),] - apply(map_taxa_select[grep("f__Ruminococcaceae|f__Lachnospiraceae|Enterococcus", rownames(map_taxa_select)),], 2, sum)
map_taxa_select[grep("^d__Bacteria p__Proteobacteria$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria p__Proteobacteria$", rownames(map_taxa_select)),] - apply(map_taxa_select[grep("Escher", rownames(map_taxa_select)),], 2, sum)
# bacteria
map_taxa_select[grep("^d__Bacteria$", rownames(map_taxa_select)),] = map_taxa_select[grep("^d__Bacteria$", rownames(map_taxa_select)),] - apply(map_taxa_select[-c(grep("^d__Bacteria$|d__Eukaryota|d__Virus|d__Archaea", rownames(map_taxa_select))),], 2, sum)

# check sums
apply(map_taxa_select, 2, sum)

## RA by subject
mpa_ra = map_taxa_select
mpa_ra = mpa_ra[order(rownames(mpa_ra), decreasing=T),]
rownames(mpa_ra)
# combine bac, euks and archaea
mpa_ra[16,] = mpa_ra[16,] + mpa_ra[17,] + mpa_ra[3,]
mpa_ra = mpa_ra[-c(17,3),]
apply(mpa_ra, 2, sum)
dim(mpa_ra)
rownames(mpa_ra)

# set taxa
stool_taxa = mpa_ra[,grep("St", rag_hv_meta$Body_site_ID)]
head(stool_taxa)
rownames(stool_taxa)
rownames(stool_taxa) = c("  Siphoviridae", "Virus/Phage", 
                         "  Escherichia coli", "Proteobacteria", 
                         "  Gemmiger qucibialis", "  Enterocloster bolteae", "  Agathobacter rectalis", "  Enterococcus faecalis", "Firmicutes", 
                         "  Parabacteroides distonasis", "  Phocaeicola dorei", "  Bacteroides fragilis", "Bacteroidota",
                         "Actinobacteriota", 
                         "Other")

# set meta
stool_meta = rag_hv_meta[grep("St", rag_hv_meta$Body_site_ID),]
paste(stool_meta$Sample, "_Counts", sep = "") == colnames(stool_taxa)

## prep data

R1 = stool_taxa[, grep("RAG_pt01", stool_meta$Subject_ID)]
colnames(R1) = stool_meta$HSCT_post_month[grep("RAG_pt01", stool_meta$Subject_ID)]
R1 = melt(t(R1))

R5 = stool_taxa[, grep("RAG_pt05", stool_meta$Subject_ID)]
colnames(R5) = stool_meta$HSCT_post_month[grep("RAG_pt05", stool_meta$Subject_ID)]
R5 = melt(t(R5))

R6 = stool_taxa[, grep("RAG_pt06", stool_meta$Subject_ID)]
colnames(R6) = stool_meta$HSCT_post_month[grep("RAG_pt06", stool_meta$Subject_ID)]
R6 = melt(t(R6))

subject_plots = as.data.frame(rbind(R1, R5, R6))
subject_plots$subject = c(rep("RAG_pt01 (child)", dim(R1)[1]),
                          rep("RAG_pt05 (adult)", dim(R5)[1]),
                          rep("RAG_pt06 (adult)", dim(R6)[1]))

# set factors
subject_plots$Var2 = factor(subject_plots$Var2,
                            levels = rev(c("  Siphoviridae", "Virus/Phage", 
                                           "  Escherichia coli", "Proteobacteria", 
                                           "  Enterocloster bolteae", "  Enterococcus faecalis", "  Agathobacter rectalis", "  Gemmiger qucibialis", "Firmicutes", 
                                           "  Parabacteroides distonasis", "  Phocaeicola dorei", "  Bacteroides fragilis", "Bacteroidota",
                                           "Actinobacteriota", 
                                           "Other")))

# add relative time variables
subject_plots$Var1[which(is.na(subject_plots$Var1) == TRUE)] = -1
subject_plots$Var1[grep("-1", subject_plots$Var1)] = "A"
subject_plots$Var1[grep("0", subject_plots$Var1)] = "B"
subject_plots$Var1[grep("1", subject_plots$Var1)] = "C"
subject_plots$Var1[grep("2", subject_plots$Var1)] = "D"
subject_plots$Var1[grep("3", subject_plots$Var1)] = "E"
subject_plots$Var1[grep("6", subject_plots$Var1)] = "F"
subject_plots$Var1[grep("A", subject_plots$Var1)] = 1
subject_plots$Var1[grep("B", subject_plots$Var1)] = 2
subject_plots$Var1[grep("C", subject_plots$Var1)] = 3
subject_plots$Var1[grep("D", subject_plots$Var1)] = 4
subject_plots$Var1[grep("E", subject_plots$Var1)] = 5
subject_plots$Var1[grep("F", subject_plots$Var1)] = 6

# plot (Figure S8)
ggplot(subject_plots,
       aes(x=as.numeric(Var1), 
           y=value, fill=Var2)) +
  #geom_area(alpha=0.9, size = 0.2, colour="black") +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  theme_bw() +
  xlab("Time (months)") +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6),
                     labels = c("B", "0", "1", "2", "3", "6")) +
  scale_fill_manual(values=rev(c("#CCCCFF", "#9933FF",
                                 #"#CCFFCC", 
                                 "#66FF66", "#006600", # greens
                                 "#CCFFFF", "#66CCFF", "#3399FF", "blue", "darkblue", # blues 
                                 "#FFCCCC", "#FF9999", "#FF3333", "darkred", #reds
                                 #"#663300", #ET
                                 "#996633",
                                 "lightgray")),
                    labels=rev(c("*Siphoviridae*", "Virus/Phage", 
                                 "*Escherichia coli*", "Proteobacteria", 
                                 "*Enterocloster bolteae*", "*Enterococcus faecalis*", "*Agathobacter rectalis*", "*Gemmiger qucibialis*", "Firmicutes", 
                                 "*Parabacteroides distonasis*", "*Phocaeicola dorei*", "*Bacteroides fragilis*", "Bacteroidota",
                                 "Actinobacteriota", 
                                 "Other"))) +
  facet_wrap(~subject, nrow=3) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 11)) +
  guides(fill=guide_legend(ncol=1)) #+
#theme(legend.position = "none")


###
### Beta diversity
###

## subset skin taxa
stool_taxa = map_ra_species[,grep("St", rag_hv_meta$Body_site_ID)]
head(stool_taxa)
# meta
stool_meta = rag_hv_meta[grep("St", rag_hv_meta$Body_site_ID),]

# remove HSCT points
stool_taxa = stool_taxa[,-c(grep("0|1|2|3|6", stool_meta$HSCT_post_month))]
stool_meta = stool_meta[-c(grep("0|1|2|3|6", stool_meta$HSCT_post_month)),]

## RAG skin microbiota: associations with subject and site
beta <- vegdist(t(stool_taxa))

# prep pcoa
pcoa <- cmdscale(beta, k=4, eig=T)
ord <- as.data.frame(pcoa$points)
names(ord) <- c("pcoa1", "pcoa2", "pcoa3", "pcoa4")
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# metadata
ord$Subject_ID = stool_meta$Subject_ID
ord$Subject_group = stool_meta$Subject_group
ord$age = stool_meta$Subject_age
ord$type = paste(ord$Subject_group, ord$age, sep = "_")

# plot pcoa (Figure 2-A)
ggplot(data=ord, 
       aes(x=pcoa1, y=pcoa2, fill=type, col=type, shape=type
       )) +
  geom_point(size=4, stroke=1.2, alpha=0.9) +
  theme_bw() +
  scale_fill_manual(values=c("blue", "red", "red"),
                    labels=c("HV (adult)", "RAG (adult)", "RAG (child)")) +
  scale_color_manual(values=c("blue", "red", "red"),
                     labels=c("HV (adult)", "RAG (adult)", "RAG (child)")) +
  scale_shape_manual(values = c(21, 21, 1)) +
  xlab("PCoA1 (18.8%)") +
  ylab("PCoA2 (16.1%)") +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank())

# PERMANOVAs
adonis2(beta ~ Subject_group, stool_meta, permutations=999)

beta <- vegdist(t(stool_taxa[,grep("adult", stool_meta$Subject_age)]))
adonis2(beta ~ Subject_group, stool_meta[grep("adult", stool_meta$Subject_age),], permutations=999)

beta <- vegdist(t(stool_taxa[,grep("rag", stool_meta$Subject_group)]))
adonis2(beta ~ Subject_age, stool_meta[grep("rag", stool_meta$Subject_group),], permutations=999)


###
### Taxa enrichment
###

# re-set order to match skin microbiome taxa
map_ra_species = map_ra_species[order(map_ra_species$avg_stool_B, decreasing = T),]
head(map_ra_species)

## MW test & t-test for significance
test = data.frame(t(map_ra_species[,1:159]), rag_hv_meta)
p_vals = as.data.frame(matrix(nrow=nrow(map_ra_species), ncol=2))
for (i in 1:nrow(map_ra_species)) {
  p_vals[i, 1] = t.test(test[intersect(intersect(grep("St", rag_hv_meta$Body_site_ID), grep("B", rag_hv_meta$Longitudinal_ABX_or_HSCT)), grep("adult", rag_hv_meta$Subject_age)), i] ~ Subject_group,
                        test[intersect(intersect(grep("St", rag_hv_meta$Body_site_ID), grep("B", rag_hv_meta$Longitudinal_ABX_or_HSCT)), grep("adult", rag_hv_meta$Subject_age)), ])$p.val
  p_vals[i, 2] = wilcox.test(test[intersect(intersect(grep("St", rag_hv_meta$Body_site_ID), grep("B", rag_hv_meta$Longitudinal_ABX_or_HSCT)), grep("adult", rag_hv_meta$Subject_age)), i] ~ Subject_group,
                             test[intersect(intersect(grep("St", rag_hv_meta$Body_site_ID), grep("B", rag_hv_meta$Longitudinal_ABX_or_HSCT)), grep("adult", rag_hv_meta$Subject_age)), ])$p.val
}
colnames(p_vals) = c("t_test", "mw_test")
rownames(p_vals) = rownames(map_ra_species)
# result
length(which(p_vals[,2] < 0.05))
p_vals[which(p_vals[,2] < 0.05),]
# result - bh correction
length(which(p.adjust(p_vals[,2], method = "BH") < 0.1))
rownames(p_vals)[which(p.adjust(p_vals[,2], method = "BH") < 0.1)]


##########################
### RAG Stool: Novelty ###
##########################

## set rag skin novelty data frame

# call RAG-patient-derived MAGs
rsn = map_unique_ra[grep("d__Bacteria", map_unique_ra$Domain),]
rsn = rsn[order(rsn$avg_all, decreasing = T), ]
#rsn = rsn[-c(grep("SMGC|hv", rownames(rsn))),]
#rsn = rsn[c(grep("rag", rownames(rsn))),]

# subset to matching sites with HV/HC
heat_plot = as.matrix(rsn[,c(grep("control", rag_hv_meta$Subject_group),
                             intersect(grep("St", rag_hv_meta$Body_site_ID), 
                                       which(is.na(rag_hv_meta$HSCT_post_month) == TRUE)))])
head(heat_plot)

# calculate average
heat_plot2 = heat_plot
colnames(heat_plot2) = rag_hv_meta$Subject_ID[c(grep("control", rag_hv_meta$Subject_group),
                                                intersect(grep("St", rag_hv_meta$Body_site_ID), 
                                                          which(is.na(rag_hv_meta$HSCT_post_month) == TRUE)))]
rownames(heat_plot2) = rownames(heat_plot)

# phylum text colors
levels(factor(rsn$Phylum))
rsn$phy_col = rep("purple", length(rsn$Phylum))
rsn$phy_col[grep("Actin", rsn$Phylum)] = "#663300"
rsn$phy_col[grep("Bactero", rsn$Phylum)] = "#FF3333"
rsn$phy_col[grep("Verr", rsn$Phylum)] = "#FFCC33"
rsn$phy_col[grep("Firm", rsn$Phylum)] = "#3399FF"
rsn$phy_col[grep("Prot", rsn$Phylum)] = "#006600"

# adjust row removal
rownames(rsn) == rownames(heat_plot2)
colnames(heat_plot2)

# remove species that are >0.1% abundant in HV/HC and controls
rsn = rsn[which(apply(heat_plot2[,c(9:14)], 1, mean)<0.001),]
heat_plot2 = heat_plot2[which(apply(heat_plot2[,c(9:14)], 1, mean)<0.001),]

rsn = rsn[which(apply(heat_plot2[,c(1:8)], 1, mean)<0.001),]
heat_plot2 = heat_plot2[which(apply(heat_plot2[,c(1:8)], 1, mean)<0.001),]

# remove species that are <0.1% abundant in RAG-skin samples
rsn = rsn[which(apply(heat_plot2[,c(15:21)], 1, mean)>0.001),]
heat_plot2 = heat_plot2[which(apply(heat_plot2[,c(15:21)], 1, mean)>0.001),]

# re-order plot
heat_plot2 = heat_plot2[,15:21]
colnames(heat_plot2)
heat_plot2 = heat_plot2[,c(1,7,2,5,6,3,4)]
colnames(heat_plot2)

heat_plot2 = heat_plot2[order(paste(rsn$Genus, rsn$Species, sep = "|")),]
rsn = rsn[order(paste(rsn$Genus, rsn$Species, sep = "|")),]

heat_plot2 = heat_plot2[order(rsn$Phylum),]
rsn = rsn[order(rsn$Phylum),]

# set color for MAG source
rsn$source_col = rep("lightgray", length(rsn$phy_col))
rsn$source_col[grep("rag", rownames(rsn))] = "black"

# set legend
max(heat_plot2)
quantile(heat_plot2)
sort(heat_plot2, decreasing = T)
mypalette = colorRampPalette(c("#FFCCFF", "#9900CC", "#660099", "#330066", "#000033", "black"))(n=370)
mypalette[1] = "white"
mycolors = seq(0, 0.37, 0.001)

# set taxa
taxa = substr(rsn$Species, start=4, stop=nchar(rsn$Species))
taxa[which(taxa == "")] = paste(substr(rsn$Genus, start=4, stop=nchar(rsn$Genus))[which(taxa == "")], "sp.")
taxa[which(taxa == " sp.")] = paste(substr(rsn$Family, start=4, stop=nchar(rsn$Family))[which(taxa == " sp.")], "sp.")

## plot heatmap (Figure 3-B)
heatmap.2(heat_plot2,
          Colv = NA,
          Rowv = NA,
          dendrogram = "none",
          col=mypalette,
          breaks=mycolors,
          cexRow = 0.8,
          cexCol = 1.1,
          scale = "none",
          margins = c(5, 20),
          density.info = c("none"),
          labRow=as.expression(lapply(taxa, function(a) bquote(italic(.(a))))),
          colRow = rsn$phy_col,
          adjCol = 0.8,
          key.title = NA,
          tracecol = NA,
          linecol= "gray")


###########################
##### ARGs: Shortbred #####
###########################

## shortbred counts table
ARG_counts = read.csv("data/ARGs/shortbred_ARGs_rag_hv.csv", row.names=1, h=T)
ARG_counts = ARG_counts[,c(order(colnames(ARG_counts)[1:159]), 160:161)]

## set metadata
rag_hv_meta = read.csv("data/metadata/rag_hv_metadata.csv", h=T, row.names=1)
rag_hv_meta = rag_hv_meta[grep("DNA", rag_hv_meta$Metagenome_type),]

# check order
substr(colnames(ARG_counts)[1:159], start = 1, stop = nchar(colnames(ARG_counts)[1:159])-6) == rag_hv_meta$Sample_name

# subset non_abx
ARG_counts = ARG_counts[,c(intersect(grep("DNA", rag_hv_meta$Metagenome_type), grep("hv_|rag_", rag_hv_meta$Subject_group)), 160:161)]
rag_hv_meta = rag_hv_meta[intersect(grep("DNA", rag_hv_meta$Metagenome_type), grep("hv_|rag_", rag_hv_meta$Subject_group)),]

# check
substr(colnames(ARG_counts)[1:127], start = 1, stop = nchar(colnames(ARG_counts)[1:127])-6) == rag_hv_meta$Sample_name

## subset ARGs with at least 30 AA seqs
ARG_counts = ARG_counts[which(ARG_counts$TotMarkerLength >= 30),]
rownames(ARG_counts) = ARG_counts$ARG_family
ARG_counts = ARG_counts[,1:127]
head(ARG_counts)
dim(ARG_counts)

## total ARGs detected by group

# RAG-skin: all sites
length(which(apply(ARG_counts[,grep("rag_skin_child", paste(rag_hv_meta$Subject_group, rag_hv_meta$Subject_age, sep="_"))], 1, sum) > 0))
length(which(apply(ARG_counts[,grep("rag_skin_adult", paste(rag_hv_meta$Subject_group, rag_hv_meta$Subject_age, sep="_"))], 1, sum) > 0))

# RAG-skin: matching sites
length(which(apply(ARG_counts[intersect(grep("rag_skin_child", paste(rag_hv_meta$Subject_group, rag_hv_meta$Subject_age, sep="_")),
                                        grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID))], 1, sum) > 0))
length(which(apply(ARG_counts[intersect(grep("rag_skin_adult", paste(rag_hv_meta$Subject_group, rag_hv_meta$Subject_age, sep="_")),
                                        grep("Ac|Vf|Ra|N", rag_hv_meta$Body_site_ID))], 1, sum) > 0))

# HV/HC-skin
length(which(apply(ARG_counts[intersect(grep("hv_skin_child", paste(rag_hv_meta$Subject_group, rag_hv_meta$Subject_age, sep="_")),
                                        grep("Ac|Vf|Ra", rag_hv_meta$Body_site_ID))], 1, sum) > 0))
length(which(apply(ARG_counts[intersect(grep("hv_skin_adult", paste(rag_hv_meta$Subject_group, rag_hv_meta$Subject_age, sep="_")),
                                        grep("Ac|Vf|Ra|N", rag_hv_meta$Body_site_ID))], 1, sum) > 0))

# RAG-stool
length(which(apply(ARG_counts[,grep("rag_stool_adult", paste(rag_hv_meta$Subject_group, rag_hv_meta$Subject_age, sep="_"))], 1, sum) > 0))
length(which(apply(ARG_counts[,intersect(grep("rag_stool_adult", paste(rag_hv_meta$Subject_group, rag_hv_meta$Subject_age, sep="_")),
                                         which(is.na(rag_hv_meta$HSCT_post_month) == TRUE))], 1, sum) > 0))

# HV-stool
length(which(apply(ARG_counts[,grep("hv_stool_adult", paste(rag_hv_meta$Subject_group, rag_hv_meta$Subject_age, sep="_"))], 1, sum) > 0))


###
### ARG enrichment-adult (skin and nares)
###

## clean table for respective LME model (adult)
ARG_counts = ARG_counts[-c(which(apply(ARG_counts[,intersect(grep("Ac|Ra|Vf|N", rag_hv_meta$Body_site_ID), 
                                                             grep("adult", rag_hv_meta$Subject_age))], 1, sum) == 0)),]
dim(ARG_counts)

# re-order
ARG_counts = ARG_counts[order(apply(ARG_counts[,grep("Ac|Ra|Vf|N", rag_hv_meta$Body_site_ID)], 1, sum), decreasing = T),]

# log transform
ARG_counts = log10(ARG_counts+1)

## LME model -- RAG vs. HV
test = data.frame(t(ARG_counts[,1:127]), rag_hv_meta)
p_vals = rep(1, nrow(ARG_counts))
for (i in 1:nrow(ARG_counts)) {
  p_vals[i] = summary(lmer(test[intersect(grep("Ac|Vf|Ra|N", test$Body_site_ID), grep("adult", test$Subject_age)),i] ~ Subject_group + (1|Body_site_ID),
                           test[intersect(grep("Ac|Vf|Ra|N", test$Body_site_ID), grep("adult", test$Subject_age)),]))$coefficients[2,5]
  }
names(p_vals) = rownames(ARG_counts)
# result
length(which(p_vals < 0.05))
p_vals[which(p_vals < 0.05)]

## plot enrichment

# load table
ARG_sub = read.csv("data/ARGs/ARGs_skin_adult_manual_edit.csv", row.names=1, h=T)
rownames(ARG_counts[which(p_vals < 0.05),]) == ARG_sub$family
ARG_sub = data.frame(ARG_counts[which(p_vals < 0.05), grep("Ac|Vf|Ra|N", rag_hv_meta$Body_site_ID)], 
                     ARG_sub[,c(2:4)])
ARG_sub = melt(ARG_sub)

# add metadata
ARG_sub$Subject_group = rep("HV", nrow(ARG_sub))
ARG_sub$Subject_group[grep(paste(rag_hv_meta$Sample[grep("rag_skin", rag_hv_meta$Subject_group)], collapse = "|"),
                        ARG_sub$variable)] = "RAG"
ARG_sub$Site_ID = rep(1, nrow(ARG_sub))

ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Ac", rag_hv_meta$Body_site_ID), grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 22
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Ac", rag_hv_meta$Body_site_ID), grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 0
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Vf", rag_hv_meta$Body_site_ID), grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 25
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Vf", rag_hv_meta$Body_site_ID), grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 6
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Ra", rag_hv_meta$Body_site_ID), grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 24
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Ra", rag_hv_meta$Body_site_ID), grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 2
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("N", rag_hv_meta$Body_site_ID), grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 23
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("N", rag_hv_meta$Body_site_ID), grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 5

# order gene layout
ARG_sub$gene = factor(ARG_sub$gene,
                      levels = c("APH(6)-Id", "ANT(2)-Ia",
                                 "murA_i", "murA_ii",
                                 "ErmC", 
                                 "mecI", "mecR1",
                                 "cmx",
                                 "fusB", "mupA", "SAT-4",
                                 "gyrA_i", "gyrB",
                                 "mtrA", "norA",
                                 "EF-Tu"
                                 ))

# plot skin ARGs - adult (Figure 4-B)
ggplot(ARG_sub[grep(paste(rag_hv_meta$Sample[intersect(grep("skin", rag_hv_meta$Subject_group),
                                                       grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                    ARG_sub$variable),],
       aes(x=gene, y=value, col=Subject_group, fill=Subject_group)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0, outlier.color = "white",
               alpha=0.4) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0),
              size=2.5, #col = "black", 
              pch=ARG_sub[grep(paste(rag_hv_meta$Sample[intersect(grep("skin", rag_hv_meta$Subject_group),
                                                                  grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                               ARG_sub$variable),]$Site_ID,
              alpha=0.7) +
  theme_bw() +
  ylim(0,4.1) +
  ylab("ARG Family log(RPKM)") +
  scale_fill_manual(values=c("red","blue")) +
  scale_color_manual(values=c("red","blue")) +
  theme(axis.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.text = element_text(size=14),
        axis.text.x = element_text(size=14, angle=45, hjust=1,
                                   color = c("#FF3333", "#FF3333", 
                                             "skyblue", "skyblue",
                                             "violet", 
                                             "orange", "orange", 
                                             "#666600",
                                             "darkgray", "darkgray", "darkgray",
                                             "darkgreen", "darkgreen",
                                             "purple", "purple", 
                                             "darkgray"
                                   )),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.ticks.x = element_line(size = 2.5, lineend = "round", 
                                    color = c("#FF3333", "#FF3333", 
                                              "skyblue", "skyblue",
                                              "violet", 
                                              "orange", "orange", 
                                              "#666600",
                                              "darkgray", "darkgray", "darkgray",
                                              "darkgreen", "darkgreen",
                                              "purple", "purple", 
                                              "darkgray"
                                    )))


###
### ARG enrichment-child (skin)
###

## shortbred counts table
ARG_counts = read.csv("data/ARGs/shortbred_ARGs_rag_hv.csv", row.names=1, h=T)
ARG_counts = ARG_counts[,c(order(colnames(ARG_counts)[1:159]), 160:161)]

## set metadata
rag_hv_meta = read.csv("data/metadata/rag_hv_metadata.csv", h=T, row.names=1)
rag_hv_meta = rag_hv_meta[grep("DNA", rag_hv_meta$Metagenome_type),]

# check order
substr(colnames(ARG_counts)[1:159], start = 1, stop = nchar(colnames(ARG_counts)[1:159])-6) == rag_hv_meta$Sample_name

# subset non_abx
ARG_counts = ARG_counts[,c(intersect(grep("DNA", rag_hv_meta$Metagenome_type), grep("hv_|rag_", rag_hv_meta$Subject_group)), 160:161)]
rag_hv_meta = rag_hv_meta[intersect(grep("DNA", rag_hv_meta$Metagenome_type), grep("hv_|rag_", rag_hv_meta$Subject_group)),]

# check
substr(colnames(ARG_counts)[1:127], start = 1, stop = nchar(colnames(ARG_counts)[1:127])-6) == rag_hv_meta$Sample_name

## subset ARGs with at least 30 AA seqs
ARG_counts = ARG_counts[which(ARG_counts$TotMarkerLength >= 30),]
rownames(ARG_counts) = ARG_counts$ARG_family
ARG_counts = ARG_counts[,1:127]
head(ARG_counts)
dim(ARG_counts)

## clean table for respective LME model (child)
ARG_counts = ARG_counts[-c(which(apply(ARG_counts[,intersect(grep("Ac|Ra|Vf|N", rag_hv_meta$Body_site_ID), 
                                                             grep("child", rag_hv_meta$Subject_age))], 1, sum) == 0)),]
dim(ARG_counts)

# re-order
ARG_counts = ARG_counts[order(apply(ARG_counts[,grep("Ac|Ra|Vf|N", rag_hv_meta$Body_site_ID)], 1, sum), decreasing = T),]

# log transform
ARG_counts = log10(ARG_counts+1)

## LME model -- RAG vs. HV
test = data.frame(t(ARG_counts[,1:127]), rag_hv_meta)
p_vals = rep(1, nrow(ARG_counts))
for (i in 1:nrow(ARG_counts)) {
  p_vals[i] = summary(lmer(test[intersect(grep("Ac|Vf|Ra|N", test$Body_site_ID), grep("child", test$Subject_age)),i] ~ Subject_group + (1|Body_site_ID),
                           test[intersect(grep("Ac|Vf|Ra|N", test$Body_site_ID), grep("child", test$Subject_age)),]))$coefficients[2,5]
}
names(p_vals) = rownames(ARG_counts)
# result
length(which(p_vals < 0.05))
p_vals[which(p_vals < 0.05)]

## plot enrichment

# load table
ARG_sub = read.csv("data/ARGs/ARGs_skin_ped_manual_edit.csv", row.names=1, h=T)
rownames(ARG_counts[which(p_vals < 0.05),]) == ARG_sub$family
ARG_sub = data.frame(ARG_counts[which(p_vals < 0.05), grep("Ac|Vf|Ra|N", rag_hv_meta$Body_site_ID)], 
                     ARG_sub[,c(2:4)])
ARG_sub = melt(ARG_sub)

# add metadata
ARG_sub$Subject_group = rep("HV", nrow(ARG_sub))
ARG_sub$Subject_group[grep(paste(rag_hv_meta$Sample[grep("rag_skin", rag_hv_meta$Subject_group)], collapse = "|"),
                        ARG_sub$variable)] = "RAG"
ARG_sub$Site_ID = rep(1, nrow(ARG_sub))

ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Ac", rag_hv_meta$Body_site_ID), grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 22
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Ac", rag_hv_meta$Body_site_ID), grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 0
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Vf", rag_hv_meta$Body_site_ID), grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 25
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Vf", rag_hv_meta$Body_site_ID), grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 6
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Ra", rag_hv_meta$Body_site_ID), grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 24
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("Ra", rag_hv_meta$Body_site_ID), grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 2
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("N", rag_hv_meta$Body_site_ID), grep("adult", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 23
ARG_sub$Site_ID[grep(paste(rag_hv_meta$Sample[intersect(grep("N", rag_hv_meta$Body_site_ID), grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                     ARG_sub$variable)] = 5

# order ARGs for plot
ARG_sub$gene = factor(ARG_sub$gene,
                      levels = c("gyrA_i", "gyrA_ii",
                                 "ErmB", "ErmF", "ErmX", "mefA", "mphC",
                                 "CfxA6", "mecI", "mecR1", 
                                 "qacA", "mupA", "CpxR"))

# set factor for groups
ARG_sub$Subject_group = factor(ARG_sub$Subject_group,
                            levels = c("RAG", "HV"))

# plot skin ARGs - adult (Figure 4-B)
ggplot(ARG_sub[grep(paste(rag_hv_meta$Sample[intersect(grep("skin", rag_hv_meta$Subject_group),
                                                       grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                    ARG_sub$variable),],
       aes(x=gene, y=value, col=Subject_group, fill=Subject_group)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0, outlier.color = "white",
               alpha=0.4) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0),
              size=2.5, #col = "black", 
              pch=ARG_sub[grep(paste(rag_hv_meta$Sample[intersect(grep("skin", rag_hv_meta$Subject_group),
                                                                  grep("child", rag_hv_meta$Subject_age))], collapse = "|"),
                               ARG_sub$variable),]$Site_ID,
              alpha=0.7) +
  theme_bw() +
  ylim(0,4.1) +
  ylab("ARG Family log(RPKM)") +
  scale_fill_manual(values=c("red","blue")) +
  scale_color_manual(values=c("red","blue")) +
  theme(axis.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.text = element_text(size=14),
        axis.text.x = element_text(size=14, angle=45, hjust=1,
                                   color = c("darkgreen", "darkgreen", 
                                             "violet", "violet", "violet", "violet", "violet",
                                             "orange", "orange", "orange", 
                                             "purple", "darkgray", "purple")),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.ticks.x = element_line(size = 2.5, lineend = "round", 
                                    color = c("darkgreen", "darkgreen", 
                                              "violet", "violet", "violet", "violet", "violet",
                                              "orange", "orange", "orange", 
                                              "purple", "darkgray", "purple")))


###
### ARG enrichment (stool)
###

## shortbred counts table
ARG_counts = read.csv("data/ARGs/shortbred_ARGs_rag_hv.csv", row.names=1, h=T)
ARG_counts = ARG_counts[,c(order(colnames(ARG_counts)[1:159]), 160:161)]

## set metadata
rag_hv_meta = read.csv("data/metadata/rag_hv_metadata.csv", h=T, row.names=1)
rag_hv_meta = rag_hv_meta[grep("DNA", rag_hv_meta$Metagenome_type),]

# check order
substr(colnames(ARG_counts)[1:159], start = 1, stop = nchar(colnames(ARG_counts)[1:159])-6) == rag_hv_meta$Sample_name

# subset non_abx
ARG_counts = ARG_counts[,c(intersect(grep("DNA", rag_hv_meta$Metagenome_type), grep("hv_|rag_", rag_hv_meta$Subject_group)), 160:161)]
rag_hv_meta = rag_hv_meta[intersect(grep("DNA", rag_hv_meta$Metagenome_type), grep("hv_|rag_", rag_hv_meta$Subject_group)),]

# check
substr(colnames(ARG_counts)[1:127], start = 1, stop = nchar(colnames(ARG_counts)[1:127])-6) == rag_hv_meta$Sample_name

## subset ARGs with at least 30 AA seqs
ARG_counts = ARG_counts[which(ARG_counts$TotMarkerLength >= 30),]
rownames(ARG_counts) = ARG_counts$ARG_family
ARG_counts = ARG_counts[,1:133]
head(ARG_counts)
dim(ARG_counts)

## clean tables for comparison (adult baseline samples)=
# ARGs
ARG_counts = ARG_counts[,intersect(intersect(grep("St", rag_hv_meta$Body_site_ID), 
                                             which(is.na(rag_hv_meta$HSCT_post_month) == TRUE)), 
                                   grep("adult", rag_hv_meta$Subject_age))]
ARG_counts = ARG_counts[which(apply(ARG_counts, 1, sum) > 0),]
dim(ARG_counts)
# metadata
rag_hv_meta = rag_hv_meta[intersect(intersect(grep("St", rag_hv_meta$Body_site_ID), 
                                              which(is.na(rag_hv_meta$HSCT_post_month) == TRUE)), 
                                    grep("adult", rag_hv_meta$Subject_age)),]

## log transform
ARG_counts = log10(ARG_counts+1)

## MW test & t-test for significance
test = data.frame(t(ARG_counts), rag_hv_meta)
p_vals = as.data.frame(matrix(nrow=nrow(ARG_counts), ncol=2))
for (i in 1:nrow(ARG_counts)) {
  p_vals[i, 1] = t.test(test[, i] ~ Subject_group, test)$p.val
  p_vals[i, 2] = wilcox.test(test[, i] ~ Subject_group, test)$p.val
}
colnames(p_vals) = c("t_test", "mw_test")
rownames(p_vals) = rownames(ARG_counts)
# result
length(which(p_vals[,2] < 0.05))
p_vals[which(p_vals[,2] < 0.05),]

###
### Plot enrichments (stool)
###

## arg subset table
ARG_sub = read.csv("data/ARGs/ARGs_stool_adult_manual_edit.csv", row.names=1, h=T)
rownames(ARG_counts[which(p_vals[,2] < 0.05),]) == ARG_sub$family
ARG_sub = data.frame(ARG_counts[which(p_vals[,2] < 0.05), intersect(intersect(grep("St", rag_hv_meta$Body_site_ID), 
                                                                              which(is.na(rag_hv_meta$HSCT_post_month) == TRUE)), 
                                                                    grep("adult", rag_hv_meta$Subject_age))], 
                     ARG_sub[,2:4])
ARG_sub = melt(ARG_sub)

# add metadata
ARG_sub$Subject_group = rep("RAG", nrow(ARG_sub))
ARG_sub$Subject_group[grep(paste(rag_hv_meta[grep("hv_stool", rag_hv_meta$Subject_group),]$Sample,
                              collapse = "|"),
                        ARG_sub$variable)] = "HV"

# order ARGs for plot
levels(factor(ARG_sub$gene))
ARG_sub$gene = factor(ARG_sub$gene,
                      levels = c("APH(2)-Ila", 
                                 "ErmB", "ErmF", "ErmG", "mefA",
                                 "sul1",
                                 "tetX",
                                 "mdtM", 
                                 "nfsA"))

# set factor for groups
ARG_sub$Subject_group = factor(ARG_sub$Subject_group,
                            levels = c("RAG", "HV"))

# plot stool ARGs (Figure 4-C)
ggplot(ARG_sub,
       aes(x=gene, y=value, col=Subject_group, fill=Subject_group)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0, outlier.color = "white",
               alpha=0.4) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0),
              size=2.5, #col = "black", 
              pch=21,
              alpha=0.7) +
  theme_bw() +
  ylim(0, 4.1) +
  ylab("ARG Family log(RPKM)") +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  theme(axis.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.text = element_text(size=14),
        axis.text.x = element_text(size=14, angle=45, hjust=1,
                                   color = c("#FF3333",
                                             "violet", "violet", "violet", "violet",
                                             "darkblue",
                                             "#660000",
                                             "purple", "darkgray")),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.ticks.x = element_line(size = 2.5, lineend = "round",
                                    color = c("#FF3333",
                                              "violet", "violet", "violet", "violet",
                                              "darkblue",
                                              "#660000",
                                              "purple", "darkgray")))

# legend
dev.off()
plot(1:50, col="white")
legend(1,49, cex = 1, legend = c("Aminoglycoside", "Fluoroquinolone", "Fosfomycin",
                                 "Macrolide", "Penam/Beta-lactam", "Phenicol", "Sulfonamide", "Tetracycline", 
                                 "Multi-drug", "Other"),
       c("#FF3333", "darkgreen", "skyblue", "violet", "orange", "#666600", "darkblue", "#660000", "purple", "darkgray")) 


###
### HV-ABX skin microbiome response (Jo et al. 2022)
###

## shortbred counts table
hv_ABX_ARG_RPKM = read.table("data/ARGs/shortbred_ARG_hv_ABX_control.txt", row.names=1, h=T)

## set metadata
hv_ABX_ARG_meta = read.table("data/metadata/hv_ABX_control_metadata.txt", h=T, row.names=1)

## subset ARGs with at least 30 AA seqs
hv_ABX_ARG_RPKM = hv_ABX_ARG_RPKM[which(hv_ABX_ARG_RPKM$TotMarkerLength >= 30),]
rownames(hv_ABX_ARG_RPKM) = hv_ABX_ARG_RPKM$ARG_family
hv_ABX_ARG_RPKM = hv_ABX_ARG_RPKM[,1:27]
hv_ABX_ARG_RPKM = hv_ABX_ARG_RPKM[which(apply(hv_ABX_ARG_RPKM, 1, sum) > 0),]
head(hv_ABX_ARG_RPKM)
dim(hv_ABX_ARG_RPKM)

# check order
hv_ABX_ARG_meta$Sample == substr(colnames(hv_ABX_ARG_RPKM), start = 1, stop = nchar(colnames(hv_ABX_ARG_RPKM))-6)

# total count
dim(hv_ABX_ARG_RPKM)

## log transform counts
hv_ABX_ARG_RPKM = log10(hv_ABX_ARG_RPKM + 1)

## LME model
test = data.frame(t(hv_ABX_ARG_RPKM), hv_ABX_ARG_meta)
p_vals = rep(1, nrow(hv_ABX_ARG_RPKM))
for (i in 1:nrow(hv_ABX_ARG_RPKM)) {
  p_vals[i] = summary(lmer(test[,i] ~ Time + (1|Site) + (1|Subject),
                           test))$coefficients[3,5]
  }
names(p_vals) = rownames(hv_ABX_ARG_RPKM)
# result
length(which(p_vals < 0.05))
p_vals[which(p_vals < 0.05)]

## prep plot for significant ARGs
ARG_sig = data.frame(t(hv_ABX_ARG_RPKM[grep("ARO_3000024|ARO_3003741|ARO_3002639|ARO_3000778|ARO_3003950", rownames(hv_ABX_ARG_RPKM)),]), 
                     hv_ABX_ARG_meta[c(3,6)])
colnames(ARG_sig)[1:5] = c("patA", "mphE", "APH3", "adeG", "msbA")
ARG_sig = melt(ARG_sig, id.vars = c("Time", "Site"))
# add shapes for site ID
ARG_sig$Site_ID = rep(22, nrow(ARG_sig))
ARG_sig$Site_ID[grep("Vf", ARG_sig$Site)] = 25
ARG_sig$Site_ID[grep("Ra", ARG_sig$Site)] = 24

## plot
ggplot(ARG_sig,
       aes(x=factor(variable, levels = c("adeG", "APH3", "mphE", "msbA", "patA")),
           y=value, color=Time, fill=Time)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0, outlier.color = "white",
               alpha=0.4) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0),
              size=2.5, #col = "black", 
              pch = ARG_sig$Site_ID,
              alpha=0.6) +
  theme_bw() +
  ylab("ARG Family log(RPKM)") +
  xlab("ARG") +
  scale_fill_manual(values = c("#330066", "purple", "violet"),
                    labels = c("Time 0", "2 weeks", "4 weeks")) +
  scale_color_manual(values = c("#330066", "purple", "violet"),
                     labels = c("Time 0", "2 weeks", "4 weeks")) +
  scale_x_discrete(labels = c("adeG", "APH(3)-Ib", "mphE", "msbA", "patA")) +
  theme(axis.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.text = element_text(size=14),
        axis.text.x = element_text(size=14, angle=45, hjust=1,
                                   color = c("purple", "#FF3333", "violet", "darkgray", "darkgreen")),
        axis.ticks.x = element_line(size = 2.5, lineend = "round",
                                    color = c("purple", "#FF3333", "violet", "darkgray", "darkgreen")),
        legend.text = element_text(size=12),
        legend.title = element_blank())


##################################
### MAGs: gtdb-tk AA alignment ###
##################################

# load alignment
tree_bac <- read.tree("data/MAGs_read_mapping/gtbtk_bac_aln.newick")
tree_bac$tip.label

# load taxa table
tree_taxa = read.csv("data/MAGs_read_mapping/rag_mags_taxa.csv", h=T, row.names=1)
head(tree_taxa)

# set tip labels
tree_bac$tip.label = paste(tree_bac$tip.label, c("a"), sep="_")
tree_taxa$user_genome = paste(tree_taxa$user_genome, c("a"), sep="_")

# re-order tips for graphics
tree_taxa1 = tree_taxa
for(i in 1:dim(tree_taxa)[1]) {
  tree_taxa1[i,] = tree_taxa[grep(tree_bac$tip.label[i], tree_taxa$user_genome, fixed = TRUE),]
}
tree_taxa1$user_genome == tree_bac$tip.label

# add shape
tree_taxa1$tip_shape = rep(24, length(tree_bac$tip.label))
tree_taxa1$tip_shape[grep("rag_skin", tree_bac$tip.label)] = 21
tree_taxa1$tip_shape[grep("rag_stool", tree_bac$tip.label)] = 1

# add colors
tree_taxa1$tip_fill = rep("gray", length(tree_bac$tip.label))
tree_taxa1$tip_fill[grep("Actin", tree_taxa1$phylum)] = "#663300"
tree_taxa1$tip_fill[grep("Firm", tree_taxa1$phylum)] = "darkblue"
tree_taxa1$tip_fill[grep("Firmicutes_A", tree_taxa1$phylum)] = "blue"
tree_taxa1$tip_fill[grep("Firmicutes_C", tree_taxa1$phylum)] = "#3399FF"
tree_taxa1$tip_fill[grep("Proteo", tree_taxa1$phylum)] = "#006600"
tree_taxa1$tip_fill[grep("Bacteroid", tree_taxa1$phylum)] = "red"
tree_taxa1$tip_fill[grep("Verru", tree_taxa1$phylum)] = "purple"
tree_taxa1$tip_fill[grep("Cyano", tree_taxa1$phylum)] = "gold"
tree_taxa1$tip_fill[grep("Campy", tree_taxa1$phylum)] = "pink"
tree_taxa1$tip_fill[grep("Desu", tree_taxa1$phylum)] = "magenta"

# add colors for outside fill
tree_taxa1$tip_col = rep("black", length(tree_bac$tip.label))
tree_taxa1$tip_col[grep("stool", tree_bac$tip.label)] = tree_taxa1$tip_fill[grep("stool", tree_taxa1$user_genome)]

# plot (Figure S2)
plot(tree_bac,
     type = "unrooted",
     show.tip.label = FALSE,
     no.margin = TRUE,
     edge.width = 1)

# add tips
tiplabels(bg = tree_taxa1$tip_fill, 
          col = tree_taxa1$tip_col,
          cex = 1.4,
          pch = tree_taxa1$tip_shape)

# add scale bar
add.scale.bar(ask = TRUE, cex = 1)

# plot legend
dev.off()
ggplot(data=data.frame(A = c(1:2),
                       B = c("A", "B")), 
       aes(x=A, y=A, shape=B)) +
  geom_point(fill="black", size = 4, alpha=0.8) +
  theme_bw() +
  scale_shape_manual(values = c(21, 1),
                     labels = c("RAG-Skin/Nares (n=111)", "RAG-Stool (n=206)")) +
  theme(legend.text = element_text(size=12),
        legend.title = element_blank())

# legend 2
dev.off()
plot(1:50, col="white")
legend(1,49, cex = 1, legend = c("Actinobacteriota", "Bacteroidota", "Campylobacterota", "Cyanobacteria", "Desulfobacterota",
                                 "Firmicutes", "Firmicutes_A", "Firmicutes_C", "Proteobacteria", "Verrumicrobiota"),
       c("#663300", "red", "pink", "gold", "magenta", "darkblue", "blue", "#3399FF", "#006600", "purple")) 


###############################
### HPV - L1 gene alignment ###
###############################

# load alignment and taxa tables
tree_hpv <- read.tree("data/HPVs/hpv_aln.newick")
tree_tips <- read.csv("data/HPVs/hpv_tree_labels.csv")

# check order
tree_tips$tip_label == tree_hpv$tip.label

# add colors
tree_tips$color[grep("alpha", tree_tips$color)] = "#999933"
tree_tips$color[grep("beta", tree_tips$color)] = "#CC9900"
tree_tips$color[grep("gamma", tree_tips$color)] = "#666666"
tree_tips$color[grep("mu", tree_tips$color)] = "black"
tree_tips$color[grep("nu", tree_tips$color)] = "white"

# plot (Figure 5-B)
plot(tree_hpv,
     type = "fan",
     show.tip.label = FALSE,
     no.margin = TRUE,
     rotate.tree = -90,
     edge.width = 1)

# add tips
tiplabels(bg = tree_tips$color, 
          col = "black", 
          cex = tree_tips$size-0.2,
          offset = 0.02,
          pch=tree_tips$shape)

# add scale bar
add.scale.bar(ask = TRUE, cex = 1.1)

# legend
ggplot(data=data.frame(test_col = c("A","B","C","D","E","F", "G"),
                       test_shape = c(1:7)), 
       aes(x=test_col, y=test_col, fill=test_col)) +
  geom_point(size=3.5, pch=c(21)) +
  theme_bw() +
  scale_fill_manual(values = c("#999933", "#CC9900", "#666666", "black", "white", "red", "red"),
                    labels = c("Alphapapillomavirus", "Betapapillomavirus",
                               "Gammapapillomavirus", "Mupapillomavirus", "Nupapillomavirus",
                               "RAG", "RAG-novel")) +
  guides(fill = guide_legend(override.aes=list(shape = c(21,21,21,21,21,21,23)))) +
  theme(legend.text = element_text(size=12))


################
### HPViewer ###
################

## load counts -- hpv
hpv_counts = read.csv("data/HPVs/hpviewer_counts.csv", h=T, row.names=1)

## clean data
hpv_list = hpv_counts$group_genome
head(hpv_counts)
dim(hpv_counts)
hpv_counts = data.frame(hpv_counts[,c(1,6)],
                        log10(as.matrix(hpv_counts[,7:132]+1)))
head(hpv_counts)

## metadata -- hpv
hpv_meta = rag_hv_meta
hpv_meta = hpv_meta[-c(grep("hv_stool|hvskin_abx", hpv_meta$Subject_group)),]
hpv_meta = hpv_meta[-c(grep("MET", hpv_meta$Sample)),]

## subset tables by duplicate site
hpv_meta$Sample == colnames(hpv_counts)[3:128]
hpv_counts = hpv_counts[, c(intersect(grep("rag_skin|hv", hpv_meta$Subject_group),
                                      grep("Ac|Ra|Vf", hpv_meta$Body_site_ID))+2)]
hpv_meta = hpv_meta[intersect(grep("rag_skin|hv", hpv_meta$Subject_group),
                              grep("Ac|Ra|Vf", hpv_meta$Body_site_ID)),]
hpv_meta$Sample == colnames(hpv_counts)

## add site
hpv_meta$meta = paste(hpv_meta$Subject_ID, hpv_meta$Body_site_ID, sep = "_")

## compute averages
hpv_counts_new = as.data.frame(matrix(nrow = nrow(hpv_counts), 
                                      ncol = length(tapply(as.numeric(hpv_counts[1,]), 
                                                           substr(hpv_meta$meta, start=1, stop=nchar(hpv_meta$meta)-1), mean))))
for (i in 1:nrow(hpv_counts)) {
  hpv_counts_new[i,] = tapply(as.numeric(hpv_counts[i,]), 
                              substr(hpv_meta$meta, start=1, stop=nchar(hpv_meta$meta)-1), mean)
}
rownames(hpv_counts_new) = hpv_list 
colnames(hpv_counts_new) = names(tapply(as.numeric(hpv_counts[1,]), 
                                        substr(hpv_meta$meta, start=1, stop=nchar(hpv_meta$meta)-1), mean))
dim(hpv_counts_new)[2]
max(hpv_counts_new)

###
### Abundances by HPV group
###

## set data
hpv_group_avg = as.data.frame(t(rbind(apply(hpv_counts_new[grep("alpha", rownames(hpv_counts_new)),], 2, mean),
                                      apply(hpv_counts_new[grep("beta", rownames(hpv_counts_new)),], 2, mean),
                                      apply(hpv_counts_new[grep("gamma", rownames(hpv_counts_new)),], 2, mean))))
colnames(hpv_group_avg) = c("alpha", "beta", "gamma")
hpv_group_avg$Subject_group = rep("HV", dim(hpv_group_avg)[1])
hpv_group_avg$Subject_group[grep("RAG", rownames(hpv_group_avg))] = "RAG"
hpv_group_avg = melt(data.frame(rownames(hpv_group_avg),
                                hpv_group_avg))

# set factor for body site (adult)
hpv_group_avg$site_pch = rep(21, dim(hpv_group_avg)[1])
hpv_group_avg$site_pch[grep("V", hpv_group_avg$rownames.hpv_group_avg.)] = 25
hpv_group_avg$site_pch[grep("A", hpv_group_avg$rownames.hpv_group_avg.)] = 22
hpv_group_avg$site_pch[grep("R", hpv_group_avg$rownames.hpv_group_avg.)] = 24

# set factor for age
hpv_group_avg$age = rep("adult", nrow(hpv_group_avg))
hpv_group_avg$age[grep(paste(rag_hv_meta$Subject_ID[grep("child", rag_hv_meta$Subject_age)], collapse = "|"),
                       hpv_group_avg$rownames.hpv_group_avg.)] = "child"

# set factor for group
hpv_group_avg$Subject_group = factor(hpv_group_avg$Subject_group,
                                  levels = c("RAG", "HV"))

# plot (Figure 5-A)
ggplot(#hpv_group_avg[grep("child", hpv_group_avg$age),],
       hpv_group_avg[grep("adult", hpv_group_avg$age),],
       aes(x=variable, y=value, fill=Subject_group, col=Subject_group)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0, outlier.color = "white",
               alpha=0.2) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.4), size=2.5, 
              #pch = hpv_group_avg$site_pch[grep("child", hpv_group_avg$age)], 
              pch = hpv_group_avg$site_pch[grep("adult", hpv_group_avg$age)], 
              alpha=0.6) +
  scale_x_discrete(labels = c("Alpha-", "Beta-", "Gamma-", "Mu-")) +
  theme_bw() +
  ylab(expression(log[10](RPKM))) +
  ylim(0,2) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        legend.position = "none")


## stats

# RAG-child enrichments
summary(lmer(value ~ Subject_group + (1|site_pch) , hpv_group_avg[intersect(grep("alpha", hpv_group_avg$variable), grep("child", hpv_group_avg$age)),]))
summary(lmer(value ~ Subject_group + (1|site_pch) , hpv_group_avg[intersect(grep("beta", hpv_group_avg$variable), grep("child", hpv_group_avg$age)),]))
summary(lmer(value ~ Subject_group + (1|site_pch) , hpv_group_avg[intersect(grep("gamma", hpv_group_avg$variable), grep("child", hpv_group_avg$age)),]))

# RAG-adult enrichments
summary(lmer(value ~ Subject_group + (1|site_pch) , hpv_group_avg[intersect(grep("alpha", hpv_group_avg$variable), grep("adult", hpv_group_avg$age)),]))
summary(lmer(value ~ Subject_group + (1|site_pch) , hpv_group_avg[intersect(grep("beta", hpv_group_avg$variable), grep("adult", hpv_group_avg$age)),]))
summary(lmer(value ~ Subject_group + (1|site_pch) , hpv_group_avg[intersect(grep("gamma", hpv_group_avg$variable), grep("adult", hpv_group_avg$age)),]))



