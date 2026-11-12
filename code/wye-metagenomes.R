############################
### Wye Pond Metagenomes ###
############################

# rblauste@umd.edu #
# last edited: 11/11/2026 #

## load packages

library(ggplot2)
library(gtools)
library(reshape2)
library(vegan)
library(gplots)
library(stringr)
library(lmerTest)


###############################################
### Microbial taxa in metagenomes (Bracken) ###
###############################################

###
### load and prep data
###

## taxa table
bracken = read.csv("R_import_data/taxa/bracken.csv", row.names=1, h=T)
bracken = bracken[,order(colnames(bracken))]
bracken$X.Classification = str_replace_all(bracken$X.Classification, "[|]", " ")

## metadata table
metadata = read.csv("R_import_data/metadata/sample-metadata.csv", row.names=1, h=T)
metadata = metadata[order(metadata$sample.id),]
metadata$sample.id == colnames(bracken)[1:24] # check that order aligns

## set taxa tables by level
#bracken_species = bracken[grep("p__", bracken$X.Classification),][-c(grep("c__|o__|f__|g__|s__", bracken[grep("p__", bracken$X.Classification), 25])),] 
bracken_x = bracken[grep("c__", bracken$X.Classification),][-c(grep("o__|f__|g__|s__", bracken[grep("c__", bracken$X.Classification), 25])),] 
bracken_x = bracken[grep("o__", bracken$X.Classification),][-c(grep("f__|g__|s__", bracken[grep("o__", bracken$X.Classification), 25])),] 
bracken_x = bracken[grep("f__", bracken$X.Classification),][-c(grep("g__|s__", bracken[grep("f__", bracken$X.Classification), 25])),] 
bracken_genus = bracken[grep("g__", bracken$X.Classification),][-c(grep("s__", bracken[grep("g__", bracken$X.Classification), 25])),] 
bracken_species = bracken[grep("s__", bracken$X.Classification),]

## subset SPECIES for plotting
bracken_species = bracken[grep("s__", bracken$X.Classification),]
rownames(bracken_species) = bracken_species$X.Classification
bracken_species = bracken_species[,1:24]

## normalize relative abundances
for(i in 1:24) {
  bracken_species[,i] = bracken_species[,i]/sum(bracken_species[,i])
}
bracken_species$avg = apply(bracken_species, 1, mean)
bracken_species$avg_water = apply(bracken_species[,grep("W", colnames(bracken_species))], 1, mean)
bracken_species$avg_sed = apply(bracken_species[,grep("S", colnames(bracken_species))], 1, mean)
bracken_species = bracken_species[order(bracken_species$avg_water, decreasing = T),]
head(bracken_species, 10)

###
### Alpha Diversity
###

## set input
colnames(bracken_species)
taxa_table = bracken_species[,1:24]
colnames(taxa_table)
dim(taxa_table)

## compute shannon and simpson diversity metrics
diversity_vec = matrix(nrow = dim(taxa_table)[2], ncol = 2)
diversity_vec = as.data.frame(diversity_vec)
for (a in 1:dim(taxa_table)[2]) {
  diversity_vec[a,1] = diversity(taxa_table[,a], index = "shannon")
  diversity_vec[a,2] = diversity(taxa_table[,a], index = "simpson")
}
colnames(diversity_vec) = c("Shannon", "Simpson")

## add alpha metric to metadata
metadata$Shannon = diversity_vec$Shannon
metadata$Simpson = diversity_vec$Simpson

## two-way ANOVA for effect of sample depth and time on alpha diversity
summary(aov(Shannon ~ sample.location.depth*time.of.day, metadata[grep("Wye_water", metadata$sample.type.pond),]))
TukeyHSD(aov(Shannon ~ sample.location.depth*time.of.day, metadata[grep("Wye_water", metadata$sample.type.pond),]))

## ANOVA for effect of sample depth on alpha diversity
summary(aov(Shannon ~ sample.location.depth, metadata[grep("Wye_water", metadata$sample.type.pond),]))
TukeyHSD(aov(Shannon ~ sample.location.depth, metadata[grep("Wye_water", metadata$sample.type.pond),]))

## ANOVA for effect of time on alpha diversity
summary(aov(Shannon ~ time.of.day, metadata[grep("Wye_water", metadata$sample.type.pond),]))
TukeyHSD(aov(Shannon ~ time.of.day, metadata[grep("Wye_water", metadata$sample.type.pond),]))

## plot alpha-diversity by pond sampling site (FIG. 1A)
ggplot(#metadata[grep("Wye_water", metadata$sample.type.pond),], 
  metadata, 
  aes(x = sample.location.depth, y = Shannon)) +
  geom_boxplot(outlier.shape = NA, size=0.8, alpha = 0.6, fill="lightgray") +
  geom_jitter(size = 3, width = 0.15, alpha = 0.6, pch=21, fill = "lightgray") +
  theme_bw() +
  ylab("Shannon Index") +
  scale_x_discrete(labels = c("Bank: surface", "Interior: surface", 
                              "Interior: 1-m depth", "Interior: 2-m depth", "Sediment")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, angle=35, hjust=1),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16))

## plot plot alpha-diversity vs. EC (FIG. 3B)
ggplot(metadata[grep("Wye_water", metadata$sample.type.pond),], 
       aes(x = Shannon, y = ECCFU)) +
  geom_point(size = 4, pch = 21, alpha = 0.8, fill = "gray") +
  theme_bw() +
  geom_smooth(method = "lm", se = F
  ) +
  xlab("Shannon index") +
  #ylab("E. coli CFU 100 ml-1") +
  ylab(expression(italic('E. coli')~plain('CFU 100 ml'^-1))) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

## correlation
cor.test(metadata[grep("Wye_water", metadata$sample.type.pond),]$Shannon,
         metadata[grep("Wye_water", metadata$sample.type.pond),]$INV,
         method = "spearman")

## plot EC levels by pond depth (FIG. 3A)
ggplot(metadata[grep("Wye_water", metadata$sample.type.pond),], 
       aes(x = factor(sample.depth.type, 
                      levels = c("surface", "depth")), y = ECCFU)) +
  geom_boxplot(outlier.shape = NA, size=0.8, alpha = 0.6, fill="lightgray") +
  geom_jitter(size = 2, width = 0.2, alpha = 0.6, pch=21, fill = "lightgray") +
  theme_bw() +
  geom_smooth(method = "lm", se = F
  ) +
  xlab("Sampling Location") +
  scale_x_discrete(labels = c("Surface", "Water Column")) +
  #ylab("E. coli CFU 100 ml-1") +
  ylab(expression(italic('E. coli')~plain('CFU 100 ml'^-1))) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank())

## stats for EC
summary(aov(ECCFU ~ sample.depth.type,
            metadata[grep("Wye_water", metadata$sample.type.pond),]))
kruskal.test(ECCFU ~ sample.depth.type,
             metadata[grep("Wye_water", metadata$sample.type.pond),])

## alpha-diversity associations with all metadata (data in TABLE 1)
alpha_cor = as.data.frame(matrix(nrow = length(colnames(metadata)[c(17:45)]), 
                                 ncol=3))
for (i in c(17:46)) {
  alpha_cor[i,1] = colnames(metadata)[i]
  alpha_cor[i,2] = cor.test(metadata$Shannon, metadata[,i], method = "spearman")$estimate
  alpha_cor[i,3] = cor.test(metadata$Shannon, metadata[,i], method = "spearman")$p.val
}
alpha_cor


###
### Beta Diversity
###

## beta-diversity measure
beta <- vegdist(t(taxa_table[,grep("Wye_water", metadata$sample.type.pond)]), 'bray', binary = T)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# add metadata
ord$sample.depth = as.character(metadata$sample.depth[grep("Wye_water", metadata$sample.type.pond)])
ord$time.of.day = as.character(metadata$time.of.day[grep("Wye_water", metadata$sample.type.pond)])
ord$sample.location = as.character(metadata$sample.location[grep("Wye_water", metadata$sample.type.pond)])
ord$sample.location.depth = as.character(metadata$sample.location.depth[grep("Wye_water", metadata$sample.type.pond)])
ord$sample.depth.broad = as.character(metadata$sample.location.depth[grep("Wye_water", metadata$sample.type.pond)])
ord$sample.depth.broad2 = "surface"
ord$sample.depth.broad2[grep("D2|D3", ord$sample.depth.broad)] = "below"
ord$ECCFU = metadata$ECCFU[grep("Wye_water", metadata$sample.type.pond)]

## plot PCoA (Fig. 1B)
ggplot(data = ord, aes(x = -pcoa1, y = -pcoa2, 
                       fill=sample.location.depth,
                       shape=time.of.day)) +
  geom_hline(yintercept = 0, lty = 2, alpha=0.7) +
  geom_vline(xintercept = 0, lty = 2, alpha=0.7) +
  geom_point(alpha = 0.7, #pch=21, 
             size=4.5) +
  theme_bw() +
  xlab('PCo1 (60.9%)') +
  ylab('PCo2 (8.8%)') +
  scale_shape_manual(values=c(22, 24, 23),
                     name = "Time of Day",
                     labels = c("9:00", "12:00", "15:00")) +
  scale_fill_manual(values=c(c("skyblue", "darkblue", "pink", "red")),
                    labels = c("Bank: surface",
                               "Interior: surface",
                               "Interior: 1-m depth",
                               "Interior: 2-m depth"),
                    name="Sample Location") +
  guides(fill = guide_legend(override.aes=list(shape = c(21)))) +
  theme(axis.text = element_text(size=14, color = "black"),
        axis.title = element_text(size=14, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color = "black")) 

## stats
adonis2(beta ~ sample.location.depth, data = ord, permutations = 999)
adonis2(beta ~ sample.depth.broad2, data = ord, permutations = 999)
adonis2(beta ~ time.of.day, data = ord, permutations = 999)
adonis2(beta ~ sample.location, data = ord, permutations = 999)
adonis2(beta ~ ECCFU, data = ord, permutations = 999)

## beta-diversity (all samples) correlations with metadata (data in TABLE 1)
beta <- vegdist(t(taxa_table[,grep("Wye_water", metadata$sample.type.pond)]), 'bray', binary = T)
# stats
adonis2(beta ~ PO4, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ NH3, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ NO3, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ ECCFU, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ ECMPN, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ TC, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ TIC, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ TN, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ TOC, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ CDOM, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ INV, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ PHYC, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)

## beta-diversity correlation with select variable (without W18, due to missing data for specific sample)
taxa_table[1,grep("Wye_water", metadata$sample.type.pond)][,-c(2)]
beta <- vegdist(t(taxa_table[,grep("Wye_water", metadata$sample.type.pond)][,-c(2)]), 'bray', binary = T)
# stats
adonis2(beta ~ Temperature, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ DO.mg.L, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ SPC.uS.cm, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ pH, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ UVNOx.mg.L, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ NTU, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ BGA.PC.RFU, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ BGA.PC.ug.L, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ Chl.RFU, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ Chl.ug.L, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ fDOM.RFU, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ fDOM.ppb, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)

## beta-diversity (NN)
taxa_table[1,grep("Wye_water", metadata$sample.type.pond)][,-c(which(is.na(metadata[grep("Wye_water", metadata$sample.type.pond),]$NN) == TRUE))]
beta <- vegdist(t(taxa_table[,grep("Wye_water", metadata$sample.type.pond)][,-c(which(is.na(metadata[grep("Wye_water", metadata$sample.type.pond),]$NN) == TRUE))]), 'bray', binary = T)

# stats
adonis2(beta ~ NN, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(which(is.na(metadata[grep("Wye_water", metadata$sample.type.pond),]$NN) == TRUE)),], permutations = 999)


###
### Taxa Abundances
###

## subset top taxa for plotting
taxa_plot = bracken_species
taxa_plot_select = rbind(apply(taxa_plot, 2, sum),
                         apply(taxa_plot[grep("_Actinobacteria", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Bacteroidetes", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Cyanobacteria", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Firmicutes", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Proteobacteria", rownames(taxa_plot)),], 2, sum),
                         taxa_plot[c(1:6),])
rownames(taxa_plot_select) = c("Other", 
                               "p__Actinobacteria",
                               "p__Bacteroidetes",
                               "p__Cyanobacteria",                               
                               "p__Firmicutes",                               
                               "p__Proteobacteria",
                               rownames(taxa_plot_select[c(7:12),]))
head(taxa_plot_select)

# clean taxa (reset abundances against broader taxonomic groups)
rownames(taxa_plot_select)
taxa_plot_select[grep("^p__Actinobacteria$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^p__Actinobacteria$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Mycolicibacterium", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^p__Cyanobacteria$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^p__Cyanobacteria$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Microcystis|Synechococcus|Sphaerospermopsis", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^p__Proteobacteria$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^p__Proteobacteria$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Alphapr|Betapr", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^Other$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^Other$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("_", rownames(taxa_plot_select)),], 2, sum)

# check sums
apply(taxa_plot_select,2,sum)

## prep plot variable labels (shortened phyla and species names)
rownames(taxa_plot_select)
rownames(taxa_plot_select) = c("Other", 
                               "p__Actinobacteria",
                               "p__Bacteroidetes",
                               "p__Cyanobacteria",                               
                               "p__Firmicutes",                               
                               "p__Proteobacteria",
                               "  Microcystis aeruginosa", "  Synechococcus sp. CBW1004", "  Limnohabitans sp. 103DPR2", 
                               "  Magnetospirillum gryphiswaldense", "  Microcystis sp. MC19", "  Polynucleobacter sp. MWH-UH24A")

## melt data
colnames(taxa_plot_select)
taxa_melt = melt(data.frame(taxa_plot_select[, 1:24],
                            rownames(taxa_plot_select)))
colnames(taxa_melt)[1] = "Taxon"

# taxa factor re-order
taxa_melt$Taxon = factor(taxa_melt$Taxon,
                         levels = c("Other", 
                                    "p__Actinobacteria", "  Mycolicibacterium", 
                                    "p__Bacteroidetes",
                                    "p__Cyanobacteria",  "  Microcystis aeruginosa", "  Synechococcus sp. CBW1004", "  Microcystis sp. MC19",                            
                                    "p__Firmicutes",                               
                                    "p__Proteobacteria", "  Limnohabitans sp. 103DPR2", "  Magnetospirillum gryphiswaldense", 
                                    "  Polynucleobacter sp. MWH-UH24A"))

# sample factor re-order
taxa_melt$variable = factor(taxa_melt$variable,
                            levels = mixedsort(levels(factor(taxa_melt$variable))))

## view barplot
head(taxa_melt)
ggplot(taxa_melt,
       aes(x=variable, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  scale_fill_manual(values=c(
    "lightgray",
    "blue",
    "yellow",
    "#003300", "#006600", "#33CC33","#99FF99", 
    "violet",
    "#330000", "#993300", "#FF6600", "#FFCC99"
  )) +
  theme(#legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    axis.text = element_text(size=14),
    axis.text.x = element_text(size=14, angle=45, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14, color = "black"))

## add metadata for plot

# set blank variables
taxa_melt$sample.depth = "X"
taxa_melt$sample.type.pond = "X"
taxa_melt$sample.location = "X"
taxa_melt$sample.water.location = "X"
taxa_melt$sample.location.depth = "X"
taxa_melt$time.of.day = "X"
taxa_melt$label.plot = "X"

# add to table
for (i in 1:nrow(taxa_melt)) {
  taxa_melt$sample.depth[grep(metadata$sample.id[i], taxa_melt$variable)] = metadata$sample.depth[i]
  taxa_melt$sample.type.pond[grep(metadata$sample.id[i], taxa_melt$variable)] = metadata$sample.type.pond[i]
  taxa_melt$sample.location[grep(metadata$sample.id[i], taxa_melt$variable)] = metadata$sample.location[i]
  taxa_melt$sample.water.location[grep(metadata$sample.id[i], taxa_melt$variable)] = metadata$sample.water.location[i]
  taxa_melt$sample.location.depth[grep(metadata$sample.id[i], taxa_melt$variable)] = metadata$sample.location.depth[i]
  taxa_melt$time.of.day[grep(metadata$sample.id[i], taxa_melt$variable)] = metadata$time.of.day[i]
  taxa_melt$label.plot[grep(metadata$sample.id[i], taxa_melt$variable)] = metadata$label.plot[i]
}

## set labels for facet
depth_names = list(
  "Bank_D1"="Bank: surface",
  "Interior_D1"="Interior: surface",
  "Interior_D2"="Interior: 1-m depth",
  "Interior_D3"="Interior: 2-m depth",
  "S"="Sediment"
)
depth_labeller <- function(variable,value){
  return(depth_names[value])
}

## set factor for sample ids
taxa_melt$variable = factor(taxa_melt$variable,
                            levels = c("W6", "W26", "W46", "W7", "W27", "W47",
                                       "W17", "W37", "W57","W8","W28","W48",
                                       "W18","W38","W58", "W9","W29","W49",
                                       "W19","W39","W59","S4","S5","S11"))
taxa_melt$label.plot = factor(taxa_melt$label.plot,
                              levels = levels(factor(taxa_melt$label.plot))[mixedorder(levels(factor(taxa_melt$label.plot)))])

# reset labels
taxa_melt$label.plot2 = substr(taxa_melt$label.plot, start = 3, stop = nchar(as.character(taxa_melt$label.plot)))
taxa_melt$label.plot2 = paste(sapply(strsplit(taxa_melt$label.plot2, "_"), function (x) x[1]),
                              substr(sapply(strsplit(taxa_melt$label.plot2, "_"), function (x) x[2]),
                                     start = 4, stop = nchar(sapply(strsplit(taxa_melt$label.plot2, "_"), function (x) x[2]))),
                              sep = "_")
taxa_melt$label.plot2 = factor(taxa_melt$label.plot2,
                               levels = levels(factor(taxa_melt$label.plot2))[mixedorder(levels(factor(taxa_melt$label.plot2)))])
levels(factor(taxa_melt$label.plot2))[mixedorder(levels(factor(taxa_melt$label.plot2)))]

## barplot (FIG. 1C)
head(taxa_melt)
ggplot(taxa_melt,
       aes(x=label.plot2, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  scale_fill_manual(values=c(
    "lightgray",
    "blue",
    "yellow",
    "#003300", "#006600", "#33CC33","#99FF99", 
    "violet",
    "#330000", "#993300", "#FF6600", "#FFCC99")) +
  facet_grid(~sample.location.depth, scales="free", space="free") +
  theme(#legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    strip.text = element_text(size=10),
    axis.text = element_text(size=14),
    axis.text.x = element_text(size=12, angle=45, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14, color = "black"))


###
### Stats for enriched species
###

## merge taxa and metadata for stats

# call table
taxa_stat = bracken_species[which(apply(bracken_species[,grep("Wye_water", metadata$sample.type.pond)], 1, mean) > 0.001), 1:24]

# order based on abundances
head(taxa_stat)
taxa_stat = taxa_stat[order(apply(taxa_stat, 1, mean), decreasing = T),]
taxa_stat = taxa_stat[c(which(as.numeric(apply(taxa_stat[,grep("surface", metadata$sample.depth.type)], 1, mean, na.rm=TRUE) + 0.001)/
                                as.numeric(apply(taxa_stat[,grep("depth", metadata$sample.depth.type)], 1, mean, na.rm=TRUE) + 0.001) > 1),
                        which(as.numeric(apply(taxa_stat[,grep("surface", metadata$sample.depth.type)], 1, mean, na.rm=TRUE) + 0.001)/
                                as.numeric(apply(taxa_stat[,grep("depth", metadata$sample.depth.type)], 1, mean, na.rm=TRUE) + 0.001) < 1)),]
taxa_stat = data.frame(t(taxa_stat), metadata)
taxa_stat = taxa_stat[4:24,] # water surface vs. depth
#taxa_stat = taxa_stat[grep("sed|depth", taxa_stat$sample.depth.type),] # depth vs. sed
#taxa_stat = taxa_stat[grep("sed|surface", taxa_stat$sample.depth.type),] # depth vs. sed

# LME stats
p.val = c(1:69)
for (i in 1:69) {
  p.val[i] = unlist(summary(lmer(taxa_stat[,i] ~ sample.depth.type + (1|time.of.day) + (1|sample.location), taxa_stat)))$coefficients10
  #p.val[i] = unlist(summary(lmer(taxa_stat[,i] ~ time.of.day + (1|sample.depth) + (1|sample.location), taxa_stat)))$coefficients10
}
p.val
which(p.adjust(p.val, method = "BH") < 0.05)
length(which(p.adjust(p.val, method = "BH") < 0.05))

# summary stats
test = matrix(nrow = 69, ncol = 4)
test = as.data.frame(test)
for (i in 1:69) {
  test[i,] = c(tapply(taxa_stat[,i], taxa_stat$sample.depth.type, mean),
               tapply(taxa_stat[,i], taxa_stat$sample.depth.type, sd))
}
rownames(test) = colnames(taxa_stat)[1:69]
colnames(test) = c("depth_mean", "surface_mean", "depth_sd", "surface_sd")
test$p.val = p.val
test$q.val = p.adjust(p.val, method = "BH")
head(test)

# export table for data in TABLE S1
#write.table(test, "Taxa_enrich2.txt", sep = "\t", quote = FALSE)


###
### Top taxa
###

# rename columns
colnames(taxa_stat)[1:69]
sapply(str_split(colnames(taxa_stat)[1:69], "__"),tail,1)
colnames(taxa_stat)[1:69] = sapply(str_split(colnames(taxa_stat)[1:69], "__"),tail,1)
colnames(taxa_stat)[1:69] = gsub("_", " ",colnames(taxa_stat)[1:69])

# plot enriched taxa by depth (FIG. 2A)
ggplot(#melt(taxa_stat[,c(which(p.adjust(p.val, method = "BH") < 0.05), 19)]),
       #melt(taxa_stat[,c(which(p.adjust(p.val, method = "BH") < 0.05), 77)]),
       melt(taxa_stat[,c(intersect(which(apply(taxa_stat[1:69], 2, mean) > 0.005),
                                   which(p.adjust(p.val, method = "BH") < 0.05)), 77)]),
       aes(y=variable, x=100*value, 
           fill=sample.depth.type, color=sample.depth.type)) +
  geom_boxplot(position = position_dodge(width=0.8), 
               alpha=0.3, outlier.size = 0, outlier.color = "white") +
  geom_jitter(position = position_dodge(width=0.8), 
              size = 0.5, alpha=0.6, pch=21) +
  xlab("Relative Abundance (%)") +
  theme_bw() +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 50, unit = "pt"))


###########################
### Functional analysis ###
###########################

###
### load data
###

## humann
humann_table = read.csv("R_import_data/annotation/humann_pathabundance_unstratified-relab-mapped_only.csv", row.names=1, h=T)
head(humann_table)

## shortbred_filt
shortbred_filt = read.csv("R_import_data/annotation/shortbred_filt.csv", row.names=1, h=T)
shortbred_filt = shortbred_filt[,c(order(colnames(shortbred_filt)[1:24]), 25:32)]

## shortbred_MAGs
shortbred_MAGs = read.csv("R_import_data/annotation/shortbred_MAGs.csv", row.names=1, h=T)

## shortbred_filt -- VFDB
shortbred_VFDB_filt = read.csv("R_import_data/annotation/shortbred_VFDB_filt.csv", row.names=1, h=T)
shortbred_VFDB_filt = shortbred_VFDB_filt[,c(order(colnames(shortbred_VFDB_filt)[1:24]), 25:28)]

## shortbred_MAGs -- VFDB
shortbred_MAGs_VFDB = read.csv("R_import_data/annotation/shortbred_MAGs_VFDB.csv", row.names=1, h=T)

## metadata -- check that sample orders align
metadata$sample.id == colnames(humann_table)[1:24]
metadata$sample.id == colnames(shortbred_filt)[1:24]
metadata$sample.id == colnames(shortbred_VFDB_filt)[1:24]


###
### Pathways (humann)
###

## order based on abundances
head(humann_table)
humann_table = humann_table[order(apply(humann_table, 1, mean), decreasing = T),]
humann_table = humann_table[c(which(as.numeric(apply(humann_table[,grep("surface", metadata$sample.depth.type)], 1, mean, na.rm=TRUE) + 0.001)/
                                      as.numeric(apply(humann_table[,grep("depth", metadata$sample.depth.type)], 1, mean, na.rm=TRUE) + 0.001) > 1),
                              which(as.numeric(apply(humann_table[,grep("surface", metadata$sample.depth.type)], 1, mean, na.rm=TRUE) + 0.001)/
                                      as.numeric(apply(humann_table[,grep("depth", metadata$sample.depth.type)], 1, mean, na.rm=TRUE) + 0.001) < 1)),]

## merge meta
colnames(humann_table) == metadata$sample.id
humann_meta = data.frame(t(humann_table), metadata)
colnames(humann_meta)
rownames(humann_meta)
humann_meta = humann_meta[4:24,] # surface water vs. water column
#humann_meta = humann_meta[grep("sed|dep", humann_meta$sample.depth.type),] # water column vs. sediments

## LME
p.val = c(1:199)
for (i in 1:199) {
  p.val[i] = unlist(summary(lmer(humann_meta[,i] ~ sample.depth.type + (1|time.of.day) + (1|sample.location), humann_meta)))$coefficients10
  #p.val[i] = unlist(summary(lmer(humann_meta[,i] ~ time.of.day + (1|sample.depth) + (1|sample.location), humann_meta)))$coefficients14
  #p.val[i] = unlist(summary(lmer(humann_meta[,i] ~ time.of.day + (1|sample.depth) + (1|sample.location), humann_meta)))$coefficients15
}
p.val
which(p.val < 0.05)
which(p.adjust(p.val, method = "BH") < 0.05)
length(which(p.adjust(p.val, method = "BH") < 0.05))
dim(humann_meta[,c(intersect(which(apply(humann_meta[1:199], 2, mean) > 0.001),
                             which(p.adjust(p.val, method = "BH") < 0.05)))])

# enrichment table
test = matrix(nrow = 199, ncol = 4)
test = as.data.frame(test)
for (i in 1:199) {
  test[i,] = c(tapply(humann_meta[,i], humann_meta$sample.depth.type, mean),
               tapply(humann_meta[,i], humann_meta$sample.depth.type, sd))
}
rownames(test) = colnames(humann_meta)[1:199]
colnames(test) = c("depth_mean", "surface_mean", "depth_sd", "surface_sd")
test$p.val = p.val
test$q.val = p.adjust(p.val, method = "BH")
head(test)

# export data for TABLE S2
#write.table(test, "Fxn_enrich.txt", sep = "\t", quote = FALSE)


## LME for 'enrich 2', removing 0's in sed
p.val = c(1:121)
for (i in 1:121) {
  p.val[i] = unlist(summary(lmer(humann_meta[,i] ~ sample.depth.type + (1|time.of.day) + (1|sample.location), humann_meta)))$coefficients10
  #p.val[i] = unlist(summary(lmer(humann_meta[,i] ~ time.of.day + (1|sample.depth) + (1|sample.location), humann_meta)))$coefficients14
  #p.val[i] = unlist(summary(lmer(humann_meta[,i] ~ time.of.day + (1|sample.depth) + (1|sample.location), humann_meta)))$coefficients15
}
p.val
which(p.val < 0.05)
which(p.adjust(p.val, method = "BH") < 0.05)
length(which(p.adjust(p.val, method = "BH") < 0.05))

# enrichment table2
dim(humann_meta)
colnames(humann_meta)
tapply(humann_meta[,1], humann_meta$sample.depth.type, mean)

test = matrix(nrow = 121, ncol = 4)
test = as.data.frame(test)
for (i in 1:121) {
  test[i,] = c(tapply(humann_meta[,i], humann_meta$sample.depth.type, mean),
               tapply(humann_meta[,i], humann_meta$sample.depth.type, sd))
}
rownames(test) = colnames(humann_meta)[1:121]
colnames(test) = c("depth_mean", "sed_mean", "depth_sd", "sed_sd")
test$p.val = p.val
test$q.val = p.adjust(p.val, method = "BH")
head(test)

# export data for TABLE S2
#write.table(test, "Fxn_enrich2.txt", sep = "\t", quote = FALSE)

## plot enriched functions by depth (FIG. 2B)
ggplot(melt(humann_meta[,c(intersect(which(apply(humann_meta[1:199], 2, mean) > 0.01),
                                     which(p.adjust(p.val, method = "BH") < 0.05)), 207)]),
       aes(y=variable, x=100*value, 
           fill=sample.depth.type, color=sample.depth.type)) +
  geom_boxplot(position = position_dodge(width=0.8), 
               alpha=0.3, outlier.size = 0, outlier.color = "white") +
  geom_jitter(position = position_dodge(width=0.8), 
              size = 1, alpha=0.6, pch=21) +
  xlab("Relative Abundance (%)") +
  theme_bw() +
  #scale_fill_manual(values = c("skyblue", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_y_discrete(limits = rev) +
  #scale_x_log10() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=14),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=14)) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt"))


###
### Beta Diversity
###

## PCoA of Wye water samples

# beta-diversity measure
beta <- vegdist(t(humann_table[,grep("Wye_water", metadata$sample.type.pond)]), 'bray', binary = T)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# add metadata
ord$sample.depth = as.character(metadata$sample.depth[grep("Wye_water", metadata$sample.type.pond)])
ord$time.of.day = as.character(metadata$time.of.day[grep("Wye_water", metadata$sample.type.pond)])
ord$sample.location = as.character(metadata$sample.location[grep("Wye_water", metadata$sample.type.pond)])
ord$sample.location.depth = as.character(metadata$sample.location.depth[grep("Wye_water", metadata$sample.type.pond)])
ord$sample.depth.broad = as.character(metadata$sample.location.depth[grep("Wye_water", metadata$sample.type.pond)])
ord$sample.depth.broad2 = "surface"
ord$sample.depth.broad2[grep("D2|D3", ord$sample.depth.broad)] = "below"
ord$ECCFU = metadata$ECCFU[grep("Wye_water", metadata$sample.type.pond)]

# plot PCoA (FIG. S2)
ggplot(data = ord, aes(x = -pcoa1, y = -pcoa2, 
                       fill=sample.location.depth,
                       shape=time.of.day)) +
  geom_hline(yintercept = 0, lty = 2, alpha=0.7) +
  geom_vline(xintercept = 0, lty = 2, alpha=0.7) +
  geom_point(alpha = 0.7, #pch=21, 
             size=4.5) +
  theme_bw() +
  xlab('PCoA1 (61.6%)') +
  ylab('PCoA2 (10.3%)') +
  scale_shape_manual(values=c(22, 24, 23),
                     name = "Time of Day",
                     labels = c("9:00", "12:00", "15:00")) +
  scale_fill_manual(values=c(c("skyblue", "darkblue", "pink", "red")),
                    labels = c("Bank: surface",
                               "Interior: surface",
                               "Interior: 1-m depth",
                               "Interior: 2-m depth"),
                    name="Sample Location") +
  guides(fill = guide_legend(override.aes=list(shape = c(21)))) +
  theme(axis.text = element_text(size=14, color = "black"),
        axis.title = element_text(size=14, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color = "black")) 

# stats
adonis2(beta ~ sample.location.depth, data = ord, permutations = 999)
adonis2(beta ~ sample.depth.broad2, data = ord, permutations = 999)
adonis2(beta ~ time.of.day, data = ord, permutations = 999)
adonis2(beta ~ sample.location, data = ord, permutations = 999)
adonis2(beta ~ ECCFU, data = ord, permutations = 999)

## beta-diversity (all samples) correlations with metadata (data in TABLE 1)
beta <- vegdist(t(humann_table[,grep("Wye_water", metadata$sample.type.pond)]), 'bray', binary = T)
# stats
adonis2(beta ~ PO4, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ NH3, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ NO3, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ ECCFU, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ ECMPN, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ TC, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ TIC, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ TN, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ TOC, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ CDOM, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ INV, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)
adonis2(beta ~ PHYC, data = metadata[grep("Wye_water", metadata$sample.type.pond),], permutations = 999)

## beta-diversity correlation with select variable (without W18, due to missing data for specific sample)
humann_table[1,grep("Wye_water", metadata$sample.type.pond)][,-c(2)]
beta <- vegdist(t(humann_table[,grep("Wye_water", metadata$sample.type.pond)][,-c(2)]), 'bray', binary = T)
# stats
adonis2(beta ~ Temperature, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ DO.mg.L, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ SPC.uS.cm, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ pH, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ UVNOx.mg.L, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ NTU, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ BGA.PC.RFU, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ BGA.PC.ug.L, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ Chl.RFU, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ Chl.ug.L, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ fDOM.RFU, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)
adonis2(beta ~ fDOM.ppb, data = metadata[grep("Wye_water", metadata$sample.type.pond),][-c(2),], permutations = 999)

## beta-diversity (NN)
taxa_table[1,grep("Wye_water", metadata$sample.type.pond)][,-c(which(is.na(metadata[grep("Wye_water", metadata$sample.type.pond),]$NN) == TRUE))]
beta <- vegdist(t(taxa_table[,grep("Wye_water", metadata$sample.type.pond)][,-c(which(is.na(metadata[grep("Wye_water", metadata$sample.type.pond),]$NN) == TRUE))]), 'bray', binary = T)


##############################################
### ARGs in metagenomes (Shortbred & CARD) ###
##############################################

## set data
ARG_heat = as.matrix(shortbred_filt[which(shortbred_filt$TotMarkerLength >= 15), 1:24])
rownames(ARG_heat) = shortbred_filt$Gene[which(shortbred_filt$TotMarkerLength >= 15)]
ARG_heat = log10(ARG_heat+1)

## merge meta
colnames(ARG_heat) == metadata$sample.id
ARG_meta = data.frame(t(ARG_heat), metadata)
colnames(ARG_meta)
rownames(ARG_meta)
ARG_meta = ARG_meta[c(4:24),-c(1,5)]
#apply(ARG_meta[,1:15], 2, sum)

## LME
summary(lmer(ARG_meta[,2] ~ sample.depth.type + (1|time.of.day) + (1|sample.location), ARG_meta))
p.val = c(1:13)
for (i in 1:13) {
  p.val[i] = unlist(summary(lmer(ARG_meta[,i] ~ sample.depth.type + (1|time.of.day) + (1|sample.location), ARG_meta)))$coefficients10
}
p.val
which(p.adjust(p.val, method = "BH") < 0.05)

## reset shortbred_filt table
shortbred_filt = read.csv("R_import_data/annotation/shortbred_filt.csv", row.names=1, h=T)
shortbred_filt = shortbred_filt[,c(order(colnames(shortbred_filt)[1:24]), 25:32)]

## set data
ARG_heat = as.matrix(shortbred_filt[which(shortbred_filt$TotMarkerLength >= 15),1:24])
rownames(ARG_heat) = shortbred_filt$Gene[which(shortbred_filt$TotMarkerLength >= 15)]
ARG_heat = log10(ARG_heat+1)
colnames(ARG_heat)

## col_order
col_order = data.frame(ID = c("W6", "W26", "W46", "W7", "W27", "W47",
                              "W17", "W37", "W57","W8","W28","W48",
                              "W18","W38","W58", "W9","W29","W49",
                              "W19","W39","W59","S4","S5","S11"),
                       ORD = (1:24))

## re-order
for (i in 1:24) {
  #col_order$ORD[i] = grep(col_order$ID[i], colnames(feature_tab))
  col_order$ORD[i] = grep(colnames(ARG_heat)[i], col_order$ID)
  #ARG_heat[,i] = ARG_heat[,grep(col_order$ID[i], colnames(ARG_heat))]
}
col_order$ID = colnames(ARG_heat)
col_order

# re-order ARG data
ARG_heat = ARG_heat[,order(colnames(ARG_heat))]
colnames(ARG_heat) == col_order$ID
ARG_heat = ARG_heat[,order(col_order$ORD)]
colnames(ARG_heat)

## set colnames
#colnames(ARG_heat) = metadata$label.plot
for (i in 1:24) {
  colnames(ARG_heat)[i] = metadata$label.plot[grep(colnames(ARG_heat)[i], metadata$sample.id)]
}
colnames(ARG_heat)

## set heat colors
min(ARG_heat)
max(ARG_heat)
quantile(ARG_heat)
mypalette <- colorRampPalette(c("white", "#00CCFF", "darkblue"))(n = 99)
mypalette <- colorRampPalette(c("white", "red", "darkred"))(n = 99)
mycolors = seq(0, 1.55, length.out = 100)

## set side colors
# location
location = rep("purple", length(metadata$sample.location))
location[grep("L5", metadata$sample.location)] = "pink"
location[grep("L11", metadata$sample.location)] = "magenta"
# time
time = rep("red", length(metadata$time.of.day))
time[grep("T1", metadata$time.of.day)] = "purple"
time[grep("T2", metadata$time.of.day)] = "pink"
time[grep("T3", metadata$time.of.day)] = "magenta"
# ARG class
ARG_class = rep("lightgray", length(shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)]))
factor(shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])
ARG_class[grep("amino", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#FFFFCC"
ARG_class[grep("disinf", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#FFFF00"
ARG_class[grep("carb", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#FFCC33"
ARG_class[grep("multidrug", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#FF6600"
ARG_class[grep("rifamycin", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#993300"
ARG_class[grep("tetracycline", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#330000"
# ARG class
levels(as.factor(shortbred_filt$Mechanism2))
ARG_mech = rep("lightgray", length(shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)]))
ARG_mech[grep("efflux", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#CCFF66"
ARG_mech[grep("inact", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#99CC33"
ARG_mech[grep("alter", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#666600"
ARG_mech[grep("protect", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#00FFFF"
ARG_mech[grep("replace", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#3399FF"
ARG_mech[grep("perm", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "darkblue"

## ARG heatmap (FIG. 4A)
dev.off()
heatmap.2(ARG_heat,
          #Rowv = FALSE,
          Colv = FALSE,
          dendrogram="row",
          key=FALSE, 
          symkey=FALSE, 
          density.info='none',
          trace = "none",
          cexCol = 1,
          srtCol = 65,
          cexRow = 1,
          labRow=as.expression(lapply(rownames(ARG_heat), function(a) bquote(italic(.(a))))),
          labCol = paste(sapply(strsplit(substr(colnames(ARG_heat), start = 3, stop = nchar(as.character(colnames(ARG_heat)))), "_"), function (x) x[1]),
                         substr(sapply(strsplit(substr(colnames(ARG_heat), start = 3, stop = nchar(as.character(colnames(ARG_heat)))), "_"), function (x) x[2]),
                                start = 4, stop = nchar(sapply(strsplit(substr(colnames(ARG_heat), start = 3, 
                                                                               stop = nchar(as.character(colnames(ARG_heat)))), "_"), function (x) x[2]))),
                         sep = "_"),
          col = mypalette, 
          breaks = mycolors,
          #RowSideColors = ARG_class,
          RowSideColors = ARG_mech,
          #ColSideColors = location,
          #ColSideColors = time,
          offsetRow = 0,
          offsetCol = 0,
          margins = c(15, 15))

# legend
dev.off()
plot(1:50, col="white")
legend(1,49, cex = 1, legend = c("Aminoglycoside", "Disinfecting agent", "Carbapenem", "Multidrug", "Rifamycin", "Tetracycline", "Other"),
       c("#FFFFCC", "#FFFF00", "#FFCC33", "#FF6600", "#993300", "#330000", "lightgray")) 

dev.off()
plot(1:50, col="white")
legend(1,49, cex = 1, legend = c("Antibiotic efflux", "Antibiotic inactivation", "Antibiotic target alteration", 
                                 "Antibiotic target protection", "Antibiotic target replacement", "Reduced permeability to antibiotic"),
       c("#CCFF66", "#99CC33", "#666600", "#00FFFF", "#3399FF", "darkblue")) 


ARG_class[grep("amino", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#FFFFCC"
ARG_class[grep("disinf", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#FFFF00"
ARG_class[grep("carb", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#FFCC33"
ARG_class[grep("multidrug", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#FF6600"
ARG_class[grep("rifamycin", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#993300"
ARG_class[grep("tetracycline", shortbred_filt$Class_subset[which(shortbred_filt$TotMarkerLength >= 15)])] = "#330000"
# ARG class
levels(as.factor(shortbred_filt$Mechanism2))
ARG_mech = rep("lightgray", length(shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)]))
ARG_mech[grep("efflux", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#CCFF66"
ARG_mech[grep("inact", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#99CC33"
ARG_mech[grep("alter", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#666600"
ARG_mech[grep("protect", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#00FFFF"
ARG_mech[grep("replace", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "#3399FF"
ARG_mech[grep("perm", shortbred_filt$Mechanism2[which(shortbred_filt$TotMarkerLength >= 15)])] = "darkblue"

# plot key for heat colors
dev.off()
hist(c(0:100), breaks = 100, col = mypalette, 
     border = mypalette, ylim = c(0,0.5),
     xlab = NULL, ylab = NULL, 
     axes = FALSE,
     main = NA)
axis(side = 1, labels = c("0.0", "0.4", "0.8", "1.2", "1.6"),
     at = c(0, 40*100/160, 80*100/160, 120*100/160, 160*100/160), las = 1, cex = 15)


#############################################
### VFs in metagenomes (Shortbred & VFDB) ###
#############################################

## set data
VFDB_heat = as.matrix(shortbred_VFDB_filt[which(shortbred_VFDB_filt$TotMarkerLength >= 15),1:24])
rownames(VFDB_heat) = shortbred_VFDB_filt$Gene2[which(shortbred_VFDB_filt$TotMarkerLength >= 15)]
VFDB_heat = log10(VFDB_heat+1)

## stats

# set table
colnames(VFDB_heat) == metadata$sample.id
VFDB_meta = data.frame(metadata, t(VFDB_heat))
colnames(VFDB_meta)

## re-order
for (i in 1:24) {
  col_order$ORD[i] = grep(colnames(VFDB_heat)[i], col_order$ID)
}
col_order$ID = colnames(VFDB_heat)
col_order

# re-order ARG data
VFDB_heat = VFDB_heat[,order(colnames(VFDB_heat))]
colnames(VFDB_heat) == col_order$ID
VFDB_heat = VFDB_heat[,order(col_order$ORD)]
colnames(VFDB_heat)

## set colnames
for (i in 1:24) {
  colnames(VFDB_heat)[i] = metadata$label.plot[grep(colnames(VFDB_heat)[i], metadata$sample.id)]
}
colnames(VFDB_heat)

## set heat colors
min(VFDB_heat)
max(VFDB_heat)
quantile(VFDB_heat)
mypalette <- colorRampPalette(c("white", "#00CCFF", "darkblue"))(n = 99)
mypalette <- colorRampPalette(c("white", "red", "darkred"))(n = 99)
mycolors = seq(0, 1.55, length.out = 100)

## plot VF heatmap (FIG. 4B)
heatmap.2(VFDB_heat,
          #Rowv = FALSE,
          Colv = FALSE,
          dendrogram="row",
          key=FALSE, 
          symkey=FALSE, 
          density.info='none',
          trace = "none",
          cexCol = 0.8,
          srtCol = 65,
          cexRow = 0.7,
          #labRow=as.expression(lapply(rownames(VFDB_heat), function(a) bquote(italic(.(a))))),
          labCol = paste(sapply(strsplit(substr(colnames(VFDB_heat), start = 3, stop = nchar(as.character(colnames(VFDB_heat)))), "_"), function (x) x[1]),
                         substr(sapply(strsplit(substr(colnames(VFDB_heat), start = 3, stop = nchar(as.character(colnames(VFDB_heat)))), "_"), function (x) x[2]),
                                start = 4, stop = nchar(sapply(strsplit(substr(colnames(VFDB_heat), start = 3, 
                                                                               stop = nchar(as.character(colnames(VFDB_heat)))), "_"), function (x) x[2]))),
                         sep = "_"),
          col = mypalette, 
          breaks = mycolors,
          #RowSideColors = ARG_class,
          #ColSideColors = location,
          #ColSideColors = time,
          offsetRow = 0,
          offsetCol = 0,
          margins = c(10,18))


##################################
### MAGs: gtdb-tk AA alignment ###
##################################

library(ape) # load here to not interfere with functions of 'vegan' package

# load alignment
tree_bac <- read.tree("R_import_data/tree/gtbtk_bac_aln.newick")
tree_tip = read.csv("R_import_data/tree/gtbtk_bac_aln_TAXA.csv", row.names=1, h=T)

# add tip colors
tree_tip$tip_fill = rep("lightgray", length(tree_tip$Phyla))
tree_tip$tip_fill[grep("Actin", tree_tip$Phyla)] = "blue"
tree_tip$tip_fill[grep("Bact", tree_tip$Phyla)] = "gold"
tree_tip$tip_fill[grep("Cyano", tree_tip$Phyla)] = "darkgreen"
tree_tip$tip_fill[grep("Proteo", tree_tip$Phyla)] = "#663300"
tree_tip$tip_fill[grep("Verru", tree_tip$Phyla)] = "purple"

# set tip labels
#tree_bac$tip.label = paste(paste(tree_bac$tip.label, tree_tip$Taxa, sep = " ("), ")", sep = "")
tree_bac$tip.label = tree_tip$Taxa2
tree_bac$tip.label = paste(tree_tip$Taxa, "              ")

# plot
dev.off()
par(mar=c(5,8,5,10))

plot(rotate(rotate(tree_bac, c(22, 22)),
            c(21, 21)),
     #type = "unrooted",
     #show.tip.label = FALSE,
     use.edge.length = FALSE,
     cex = 0.6,
     tip.color = tree_tip$tip_fill,
     #no.margin = TRUE,
     label.offset = 1.2,
     align.tip.label = TRUE,
     edge.width = 1)

# add tips
tiplabels(bg = tree_tip$tip_fill, 
          col = "black",
          #cex = 1.5,
          cex = (tree_tip$AMRF_50 + 1)/3,
          offset = 0.3,
          pch = 21)

# add scale bar
add.scale.bar(ask = TRUE, cex = 0.7)

# legend
dev.off()
plot(1:50, col="white")
legend(1,49, cex = 1, legend = c("Actinobacteriota", "Bacteroidota", "Cyanobacteria", "Proteobacteria", "Verrucomicrobiota"),
       c("blue", "gold", "darkgreen", "#663300", "purple")) 

