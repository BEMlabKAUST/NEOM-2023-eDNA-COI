#!/usr/bin/Rscript
# last updated by WRF 2025-06-30

library(dplyr)
library(ggplot2)
library(readxl)
library(vegan)
library(dendextend)



################################################################################
# function to reduce taxonomy table to best single name
# Genus species, if possible, otherwise highest taxon level
get_top_taxonomy = function(x){
  x[,"Kingdom"] = ifelse( x[,"Kingdom"]=="Unclassified","Unknown",x[,"Kingdom"])
  x[,"Phylum"] = ifelse( x[,"Phylum"]=="Unclassified",x[,"Kingdom"],x[,"Phylum"])
  x[,"Class"] = ifelse( x[,"Class"]=="Unclassified",x[,"Phylum"],x[,"Class"])
  x[,"Order"] = ifelse( x[,"Order"]=="Unclassified",x[,"Class"],x[,"Order"])
  x[,"Family"] = ifelse( x[,"Family"]=="Unclassified",x[,"Order"],x[,"Family"])
  genus_sp = ifelse( x[,"Genus"]=="Unclassified",x[,"Family"],TRUE)
  genus_sp[genus_sp==TRUE] = ifelse( (x[,"Species"]=="Unclassified")[genus_sp==TRUE], 
                                     paste(x[,"Genus"],"sp")[genus_sp==TRUE], 
                                     x[,"Species"][genus_sp==TRUE] )
  return(genus_sp)
}

################################################################################




# read two files, one taxonomy, one clustered OTUs
taxa.unclust <- read.table("~/project/NEOM_COI/ASVs_taxonomy-no-contam-final-clean.tsv", header=TRUE, sep="\t")
asv_table.clust <- read.table("~/project/NEOM_COI/ASVs_counts-no-contam-final-clean.clusters.tsv", header=TRUE, sep="\t")
# rename headers to short names
# starting from   NEOMeDNAFiltersite12R5_NEOMeDNAFiltersite12R5_N718.N503_L001_R1_001.fastq.cut.gz
# change to       NEOMeDNA12R5
neom_sites_shortnames = gsub("Filtersite","", sapply( strsplit(colnames(asv_table.clust)[-1],"_"), FUN = `[`, 1 ) )
colnames(asv_table.clust) <- c("asv_id", neom_sites_shortnames )
head(asv_table.clust)

# reorder NEOM1 - NEOM12 by Gulf of Aqaba, Inshore, Offshore
site_order_for_paper = c(2,3,7,1, 8,4,9,10, 12,11,5,6)

# this is only for reference, never called
# site 9 has 5 samples since 1 was removed
# site 11 has only 4 samples
# site 12 has 5 samples
base_column_order = c(  9:14, 15:20, 21:26, 27:32, 33:38, 39:44, 
                        45:50, 51:56, 57:61, 62:67, 68:71, 72:76 )
# this is the site order
column_order_for_paper = c( 15:20, 21:26, 45:50,  9:14,
                            51:56, 27:32, 57:61, 62:67,
                            72:76, 68:71, 33:38, 39:44 )
# these have the order of site_order_for_paper , as 2,3,7,1,etc
site_x = rep( seq(1,12,1)[site_order_for_paper], rep(c(6,5,6,4,5),c(8,1,1,1,1))[site_order_for_paper] )
# this has the number of reps of site_x , but ordered 1 to 12
base_x_points = rep( seq(1,12,1), rep(c(6,5,6,4,5),c(8,1,1,1,1))[site_order_for_paper])

# this spaces the bars by site, for visual clarity
bar_spacing_for_paper = c(0.1,ifelse( diff(site_x)==0, 0, 1)+0.1)

# in order of 1 to 13
site_names = c("MAAR", "MQNA", "SMJW", "JUQS", 
               "YUBA", "JZSL", "TINW", "UMHS", 
               "ASSA", "SHAR", "SHSH", "SFIE", "SHJU")[site_order_for_paper]
neom_site_colors = c( "#3BBFDBff", "#0A4470ff", "#4876B8ff", "#10EB88ff", "#EBC849ff", "#EBB449ff", 
                      "#6CB1E6ff", "#BDFA84ff", "#069926ff", "#08701Eff", "#EBDB49ff", "#DCEB49ff" )
neom_site_color_13 = c("#EB9D49ff")
neom_site_colors.lite = gsub("ff","aa",neom_site_colors)



# remove sample NE09eDNAR6, since it has 3401 reads
total_reads_by_sample.all = colSums(asv_table.clust[,-1])
has_enough_reads = which(total_reads_by_sample.all > 20000)+1
asv_table.clust.has_enough = asv_table.clust[,has_enough_reads]
rownames(asv_table.clust.has_enough) <- asv_table.clust$asv_id

# keep unrarefied version, for some numbers in the paper
asv_table.clust.w_taxa.unrare = right_join(taxa.unclust, asv_table.clust , by = "asv_id" )

# for estimate of rough number of reads per sample
png(file="~/project/NEOM_COI/neom_read_counts_by_sample.png", height=800, width=800, res=200 )
hist(total_reads_by_sample.all, breaks=20, col=rep(c("#888888","#f2bf5fff"),c(1,20)), 
     xlab="Total reads per sample (69 samples)" , ylab="Number of samples",
     main="")
dev.off()

# rarefy, for one random iteration
asv_table.clust.rare = apply( asv_table.clust.has_enough, 2, function(x) rrarefy(x, 20000))
# remove OTUs that now occur 0 times
asv_table.clust.rare.has_no0 = which( rowSums(asv_table.clust.rare) > 0 )
asv_table.clust.rare.w_no0 = data.frame( asv_id=rownames(asv_table.clust.has_enough), 
                                         asv_table.clust.rare )[asv_table.clust.rare.has_no0,]

rarefied_otu_table_file = "~/project/NEOM_COI/clustered_otu_table_w_taxonomy.rare.Rds"
#rarefied_otu_table_file = "~/project/NEOM_COI/clustered_otu_table_w_taxonomy.rare.tsv"
if (file.exists(rarefied_otu_table_file)) {
  asv_table.clust.w_taxa = readRDS(rarefied_otu_table_file)
  #asv_table.clust.w_taxa = read.table(rarefied_otu_table_file, sep="\t", header=TRUE )
} else {
  asv_table.clust.w_taxa = right_join(taxa.unclust, asv_table.clust.rare.w_no0 , by = "asv_id" )
  saveRDS(asv_table.clust.w_taxa, rarefied_otu_table_file )
  #write.table(asv_table.clust.w_taxa, rarefied_otu_table_file, sep="\t", quote = FALSE , row.names = FALSE )
}

dim(asv_table.clust.w_taxa)
#[1] 5534   76

# USE THIS TABLE FOR THE REST OF ANALYSIS
species_totals = rowSums(asv_table.clust.w_taxa[,9:76])
asv_table.clust.reorder = arrange(asv_table.clust.w_taxa, desc(species_totals) )[,c(1:8,column_order_for_paper)]
# reassign 2 taxa that somehow hit "Coleoptera"
# https://www.ncbi.nlm.nih.gov/nucleotide/AY955056.1
asv_table.clust.reorder[which(asv_table.clust.reorder$asv_id=="ASV_1100"),3:8] = 
  c("Chlorophyta","Mamiellophyceae","Mamiellales","Mamiellaceae","Micromonas","Micromonas pusilla")
#https://www.ncbi.nlm.nih.gov/nucleotide/AY955049.1
asv_table.clust.reorder[which(asv_table.clust.reorder$asv_id=="ASV_3073"),3:8] = 
  c("Chlorophyta","Mamiellophyceae","Mamiellales","Mamiellaceae","Micromonas","Micromonas pusilla")
# refine one arthropod to Diptera
asv_table.clust.reorder[which(asv_table.clust.reorder$asv_id=="ASV_6913"),4:5] = 
  c("Insecta","Diptera")
saveRDS(asv_table.clust.reorder, "~/project/NEOM_COI/clustered_otu_table_w_taxonomy.rare.sorted.Rds")


total_reads_by_site = colSums(asv_table.clust.reorder[,9:76])


################################################################################

# make presence absence table, and species counts
asv_table.clust.presence = apply( asv_table.clust.reorder[,c(-1:-8)] , 2 , FUN= function(x){as.integer(as.logical(as.integer(x)))} )
n_taxa_by_site = colSums(asv_table.clust.presence)


species_totals.reorder = rowSums(asv_table.clust.reorder[,9:76])
min(which(cumsum(species_totals.reorder)>0.9*sum(species_totals.reorder)))

# rarefaction curve of taxon abundance
#pdf(file="~/project/NEOM_COI/neom_top_250_abund_curve.pdf", height=4, width=4, title="NEOM Sites Top 250 taxa")
png(file="~/project/NEOM_COI/neom_top_250_abund_curve.png", height=800, width=800, res=200 )
plot( cumsum(species_totals.reorder)/sum(species_totals.reorder), type='h', xlim=c(0,250) , ylim=c(0,1), 
      xlab="N species", ylab="Cumulative fraction of total reads", col="#f2bf5fff")
segments(0,0.9,min(which(cumsum(species_totals.reorder)>0.9*sum(species_totals.reorder))),0.9, lwd=2, col="#00000066")
dev.off()



coi_alpha_diversity.raw = diversity(asv_table.clust.reorder[,seq(9,ncol(asv_table.clust.reorder))], MARGIN = 2 )
coi_alpha_diversity.norm = diversity(asv_table.clust.presence, MARGIN = 2 )
coi_alpha_diversity.even = coi_alpha_diversity.raw / coi_alpha_diversity.norm

# counts of OTUs and alpha diversity
# FOR SUPPLEMENT
pdf(file="~/project/NEOM_COI/neom_COI_species_by_site_v4.pdf", height=6, width=5, title="COI diversity in NEOM Sites")
#png(file="~/project/NEOM_COI/neom_COI_species_by_site_v4.png", height=1200, width=1000, res=200)
par(mar=c(3.7,4.5,2.5,1), mfrow=c(2,1))
plot( base_x_points, n_taxa_by_site, axes=FALSE, ylim=c(0,850), 
      xlab="", ylab="Number of unique OTUs (rarefied)", 
      bg=neom_site_colors.lite[site_x], pch=21, cex=3 )
segments(c(0.5),seq(0,800,200),c(12.5),seq(0,800,200), col="#00000022" )
axis(1, at=1:12, labels=site_names, tick=FALSE, las=2 ,line = -0.8)
axis(2, at=seq(0,800,200), labels=c(0,NA,400,NA,800) )
mtext("A",adj=0, cex=2)
plot( base_x_points, coi_alpha_diversity.raw, axes=FALSE, ylim=c(0,7), 
      xlab="", ylab="Alpha diversity (Shannon)", 
      bg=neom_site_colors.lite[site_x], pch=21, cex=3 )
segments(c(0.5),seq(0,7,1),c(12.5),seq(0,7,1), col="#00000022" )
axis(1, at=1:12, labels=site_names, tick=FALSE, las=2 )
axis(2, at=seq(0,7,1) )
mtext("B",adj=0, cex=2)
dev.off()

################################################################################

# colors for top 14 phyla annotations
#  [1] "Unclassified"     "Arthropoda"       "Haptophyta"       "Cnidaria"         "Bacillariophyta"  "Porifera"        
#  [7] "Rhodophyta"       "Ochrophyta"       "Chordata"         "Oomycota"         "Mollusca"         "Annelida"        
# [13] "Chlorophyta"      "Heterokontophyta"
phyla_colors = c("#999999ff", "#f2bf5fff", "#f57acfff", "#d0f076ff", "#aaf093ff", "#f78585ff",
                 "#e65232ff", "#93f3f0ff", "#e5ea42ff", "#e8eac8ff", "#61abf0ff", "#8c5f27ff",
                 "#2a8043ff", "#bb2972ff", "#b2c3e2ff" )

# this is number of OTUs in each phyla, for Supplement
all_phyla = table(asv_table.clust.reorder$Phylum)
phyla_table_by_sample = sapply(asv_table.clust.reorder[9:76], FUN=function(x){
  tab <- table(factor(asv_table.clust.reorder$Phylum[which(as.logical(x))], levels = names(all_phyla)))
  as.integer(tab)
})
rownames(phyla_table_by_sample) <- names(all_phyla)

# total fraction annotated
( sum(all_phyla) - all_phyla[["Unclassified"]] ) / sum(all_phyla)
( sum(all_phyla) - all_phyla[["Unclassified"]] )

unclass_counts = c( ( sum(all_phyla) - table(asv_table.clust.reorder$Kingdom)[["Unclassified"]] ) / sum(all_phyla),
  ( sum(all_phyla) - table(asv_table.clust.reorder$Phylum)[["Unclassified"]] ) / sum(all_phyla),
( sum(all_phyla) - table(asv_table.clust.reorder$Class)[["Unclassified"]] ) / sum(all_phyla),
( sum(all_phyla) - table(asv_table.clust.reorder$Order)[["Unclassified"]] ) / sum(all_phyla),
( sum(all_phyla) - table(asv_table.clust.reorder$Family)[["Unclassified"]] ) / sum(all_phyla),
( sum(all_phyla) - table(asv_table.clust.reorder$Genus)[["Unclassified"]] ) / sum(all_phyla) )

png(file="~/project/NEOM_COI/neom_COI_counts_unclassified_taxa_v1.png", height=600, width=1000, res=200)
par(mar=c(4.5,4.5,1,1))
b = barplot(unclass_counts, ylim=c(0,1), ylab="Fraction annotated", 
            names.arg = names(asv_table.clust.reorder)[2:7], las=2, col="#f2bf5fff")
text(b, unclass_counts+0.05, round(unclass_counts,2) )
dev.off()

#
phyla_table_order = match( names(sort(rowSums(phyla_table_by_sample), decreasing = TRUE))[1:14], names(all_phyla) )
phyla_table.sorted = rbind(phyla_table_by_sample[phyla_table_order,], Others=colSums(phyla_table_by_sample[-(phyla_table_order),]) )
phyla_table.sorted.norm = apply( phyla_table.sorted, 2, FUN=function(x){x/sum(x)*100 } )

# this is number of reads for each phylum
phyla_read_counts = aggregate( asv_table.clust.reorder[,9:76], by=list(asv_table.clust.reorder$Phylum) , FUN = sum )
dim(phyla_read_counts)
phyla_read_counts.sorted = rbind(phyla_read_counts[phyla_table_order,-1], Others=colSums(phyla_read_counts[-1*(phyla_table_order),-1]) )
phyla_read_counts.norm = apply( phyla_read_counts.sorted, 2, FUN=function(x){x/sum(x)*100 } )

# get phyla counts for each site
phyla_read_counts.sites = t(aggregate( t(phyla_read_counts[,-1]), by=list(site_x) , FUN = sum ))[-1,]
phyla_read_counts.sites_sorted = rbind(phyla_read_counts.sites[phyla_table_order,], Others=colSums(phyla_read_counts.sites[-(phyla_table_order),]) )
phyla_read_counts.sites.norm = apply( phyla_read_counts.sites_sorted, 2, FUN=function(x){x/sum(x)*100 } )
rownames(phyla_read_counts.sites) <- phyla_read_counts$Group.1

# number of protist taxa
protist_phyla = c("Amoebozoa", "Apicomplexa", "Ascomycota", "Bacillariophyta", "Basidiomycota",
                  "Cercozoa", "Cryptophyta", "Discosea", "Evosea",  "Chlorophyta", 
                  "Gastrotricha", "Haptophyta", "Heterolobosea", "Mucoromycota", 
                  "Oomycota" )
is_protist = which(rownames(phyla_table_by_sample) %in% protist_phyla )
# count of OTUs
sort( colSums(phyla_table_by_sample[is_protist,] / colSums(phyla_table_by_sample) ) * 100 )
# count of reads
sort( colSums(phyla_read_counts[is_protist,-1] / colSums(phyla_read_counts[,-1]) ) * 100 )
sort( colSums(phyla_read_counts[which(rownames(phyla_table_by_sample)=="Chlorophyta"),-1] / colSums(phyla_read_counts[,-1]) ) * 100 )


# count of Cnidaria
sort(phyla_read_counts.norm[2,])
# count of Cnidaria
sort(phyla_read_counts.norm[3,])



# barplot of phyla aggregated by samples within a site
# NOT HELPFUL
# png(file="~/project/NEOM_COI/neom_COI_phyla_barplot_agg_v3.png", height=1200, width=1600, res=200)
# par(mar=c(6,4,1,1), xpd=TRUE)
b = barplot( as.matrix(phyla_read_counts.sites.norm), col=phyla_colors, axes=TRUE, ylim=c(0,100),
         las=1, cex.names=1.1 , ylab="Number of reads (rarefied)", names.arg = site_names, las=2 )
points(b, rep(-1,12), pch=21, bg=neom_site_colors[site_order_for_paper] )
# dev.off()

# barplot of number of taxa for each phylum, by sample
#pdf(file="~/project/NEOM_COI/neom_COI_phyla_barplot_norm_v3.pdf", height=6, width=8, title="Phyla diversity in NEOM Sites")
png(file="~/project/NEOM_COI/neom_COI_phyla_barplot_norm_v3.png", height=1200, width=1600, res=200)
par(mar=c(6,4,1,1), xpd=TRUE)
b = barplot( phyla_table.sorted.norm, col=phyla_colors, axes=FALSE, ylim=c(0,125),
             las=2, cex.names=0.5 , ylab="Percentage of OTUs (rarefied)", border="#aaaaaaff",
             space= bar_spacing_for_paper , names.arg = gsub("Filtersite","",colnames(phyla_table.sorted.norm)) )
axis(2, at=seq(0,100,20))
legend("top", legend=rownames(phyla_table.sorted), ncol=5, pch=15, col=phyla_colors, text.width = NA )
points( b, rep(-2,68), pch=21, bg=neom_site_colors[site_x] )
dev.off()


# for PCA or NMDS
asv_table.clust.short = asv_table.clust.reorder[,9:76]

# NMDS plot, standalone, also part of Figure 4
sites_to_remove = -1 * which( names(asv_table.clust.short) %in% c("NEOMeDNA11R1","NEOMeDNA11R4") )
coi_nmds = metaMDS( t(asv_table.clust.short[,sites_to_remove]), distance = "hellinger")
pdf(file="~/project/NEOM_COI/neom_coi_nmds_v2.pdf", height=5, width=5, title="COI in NEOM Sites")
plot(coi_nmds, type='n', xlim=c(-0.6,0.6),  ylim=c(-0.6,0.6) )
ordiellipse(coi_nmds,  groups = rep(c(1,2,3),c(24,23,21))[sites_to_remove], col=neom_site_colors[c(1,9,5)], kind = "ehull", draw="polygon", lty=0, alpha=0.2)
points(coi_nmds, display = "sites", cex = 2, pch=21 , bg=neom_site_colors.lite[site_x[sites_to_remove]] )
dev.off()


is_mollusca = which( asv_table.clust.reorder$Phylum=="Mollusca" & asv_table.clust.reorder$Genus!="Tridacna" )
is_tridacna = which(asv_table.clust.reorder$Genus=="Tridacna")

# for ANOVA test
has_tridacna = ifelse(site_x %in% c(3,4,8),1,0)
mollusc_sums = colSums(asv_table.clust.reorder[is_mollusca,9:76])
mollusc_anova = aov(mollusc_sums ~ has_tridacna + site_x  )
summary(mollusc_anova)
# N species of mollusc
sort(apply( asv_table.clust.reorder[is_mollusca,9:76], MARGIN = 2, FUN=function(x){sum(as.integer(as.logical(x)))}))
# percent of mollusc reads
colSums(asv_table.clust.reorder[is_mollusca,9:76]) / total_reads_by_site * 100




# paper Figure 4
# combines barplot of read counts by phylum with NMDS plot
#pdf(file="~/project/NEOM_COI/figure4_neom_COI_phyla_barplot_reads_v5.pdf", height=5, width=8, title="Phyla diversity in NEOM Sites")
png(file="~/project/NEOM_COI/figure4_neom_COI_phyla_barplot_reads_v5.png", height=1500, width=2400, res=300 )
par(mar=c(4,4.2,2,0.5), xpd=TRUE)
layout( matrix(c(1,3,2,4),ncol=2), widths=c(3,2), heights = c(2,1) )
b = barplot( as.matrix(phyla_read_counts.norm), col=phyla_colors, axes=FALSE, ylim=c(0,100),
             las=2, cex.names=0.6 , ylab="Percent of rarefied reads", border="#aaaaaaff",
             space= bar_spacing_for_paper , axisnames = FALSE , cex.lab=1.2 )
axis(2, at=seq(0,100,25), labels=seq(0,100,25) )
mtext(text = site_names, side=1, at=b[c(1,which(diff(bar_spacing_for_paper)==-1))+2], adj=1, line=1 , las=2)
points( b, rep(-2,68), pch=21, bg=neom_site_colors.lite[site_x] )
mtext("A", side=3, at=-18, adj = 0, cex=2, font=2, line=0)
mtext("C", side=1, at=-18, adj = 0, cex=2, font=2, line=2)
par(mar=c(4,2,2,1), xpd=FALSE)
plot(coi_nmds, type='n', xlim=c(-0.6,0.6),  ylim=c(-0.6,0.6) , xlab="", ylab="")
ordiellipse(coi_nmds,  groups = rep(c(1,2,3),c(24,23,21))[sites_to_remove], col=neom_site_colors[c(1,9,5)], kind = "ehull", draw="polygon", lty=0, alpha=0.2)
points(coi_nmds, display = "sites", cex = 2, pch=21 , bg=neom_site_colors.lite[site_x[sites_to_remove]] )
mtext("B", side=3, at=-0.9, adj = 0, cex=2, font=2, line=0)
text(c(-0.4,0.3,0.4),c(0.57,0.4,-0.5),c("Gulf of Aqaba","Nearshore NRS","Offshore NRS"), col=neom_site_colors[c(1,9,5)], cex=1.1)
par(mar=c(1,4.2,1,0.5) , xpd=TRUE )
b = barplot( as.matrix( rbind( colSums(asv_table.clust.reorder[is_mollusca,9:76]) / total_reads_by_site * 100, 
                               colSums(asv_table.clust.reorder[is_tridacna,9:76]) / total_reads_by_site * 100) ),
             col=c("#86bff4ff","#5e688eff"), axes=FALSE , ylab="Percent of rarefied reads", ylim=c(0,25), las=2,
             axisnames = FALSE , space = bar_spacing_for_paper )
axis(2, at=seq(0,25,5), labels=seq(0,25,5) )
legend("topright", col=rev(c("#86bff4ff","#5e688eff")), legend=rev(c("All other molluscs", "Tridacna maxima")), 
       pch=15, pt.cex=2, text.font=c(3,1) )
points( b , rep(-1.2,68), pch=21, bg=neom_site_colors.lite[site_x] )
par(mar=c(1,1,0,0) )
plot(0,0,type='n', xlab=NA, ylab=NA, axes=FALSE )
legend("left", legend=rownames(phyla_table.sorted), ncol=2, pch=15, col=phyla_colors, text.width = NA , pt.cex=2 )
dev.off()


################################################################################

# read benthic data, aggregated
benthic_cover_data = read_xlsx(path = "~/project/NEOM/Benthic cover Neom 2023.xlsx", )[,c(1:4,9:16,5:8,17)]
# reassign 3 categories
# correct name of Hard Coral to lowercase
# merge the few Sediment counts into Sand
# change name of Pavement to Bare rock
benthic_cover_data$Category[benthic_cover_data$Category=="Hard Coral"] = "Hard coral"
benthic_cover_data$Category[benthic_cover_data$Category=="Sediment"] = "Sand"
benthic_cover_data$Category[benthic_cover_data$Category=="Pavement"] = "Bare rock"
# remove Shadow , Tape , and two lines for totals
benth_cats_removed = c(-67,-74,-85,-86)
benthic_cover_data.nm.t = data.frame(t(benthic_cover_data[benth_cats_removed,(site_order_for_paper+3)]))
benthic_cover_data.nm.t[is.na(benthic_cover_data.nm.t)] = 0
dim(benthic_cover_data.nm.t)
benthic_cover_data.summed = aggregate(t(benthic_cover_data.nm.t), by=list(benthic_cover_data$Category[benth_cats_removed]), FUN = sum )
benthic_cover_data.summed.norm = apply( benthic_cover_data.summed[,-1], 2, FUN=function(x){x/sum(x)*100 } )

benthic_cover_data.summed.total = rowSums(benthic_cover_data.summed[,-1])
benthic_cover_data.summed.total
benthic_cover_data.summed.total/sum(benthic_cover_data.summed.total[c(1,2,4,5,11)])

benthic_bar_colors = c("#79b387ff", "#b59bf9ff", "#19752fff", "#94566fff", "#e65232ff", "#5e688eff", "#888888ff", 
                       "#8c8c64ff", "#e8eac8ff", "#f57acfff", "#2a4921ff" )
category_order = c(5,4,10,6,1,3,11,8,9,7,2)



# individual transect counts, for PCA
# make into site order
# except that Site 3 had only 2 transects
site_order_for_paper.transects = c(7:9, 10:11, 21:23, 4:6, 24:26, 12:14, 27:29, 30:32, 36:38, 33:35, 15:17, 18:20)
benthic_cover_data.transects = read_xlsx(path = "~/project/NEOM/Benthic cover NEOM23.xlsx", sheet = 1 )[c(-67,-74,-85),c(1:3,4:6,19:41,7:15)][,c(1:3,site_order_for_paper.transects)]
benthic_cover_data.transects$`Funtional group`[benthic_cover_data.transects$`Funtional group`=="Hard Coral"] = "Hard coral"
benthic_cover_data.transects$`Funtional group`[benthic_cover_data.transects$`Funtional group`=="Sediment"] = "Sand"
benthic_cover_data.transects$`Funtional group`[benthic_cover_data.transects$`Funtional group`=="Pavement"] = "Bare rock"

site_x.transects = rep(site_order_for_paper, rep(c(3,2,3),c(1,1,10)) )
benthic_cover_data.transects.summed = aggregate(benthic_cover_data.transects[,4:38], by=list(benthic_cover_data.transects$`Funtional group`), FUN = sum )

benthic_pca = pca( t(benthic_cover_data.transects.summed[,-1]) )
benthic_constraints.transect = scores(benthic_pca, display = "species")
# standalone PCA, main version in Figure 2
plot(benthic_pca, type='n')
points(benthic_pca, display = "sites", cex = 2, pch=16 , col=neom_site_colors.lite[site_x.transects] )
benthic_constraints = scores(benthic_pca, display = "species")
arrows(c(0),c(0), benthic_constraints.transect[,1], benthic_constraints.transect[,2], length = 0.1, col="#0000aa44" )
text(benthic_pca, display = "species", cex=0.8, col="#0000aaff", labels = c(benthic_cover_data.transects.summed$Group.1) )

benthic_cat_scores = scores(benthic_pca, display = "species")
magnitudes_PC1PC2 <- sqrt(rowSums(benthic_cat_scores[, 1:2]^2))
relative_magnitudes = magnitudes_PC1PC2 / sum(magnitudes_PC1PC2) * 100
sum(relative_magnitudes[c(2,5,9)])




benthic_category_names = benthic_cover_data.summed$Group.1[category_order]

benthic_cover_data.summed.rn = benthic_cover_data.summed
colnames(benthic_cover_data.summed.rn) <- c("Group",site_names)
benthic_cover_data.summed.dist = dist( t(benthic_cover_data.summed.rn[,-1]) )
benthic_cover_data.summed.hc = hclust(benthic_cover_data.summed.dist)


# for paper Figure 3
# dendrogram of benthic data, barplot, and PCA
pdf(file="~/project/NEOM_COI/figure3_dendro_benthicbars_pca_v5.pdf", height=4, width=8, title="Benthic cover in NEOM Sites")
#png(file="~/project/NEOM_COI/figure3_dendro_benthicbars_pca_v5.png", height=1200, width=2400, res=300 )
layout( matrix(c(1,1,2,2,3,4), ncol=3), widths=c(0.9,2.5,2), heights = c(2,1) )
par(mar=c(3,1,1,4), xpd=TRUE )
benthic_cover_data.summed.hc.d = as.dendrogram( benthic_cover_data.summed.hc )
plot( benthic_cover_data.summed.hc.d, horiz=TRUE, axes=FALSE,  )
mtext("A", side=3, adj = 0, cex=2, font=2, line=-2)
points(1:12,1:12, pch=21, cex=2, bg=neom_site_colors[site_order_for_paper][benthic_cover_data.summed.hc$order] )
par(mar=c(3,0,1,1), xpd=FALSE )
barplot( as.matrix(benthic_cover_data.summed.norm[category_order,][,benthic_cover_data.summed.hc$order]) , 
         las=1, axisnames = FALSE, xlim=c(0,110), axes=FALSE,
         col=benthic_bar_colors[category_order] , xlab="Benthic cover %", #names.arg = c(1:13)[benthic_cover_data.summed.hc$order],
         border="#00000022", horiz=TRUE )
mtext("B", side=3, adj = 1, cex=2, font=2, line=-2)
par(mar=c(2,2,1,1) )
benthic_constraints.transect = scores(benthic_pca, display = "species")
plot(benthic_pca, type='n', ylim=c(-18,18))
points(benthic_pca, display = "sites", cex = 2, pch=21 , bg=neom_site_colors.lite[site_x.transects] )
arrows(c(0),c(0), benthic_constraints.transect[,1], benthic_constraints.transect[,2], length = 0.1, col="#0000aa66" )
text(benthic_pca, display = "species", cex=1, col="#0000aaff", labels = gsub(" ","\n",ifelse(benthic_cover_data.transects.summed$Group.1%in%c("Hard coral","Bare rock","Sand"),benthic_cover_data.transects.summed$Group.1,"") ) )
legend("topleft", legend=c("Gulf of Aqaba","Nearshore NRS","Offshore NRS"),
       pch=16, col = neom_site_colors.lite[c(1,9,5)] , cex=0.9, pt.cex=1.1 )
par(mar=c(1,1,1,1))
plot(0,0,type='n', axes=FALSE, xlab=NA, ylab=NA )
legend("left", legend=benthic_category_names , cex=1.1, col=benthic_bar_colors[category_order], ncol=2, pch=15, pt.cex=2 )
dev.off()



################################################################################
# relationships between benthos and taxa


# for Supplemental Figures
is_chlorophyta = asv_table.clust.reorder$Phylum=="Chlorophyta" & asv_table.clust.reorder$Genus!="Micromonas"
is_micromonas = asv_table.clust.reorder$Genus=="Micromonas"


#pdf(file="~/project/NEOM_COI/micromonas_cover_v2.pdf", height=5, width=6, title="Micromonas vs benthic cover")
png(file="~/project/NEOM_COI/micromonas_cover_v2.png", height=1000, width=1200, res=200 )
par(mar=c(4,5,2,1), xpd=TRUE)
layout( matrix(c(1,2,3),nrow=3), heights = c(3,2,2) )
b = barplot( as.matrix( rbind( colSums(asv_table.clust.reorder[is_chlorophyta,9:76]) / total_reads_by_site * 100, 
                               colSums(asv_table.clust.reorder[is_micromonas,9:76]) / total_reads_by_site * 100) ),
             las=2 , cex.names = 0.5, ylim=c(-2,20),
             col=c("#bbbbbbff","#28b74aaa"), axes=TRUE, ylab="Micromonas\nOTU reads (% total)", 
             axisnames = FALSE, space = bar_spacing_for_paper )
mtext(text = site_names, side=1, at=b[c(1,which(diff(bar_spacing_for_paper)==-1))+2], adj=1, line=1 , las=2)
points( b, rep(-1,68), pch=21, bg=neom_site_colors[site_x] )
par(mar=c(1.5,5,1,1))
b = barplot( unlist(benthic_cover_data.nm.t[,c(4)])/ rowSums(benthic_cover_data.nm.t)*100 , ylim=c(0,55),
             col="#b59bf9ff" , ylab="Bare rock\n(% total cover)", axisnames = FALSE , axes=TRUE )
points( b, rep(-4,12), cex=2,pch=21, bg=neom_site_colors[site_order_for_paper] )
b = barplot( as.matrix( rbind(unlist(benthic_cover_data.nm.t[,c(18)])/rowSums(benthic_cover_data.nm.t)*100,
                              unlist(benthic_cover_data.nm.t[,c(35)])/rowSums(benthic_cover_data.nm.t)*100)) , 
             ylim=c(0,15), ylab="Encrusting algae\n(% total cover)", 
             col=c("#8c8c64ff","#8cac64ff") , names.arg = rep("",12) , axes=FALSE )
legend("topright", legend=c("Encrusting algae","Liagora"), col=c("#8c8c64ff","#8cac64ff"), pch=15, pt.cex=2, text.font=c(1,3) )
axis(2, at=c(0,5,10), labels=c(0,5,10) )
dev.off()




pdf(file="~/project/NEOM_COI/tridacna_counts_barplot_v2.pdf", height=5, width=6, title="Tridacna and other mollusc counts")
#png(file="~/project/NEOM_COI/tridacna_counts_barplot_v2.png", height=1000, width=1200, res=200 )
par( mar=c(4,4,1,1), xpd=TRUE)
layout(matrix(c(1,2),nrow=2), heights=c(2,1)) 
b = barplot( as.matrix( rbind( colSums(asv_table.clust.reorder[is_mollusca,9:76]) / total_reads_by_site * 100, 
                               colSums(asv_table.clust.reorder[is_tridacna,9:76]) / total_reads_by_site * 100) ),
             col=c("#bbbbbbff","#5e688eff"), axes=FALSE , ylab="Percent of rarefied reads", ylim=c(0,25), las=2,
             axisnames = FALSE , space = bar_spacing_for_paper )
mtext(text = site_names, side=1, at=b[c(1,which(diff(bar_spacing_for_paper)==-1))+2], adj=1, line=1 , las=2)
axis(2, at=seq(0,25,5), labels=seq(0,25,5) )
legend("topright", col=rev(c("#bbbbbbff","#5e688eff")), legend=rev(c("All other molluscs", "Tridacna maxima")), 
       pch=15, pt.cex=2, text.font=c(3,1) )
points( b , rep(-1.2,68), pch=21, bg=neom_site_colors.lite[site_x] )
par( mar=c(1,4,1,1) )
b = barplot( unlist(benthic_cover_data.nm.t[,22])/rowSums(benthic_cover_data.nm.t)*100 , ylim=c(0,1.2),
             col="#5e688eff" ,ylab="Tridacna (% cover)", axisnames = FALSE , axes=TRUE, las=2 )
dev.off()




# for all other taxa that are also found in benthos

coral_families = c("Xeniidae", "Nephtheidae")
coral_genera = c("Pocillopora", "Galaxea" , "Goniopora" , "Pavona" , "Porites",
                 "Psammocora" , "Sarcophyton" , "Tubipora" )
algae_genera = c( "Dictyota" , "Liagora" , "Lobophora" )

figure_letters = strsplit("ABCDEFGHIJKLMNOPQRSTUV","")[[1]]

pdf(file="~/project/NEOM_COI/benthic_coral_reads_v1.pdf", height=11, width=9, title="Benthic coral cover vs read counts")
#png(file="~/project/NEOM_COI/benthic_coral_reads_v1.p1.png", height=8, width=8, res=200 )
par(mfcol=c(6,2) , mar=c(4,5,3,1), xpd=TRUE )
for (i in 1:length(coral_families)){
  is_coral = asv_table.clust.reorder$Family==coral_families[i]
  coi_barcounts = colSums( asv_table.clust.reorder[is_coral,9:76]) / total_reads_by_site * 100
  benth_barcounts = unlist(benthic_cover_data.nm.t[,match( coral_families[i] , benthic_cover_data$Name[benth_cats_removed] )])/rowSums(benthic_cover_data.nm.t)*100
  #benth_barcounts.reps = c(benth_barcounts)[match(1:12, site_order_for_paper)][site_x]
  b = barplot( coi_barcounts ,  cex.names = 0.3, axisnames = FALSE ,
               space = bar_spacing_for_paper, main=coral_families[i], 
               col="#f57acfff", axes=TRUE, ylab="OTU reads (% total)" )
  mtext(figure_letters[i], side=3, adj=0, cex=2, line=1, at=-18 )
  mtext(text = site_names, side=1, at=b[c(1,which(diff(bar_spacing_for_paper)==-1))+2], adj=1, line=1 , las=2)
  points( b , rep(-1*max(barcounts)/30,68), pch=21, bg=neom_site_colors.lite[site_x] )
  b = barplot(  benth_barcounts, 
               col="#f57acfff" , ylab=paste(coral_families[i],"(% total cover)"), # main=coral_families[i], 
               axisnames = FALSE )
}
barcounts = colSums( asv_table.clust.reorder[is_tridacna,9:76]) / total_reads_by_site * 100
b = barplot( barcounts ,  cex.names = 0.3, axisnames = FALSE ,
             space = bar_spacing_for_paper, main="Tridacna", 
             col="#5e688eff", axes=TRUE, ylab="OTU reads (% total)" )
points( b , rep(-1*max(barcounts)/30,68), pch=21, bg=neom_site_colors.lite[site_x] )
mtext(figure_letters[3], side=3, adj=0, cex=2, line=1, at=-18 )
mtext(text = site_names, side=1, at=b[c(1,which(diff(bar_spacing_for_paper)==-1))+2], adj=1, line=1 , las=2)
b = barplot( unlist(benthic_cover_data.nm.t[,match( "Tridacna" , benthic_cover_data$Name[benth_cats_removed] )])/rowSums(benthic_cover_data.nm.t)*100 , 
             col="#5e688eff" , ylab="Tridacna (% total cover)", #  main="Tridacna", 
             axisnames = FALSE )
for (i in 1:length(coral_genera)){
  target_coral = coral_genera[i]
  is_coral = asv_table.clust.reorder$Genus==coral_genera[i]
  coi_barcounts = colSums( asv_table.clust.reorder[is_coral,9:76]) / total_reads_by_site * 100
  b = barplot( coi_barcounts ,  cex.names = 0.3, axisnames = FALSE ,
               space = bar_spacing_for_paper, main=coral_genera[i], 
               col="#e65232ff", axes=TRUE, ylab="OTU reads (% total)" )
  points( b , rep(-1*max(barcounts)/30,68), pch=21, bg=neom_site_colors.lite[site_x] )
  mtext(figure_letters[i+3], side=3, adj=0, cex=2, line=1, at=-18 )
  mtext(text = site_names, side=1, at=b[c(1,which(diff(bar_spacing_for_paper)==-1))+2], adj=1, line=1 , las=2)
  target_coral = ifelse(target_coral %in% c("Porites"), "Porites massive", target_coral)
  match_pos = match( target_coral , benthic_cover_data$Name[benth_cats_removed] )
  benth_barcounts = unlist(benthic_cover_data.nm.t[,match_pos])/rowSums(benthic_cover_data.nm.t)*100
  # benth_barcounts.reps = c(benth_barcounts)[match(1:12, site_order_for_paper)][site_x]
  b = barplot( benth_barcounts , 
               col="#e65232ff" ,ylab=paste(coral_genera[i],"(% total cover)"), # main=target_coral,
               axisnames = FALSE )
  # plot(benth_barcounts.reps, coi_barcounts,  main=coral_genera[i], 
  #      ylab="(% total reads)", xlab="benthic cover (%)",
  #      pch=21, bg=neom_site_colors.lite[site_x], cex=2 )
}
for (i in 1:length(algae_genera)){
  is_algae = asv_table.clust.reorder$Genus==algae_genera[i]
  barcounts = colSums( asv_table.clust.reorder[is_coral,9:76]) / total_reads_by_site * 100
  b = barplot( barcounts ,  cex.names = 0.3, axisnames = FALSE ,
               space = bar_spacing_for_paper, main=algae_genera[i], 
               col="#79b387ff", axes=TRUE, ylab="OTU reads (% total)" )
  points( b , rep(-1*max(barcounts)/30,68), pch=21, bg=neom_site_colors.lite[site_x] )
  mtext(figure_letters[i+11], side=3, adj=0, cex=2, line=1, at=-18 )
  mtext(text = site_names, side=1, at=b[c(1,which(diff(bar_spacing_for_paper)==-1))+2], adj=1, line=1 , las=2)
  b = barplot( unlist(benthic_cover_data.nm.t[,match( algae_genera[i] , benthic_cover_data$Name[benth_cats_removed] )])/rowSums(benthic_cover_data.nm.t)*100 , 
               col="#79b387ff" ,ylab=paste(algae_genera[i],"(% total cover)"), # main=algae_genera[i], 
               axisnames = FALSE )
}
dev.off()


################################################################################

# take only phyla with greater than 100 reads total
phyla_read_counts.big = phyla_read_counts[rowSums(phyla_read_counts[,-1])>100,]
# flatten data to give same Y data point for each X
benthic_cover_data.summed_reps = benthic_cover_data.summed[,-1][,match(1:12, site_order_for_paper)][,site_x]

# calculate all ANOVA
all_p_values = c()
all_lm_corrs = c()
for (i in 1:length(phyla_read_counts.big$Group.1)){
  pvals = c()
  lm_corrs = c()
  for (j in 1:length(benthic_cover_data.summed$Group)){
    is_target_phylum = asv_table.clust.reorder$Phylum==phyla_read_counts.big$Group.1[i]
    counts_for_target_phylum = colSums(asv_table.clust.reorder[is_target_phylum,9:76])
    current_phy_count = c(counts_for_target_phylum)
    print( paste(phyla_read_counts.big$Group.1[i], "vs", benthic_cover_data.summed$Group[j] ) )
    bc = unlist( benthic_cover_data.summed_reps[j,])
    site_groups = factor(site_x)
    aov_model <- aov(current_phy_count ~ bc + site_groups )
    #summary(aov_model)
    summary_model = summary(aov_model)
    summary_model[[1]]$`Pr(>F)`[1]
    pvals = c(pvals, summary_model[[1]]$`Pr(>F)`[1])
    lm_model <- lm(current_phy_count ~ bc)
    lm_corrs = c(lm_corrs, lm_model$coefficients["bc"])
  }
  all_p_values = rbind(all_p_values, pvals)
  all_lm_corrs = rbind(all_lm_corrs, lm_corrs)
}


# for Figure 5
#pdf(file="~/project/NEOM_COI/figure5_all_phyla_vs_benthic_anova_v2.pdf", height=6, width=8, title="Benthic cover vs Phyla")
png(file="~/project/NEOM_COI/figure5_all_phyla_vs_benthic_anova_v2.png", height=1800, width=2400, res=300)
par(mar=c(8,8,1,2), xpd=TRUE)
layout(matrix(c(1,1,1,1,1,1,2,4,6,3,5,7),ncol=4), widths=c(2,2,2,2))
plot(0,0, type='n', xlim=c(1,11), ylim=c(1,25), axes=FALSE, xlab="", ylab="" )
axis(1, at=seq(1,11,1), labels=benthic_cover_data.summed$Group[category_order], las=2, tick=FALSE)
axis(2, at=seq(1,25,1), labels=rev(phyla_read_counts.big$Group.1), las=2, tick=FALSE)
for (i in 1:25){
  point_size_log = -1*log10(all_p_values[i,])
  # cutoff at 0.01 / 275
  point_size_adj = ifelse(point_size_log>4.4, sqrt(point_size_log-4.4), 0.1 )[category_order]
  points( 1:11, rep(26-i,11), cex= point_size_adj, pch=21, col=benthic_bar_colors[category_order], bg=gsub("ff","aa", benthic_bar_colors)[category_order] )
  corr_direction = ifelse(point_size_log>2, ifelse(all_lm_corrs[i,] >0, "+", "–"), "")
  #corr_direction = ifelse(point_size_log>2, ifelse(all_lm_corrs[i,] >0, benthic_text_colors[category_order], "black"), benthic_text_colors[category_order])
  #corr_direction
  points( 1:11, rep(26-i,11), cex=point_size_adj*0.6, pch=corr_direction[category_order], col="#000000bb")
}
mtext("A", side=3, at=-3.2, adj = 0, cex=2, font=2, line=-2)
par(mar=c(1,1,3,1))
plot(0,0,type='n',axes=FALSE, xlab=NA, ylab=NA )
legend("top", legend=c("NS", expression( ~ 10^-5),  expression( ~ 10^-10),  expression( ~ 10^-20)), pt.cex=c(0.1, sqrt(0.6), sqrt(5.6), sqrt(15.6)), pch=21, 
       col=benthic_bar_colors[5], pt.bg=gsub("ff","aa", benthic_bar_colors)[5], cex=1.4 )
par(mar=c(4,4,3,1))
plot( unlist( benthic_cover_data.summed_reps[benthic_cover_data.summed$Group=="Algae",]) / colSums(benthic_cover_data.summed_reps) * 100,
      c(colSums(asv_table.clust.reorder[asv_table.clust.reorder$Phylum=="Chlorophyta",9:76]) / total_reads_by_site * 100 ), 
      ylab="Chlorophyta (% total reads)", xlab="Algae benthic cover (%)", xlim=c(0,15), ylim=c(0,25),
      pch=21, bg=neom_site_colors.lite[site_x] )
mtext("B", side=3, adj = 0, cex=2, font=2, line=0, at=-3.5 )
plot( unlist( benthic_cover_data.summed_reps[benthic_cover_data.summed$Group=="Bare rock",]) / colSums(benthic_cover_data.summed_reps) * 100,
      c(colSums(asv_table.clust.reorder[asv_table.clust.reorder$Phylum=="Chlorophyta",9:76]) / total_reads_by_site * 100 ), 
      ylab="Chlorophyta (% total reads)", xlab="Bare rock benthic cover (%)",  xlim=c(0,60), ylim=c(0,25),
      pch=21, bg=neom_site_colors.lite[site_x] )
mtext("C", side=3, adj = 0, cex=2, font=2, line=0.5, at=-14 )
plot( unlist( benthic_cover_data.summed_reps[benthic_cover_data.summed$Group=="Sand",])/ colSums(benthic_cover_data.summed_reps) * 100,
      c(colSums(asv_table.clust.reorder[asv_table.clust.reorder$Phylum=="Chlorophyta",9:76]) / total_reads_by_site * 100 ), 
      ylab="Chlorophyta (% total reads)", xlab="Sand benthic cover (%)", xlim=c(0,60), ylim=c(0,25),
      pch=21, bg=neom_site_colors.lite[site_x] )
mtext("D", side=3, adj = 0, cex=2, font=2, line=0.5, at=-14 )
plot( unlist( benthic_cover_data.summed_reps[benthic_cover_data.summed$Group=="Soft coral",])/ colSums(benthic_cover_data.summed_reps) * 100,
      c(colSums(asv_table.clust.reorder[asv_table.clust.reorder$Phylum=="Ochrophyta",9:76]) / total_reads_by_site * 100 ), 
      ylab="Ochrophyta (% total reads)", xlab="Soft coral cover (%)", xlim=c(0,20), ylim=c(0,5),
      pch=21, bg=neom_site_colors.lite[site_x] )
mtext("E", side=3, adj = 0, cex=2, font=2, line=0.5, at=-4.5 )
plot( unlist( benthic_cover_data.summed_reps[benthic_cover_data.summed$Group=="Invertebrates",])/ colSums(benthic_cover_data.summed_reps) * 100,
      c(colSums(asv_table.clust.reorder[asv_table.clust.reorder$Phylum=="Ochrophyta",9:76]) / total_reads_by_site * 100 ), 
      ylab="Ochrophyta (% total reads)", xlab="Invertebrate benthic cover (%)", xlim=c(0,3), ylim=c(0,5),
      pch=21, bg=neom_site_colors.lite[site_x] )
mtext("F", side=3, adj = 0, cex=2, font=2, line=0.5, at=-0.7 )
dev.off()

################################################################################

# try constrained ordination
# does not really work, but include in supplement anyway

top250 = asv_table.clust.reorder[1:250,]
top250_merged = aggregate( t(top250[,9:76]), by=list(site_x), sum)[,-1]
benthic_cover.12only = t(benthic_cover_data.summed[,2:13])
colnames(benthic_cover.12only) <- benthic_cover_data.summed$Group.1
top250_benthic_cca = cca(formula = top250_merged ~ . , data = data.frame(benthic_cover.12only) )
#pdf(file="~/project/NEOM_COI/benthic_top250_cca_v1.pdf", height=5, width=7, title="Benthic cover vs Top 250 species")
png(file="~/project/NEOM_COI/benthic_top250_cca_v1.png", height=1000, width=1400, res=200)
par(mar=c(4,4,1,1) )
plot(top250_benthic_cca, type='n')
points(top250_benthic_cca, display="sites", c(1:12), col=neom_site_colors.lite, cex=2, pch=16 )
points(top250_benthic_cca, display = "species", cex=0.5, pch=2, col="#00000066" )
#text(top250_benthic_cca, display="species", labels=top250$asv_id )
text(top250_benthic_cca, display = "cn", col="#0000ffcc"  )
#text(top250_benthic_cca, display="sites", c(1:12))
dev.off()

################################################################################

# make mantel test with bray curtis distance for both datasets
asv_table.clust.reorder.sites = t(aggregate( t(asv_table.clust.reorder[,9:76]), by=list(site_x) , FUN = sum ))[-1,]
asv_table.clust.reorder.sites.no0 = asv_table.clust.reorder.sites[which(rowSums(asv_table.clust.reorder.sites)>1),]
asv_table.clust.reorder.sites.dist <- vegdist( t(asv_table.clust.reorder.sites.no0), method = "bray")
benthic_cover_data.summed.dist <- vegdist( t(benthic_cover_data.summed[,-1]), method = "bray")
mantel_result_bray <- mantel(asv_table.clust.reorder.sites.dist, benthic_cover_data.summed.dist, permutations = 999)
mantel_result_bray
summary(mantel_result_bray)
# Mantel statistic r: -0.4155 
# Significance: 0.992 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
#   0.314 0.412 0.472 0.564 
# Permutation: free
# Number of permutations: 999

# Mantel test for each parameter individually
mantel_result.solo_r = c()
mantel_result.solo_p = c()
for (i in 1:11){
  benthic_cover_data.summed.solo = t(benthic_cover_data.summed[i,-1])
  benthic_cover_data.summed.is_not_zero = which(benthic_cover_data.summed.solo>0)
  benthic_cover_data.summed.dist.solo = vegdist( t(benthic_cover_data.summed[i,-1][,benthic_cover_data.summed.is_not_zero]), method = "bray")
  asv_table.clust.reorder.sites.dist.solo <- vegdist( t(asv_table.clust.reorder.sites.no0[,benthic_cover_data.summed.is_not_zero]), method = "bray")
  mantel_result_bray.solo = mantel(asv_table.clust.reorder.sites.dist.solo, benthic_cover_data.summed.dist.solo, permutations = 999)
  mantel_result.solo_r = c(mantel_result.solo_r, mantel_result_bray.solo$statistic )
  mantel_result.solo_p = c(mantel_result.solo_p, mantel_result_bray.solo$signif)
}
mantel_result.solo_r
mantel_result.solo_p


# for PERMANOVA
asv_table.clust.reorder.dist <- vegdist( t(asv_table.clust.reorder[,9:76]), method = "bray")
site_n_group = data.frame( site=names(asv_table.clust.short), groups = rep(c(1,2,3),c(24,23,21)))
permanova_bray = adonis2( asv_table.clust.reorder.dist ~ groups, data = site_n_group )
permanova_bray
# adonis2(formula = asv_table.clust.reorder.dist ~ groups, data = site_n_group)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     1   0.8738 0.11615 8.6732  0.001 ***
# Residual 66   6.6496 0.88385                  
# Total    67   7.5234 1.00000   

#