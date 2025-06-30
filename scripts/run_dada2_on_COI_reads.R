#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: run_dada2_on_16S_reads.R <input_folder> <output_folder> <log_file>")
}

library(dada2)
library(ggplot2)

MAKE_QC_FIGS <- TRUE
N_THREADS <- 6

input_folder <- normalizePath(args[1])
output_folder <- normalizePath(args[2])
input_folder <- normalizePath("~/project/NEOM_COI/COIrawdata_1/240118_M06860_0205_000000000-LB4YC/Lane1/version_01/")
output_folder <- normalizePath("~/project/NEOM_COI/output")

# Ensure output directory exists
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

log_file <- args[3]

# Start logging
sink(log_file, append = TRUE)
cat("Starting R script at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input folder:", input_folder, "\n")
cat("Output folder:", output_folder, "\n")

# Get full filenames of reads
reads_forward.full <- list.files(input_folder, "R1_001.fastq.cut.gz", full.names = TRUE)
reads_reverse.full <- list.files(input_folder, "R2_001.fastq.cut.gz", full.names = TRUE)

if (length(reads_forward.full) == 0 || length(reads_reverse.full) == 0) {
  stop("Error: No FASTQ files found in the input folder.")
}
if (length(reads_forward.full) != length(reads_reverse.full) ) {
  stop( paste("Error: Different numbers of FORWARD",length(reads_forward.full),"and REVERSE",length(reads_reverse.full),"FASTQ files.") )
}

reads_forward <- basename(reads_forward.full)
reads_reverse <- basename(reads_reverse.full)

# Define output paths
filtered_reads_dir <- file.path(output_folder, "filtered_reads")
dir.create(filtered_reads_dir, showWarnings = FALSE)
reads_forward.filt <- file.path(filtered_reads_dir, reads_forward)
reads_reverse.filt <- file.path(filtered_reads_dir, reads_reverse)

# Generate quality control plots
qc_fig_dir <- file.path(output_folder, "qc_figures")
dir.create(qc_fig_dir, showWarnings = FALSE)
if (MAKE_QC_FIGS) {
  qual_fwd_file <- file.path(qc_fig_dir, "qual_forward_2plot.pdf")
  qual_rev_file <- file.path(qc_fig_dir, "qual_reverse_2plot.pdf")
  
  if (!file.exists(qual_fwd_file)) {
    cat("Plotting forward read quality:", qual_fwd_file, "\n")
    qualp_fwd <- plotQualityProfile(reads_forward.full[2:3])
    ggsave(qual_fwd_file, qualp_fwd, device="pdf", width=7, height=4)
  }
  
  if (!file.exists(qual_rev_file)) {
    cat("Plotting reverse read quality:", qual_rev_file, "\n")
    qualp_rev <- plotQualityProfile(reads_reverse.full[1:2])
    ggsave(qual_rev_file, qualp_rev, device="pdf", width=7, height=4)
  }
}

# Run filter and trim
filter_stats_file <- file.path(output_folder, "filter_stats.txt")
if (file.exists(filter_stats_file)) {
  cat("Reading filter stats from:", filter_stats_file, "\n")
  filter_stats_out <- read.table(filter_stats_file, sep = "\t")
} else {
  cat("Filtering reads:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  filter_stats_out <- filterAndTrim(
    reads_forward.full, reads_forward.filt,
    reads_reverse.full, reads_reverse.filt,
    truncLen = c(270,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
    compress = TRUE, multithread = N_THREADS
  )
  cat( "Writing", nrow(filter_stats_out), "rows of filter stats to:", filter_stats_file, "\n")
  write.table(filter_stats_out, file = filter_stats_file, sep = "\t", quote = FALSE)
}

# Learn error profiles
errF_file <- file.path(output_folder, "err_profile_F_reads.rds")
errR_file <- file.path(output_folder, "err_profile_R_reads.rds")
if (file.exists(errF_file)) {
  cat("Reading error profiles for -FORWARD- reads from:", errF_file, "\n")
  errF <- readRDS(errF_file)
} else {
  cat( "learning error profiles for -FORWARD- reads", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  errF <- learnErrors(reads_forward.filt, multithread=N_THREADS )
  saveRDS(errF, errF_file)
}

if (file.exists(errR_file)) {
  cat("Reading error profiles for -REVERSE- reads from:", errR_file, "\n")
  errR <- readRDS(errR_file)
} else {
  cat( "learning error profiles for -REVERSE- reads", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  errR <- learnErrors(reads_reverse.filt, multithread=N_THREADS )
  saveRDS(errR, errR_file)
}



# Run DADA2
fw_file <- file.path(output_folder, "dada_F_reads.rds")
rv_file <- file.path(output_folder, "dada_R_reads.rds")
if (file.exists(fw_file)) {
  cat("Reading DADA forward reads from:", fw_file, "\n")
  dada_fw <- readRDS(fw_file)
} else {
  cat("Running DADA2 on forward reads:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  dada_fw <- dada(reads_forward.filt, err=errF, multithread=N_THREADS )
  saveRDS(dada_fw, fw_file)
}

if (file.exists(rv_file)) {
  cat("Reading DADA reverse reads from:", rv_file, "\n")
  dada_rv <- readRDS(rv_file)
} else {
  cat("Running DADA2 on reverse reads:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  dada_rv <- dada(reads_reverse.filt, err=errR, multithread=N_THREADS )
  saveRDS(dada_rv, rv_file)
}

# Merge paired reads
merged_file <- file.path(output_folder, "dada_merged_reads.rds")
if (file.exists(merged_file)) {
  cat("Reading merged reads reads from:", merged_file, "\n")
  merged_pairs <- readRDS(merged_file)
} else {
  cat("Merging read pairs:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  merged_pairs <- mergePairs(dada_fw, reads_forward.filt, dada_rv, reads_reverse.filt, verbose=TRUE)
  saveRDS(merged_pairs, merged_file)
}

# Generate ASV table
asv_file <- file.path(output_folder, "asv_table_raw.rds")
if (file.exists(asv_file)) {
  cat("Reading ASV table from:", asv_file, "\n")
  asv_table <- readRDS(asv_file)
} else {
  cat("Making ASV table:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  asv_table <- makeSequenceTable(merged_pairs)
  saveRDS(asv_table, asv_file)
}

N_SITES = nrow(asv_table)
N_SPECIES = ncol(asv_table)

# compile final read tracking results
# and write to table
getN <- function(x) sum(getUniques(x))
read_tracking_file <- file.path(output_folder, paste0("read_tracking_table_raw_",N_SITES,"x",N_SPECIES,".tsv") )
if (!file.exists(read_tracking_file)) {
  cat("Making table of read tracking:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  read_tracking <- cbind( gsub(".fastq.cut.gz","",reads_forward),
                       filter_stats_out, 
                       filter_stats_out[,2] / filter_stats_out[,1] *100 ,
                       sapply(dada_fw, getN), 
                       sapply(dada_rv, getN), 
                       sapply(merged_pairs, getN),
                       sapply(merged_pairs, getN)/ filter_stats_out[,1] *100 )
  colnames(read_tracking) <- c("sample_name", "input", "filtered", "filtered_pct", "denoisedF", "denoisedR", "merged", "merged_pct")
  write.table( read_tracking , file=read_tracking_file, row.names = FALSE, sep="\t", quote = FALSE)
}


# possibly for checking translation codes
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
seq_lengths = nchar(getSequences(asv_table))
seq_lengths.raw = table(nchar(getSequences(asv_table)))
is_right_length = seq_lengths==313 | seq_lengths==310 | seq_lengths==316
is_right_length.tab = table(is_right_length)
is_long_coi = seq_lengths==349 | seq_lengths==352 # these may be junk or bacteria
asv_table.length_filt = asv_table[,which(is_right_length)]

seq_length_plot_file <- file.path(qc_fig_dir, "coi_seq_length_plot.pdf")
if (MAKE_QC_FIGS){
  base_colors = rep("#592f2fff", length( seq_lengths.raw ))
  base_colors[which(names(seq_lengths.raw) %in% c(310,313,316))] = "#408ad5ff"
  pdf(file=seq_length_plot_file, width=6, height=5, title="Unfiltered ASVs by length")
  plot(as.integer(names(seq_lengths.raw)), as.integer(seq_lengths.raw), type='n', 
       xlab="Merged ASV sequence length", ylab="Number of ASVs", xlim=c(300,400),
       main=paste(sum(seq_lengths.raw), "total ASVs;", is_right_length.tab["TRUE"], "are expected length" ),
       cex.lab=1.3, cex.axis=1.3 )
  segments(as.integer(names(seq_lengths.raw)), rep(0,length(seq_lengths.raw)),
           as.integer(names(seq_lengths.raw)), as.integer(seq_lengths.raw), col=base_colors, lwd=4)
  dev.off()
}



# prep files for clustering
# make headers as unique numbers, 1000000 + the column number, for easy sorting
BUFFER_NUMBER = 10^(round(log10(ncol(asv_table.length_filt)))+1)
seq_headers = paste0("ASV", (1:ncol(asv_table.length_filt))+BUFFER_NUMBER )
asv_table_file <- file.path(output_folder, "asv_table_raw.tab" )
if (!file.exists(asv_table_file)){
  cat("Writing ASV table of", length(seq_headers),"samples for", N_SITES, "sites to:", asv_table_file, "\n")
  write.table(t(asv_table.length_filt), file=asv_table_file, 
              quote=FALSE, row.names = seq_headers, sep="\t", col.names =NA )
}
seq_names = getSequences(asv_table.length_filt)
seq_fasta_headers = paste0(">ASV", (1:ncol(asv_table.length_filt))+BUFFER_NUMBER )
asv_fasta_file <- file.path(output_folder, "asv_table_raw.fasta" )
if (!file.exists(asv_fasta_file)){
  cat("Writing", length(seq_names), "ASV sequences to:", asv_fasta_file, "\n")
  write(rbind(seq_fasta_headers,seq_names), asv_fasta_file )
}


cat("Processing completed successfully", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()
#

# values taken from vsearch
N_CLUSTER = 12720
N_SINGLETONS = 9238
barplot( matrix( c(0, ncol(asv_table), 0, ncol(asv_table.length_filt), N_SINGLETONS, N_CLUSTER-N_SINGLETONS), ncol=3 ),
         names.arg = c("Total\nASVs", "Correct\nlength ASVs", "Clusters/\nSingletons"),
         col=c("#85c4f3ff", "#284ab7ff"))



#