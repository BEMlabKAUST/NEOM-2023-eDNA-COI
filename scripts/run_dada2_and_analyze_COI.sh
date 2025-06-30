#!/bin/bash

# Usage: ./run_analysis.sh <input_folder> <output_folder>

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_folder> <output_folder>" >&2
    exit 1
fi

# USERS MUST SET PATHS TO SOME PROGRAMS

# https://github.com/marcelm/cutadapt/
# https://cutadapt.readthedocs.io/en/stable/
CUTADAPT_BIN=/usr/bin/cutadapt

# https://github.com/torognes/vsearch
VSEARCH_BIN=~/vsearch-2.29.1-linux-x86_64/bin/vsearch

#
MERGE_CLUSTERS_PY=~/git/lavaLampPlot/merge_uc_clusters.py

INPUT_FOLDER=$1
OUTPUT_FOLDER=$2
START_TIMESTAMP="$(date +"%Y-%m-%d_%H-%M-%S")"
LOG_FILE="$OUTPUT_FOLDER/run_analysis_$START_TIMESTAMP.log"

CUTADAPT_LOG_FILE="$OUTPUT_FOLDER/run_cutadapt.log"

# Ensure output folder exists
mkdir -p "$OUTPUT_FOLDER"

# Start logging
echo "Starting analysis at $START_TIMESTAMP" > "$LOG_FILE"
echo "Input Folder: $INPUT_FOLDER" >> "$LOG_FILE"
echo "Output Folder: $OUTPUT_FOLDER" >> "$LOG_FILE"

# remove COI primers with cutadapt
for FILE in "$INPUT_FOLDER"/*_R1_001.fastq.gz
do
  BASE="${FILE%_R1_001.fastq.gz}"
  echo "Processing: ${BASE}_R1_001.fastq.gz and ${BASE}_R2_001.fastq.gz" >> "$LOG_FILE"
  $CUTADAPT_BIN -e 0.07 --discard-untrimmed -g ^GGWACWGGWTGAACWGTWTAYCCYCC -G ^TAIACYTCIGGRTGICCRAARAAYCA -o ${BASE}_R1_001.fastq.cut.gz -p ${BASE}_R2_001.fastq.cut.gz ${BASE}_R1_001.fastq.gz ${BASE}_R2_001.fastq.gz >> "$CUTADAPT_LOG_FILE" 2>&1
done 

# Run the R script with user-specified folders
R_DADA_LOG="$OUTPUT_FOLDER/run_DADA2_$START_TIMESTAMP.log"
if Rscript ~/project/NEOM/run_dada2_on_COI_reads.R "$INPUT_FOLDER" "$OUTPUT_FOLDER" "$R_DADA_LOG" >> "$LOG_FILE" 2>&1; then
    cat "$R_DADA_LOG" >> "$LOG_FILE"
    echo "DADA2 pipeline completed successfully at $(date)" >> "$LOG_FILE"
else
    cat "$R_DADA_LOG" >> "$LOG_FILE"
    echo "Error: DADA2 pipeline failed. Check the log file: $LOG_FILE" >&2
    exit 1
fi

# Cluster ASVs using vsearch
ASV_FASTA="$OUTPUT_FOLDER/asv_table_raw.fasta"
CENTROIDS_FASTA="$OUTPUT_FOLDER/asv_table_raw.centroids.fas"
CLUSTERS_UC="$OUTPUT_FOLDER/asv_table_raw.clusters.uc"

if [ -f "$ASV_FASTA" ]; then
    echo "Running vsearch clustering on ASV table..." >> "$LOG_FILE"
    echo "$VSEARCH_BIN --cluster_fast \"$ASV_FASTA\" --id 0.97 --centroids \"$CENTROIDS_FASTA\" --uc \"$CLUSTERS_UC\"" >> "$LOG_FILE"
    $VSEARCH_BIN --cluster_fast "$ASV_FASTA" --id 0.97 \
                 --centroids "$CENTROIDS_FASTA" --uc "$CLUSTERS_UC" >> "$LOG_FILE" 2>&1
    echo "vsearch clustering completed." >> "$LOG_FILE"
else
    echo "Error: ASV FASTA file not found. Cannot do vsearch clustering." >> "$LOG_FILE"
    exit 1
fi

# Merge UC clusters using Python script
MERGED_CLUSTERS_TAB="$OUTPUT_FOLDER/asv_table_clust.clusters.tab"
if [ -f "$CLUSTERS_UC" ]; then
    echo "Running merge_uc_clusters.py to process vsearch results..." >> "$LOG_FILE"
    $MERGE_CLUSTERS_PY "$CLUSTERS_UC" "$OUTPUT_FOLDER/asv_table_raw.tab" > "$MERGED_CLUSTERS_TAB" 2>> "$LOG_FILE"
    echo "UC clustering merge completed." >> "$LOG_FILE"
else
    echo "Error: UC clusters file not found. Cannot run merge step." >> "$LOG_FILE"
    exit 1
fi

# Run the R script with user-specified folders
#if Rscript ~/project/NEOM/run_taxa_on_16S_reads.R "$MERGED_CLUSTERS_TAB" >> "$LOG_FILE" 2>&1; then
#    echo "Analysis completed successfully at $(date)" >> "$LOG_FILE"
#else
#    echo "Error: Analysis failed. Check the log file: $LOG_FILE" >&2
#    exit 1
#fi

