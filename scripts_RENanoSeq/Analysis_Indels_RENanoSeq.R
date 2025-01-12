# ANALYSIS OF SOMATIC INDELS IN CHOLANGIOCYTE ORGANOIDS (Restriction Enzyme NanoSeq)
# Adrian Baez-Ortega, Wellcome Sanger Institute, 2023-24

# This analysis is part of the publication:
# Petrus-Reurer, Baez-Ortega, Lei et al., (2025). 
# Immune and Mutational Profile of Gene-Edited ‘Universal’ Low-Immunogenic Human Primary Cholangiocyte Organoids


# *NB*. This script requires access to a FASTA file for the reference human genome (GRCh37)*
GENOME.PATH = "/path/to/Homo_sapiens.GRCh37.XX.genome.fna"


# Internal sample IDs
SAMPLES = c("325_WT_P0"  = "PD56206b_ds0001",  "325_WT_P1"   = "PD56206h_ds0001", 
            "325_DKO_P1" = "PD56206i_ds0001",  "299_WT_Prim" = "PD61632b_ds0001",
            "299_WT_P0"  = "PD61632c_ds0001",  "299_WT_P1"   = "PD61632e_ds0001",
            "299_WT_P5"  = "PD61632g_ds0001",  "299_DKO_P1"  = "PD61632d_ds0001",
            "299_DKO_P5" = "PD61632f_ds0001",  "808_WT_Prim" = "PD61633b_ds0001",
            "808_WT_P0"  = "PD61633c_ds0001",  "808_WT_P1"   = "PD61633d_ds0001",
            "808_DKO_P1" = "PD61633e_ds0001")

# Load required packages: Biostrings, dndscv, sigfit
library(dndscv)
library(sigfit)
library(Biostrings)


# Source function to create indel spectra
# (Indelwald.R by Maximilian Stammnitz: https://github.com/MaximilianStammnitz/Indelwald)
source("scripts_aux/Indelwald.R")


# (1) Load reference genome (GRCh37) and NanoSeq indels
genome = readDNAStringSet(GENOME.PATH, format="fasta", use.names=TRUE)
names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)

indels.nanoseq = sapply(SAMPLES, function(id) {
    v = read.table(gzfile(paste0("data_RENanoSeq/", id, ".botseq.indel.vcf.gz")), header=F, sep="\t", as.is=T,
                   col.names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"))
    v[grep("PASS", v$FILTER), ]
}, simplify=F)


# (2) Generate and plot indel spectra
counts.nanoseq = t(sapply(indels.nanoseq, function(id) {
    indel.spectrum(id[, c("CHROM", "POS", "REF", "ALT")], genome)
}))

plot_spectrum(counts.nanoseq, paste0(Sys.Date(), "_Spectra_Indels_RENanoSeq.pdf"))

