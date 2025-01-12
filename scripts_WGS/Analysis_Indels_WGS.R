# ANALYSIS OF SOMATIC INDELS IN CHOLANGIOCYTE ORGANOIDS (WGS)
# Adrian Baez-Ortega, Wellcome Sanger Institute, 2023-24

# This analysis is part of the publication:
# Petrus-Reurer, Baez-Ortega, Lei et al., (2025). 
# Immune and Mutational Profile of Gene-Edited ‘Universal’ Low-Immunogenic Human Primary Cholangiocyte Organoids


# *NB*. This script requires access to a FASTA file for the reference human genome (GRCh37)*
GENOME.PATH = "/path/to/Homo_sapiens.GRCh37.XX.genome.fna"


# Internal sample IDs and mutation categories
SAMPLES = c("325_DKO_P1" = "PD56206i",
            "299_DKO_P1" = "PD61632d",
            "808_DKO_P1" = "PD61633e")
MUT.TYPES = c("Ins", "Del")


# Load required packages: Biostrings, dndscv, sigfit, stringr
library(dndscv)
library(sigfit)
library(stringr)
library(Biostrings)


# Source function to create indel spectra
# (Indelwald.R by Maximilian Stammnitz: https://github.com/MaximilianStammnitz/Indelwald)
source("scripts_aux/Indelwald.R")


# (1) Load reference genome, Cancer Gene Census table, WGS mutation calls
genome = readDNAStringSet(GENOME.PATH, format="fasta", use.names=TRUE)
names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)

cgc = read.table("data_WGS/COSMIC_CancerGeneCensus_v98_GRCh37.tsv", sep="\t", header=T, quote="", as.is=T)
cgc = cgc[grepl("(Mis)|F|N|S|D", cgc$MUTATION_TYPES), ]

vars = sapply(SAMPLES, function(id) {
    read.table(gzfile(paste0("data_WGS/", id, ".pindel.annot.vcf.gz")), sep="\t", header=F, as.is=T,
                   col.names=c("Chr", "Start", "", "Ref", "Alt", "", "Filter", "Info", "Format", "Nrm", "Tum"))
}, simplify=F)


# (2) Retrieve variant information
for (i in 1:length(vars)) {
    # Coverage and VAF
    vars[[i]]$NR = sapply(vars[[i]]$Tum, function(x) {
        sum(as.numeric(str_split_fixed(x, ":", 18)[, 14:15]))
    })
    vars[[i]]$NV = as.numeric(str_split_fixed(vars[[i]]$Tum, ":", 18)[, 14])
    vars[[i]]$VAF = as.numeric(str_split_fixed(vars[[i]]$Tum, ":", 18)[, 18])
    # Variant type
    vars[[i]]$Type = ifelse(nchar(vars[[i]]$Ref) > nchar(vars[[i]]$Alt), "Del", "Ins")
}


# (3) Plot mutational spectra
spectra = t(sapply(vars, function(v) {
    indel.spectrum(structure(as.matrix(v[, c(1:2,4:5)]),
                             dimnames=list(NULL, c("CHROM","POS","REF","ALT"))), genome)
}))
rownames(spectra) = names(vars)
plot_spectrum(spectra, paste0(Sys.Date(), "_Spectra_Indels_WGS.pdf"))


# (4) Plot coverage and VAF histograms
cairo_pdf(paste0(Sys.Date(), "_Coverage-VAF_Indels_WGS.pdf"), 12, 6, onefile=T)
par(mfrow=c(1,2), mgp=c(2,0.6,0), tck=-0.025, cex.lab=1.2, cex.main=1.3)
for (i in 1:length(vars)) {
    plot(density(vars[[i]]$NR, adjust=0.5),
         col="dodgerblue4", lwd=4, xlim=c(0, 100), xlab="Coverage at variant site",
         main=paste0(names(vars)[i], " (", nrow(vars[[i]]), " PASS indels)\nSequencing coverage"))
    abline(v=median(vars[[i]]$NR), lty=2, lwd=2, col="grey40")
    plot(density(vars[[i]]$VAF, adjust=0.5),
         col="dodgerblue4", lwd=4, xlim=c(0, 1), xlab="Variant allele fraction",
         main=paste0(names(vars)[i], " (", nrow(vars[[i]]), " PASS indels)\nVariant allele fraction"))
    abline(v=median(vars[[i]]$VAF), lty=2, lwd=2, col="grey40")
}
dev.off()


# (5) Annotate indels with dNdScv
# Create mutation table
mut.table = NULL
for (i in 1:length(vars)) {
    mut.table = rbind(mut.table, cbind("Sample"=names(vars)[i], vars[[i]][, c("Chr", "Start", "Ref", "Alt")]))
}
# Add one substitution to avoid dndscv error
mut.table = rbind(mut.table,
                  data.frame("Sample"="325_DKO_P1", "Chr"=4, "Start"=76692246, "Ref"="T", "Alt"="C",
                             stringsAsFactors=F))

# Run dndscv
annotmuts = dndscv(mutations=mut.table)$annotmuts

# Output genes with nonsynonymous mutations
for (id in names(vars)) {
    genes = sort(unique(annotmuts$gene[annotmuts$sampleID == id & annotmuts$impact != "Synonymous"]))
    if (length(genes) == 0) genes = "None"
    cat(id,
        "\nMutated (indel) genes:", paste(genes, collapse=", "),
        "\nMutated genes in CGC:",
        if (!any(genes %in% cgc$GENE_SYMBOL)) {
            "None"
        } else {
            paste(genes[genes %in% cgc$GENE_SYMBOL], collapse=", ")
        }, "\n\n",
        append=T, file=paste0(Sys.Date(), "_Annotation_Indels_WGS.txt"))
}

