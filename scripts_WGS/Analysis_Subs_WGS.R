# ANALYSIS OF SOMATIC SUBSTITUTIONS IN CHOLANGIOCYTE ORGANOIDS (WGS)
# Adrian Baez-Ortega, Wellcome Sanger Institute, 2023-24

# This analysis is part of the publication:
# Petrus-Reurer, Baez-Ortega, Lei et al., (2025). 
# Immune and Mutational Profile of Gene-Edited ‘Universal’ Low-Immunogenic Human Primary Cholangiocyte Organoids


# *NB*. This script requires access to a FASTA file for the reference human genome (GRCh37)*
GENOME.PATH = "/path/to/Homo_sapiens.GRCh37.XX.genome.fna"


# Internal sample IDs and mutation type equivalences
SAMPLES = c("325_DKO_P1" = "PD56206i",
            "299_DKO_P1" = "PD61632d",
            "808_DKO_P1" = "PD61633e")
MUT.TYPES = c("G>T" = "C>A", "G>C" = "C>G", "G>A" = "C>T",
              "A>T" = "T>A", "A>G" = "T>C", "A>C" = "T>G")


# Load required packages: Biostrings, dndscv, sigfit
library(dndscv)
library(sigfit)
library(Biostrings)


# Function: reverse complement (for string vectors)
rev.comp = function(nucleotide.list) {
    sapply(nucleotide.list, function(nucleotides) {
        paste(
            rev(sapply(strsplit(nucleotides, "")[[1]], function(nuc) {
                if (nuc == "A") "T"
                else if (nuc == "C") "G"
                else if (nuc == "G") "C"
                else if (nuc == "T") "A"
            })),
            collapse="")
    }, USE.NAMES=F)
}


# (1) Load reference genome, Cancer Gene Census table, WGS mutation calls
genome = readDNAStringSet(GENOME.PATH, format="fasta", use.names=TRUE)
names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)

cgc = read.table("data_WGS/COSMIC_CancerGeneCensus_v98_GRCh37.tsv", sep="\t", header=T, quote="", as.is=T)
cgc = cgc[grepl("(Mis)|F|N|S|D", cgc$MUTATION_TYPES), ]

vars = sapply(SAMPLES, function(id) {
    read.table(paste0("data_WGS/", id, "_filtered.txt"), sep="\t", header=T, as.is=T)
}, simplify=F)


# (2) Retrieve variant information
for (i in 1:length(vars)) {
    # Trinucleotide context
    vars[[i]]$Context = as.character(padAndClip(genome[vars[[i]]$Chr],
                                                IRanges(vars[[i]]$Start - 1, vars[[i]]$Start + 1),
                                                Lpadding.letter=".", Rpadding.letter="."))
    stopifnot(identical(vars[[i]]$Ref, substr(vars[[i]]$Context, 2, 2)))
    # Coverage and VAF
    vars[[i]]$NR = vars[[i]]$depth
    vars[[i]]$VAF = vars[[i]]$alt.freq
    # Variant type
    vars[[i]]$Type = paste0(vars[[i]]$Ref, ">", vars[[i]]$Alt)
    for (j in 1:length(MUT.TYPES)) {
        vars[[i]]$Type[vars[[i]]$Type == names(MUT.TYPES)[j]] = MUT.TYPES[j]
    }
}


# (3) Plot mutational spectra
mut.table = NULL
for (i in 1:length(vars)) {
    mut.table = rbind(mut.table, cbind(names(vars)[i], vars[[i]][, c("Ref", "Alt", "Context")]))
}
plot_spectrum(build_catalogues(mut.table), paste0(Sys.Date(), "_Spectra_Subs_WGS.pdf"))


# (4) Plot coverage and VAF histograms
cairo_pdf(paste0(Sys.Date(), "_Coverage-VAF_Subs_WGS.pdf"), 12, 6, onefile=T)
par(mfrow=c(1,2), mgp=c(2,0.6,0), tck=-0.025, cex.lab=1.2, cex.main=1.3)
for (i in 1:length(vars)) {
    plot(density(vars[[i]]$NR, adjust=0.5),
         col="dodgerblue4", lwd=4, xlim=c(0, 100), xlab="Coverage at variant site",
         main=paste0(names(vars)[i], " (", nrow(vars[[i]]), " mutations)\nSequencing coverage"))
    abline(v=median(vars[[i]]$NR), lty=2, lwd=2, col="grey40")
    plot(density(vars[[i]]$VAF, adjust=0.5),
         col="dodgerblue4", lwd=4, xlim=c(0, 1), xlab="Variant allele fraction",
         main=paste0(names(vars)[i], " (", nrow(vars[[i]]), " mutations)\nVariant allele fraction"))
    abline(v=median(vars[[i]]$VAF), lty=2, lwd=2, col="grey40")
}
dev.off()


# (5) Annotate substitutions with dNdScv
# Create mutation table
mut.table = NULL
for (i in 1:length(vars)) {
    mut.table = rbind(mut.table, cbind(names(vars)[i], vars[[i]][, c("Chr", "Start", "Ref", "Alt")]))
}

# Run dndscv
annotmuts = dndscv(mutations=mut.table, outp=1)$annotmuts

# Output genes with nonsynonymous mutations
for (id in names(vars)) {
    genes = sort(unique(annotmuts$gene[annotmuts$sampleID == id & annotmuts$impact != "Synonymous"]))
    if (length(genes) == 0) genes = "None"
    cat(id,
        "\nMutated (nonsyn) genes:", paste(genes, collapse=", "),
        "\nMutated genes in CGC:",
        if (!any(genes %in% cgc$GENE_SYMBOL)) {
            "None"
        } else {
            paste(genes[genes %in% cgc$GENE_SYMBOL], collapse=", ")
        }, "\n\n",
        append=T, file=paste0(Sys.Date(), "_Annotation_Subs_WGS.txt"))
}

