# ANALYSIS OF OFF-TARGET SITES IN CHOLANGIOCYTE ORGANOIDS (WGS)
# Adrian Baez-Ortega, Wellcome Sanger Institute, 2023-24

# This analysis is part of the publication:
# Petrus-Reurer, Baez-Ortega, Lei et al., (2025). 
# Immune and Mutational Profile of Gene-Edited ‘Universal’ Low-Immunogenic Human Primary Cholangiocyte Organoids


# *NB*. This script requires access to a GTF file for the human genome (Ensembl GRCh37)*
GTF.PATH = "/path/to/Homo_sapiens.GRCh37.XX.chr.gtf.gz"


# Internal sample IDs and mutation type equivalences
SAMPLES = c("325_DKO_P1" = "PD56206i",
            "299_DKO_P1" = "PD61632d",
            "808_DKO_P1" = "PD61633e")


# Load required packages: GenomicRanges, rtracklayer, stringr
library(stringr)
library(rtracklayer)
library(GenomicRanges)


# (1) Load data objects
# CRISPR off-target site predictions (from Cas-OFFinder using ≤9 mismatches, bulge sizes ≤2)
offtarget = read.table("data_WGS/CasOFFinder_Offtarget_B2M-RFX5_m9-bs2-bs2.txt",
                       sep="\t", header=T, comment.char="", as.is=T)
offtarget$Chromosome = gsub("chr", "", offtarget$Chromosome)

# Merge overlapping off-target regions (100-bp window)
offtarget$Start = offtarget$Position - 50
offtarget$End = offtarget$Position + 50
offtarget.gr = makeGRangesFromDataFrame(offtarget)
offtarget.gr = union(offtarget.gr, offtarget.gr)

# Ensembl GTF hg19 annotation
gtf.gr = import(GTF.PATH)
cds.gr = gtf.gr[gtf.gr@elementMetadata$type == "CDS"]

# Cancer Gene Census (genes with causative SNVs/indels or large deletions)
cgc = read.table("data_WGS/COSMIC_CancerGeneCensus_v98_GRCh37.tsv", sep="\t", header=T, as.is=T, quote="")
cgc = cgc[grepl("(Mis)|F|N|S|D", cgc$MUTATION_TYPES), ]
cgc.gr = makeGRangesFromDataFrame(cgc, seqnames.field="CHROMOSOME",
                                  start.field="GENOME_START", end.field="GENOME_STOP")
# Mutation calls
subs = sapply(SAMPLES, function(id) {
    read.table(paste0("data_WGS/", id, "_filtered.txt"), sep="\t", header=T, as.is=T)
}, simplify=F)

indels = sapply(SAMPLES, function(id) {
    x = read.table(gzfile(paste0("data_WGS//", id, ".pindel.annot.vcf.gz")), sep="\t", header=F, as.is=T,
                   col.names=c("Chr", "Start", "", "Ref", "Alt", "", "Filter", "Info", "Format", "Nrm", "Tum"))
    x = x[x$Filter == "PASS", ]
    x$End = x$Start + abs(nchar(x$Ref) - nchar(x$Alt))
    x$Func.refGene = as.character(sapply(str_split(x$Info, ";"), function(y) {
        for (i in 1:length(y)) {
            if (substr(y[i],1,2) == "VC") return(substr(y[i],4,nchar(y[i])))
        }
    }))
    x$Gene.refGene = as.character(sapply(str_split(x$Info, ";"), function(y) {
        for (i in 1:length(y)) {
            if (substr(y[i],1,2) == "VD") return(str_split(substr(y[i],4,nchar(y[i])), "\\|")[[1]][1])
        }
    }))
    x
}, simplify=F)

# mpileup calls at predicted CRISPR off-target sites
offtg.calls = NULL
for (id in names(SAMPLES)) {
    x = read.table(gzfile(paste0("data_WGS/Offtarget_", id, ".vcf.gz")), sep="\t", header=F, as.is=T,
                   col.names=c("CHROM", "POS", ".", "REF", "ALT", ".", ".", "INFO", ".", "TUM"))
    x = x[!duplicated(paste0(x$CHROM, ":", x$POS, ":", x$REF, ":", x$ALT)), ]
    offtg.calls = rbind(offtg.calls,
                        cbind("SAMPLE"=id, x[order(x$CHROM, x$POS), c("CHROM","POS","REF","ALT","INFO","TUM")]))
}

# Merge mutations for each sample
muts = sapply(names(subs), function(id) {
    m = rbind(subs[[id]][, c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")],
              indels[[id]][, c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene")])
    m[!(m$Gene.refGene %in% c("RFX5", "B2M")), ]  # Exclude target genes
}, simplify=F)
muts.gr = sapply(muts, function(m) makeGRangesFromDataFrame(m))


# (2) Overlap off-target sites with exons, CGC genes, mutations
outfile = paste0(Sys.Date(), "_Offtarget_Overlaps_WGS.txt")
cat("CDS overlaps:\n",
    sum(countOverlaps(offtarget.gr, cds.gr) > 0),
    " (", round(mean(countOverlaps(offtarget.gr, cds.gr) > 0)*100, 2),
    "%) off-target sites overlapping CDS regions\n\n",
    sep="", file=outfile)
cat("Cancer Gene Census overlaps:\n",
    sum((countOverlaps(offtarget.gr, cds.gr) > 0) & (countOverlaps(offtarget.gr, cgc.gr) > 0)),
    " (", round(mean((countOverlaps(offtarget.gr, cds.gr) > 0) & (countOverlaps(offtarget.gr, cgc.gr) > 0))*100, 2),
    "%) off-target sites overlapping CDS of cancer genes\n",
    sep="", file=outfile, append=T)
for (i in 1:length(muts.gr)) {
    overlap.idx = countOverlaps(muts.gr[[i]], offtarget.gr) > 0
    coding.idx = muts[[i]]$Func.refGene %in% c("ess_splice", "splice_region", "frameshift", "inframe", "exonic")
    cat("\n-------------------------", names(muts.gr)[i], "-------------------------\n\n", file=outfile, append=T)
    cat("Mutations overlapping off-target sites: ",
        sum(overlap.idx),
        " (", round(mean(overlap.idx)*100, 2), "%)\n",
        "Mutations overlapping off-target sites and CDS regions: ",
        sum(countOverlaps(muts.gr[[i]], subsetByOverlaps(offtarget.gr, cds.gr)) > 0),
        " (", round(mean(countOverlaps(muts.gr[[i]], subsetByOverlaps(offtarget.gr, cds.gr)) > 0)*100, 2), "%)\n",
        "Mutations overlapping off-target sites and CGC genes: ",
        sum(unlist(str_split(muts[[i]]$Gene.refGene[overlap.idx], ",")) %in% cgc$GENE_SYMBOL),
        "\n\nDetected off-target mutations:\n", 
        sep="", file=outfile, append=T)
    write.table(muts[[i]][overlap.idx, ], sep="\t", quote=F, row.names=F, file=outfile, append=T)
}


# (3) Examine mpileup calls at predicted off-target sites for clusters of variants
# Search for runs of 5 variants within 20 bp
offtg.calls.cluster = offtg.calls[abs(diff(offtg.calls$POS, lag=4)) <= 20, ]
offtg.calls.cluster.gr = makeGRangesFromDataFrame(offtg.calls.cluster, start.field="POS", end.field="POS")

# Check for overlaps with coding sequences
offtg.calls.cluster = cbind("OVERLAPS_CDS" = countOverlaps(offtg.calls.cluster.gr, cds.gr) > 0,
                            offtg.calls.cluster)
write.table(offtg.calls.cluster, sep="\t", quote=F, row.names=F,
            file=paste0(Sys.Date(), "_Offtarget_Clusters_WGS.tsv"))

