# ANALYSIS OF SOMATIC SUBSTITUTIONS IN CHOLANGIOCYTE ORGANOIDS (Restriction Enzyme NanoSeq)
# Adrian Baez-Ortega, Wellcome Sanger Institute, 2023-24

# This analysis is part of the publication:
# Petrus-Reurer, Baez-Ortega, Lei et al., (2025). 
# Immune and Mutational Profile of Gene-Edited ‘Universal’ Low-Immunogenic Human Primary Cholangiocyte Organoids


# Human genome sizes (total & coding) and internal sample IDs
GENOME.SIZE = 6.2e9     # total human genome size
CODING.SIZE = 69553910  # coding genome size from dNdScv [data("refcds_hg19"); sum(width(gr_genes))*2]
SAMPLES = c("325_WT_P0"  = "PD56206b_ds0001",  "325_WT_P1"   = "PD56206h_ds0001", 
            "325_DKO_P1" = "PD56206i_ds0001",  "299_WT_Prim" = "PD61632b_ds0001",
            "299_WT_P0"  = "PD61632c_ds0001",  "299_WT_P1"   = "PD61632e_ds0001",
            "299_WT_P5"  = "PD61632g_ds0001",  "299_DKO_P1"  = "PD61632d_ds0001",
            "299_DKO_P5" = "PD61632f_ds0001",  "808_WT_Prim" = "PD61633b_ds0001",
            "808_WT_P0"  = "PD61633c_ds0001",  "808_WT_P1"   = "PD61633d_ds0001",
            "808_DKO_P1" = "PD61633e_ds0001")


# Load required packages: sigfit
library(sigfit)


# (1) Load NanoSeq mutations, spectra and burdens, and sample dates
subs.nanoseq = sapply(SAMPLES, function(id) {
    v = read.table(gzfile(paste0("data_RENanoSeq/", id, ".botseq.muts.vcf.gz")), header=F, sep="\t", as.is=T,
                   col.names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
    v[grep("PASS", v$FILTER), ]
}, simplify=F)

counts.nanoseq = t(sapply(SAMPLES, function(id) {
    read.table(list.files("data_RENanoSeq/", paste0(id, "_vs_PD.*.trint_subs_obs_corrected.tsv"), full.names=T),
               sep="\t", header=T, as.is=T)[, "trint_onto_genome"]
}))

burdens.nanoseq = as.data.frame(t(sapply(SAMPLES, function(id) {
    as.matrix(read.table(list.files("data_RENanoSeq/", paste0(id, "_vs_PD.*.mut_burden.tsv"), full.names=T),
                         sep="\t", header=T, as.is=T))["corrected", ]
})))

sample.dates = read.table("data_RENanoSeq/Sample_Passages_Dates.txt", sep="\t", header=T, as.is=T)
sample.dates = sample.dates[match(rownames(burdens.nanoseq), sample.dates$sample), ]


# (2) Calculate mutation burdens and rates per diploid genome
burdens.nanoseq$burden_genome = burdens.nanoseq$burden * GENOME.SIZE
burdens.nanoseq$burden_genome_lci = burdens.nanoseq$burden_lci * GENOME.SIZE
burdens.nanoseq$burden_genome_uci = burdens.nanoseq$burden_uci * GENOME.SIZE
burdens.nanoseq$burden_coding = burdens.nanoseq$burden * CODING.SIZE
burdens.nanoseq$days_since_p0 = sample.dates$days_since_p0
burdens.nanoseq$burden_genome_p0 = burdens.nanoseq[sample.dates$p0_sample, "burden_genome"]
burdens.nanoseq$mut_rate = ifelse(burdens.nanoseq$days_since_p0 == 0, 0,
                                  (burdens.nanoseq$burden_genome - burdens.nanoseq$burden_genome_p0) /
                                      burdens.nanoseq$days_since_p0)
burdens.nanoseq$sample = rownames(burdens.nanoseq)
write.table(burdens.nanoseq[, c(ncol(burdens.nanoseq), 1:(ncol(burdens.nanoseq)-1))],
            file=paste0(Sys.Date(), "_Burdens_Subs_RENanoSeq.tsv"), sep="\t", quote=F, row.names=F)


# (3) Plot mutation burdens and spectra
# Burdens
cairo_pdf(paste0(Sys.Date(), "_Burdens_Subs_RENanoSeq.pdf"), 16, 8, onefile=T)
par(mar=c(6.5,6.5,6,0))
b = barplot(burdens.nanoseq$burden_genome,
            col=rep(c("steelblue3", "steelblue4", "steelblue3"), c(3, 6, 4)), border=NA, las=2,
            cex.lab=1.3, cex.names=1e-9, ylim=c(0, 4000), ylab="Mutations per diploid genome\n")
text(x=b, y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=paste0(burdens.nanoseq$sample, "\n(", burdens.nanoseq$days_since_p0, " days)"),
     srt=45, adj=1, xpd=TRUE, cex=1.2)
text(b, burdens.nanoseq$burden_genome_uci + 230, xpd=NA, font=2,
     labels=paste0("Total: ", round(burdens.nanoseq$burden_genome), "\n\n"))
text(b, burdens.nanoseq$burden_genome_uci + 230, xpd=NA, 
     labels=paste0("\n(", round(burdens.nanoseq$burden_genome_lci), "–",
                   round(burdens.nanoseq$burden_genome_uci), ")",
                   "\nCoding: ", round(burdens.nanoseq$burden_coding)))
segments(x0=b, y0=burdens.nanoseq$burden_genome_uci, y1=burdens.nanoseq$burden_genome_lci, lwd=2, xpd=NA)
dev.off()

# Spectra
plot_spectrum(counts.nanoseq, paste0(Sys.Date(), "_Spectra_Subs_RENanoSeq.pdf"))


# (4) Extract mutational signatures
# Use FitExt with fixed SBS1, SBS5 + 2 de novo signatures
data("cosmic_signatures_v3.2")
sig.fitext.2.2 = fit_extract_signatures(counts.nanoseq,
                                        signatures=cosmic_signatures_v3.2[c("SBS1", "SBS5"), ],
                                        num_extra_sigs=2, iter=15000, warmup=5000, seed=0xC0FFEE)
plot_all(sig.fitext.2.2, out_path=paste0(Sys.Date(), "_Signatures_FitExt_2+2"))

# Re-fit signatures to all samples,
# using only those signatures with credible exposures (lower_95 > 0.01) for each sample 
sigs = retrieve_pars(sig.fitext.2.2, "signatures")$mean
rownames(sigs) = c("SBS1", "SBS5", "Signature A", "Signature B")
expos.lwr95 = retrieve_pars(sig.fitext.2.2, "exposures")$lower_95
sig.idx = expos.lwr95 > 0.01
sig.idx[, 1:2] = TRUE  # Force SBS1 & SBS5
sig.fit.4 = sapply(rownames(counts.nanoseq), function(smp) {
    cat(smp, "\n")
    fit_signatures(counts.nanoseq[smp,,drop=F],
                   sigs[sig.idx[smp,], ],
                   iter=4000, chains=3, cores=3, seed=0xC0FFEE)
}, simplify=F)

# Retrieve final signature exposures
exposures = list(
    "mean" = t(sapply(sig.fit.4, function(fit) {
        x = retrieve_pars(fit, "exposures")$mean
        unlist(ifelse(rownames(sigs) %in% names(x), x, 0))
    })),
    "lower_95" = t(sapply(sig.fit.4, function(fit) {
        x = retrieve_pars(fit, "exposures")$lower_95
        unlist(ifelse(rownames(sigs) %in% names(x), x, 0))
    })),
    "upper_95" = t(sapply(sig.fit.4, function(fit) {
        x = retrieve_pars(fit, "exposures")$upper_95
        unlist(ifelse(rownames(sigs) %in% names(x), x, 0))
    }))
)
colnames(exposures$mean) = rownames(sigs)

# Plot signature exposures
cols = c("#228B22", "#FF8C00", "#9932CC", "firebrick3")
dir.create(paste0(Sys.Date(), "_Signatures_Refit"))
plot_spectrum(sigs, paste0(Sys.Date(), "_Signatures_Refit/Signatures_Refit_", Sys.Date(), ".pdf"))
plot_exposures(counts=counts.nanoseq, exposures=exposures,
               pdf_path=paste0(Sys.Date(), "_Signatures_Refit/Exposures_Refit_", Sys.Date(), ".pdf"),
               legend_pos="topright", sig_color_palette=cols, signature_names=rownames(sigs))

# Plot signature-specific mutation burdens
cairo_pdf(paste0(Sys.Date(), "_Signatures_Refit/Burdens_Refit_", Sys.Date(), ".pdf"), 18, 5, onefile=T)
sig.burdens.1 = list("mle" = burdens.nanoseq$burden_genome * exposures$mean,
                     "lci" = burdens.nanoseq$burden_genome_lci * exposures$mean,
                     "uci" = burdens.nanoseq$burden_genome_uci * exposures$mean)
sig.burdens.2 = list("mle" = burdens.nanoseq$burden_genome * exposures$mean,
                     "lci" = burdens.nanoseq$burden_genome_lci * exposures$lower_95,
                     "uci" = burdens.nanoseq$burden_genome_uci * exposures$upper_95)
b = barplot(t(sig.burdens.1$mle), beside=T, ylim=c(0, 2000), col=cols, border=cols, 
            main="Signature-specific burdens per genome, using burden CIs")
segments(x0=as.numeric(b), y0=as.numeric(t(sig.burdens.1$lci)), y1=as.numeric(t(sig.burdens.1$uci)),
         col="darkgrey", lwd=1.75, xpd=NA)
legend("topleft", legend=rownames(sigs), bty="n", fill=cols, border=NA, inset=c(0.02, 0), cex=1.1)
b = barplot(t(sig.burdens.2$mle), beside=T, ylim=c(0, 2000), col=cols, border=cols, 
            main="Signature-specific burdens per genome, using burden and exposure CIs")
segments(x0=as.numeric(b), y0=as.numeric(t(sig.burdens.2$lci)), y1=as.numeric(t(sig.burdens.2$uci)),
         col="darkgrey", lwd=1.75, xpd=NA)
legend("topleft", legend=rownames(sigs), bty="n", fill=cols, border=NA, inset=c(0.02, 0), cex=1.1)
dev.off()

