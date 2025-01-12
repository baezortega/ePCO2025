# ANALYSIS OF SOMATIC SUBSTITUTIONS IN CHOLANGIOCYTE ORGANOIDS (Targeted NanoSeq)
# Adrian Baez-Ortega, Wellcome Sanger Institute, 2023-24

# This analysis is part of the publication:
# Petrus-Reurer, Baez-Ortega, Lei et al., (2025). 
# Immune and Mutational Profile of Gene-Edited ‘Universal’ Low-Immunogenic Human Primary Cholangiocyte Organoids


# Internal sample IDs
SAMPLES = c("299_WT_P0"  = "PD61632c_tds0001",  "299_WT_P1"  = "PD61632e_tds0001",
            "299_DKO_P1" = "PD61632d_tds0001",  "808_WT_P0"  = "PD61633c_tds0001",
            "808_WT_P1"  = "PD61633d_tds0001",  "808_DKO_P1" = "PD61633e_tds0001")


# Load required packages: dndscv, grid, lattice
library(grid)
library(dndscv)
library(lattice)


# (1) Load NanoSeq mutations
muts = read.table("data_TNanoSeq/plate_043_SQPP-25885-G_targeted_2_decontamination_final_muts_dnds_annotated.tsv",
                  header=T, sep="\t", as.is=T)
muts$sampleID2 = names(SAMPLES)[match(muts$sampleID, SAMPLES)]
stopifnot(all(muts$sampleID %in% SAMPLES))


# (2) Load mean gene coverage table; exclude genes with duplex coverage = 0
gene.cov = read.csv("data_TNanoSeq/plate_043_SQPP-25885-G_targeted_2_decontamination_mean_gene_cov.csv",
                    as.is=T, row.names=1)
gene.cov = gene.cov[rowSums(gene.cov) > 0, ]


# (3) Run dNdScv on deduplicated mutations from each ePCO line
muts.dedup = unique(cbind(substr(muts$sampleID2, 1, 3), muts[, c("chr", "pos", "ref", "mut")]))
dnds.out = dndscv(muts.dedup, gene_list=rownames(gene.cov), dc=rowSums(gene.cov),
                  max_muts_per_gene_per_sample=Inf)

dnds.fname = paste0(Sys.Date(), "_Targeted_dNdS.txt")
cat("GLOBAL DN/DS TABLE\n\n", file=dnds.fname)
write.table(dnds.out$globaldnds, sep="\t", row.names=F, quote=F, append=T, file=dnds.fname)
cat("\n\nDN/DS PER GENE (TOP 20)\n\n", append=T, file=dnds.fname)
write.table(dnds.out$sel_cv[1:20, ], sep="\t", row.names=F, quote=F, append=T, file=dnds.fname)


# (4) For each mutated gene, calculate the fraction of mutant cells in each sample:
#       sum(2*duplex_vaf) across all mutations in that gene and sample
gene.cell.fraction = t(sapply(sort(unique(na.omit(muts$gene))), function(gene) {
    sapply(rev(SAMPLES), function(id) {
        idx = (muts$gene == gene) %in% T & muts$sampleID == id
        if (!any(idx)) {
            0
        } else {
            sum(2 * muts$duplex_vaf[idx])
        }
    })
}))


# (5) Plot mutant cell fractions for all mutated genes
cairo_pdf(paste0(Sys.Date(), "_CellFractions_Targeted.pdf"), 25, 4)
levelplot(gene.cell.fraction,
          xlab="", ylab="", scales=list(x=list(rot=45, tck=0), y=list(tck=0)),
          col.regions=colorRampPalette(c("white", "orange", "firebrick")), 
          main=list("\nFraction of mutant cells per gene per sample", cex=1.5),
          panel=function(x, y, z, ...) {  # to label values in each cell
              grid.rect(gp=gpar(col=NA, fill="grey95"))
              panel.levelplot(x, y, z, ...)
              panel.text(x, y, ifelse(z == -1, "-", as.character(round(z, 2))), cex=0.5)
              panel.abline(h=3.5)
          })
dev.off()


# (6) Plot only genes with >2% total mutant cells across all samples
cairo_pdf(paste0(Sys.Date(), "_CellFractions_Targeted_2pc.pdf"), 10, 4)
gene.cf.2 = gene.cell.fraction[rowSums(gene.cell.fraction) > 0.02, ] * 100
colnames(gene.cf.2) = c(paste("Line 3", c(" DKO (P1)", "   WT (P1)", "   WT (P0)")),
                        paste("Line 2", c(" DKO (P1)", "   WT (P1)", "   WT (P0)")))
levelplot(gene.cf.2,
          xlab="", ylab="", scales=list(x=list(rot=45, tck=0), y=list(tck=0)),
          col.regions=colorRampPalette(c("white", "orange", "firebrick")), 
          main=list("\nPercentage of mutant cells per gene per sample (genes with >2% across all samples)", cex=1),
          panel=function(x, y, z, ...) {  # to label values in each cell
              grid.rect(gp=gpar(col=NA, fill="grey95"))
              panel.levelplot(x, y, z, ...)
              panel.text(x, y, ifelse(z == -1, "-", as.character(round(z, 1))), cex=0.55)
              panel.abline(h=3.5)
          })
dev.off()

