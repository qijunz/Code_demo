### Purpose: QTL mapping for 16S traits in DO mice
### Created: 2024-02-16

# load required libraries

source("/Users/rootqz/R/QZ_functions.R")
options(stringsAsFactors = FALSE)

setwd("/Users/rootqz/Desktop/ReyLab/project/DO_diet/")

# load geno prob and pmap
probs <- readRDS("data/final_genoprobs_1176.rds")
kinship <- readRDS("data/K_1176_DietDO.rds")
pmap <- readRDS("data/gigamuga_map_test_v2.RDS")
markers <- fread("data/marker_with_IDs_for2nd_map.csv")
markers <- markers %>%
    select(-V1, -probe) %>%
    distinct(remarker, .keep_all = TRUE)

# load phenotype data
phe_all <- fread("data/Warren_MF_DO_Study_SCFA_Rel_Abun_Attie_Phenotype_Metadata_2_16_2024_v2.csv") %>% as.data.frame()

# set covar
phe_covar <- phe_all %>%
    mutate(Diet = as.factor(Diet),
           Sex = as.factor(Sex),
           Wave = as.factor(Wave),
           DNA_Batch = as.factor(DNA_Batch),
           Seq_Run = as.factor(Seq_Run)) %>%
    column_to_rownames("Sample_ID")

covar_add <- model.matrix(~Sex*Diet+Wave+DNA_Batch+Seq_Run, data=phe_covar)[,-1]
covar_int <- t(t(model.matrix(~Diet, data=phe_covar)[,-1]))

traits_16s <- names(phe_all)[38:119]

## get QTL peaks
qtl_scan_add <- readRDS("qtl_scan/16S_20240216/qtl_scan_16s_add_20240216.rds")
qtl_scan_int <- readRDS("qtl_scan/16S_20240216/qtl_scan_16s_int_20240216.rds")

qtl_scan_HF <- readRDS("qtl_scan/16S_20240216/qtl_scan_16s_HF_20240216.rds")
qtl_scan_HC <- readRDS("qtl_scan/16S_20240216/qtl_scan_16s_HC_20240216.rds")

# diet as interactive covar, LOD > 10
qtl_peak <- as.data.frame(matrix(ncol = 7, nrow = 0)) 
colnames(qtl_peak)<- c("lodindex", "lodcolumn", "chr", "pos", "lod", "ci_lo", "ci_hi")

for (trait in traits_16s) {
    
    # get peaks
    peak_all <- find_peaks(qtl_scan_int[, trait, drop=F], pmap, threshold = 10, prob=0.95, peakdrop = 2)
    if (nrow(peak_all) > 0 ) {
        qtl_peak <<- rbind(qtl_peak, peak_all)
    }
}

qtl_peak$lodindex <- NULL                      
colnames(qtl_peak)<- c("pheno", "chr", "peak_mbp", "lod", "ci_lo", "ci_hi")

qtl_peak_int <- qtl_peak %>%
    left_join(markers %>% rename_at("bp_mm10", ~ "peak_mbp"), by = c("chr", "peak_mbp")) %>%
    select(pheno, chr, peak_mbp, lod, ci_lo, ci_hi, remarker, marker, snp) %>%
    filter(lod > 7) %>%
    arrange(desc(lod)) %>%
    mutate(add_covar = "Sex*Diet+Wave+DNA_Batch+Seq_Run",
           int_covar = "Diet")

fwrite(qtl_peak_int, file = "qtl_scan/16S_20240216/qtl_peak_16s_int_LOD10.csv", sep = ",")

# diet as additive covar, LOD > 7
qtl_peak <- as.data.frame(matrix(ncol = 7, nrow = 0)) 
colnames(qtl_peak)<- c("lodindex", "lodcolumn", "chr", "pos", "lod", "ci_lo", "ci_hi")

for (trait in traits_16s) {
    
    # get peaks
    peak_all <- find_peaks(qtl_scan_add[, trait, drop=F], pmap, threshold = 7, prob=0.95, peakdrop = 2)
    if (nrow(peak_all) > 0 ) {
        qtl_peak <<- rbind(qtl_peak, peak_all)
    }
}

qtl_peak$lodindex <- NULL                      
colnames(qtl_peak)<- c("pheno", "chr", "peak_mbp", "lod", "ci_lo", "ci_hi")

qtl_peak_add <- qtl_peak %>%
    left_join(markers %>% rename_at("bp_mm10", ~ "peak_mbp"), by = c("chr", "peak_mbp")) %>%
    select(pheno, chr, peak_mbp, lod, ci_lo, ci_hi, remarker, marker, snp) %>%
    filter(lod > 7) %>%
    arrange(desc(lod)) %>%
    mutate(add_covar = "Sex*Diet+Wave+DNA_Batch+Seq_Run",
           int_covar = "NA")

fwrite(qtl_peak_add, file = "qtl_scan/16S_20240216/qtl_peak_16s_add_LOD7.csv", sep = ",")

## Manhattan plot 
for (trait in traits_16s) {
    
    ymax <- max(max(qtl_scan_add[, trait, drop=F],
                    qtl_scan_int[, trait, drop=F]))*1.02
    
    ymax_diet <- max(max(qtl_scan_HF[, trait, drop=F]), 
                     max(qtl_scan_HC[, trait, drop=F]))*1.02
    
    pdf(paste0("qtl_scan/16S_20240216/qtl_16s/", trait, ".pdf"), width = 12, height = 16)
    
    par(mfcol=c(4,1), mar = c(2.5, 2.5, 2.5, 2.5), cex=1.2)
    plot_scan1(qtl_scan_int[, trait, drop=F], pmap, ylim = c(0, ymax), xlab ="", col = "#1F2544")
    plot_scan1(qtl_scan_add[, trait, drop=F], pmap, ylim = c(0, ymax), xlab ="", col = "#D04848", add = T)
    title(trait)
    
    lod_diff <- qtl_scan_int[, trait, drop=F]-qtl_scan_add[, trait, drop=F]
    plot_scan1(lod_diff, pmap, ylim = c(0, max(lod_diff)*1.1), xlab ="", col = "#12372A")
    title("LOD difference")
    
    plot_scan1(qtl_scan_HF[, trait, drop=F], ylim = c(0, ymax_diet), pmap, xlab ="", col = "#1F2544")
    title("HF mice")
    
    plot_scan1(qtl_scan_HC[, trait, drop=F], ylim = c(0, ymax_diet), pmap, xlab ="", col = "#1F2544")
    title("HC mice")
    
    dev.off()
}

## SNP association
# load gene and snp database
query_gene <- create_gene_query_func("~/Desktop/DO/data/mouse_genes_mgi.sqlite")
query_variants <- create_variant_query_func("~/Desktop/DO/data/cc_variants.sqlite")

phe_16s <- phe_all %>%
    select(Sample_ID, Diet, Sex, any_of(traits_16s)) %>%
    column_to_rownames("Sample_ID") %>%
    mutate(across(Shannon:`Unknown Ruminococcaceae (Family)`, rankZ)) %>%
    select(-Diet, -Sex)

for (i in c(1:nrow(qtl_peak_int))) {
    
    this_peak <- qtl_peak_int[i,,drop=F]
    
    chr <- this_peak$chr
    phename <- this_peak$pheno
    marker <- this_peak$marker
    x_min <- this_peak$peak_mbp-5
    x_max <- this_peak$peak_mbp+5
    
    cat(paste0(phename, " Chr", chr, "...\n"))
    
    phe <- phe_16s %>% 
        select(one_of(phename)) %>%
        na.omit()
    
    qtl_snp <- scan1snps(probs, pmap, phe, kinship[[chr]], 
                         addcovar = covar_add, intcovar = covar_int, 
                         chr = chr, start = x_min, end = x_max, query_func=query_variants, keep_all_snps=TRUE)
    snpinfo <- index_snps(pmap, qtl_snp$snpinfo)                
    snpprobs <- genoprob_to_snpprob(probs, snpinfo)         
    scan_snppr <- scan1(snpprobs, phe, kinship[[chr]], addcovar = covar_add, intcovar = covar_int)
    qtl_snp$lod <- scan_snppr
    
    saveRDS(qtl_snp, file = paste0("qtl_scan/16S_20240216/snp_scan/", phename, "-Chr", chr, "-", round(this_peak$peak_mbp), "Mbp-snp.rds"))
    
    # plot snp association and gene
    qtl_gene <- query_gene(chr = chr, x_min, x_max) %>% filter(bioType == "protein coding gene")
    
    pdf(paste0("qtl_scan/16S_20240216/snp_gene/", phename, "-Chr", chr, "-", round(this_peak$peak_mbp), "Mbp-snp.pdf"), width = 12, height = 6)
    par(cex=1.2)
    plot_snpasso(scan_snppr, qtl_snp$snpinfo, genes=qtl_gene, colors = "black", top_panel_prop = 0.5)
    dev.off()
    
    # plot alleles of lead snp
    lead_snp <- rownames(scan_snppr)[which(scan_snppr[,1] == max(scan_snppr))]
    lead_snp_allele <- as.data.frame(snpprobs[[1]][,,lead_snp]) %>%
        mutate(allele = case_when(A > 0.9 ~ "AA", 
                                  B > 0.9 ~ "BB",
                                  abs(A-B) < 0.1 ~ "AB",
                                  .default = NA))
    
    allele_boxplot <- lead_snp_allele %>%
        rownames_to_column("Sample_ID") %>%
        left_join(phe_all, by ="Sample_ID") %>%
        filter(!is.na(allele)) %>%
        filter(!is.na(Diet)) %>%
        mutate() %>%
        ggplot(aes(x = allele, y = get(phename), color = allele)) +
        geom_boxplot(outlier.shape = NA) + 
        geom_point(position=position_jitterdodge(), size = 1.2) +
        theme_classic() +
        theme(axis.title = element_text(size = 15),
              axis.text = element_text(size = 15),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              strip.text = element_text(size = 15),
              legend.position = "none",
              title = element_text(size=15)) +
        facet_grid(. ~ Diet) +
        scale_color_manual(values = c("#265073", "#2D9596", "#9AD0C2")) +
        scale_y_continuous(limits = c(0, NA)) +
        ggtitle(lead_snp) +
        ylab(phename)
    
    pdf(paste0("qtl_scan/16S_20240216/snp_allele/", phename, "-Chr", chr, "-", round(this_peak$peak_mbp), "Mbp-snp.pdf"), width = 5, height = 4)
    print(allele_boxplot)
    dev.off()
    
}


## allele effects, using function scan1coef()
qtl_coef_int <- data.frame(matrix(nrow = nrow(qtl_peak_int), ncol = 8))
names(qtl_coef_int) <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
qtl_coef_int <- cbind(qtl_peak_int, qtl_coef_int)

for (i in c(1:nrow(qtl_peak_int))) {
    
    this_peak <- qtl_peak_int[i,,drop=F]
    
    chr <- this_peak$chr
    phename <- this_peak$pheno
    marker <- this_peak$marker
    x_min <- this_peak$peak_mbp-5
    x_max <- this_peak$peak_mbp+5
    
    cat(paste0(phename, " Chr", chr, "...\n"))
    
    qtl_coef <- readRDS(paste0("qtl_scan/16S_20240216/coef_16s/", phename, "-Chr", chr, "-", round(this_peak$peak_mbp), "Mbp-coef.rds"))
    qtl_coef_zoomin <- qtl_coef %>%
        as.data.frame() %>%
        select(A:H) %>%
        rownames_to_column("marker") %>%
        left_join(markers, by = "marker") %>%
        filter(bp_mm10 > x_min, bp_mm10 < x_max)
    
    ylim <- max(abs(qtl_coef_zoomin[,2:9]))*1.05
    
    pdf(paste0("qtl_scan/16S_20240216/qtl_coef/", phename, "-Chr", chr, "-", round(this_peak$peak_mbp), "Mbp-coef.pdf"), width = 8, height = 5)
    plot_coefCC(qtl_coef, scan1_output = qtl_scan_int[, phename, drop=F], pmap, bgcolor="white", col = CCcolor, xlim=c(x_min, x_max), ylim=c(-ylim, ylim))
    dev.off()
    
    # add coef to QTL table
    qtl_coef_int[i,12:19] <- qtl_coef[which(rownames(qtl_coef) == marker),c(1:8)]
}

fwrite(qtl_coef_int, file = "qtl_scan/16S_20240216/qtl_coef_16s_int_LOD10.csv", sep = ",")
