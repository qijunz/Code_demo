### Purpose: QTL mapping for 16S traits in DO mice
### Created: 2024-02-16

# load required libraries
library(qtl2)
library(dplyr)
library(tibble)
library(data.table)

options(stringsAsFactors = FALSE)

# load geno prob and pmap
probs <- readRDS("final_genoprobs_1176.rds")
kinship <- readRDS("K_1176_DietDO.rds")
pmap <- readRDS("gigamuga_map_test_v2.RDS")

# rankZ function
rankZ = function(x) {
    x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
    return(qnorm(x))
}

# load phenotype data
phe_all <- fread("Warren_MF_DO_Study_SCFA_Rel_Abun_Attie_Phenotype_Metadata_2_16_2024_v2.csv") %>% as.data.frame()

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

## QTL for all mice
phe_16s <- phe_all %>%
    select(Sample_ID, Diet, Sex, any_of(traits_16s)) %>%
    column_to_rownames("Sample_ID") %>%
    mutate(across(Shannon:`Unknown Ruminococcaceae (Family)`, rankZ)) %>%
    select(-Diet, -Sex)

## diet as interactive covariate
qtl_scan_list <- list()
for (trait in traits_16s) {
    cat ("On phenotype ", trait, "\n")
    
    phe <- phe_16s %>% select(one_of(trait)) %>% na.omit()
    qtl_scan <- scan1(probs, phe, kinship, addcovar = covar_add, intcovar = covar_int, cores = 0)
    
    qtl_scan_list[[trait]] <- qtl_scan
}
qtl_scan_all <- do.call("cbind", qtl_scan_list)
saveRDS(qtl_scan_all, file = "qtl_scan_16s_int_20240216.rds")

## diet as additive covariate
qtl_scan_list <- list()
for (trait in traits_16s) {
    cat ("On phenotype ", trait, "\n")
    
    phe <- phe_16s %>% select(one_of(trait)) %>% na.omit()
    qtl_scan <- scan1(probs, phe, kinship, addcovar = covar_add, cores = 0)
    
    qtl_scan_list[[trait]] <- qtl_scan
}
qtl_scan_all <- do.call("cbind", qtl_scan_list)
saveRDS(qtl_scan_all, file = "qtl_scan_16s_add_20240216.rds")

