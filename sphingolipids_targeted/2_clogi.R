#########################################################
##### Conditional Logistic Models                   #####
#########################################################


##### 1. Setup #####

## Packages

pkg <- c("tidyverse", "data.table", "here", "survival", "gtools")
for (p in pkg) {
        if (require(p, character.only = T)) {
                print(paste0(p, " loaded successfully"))
        } else {
                install.packages(p)
                require(p, character.only = T)
                print(paste0(p, " downloaded and loaded successfully"))
        }
}


## Paths

dat_dir <- "data/"
res_dir <- "results/"


## Filenames

dat_fname <- here(dat_dir, "eda_updated.RData")


## Hard-coded numbers

missing_thld <- 0.10
sig_digits <- 3
alpha_thld <- 0.05
z_95 <- qnorm(p = 0.975)



##### 2. Load data #####

load(dat_fname, verbose = T)


## Standardize the metabolites

mets_pheno[, (sph_list) := lapply(.SD, function(x) scale(x, center = T, scale = T)), .SDcols = sph_list]

unmatched_stratum <- mets_pheno[is.na(age), unique(stratum)]
mets_pheno_trim <- mets_pheno[stratum != unmatched_stratum, ]



##### 3. Conditional logistic models #####

clogi_res_fnc <- function(outc, mets_list, dat, add_covar = "", 
                          first_cols = c("metabolite", "beta", "pval", "fdr_bh", "or", "lower95", "upper95", 
                                         "totaln", "casen", "pct_na_met", "cv_qcinjs", "cv_qcexts", "skewness_met"), 
                          na_thld = missing_thld) {
        
        clogi_mdl_tmp <- function(met) {
                
                fit <- eval(parse(text = str_c("clogit(", outc, " ~ ", met, " + strata(stratum)", 
                                               add_covar, ", data = ", dat, ")")))
                
                coef <- summary(fit)$coefficients
                coef <- coef[rownames(coef) == met, ]
                
                confint <- summary(fit)$conf.int
                confint <- confint[rownames(confint) == met, -c(1:2)]
                
                totaln <- summary(fit)$n
                casen <- summary(fit)$nevent
                
                output <- t(c(outc, met, coef, confint, totaln, casen)) %>% as.data.table()
                colnames(output) <- c("outc", "metabolite", "beta", "or", "se", "zval", "pval", "lower95", "upper95", "totaln", "casen")
                output[, (colnames(output)[-c(1:2)]) := lapply(.SD, as.numeric), .SDcols = colnames(output)[-c(1:2)]]
                
                return(output)
        }
        
        res <- lapply(mets_list, function(x) {clogi_mdl_tmp(met = x)}) %>% rbindlist()
        res <- merge(res, 
                     mets_info[, .(metabolite, pct_na_met, cv_qcinjs, cv_qcexts, skewness_met)], 
                     by = "metabolite")
        res[pct_na_met <= na_thld, fdr_bh := p.adjust(pval, method = "BH")] # Only calculate FDR for metabolites with <10% missing
        
        setcolorder(res, first_cols)
        res <- res[order(pval)]
        return(res)
}


## 3.1. Crude: only accounting for matched strata

length(sph_list) # 78

res_crd <- clogi_res_fnc(outc = "asthma_sev", mets_list = sph_list, dat = "mets_pheno")
res_crd[1:10, ] %>% print(digits = sig_digits) # Lowest p = 1.55e-06
res_crd[pval < alpha_thld, ] %>% nrow() # 28
res_crd[fdr_bh < alpha_thld, ] %>% nrow() # 21


## Crude in subjects with BMI (not adjusted for BMI)

mets_pheno[!is.na(BMI), .N, stratum][, .N, N]

res_crd_wbmi <- clogi_res_fnc(outc = "asthma_sev", mets_list = sph_list, dat = "mets_pheno[!is.na(BMI), ]")
res_crd_wbmi[1:10, ] %>% print(digits = sig_digits) # Lowest p = 0.00103
res_crd_wbmi[pval < alpha_thld, ] %>% nrow() # 14


## 3.2. Adjusting for BMI, will reduce sample size to 293 pairs

res_bmi <- clogi_res_fnc(outc = "asthma_sev", mets_list = sph_list, dat = "mets_pheno", 
                         add_covar = " + BMI")
res_bmi[1:10, ] %>% print(digits = sig_digits) # Lowest p = 0.0105
res_bmi[pval < alpha_thld, ] %>% nrow() # 9
res_bmi[fdr_bh < alpha_thld, ] %>% nrow() # 0


## 3.2. Adjusting for sample redness & milkiness

res_red_mlk <- clogi_res_fnc(outc = "asthma_sev", mets_list = sph_list, dat = "mets_pheno", 
                             add_covar = " + redness + milkiness")
res_red_mlk[1:10, ] %>% print(digits = sig_digits) # Lowest p = 2.62e-06
res_red_mlk[pval < alpha_thld, ] %>% nrow() # 28
res_red_mlk[fdr_bh < alpha_thld, ] %>% nrow() # 17


## Exclude samples red/milky

res_normal <- clogi_res_fnc(outc = "asthma_sev", mets_list = sph_list, 
                          dat = "mets_pheno[redness == 'normal' & milkiness == 'normal', ]")
res_normal[1:10, ] %>% print(digits = sig_digits) # Lowest p = 2.62e-06
res_normal[pval < alpha_thld, ] %>% nrow() # 29
res_normal[fdr_bh < alpha_thld, ] %>% nrow() # 16


## 3.3. Adjusting for other covariates

### food_allergy_diagnoses

# ggplot(mets_pheno_trim, aes(x = Group, y = food_allergy_diagnoses)) + 
#         geom_boxplot(colour = "grey50") + 
#         geom_jitter(aes(color = Group), alpha = 0.5)
res_fa <- clogi_res_fnc(outc = "asthma_sev", mets_list = sph_list, dat = "mets_pheno_trim", 
                        add_covar = " + food_allergy_diagnoses")
res_fa[1:10, ] %>% print(digits = sig_digits) # Lowest p = 2.15e-06


### unspecified_allergy_diagnoses

# ggplot(mets_pheno_trim, aes(x = Group, y = unspecified_allergy_diagnoses)) + 
#         geom_boxplot(colour = "grey50") + 
#         geom_jitter(aes(color = Group), alpha = 0.5)
res_alg <- clogi_res_fnc(outc = "asthma_sev", mets_list = sph_list, dat = "mets_pheno_trim", 
                         add_covar = " + unspecified_allergy_diagnoses")
res_alg[1:10, ] %>% print(digits = sig_digits) # Lowest p = 1.80e-06


### allergic_rhinitis

mets_pheno_trim[, .N, .(Group, allergic_rhinitis)][order(Group, allergic_rhinitis)]
#        Group allergic_rhinitis   N
# 1: Asthmatic             FALSE 233
# 2: Asthmatic              TRUE 336
# 3:   Control             FALSE 513
# 4:   Control              TRUE  56
res_arh <- clogi_res_fnc(outc = "asthma_sev", mets_list = sph_list, dat = "mets_pheno_trim", 
                         add_covar = " + allergic_rhinitis")
res_arh[1:10, ] %>% print(digits = sig_digits) # Lowest p = 5.39e-05



##### 4. number of ICS/OCS and sphingolipids #####

linear_res_fnc <- function(expo, mets_list, dat, add_covar = "", 
                           first_cols = c("metabolite", "beta", "pval", "fdr_bh", "totaln", 
                                          "pct_na_met", "cv_qcinjs", "cv_qcexts", "skewness_met"), 
                           na_thld = missing_thld) {
        
        linear_mdl_tmp <- function(met) {
                
                fit <- eval(parse(text = str_c("lm(", met, " ~ ", expo, " + Age + Gender + Race", 
                                               add_covar, ", data = ", dat, ")")))
                
                coef <- summary(fit)$coefficients
                coef <- coef[rownames(coef) == expo, ]

                totaln <- summary(fit)$fstatistic[2] + summary(fit)$fstatistic[3] + 2
                
                output <- t(c(expo, met, coef, totaln)) %>% as.data.table()
                colnames(output) <- c("expo", "metabolite", "beta", "se", "tval", "pval", "totaln")
                output[, (colnames(output)[-c(1:2)]) := lapply(.SD, as.numeric), .SDcols = colnames(output)[-c(1:2)]]
                
                return(output)
        }
        
        res <- lapply(mets_list, function(x) {linear_mdl_tmp(met = x)}) %>% rbindlist()
        res <- merge(res, 
                     mets_info[, .(metabolite, pct_na_met, cv_qcinjs, cv_qcexts, skewness_met)], 
                     by = "metabolite")
        res[pct_na_met <= na_thld, fdr_bh := p.adjust(pval, method = "BH")] # Only calculate FDR for metabolites with <10% missing
        
        setcolorder(res, first_cols)
        res <- res[order(pval)]
        return(res)
}

mets_pheno_asth <- mets_pheno[Group == 'Asthmatic', ]

ggplot(mets_pheno_trim[Group == "Asthmatic", ]) + 
        geom_histogram(aes(inhaled_steriod_prescriptions))
res_mets_ics <- linear_res_fnc(expo = "inhaled_steriod_prescriptions", mets_list = sph_list, 
                               dat = "mets_pheno_asth")
res_mets_ics_bmi <- linear_res_fnc(expo = "inhaled_steriod_prescriptions", mets_list = sph_list, 
                                   dat = "mets_pheno_asth", add_covar = " + BMI")

ggplot(mets_pheno_trim[Group == "Asthmatic", ]) + 
        geom_histogram(aes(oral_steroid_prescriptions))
res_mets_ocs <- linear_res_fnc(expo = "oral_steroid_prescriptions", mets_list = sph_list, 
                               dat = "mets_pheno_asth")
res_mets_ocs_bmi <- linear_res_fnc(expo = "oral_steroid_prescriptions", mets_list = sph_list, 
                                   dat = "mets_pheno_asth", add_covar = " + BMI")

ggplot(mets_pheno_trim[Group == "Asthmatic", ]) + 
        geom_histogram(aes(bronchodialator_prescriptions))
res_mets_bd <- linear_res_fnc(expo = "bronchodialator_prescriptions", mets_list = sph_list, 
                              dat = "mets_pheno_asth")
res_mets_bd_bmi <- linear_res_fnc(expo = "bronchodialator_prescriptions", mets_list = sph_list, 
                                  dat = "mets_pheno_asth", add_covar = " + BMI")


##### 5. Save results #####

save(list = c(grep("^res_", ls(), value = T), grep("^mets_", ls(), value = T), grep("_list$", ls(), value = T)), 
     file = here(res_dir, "res_clogi.RData"))

sink("scripts/2_clogi_sessinfo.txt")
sessionInfo()
sink()
