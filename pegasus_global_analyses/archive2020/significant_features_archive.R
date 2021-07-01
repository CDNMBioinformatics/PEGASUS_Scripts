
# 3. HILIC positive correlation structure
## 3.1. Annotated features with p<0.05

```{r, eval = F}

load(hilic_pos_dat_fname, verbose = T)

sig_kn_hilic_pos <- res_hilic_pos[pval < alpha_thld & tentative_annotation != "Unknown" & 
                                          subgroup == "overall", feature]
med_annot <- res_hilic_pos[pval < alpha_thld & tentative_annotation != "Unknown" & 
                                   subgroup == "overall", .(feature, med_feat)] %>% as.data.frame()
rownames(med_annot) <- med_annot$feature
med_annot$feature <- NULL
med_annot$med_feat <- as.factor(med_annot$med_feat)
hmap_kn_hilic_pos <- cor(hilic_pos_dat[, ..sig_kn_hilic_pos]) %>% 
        pheatmap(annotation_col = med_annot, cutree_cols = 3, cutree_rows = 3, 
                 main = "Correlation among HILIC-pos annotated features with p<0.05")

save_pheatmap_png(hmap_kn_hilic_pos, here(fig_dir, "corr_sig_kn_hilic_pos.png"))


## Long format

cor_long_kn_hilic_pos <- cor(hilic_pos_dat[, ..sig_kn_hilic_pos]) %>% 
        cor_gather() %>% as.data.table() %>% setnames(c("feature", "feature2", "cor"))
cor_long_kn_hilic_pos <- merge(cor_long_kn_hilic_pos, 
                               hilic_pos_mets_info[, .(feature, tentative_annotation)], by = "feature")
setcolorder(cor_long_kn_hilic_pos, c("feature", "tentative_annotation"))
datatable(cor_long_kn_hilic_pos[feature != feature2][order(-cor)], filter = "top") %>% 
        formatSignif(columns = c("cor"), digits = 3)

```

## 3.2. Features in high correlation

```{r, eval = F}

## Cluster 1

datatable(res_hilic_pos[feature %in% c("F1500", "F1120", "F162", "F1426", "F1323", "F467", "F375"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F381", "F475"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F1258", "F461", "F365"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F299", "F526", "F784"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F225", "F544"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)


## Cluster 2

datatable(res_hilic_pos[feature %in% c("F522", "F137", "F2750", "F1439", "F574", "F195", "F796", "F1467", "F1461"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)


## Cluster 3

datatable(res_hilic_pos[feature %in% c("F1537", "F1663", "F1355", "F1279", "F1827", "F990", "F1444", "F1489"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F423", "F339", "F340"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F745", "F756"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F1466", "F883"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F1301", "F1372"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F504", "F311"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F437", "F351", "F443"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_hilic_pos[feature %in% c("F1307", "F1501", "F1511"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

```

## 3.3. Annotated features vs unknown features with pval<0.05

```{r, eval = F}

sig_unkn_hilic_pos <- res_overall_hilic_pos[pval < alpha_thld & tentative_annotation == "Unknown", feature]
hmap_kn_vs_unkn_hilic_pos <- cor(hilic_pos_dat[, ..sig_kn_hilic_pos], hilic_pos_dat[, ..sig_unkn_hilic_pos]) %>% 
        pheatmap(annotation_row = med_annot, cutree_cols = 3, cutree_rows = 3, 
                 main = "Correlation among HILIC-pos annotated vs unknown features with p<0.05")

save_pheatmap_png(hmap_kn_vs_unkn_hilic_pos, here(fig_dir, "corr_sig_kn_vs_unkn_hilic_pos.png"), width = 60)


## Long format

cor_long_kn_vs_unkn_hilic_pos <- cor(hilic_pos_dat[, ..sig_kn_hilic_pos], hilic_pos_dat[, ..sig_unkn_hilic_pos]) %>% 
        cor_gather() %>% as.data.table() %>% setnames(c("feature", "unknown_feature", "cor"))
cor_long_kn_vs_unkn_hilic_pos <- merge(cor_long_kn_vs_unkn_hilic_pos, 
                                       hilic_pos_mets_info[, .(feature, tentative_annotation)], by = "feature")
setcolorder(cor_long_kn_vs_unkn_hilic_pos, c("feature", "tentative_annotation"))
datatable(cor_long_kn_vs_unkn_hilic_pos[abs(cor) >= 0.5][order(-cor)], filter = "top") %>% 
        formatSignif(columns = c("cor"), digits = 3)

fwrite(cor_long_kn_vs_unkn_hilic_pos[abs(cor) >= 0.5][order(-cor)], 
       file = here(res_dir, "hilic_pos_sig_unknown_feature_high_cor_w_annotated_feature.csv"))
fwrite(res_overall_hilic_pos[pval < alpha_thld & tentative_annotation == "Unknown", ], 
       file = here(res_dir, "hilic_pos_sig_unknown_feature.csv"))
fwrite(hilic_pos_mets_info, file = here(dat_dir, "hilic_pos_feature_info_qc.csv"))

```

# 4. Lipidomics positive correlation structure
## 4.1. Annotated features with p<0.05

```{r, eval = F}

load(lipid_pos_dat_fname, verbose = T)

sig_kn_lipid_pos <- res_lipid_pos[pval < alpha_thld & tentative_annotation != "Unknown" & 
                                          subgroup == "overall", feature]
hmap_kn_lipid_pos <- cor(lipid_pos_dat[, ..sig_kn_lipid_pos]) %>% 
        pheatmap(cutree_cols = 4, cutree_rows = 4, 
                 main = "Correlation among lipid-pos annotated \n features with p<0.05")

save_pheatmap_png(hmap_kn_lipid_pos, here(fig_dir, "corr_sig_kn_lipid_pos.png"), width = 6, height = 6)


## Long format

cor_long_kn_lipid_pos <- cor(lipid_pos_dat[, ..sig_kn_lipid_pos]) %>% 
        cor_gather() %>% as.data.table() %>% setnames(c("feature", "feature2", "cor"))
cor_long_kn_lipid_pos <- merge(cor_long_kn_lipid_pos, 
                               lipid_pos_mets_info[, .(feature, tentative_annotation)], by = "feature")
setcolorder(cor_long_kn_lipid_pos, c("feature", "tentative_annotation"))
datatable(cor_long_kn_lipid_pos[feature != feature2][order(-cor)], filter = "top") %>% 
        formatSignif(columns = c("cor"), digits = 3)

```

## 4.2. Features in high correlation

```{r, eval = F}

## Cluster 1

datatable(res_lipid_pos[feature %in% c("F2720", "F2723"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

```

## 4.3. Annotated features vs unknown features with pval<0.05

```{r, eval = F}

sig_unkn_lipid_pos <- res_overall_lipid_pos[pval < alpha_thld & tentative_annotation == "Unknown", feature]

hmap_kn_vs_unkn_lipid_pos <- cor(lipid_pos_dat[, ..sig_kn_lipid_pos], lipid_pos_dat[, ..sig_unkn_lipid_pos]) %>% 
        pheatmap(cutree_cols = 6, cutree_rows = 4, 
                 main = "Correlation among lipid-pos annotated vs unknown features with p<0.05")

save_pheatmap_png(hmap_kn_vs_unkn_lipid_pos, here(fig_dir, "corr_sig_kn_vs_unkn_lipid_pos.png"), 
                  width = 20, height = 6)


## Long format

cor_long_kn_vs_unkn_lipid_pos <- cor(lipid_pos_dat[, ..sig_kn_lipid_pos], lipid_pos_dat[, ..sig_unkn_lipid_pos]) %>% 
        cor_gather() %>% as.data.table() %>% setnames(c("feature", "unknown_feature", "cor"))
cor_long_kn_vs_unkn_lipid_pos <- merge(cor_long_kn_vs_unkn_lipid_pos, 
                                       lipid_pos_mets_info[, .(feature, tentative_annotation)], by = "feature")
setcolorder(cor_long_kn_vs_unkn_lipid_pos, c("feature", "tentative_annotation"))
datatable(cor_long_kn_vs_unkn_lipid_pos[abs(cor) >= 0.5][order(-cor)], filter = "top") %>% 
        formatSignif(columns = c("cor"), digits = 3)

fwrite(cor_long_kn_vs_unkn_lipid_pos[abs(cor) >= 0.5][order(-cor)], 
       file = here(res_dir, "lipid_pos_sig_unknown_feature_high_cor_w_annotated_feature.csv"))
fwrite(res_overall_lipid_pos[pval < alpha_thld & tentative_annotation == "Unknown", ], 
       file = here(res_dir, "lipid_pos_sig_unknown_feature.csv"))
fwrite(lipid_pos_mets_info, file = here(dat_dir, "lipid_pos_feature_info_qc.csv"))

```

# 5. Lipidomics negative correlation structure
## 5.1. Annotated features with p<0.05

```{r, eval = F}

load(lipid_neg_dat_fname, verbose = T)

sig_kn_lipid_neg <- res_lipid_neg[pval < alpha_thld & tentative_annotation != "Unknown" & 
                                          subgroup == "overall", feature]
hmap_kn_lipid_neg <- cor(lipid_neg_dat[, ..sig_kn_lipid_neg]) %>% 
        pheatmap(cutree_cols = 3, cutree_rows = 3, 
                 main = "Correlation among lipid-neg annotated \n features with p<0.05")

save_pheatmap_png(hmap_kn_lipid_neg, here(fig_dir, "corr_sig_kn_lipid_neg.png"), width = 8, height = 8)


## Long format

cor_long_kn_lipid_neg <- cor(lipid_neg_dat[, ..sig_kn_lipid_neg]) %>% 
        cor_gather() %>% as.data.table() %>% setnames(c("feature", "feature2", "cor"))
cor_long_kn_lipid_neg <- merge(cor_long_kn_lipid_neg, 
                               lipid_neg_mets_info[, .(feature, tentative_annotation)], by = "feature")
setcolorder(cor_long_kn_lipid_neg, c("feature", "tentative_annotation"))
datatable(cor_long_kn_lipid_neg[feature != feature2][order(-cor)], filter = "top") %>% 
        formatSignif(columns = c("cor"), digits = 3)

```

## 5.2. Features in high correlation

```{r, eval = F}

## Cluster 1

datatable(res_lipid_neg[feature %in% c("F2573", "F2494"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

## Cluster 2

datatable(res_lipid_neg[feature %in% c("F265", "F370"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_lipid_neg[feature %in% c("F1995", "F2533"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

datatable(res_lipid_neg[feature %in% c("F389", "F747"), ]) %>% 
        formatSignif(columns = num_cols, digits = 3)

```

## 5.3. Annotated features vs unknown features with pval<0.05

```{r, eval = F}

sig_unkn_lipid_neg <- res_overall_lipid_neg[pval < alpha_thld & tentative_annotation == "Unknown", feature]

hmap_kn_vs_unkn_lipid_neg <- cor(lipid_neg_dat[, ..sig_kn_lipid_neg], lipid_neg_dat[, ..sig_unkn_lipid_neg]) %>% 
        pheatmap(cutree_cols = 3, cutree_rows = 3, 
                 main = "Correlation among lipid-neg annotated vs unknown features with p<0.05")

save_pheatmap_png(hmap_kn_vs_unkn_lipid_neg, here(fig_dir, "corr_sig_kn_vs_unkn_lipid_neg.png"), 
                  width = 12, height = 8)


## Long format

cor_long_kn_vs_unkn_lipid_neg <- cor(lipid_neg_dat[, ..sig_kn_lipid_neg], lipid_neg_dat[, ..sig_unkn_lipid_neg]) %>% 
        cor_gather() %>% as.data.table() %>% setnames(c("feature", "unknown_feature", "cor"))
cor_long_kn_vs_unkn_lipid_neg <- merge(cor_long_kn_vs_unkn_lipid_neg, 
                                       lipid_neg_mets_info[, .(feature, tentative_annotation)], by = "feature")
setcolorder(cor_long_kn_vs_unkn_lipid_neg, c("feature", "tentative_annotation"))
datatable(cor_long_kn_vs_unkn_lipid_neg[abs(cor) >= 0.5][order(-cor)], filter = "top") %>% 
        formatSignif(columns = c("cor"), digits = 3)

fwrite(cor_long_kn_vs_unkn_lipid_neg[abs(cor) >= 0.5][order(-cor)], 
       file = here(res_dir, "lipid_neg_sig_unknown_feature_high_cor_w_annotated_feature.csv"))
fwrite(res_overall_lipid_neg[pval < alpha_thld & tentative_annotation == "Unknown", ], 
       file = here(res_dir, "lipid_neg_sig_unknown_feature.csv"))
fwrite(lipid_neg_mets_info, file = here(dat_dir, "lipid_neg_feature_info_qc.csv"))

```
