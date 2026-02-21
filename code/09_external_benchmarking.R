library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(GSVA)          # gsva/ssgsea
library(GSEABase)      # GeneSet
library(broom)
library(effsize)       # cohen.d
library(rcompanion)    
library(ggplot2)
library(patchwork)

gene_sets <- list(
  Core21 = unique(core21_genes),
  SenSkin = unique(senskin_genes)
)
# ssGSEA
score_ssgsea <- function(pb_counts,meta,site_label,
                                  core_genes, senskin_genes,
                                  col_donor, col_age,min_core = 5, min_senskin = 10) {
  
  
  core_use  <- intersect(core_genes, rownames(pb_counts))
  sens_use  <- intersect(senskin_genes, rownames(pb_counts))
  
  if (length(core_use) < min_core) warning(site_label, ": Core21 overlap too small (", length(core_use), "). Check gene symbols/IDs.")
  if (length(sens_use) < min_senskin) warning(site_label, ": SenSkin overlap too small (", length(sens_use), "). Check gene symbols/IDs.")
  
  gene_sets <- list(Core21 = core_use, SenSkin = sens_use)
  
  expr_mat <- as.matrix(pb_counts)
  ssgsea <- GSVA::gsva(
    expr = expr_mat,
    gset.idx.list = gene_sets,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = TRUE,
    min.sz = 5,
    max.sz = 5000,
    verbose = FALSE
  )  # sets x donors
  ssgsea_z <- t(scale(t(ssgsea)))
  donor_meta <- meta %>%
    dplyr::select(all_of(c(col_donor, col_age))) %>%
    distinct() %>%
    mutate(donor = .data[[col_donor]]) %>%
    filter(donor %in% colnames(expr_mat)) %>%
    arrange(match(donor, colnames(expr_mat)))
  out <- donor_meta %>%
    transmute(
      donor = .data[[col_donor]],
      age   = as.numeric(.data[[col_age]]),
      site  = site_label,
      core21_ssgsea  = as.numeric(ssgsea["Core21", donor]),
      senskin_ssgsea = as.numeric(ssgsea["SenSkin", donor]),
      core21_ssgsea_z  = as.numeric(ssgsea_z["Core21", donor]),
      senskin_ssgsea_z = as.numeric(ssgsea_z["SenSkin", donor])
    )
  
  return(out)
}


# eyelid meta: sample_id, age_years, sample_group
# protected meta: subj, age_years, sample_group
meta_protected <- read_tsv('PS/Analysis/public_data/GSE130973/sample_info.txt')
meta_protected2 <- meta_protected %>%
  rename(sample_id = subj) %>%
  rename(sample_group = age) %>%
  mutate(sample_group = ifelse(str_to_upper(sample_group) %in% c("Y","YOUNG"), "Y", "O"))
meta_eyelid <- read_tsv('PS/Analysis/public_data/Zou_Liu_DevC_2021/Analysis/sample_info.txt')
meta_eyelid2 <- meta_eyelid

pb_eyelid <- read.table('PS/Analysis/public_data/Zou_Liu_DevC_2021/Analysis/1.pseudobulk/FB/DaMiRseq_adjust/FB_SV_adjusted_expression.txt',sep = '\t')
pb_protected <- read.table('PS/Analysis/public_data/GSE130973/1.pseudobulk/FB/DaMiRseq_adjust/FB_SV_adjusted_expression.txt',sep = '\t')

sc_eyelid <- score_ssgsea(pb_eyelid, core_genes = core21_genes,
                          senskin_genes = senskin_genes,
                          site_label = "eyelid",meta=meta_eyelid2,
                          col_donor = 'sample_id', col_age = 'age_years'
                          )   # signature x donor
sc_prot   <- score_ssgsea(pb_protected, core_genes = core21_genes,
                          senskin_genes = senskin_genes,
                          site_label = "protected",meta=meta_protected2,
                          col_donor = 'sample_id', col_age = 'age_years')
sc_eyelid <- sc_eyelid %>% mutate(sample_group = ifelse(age <= 24, "Y",
                                                        ifelse(age >= 70, "O","M")))
sc_prot <- sc_prot %>% mutate(sample_group = ifelse(age <= 30, "Y", "O"))

df <- bind_rows(sc_eyelid, sc_prot)
df$site <- factor(df$site)
df$site <- relevel(df$site, ref = "protected") 


df_all <- df %>%
  pivot_longer(
    cols = c(core21_ssgsea, senskin_ssgsea),
    names_to = "signature",
    names_prefix = "_ssgsea",  
    names_transform = list(signature = ~str_replace(., "_ssgsea", "")),  
    values_to = "score"
  ) %>%
  mutate(signature = str_to_title(signature))  
fit_A <- df_all %>%
  group_by(signature) %>%
  do(broom::tidy(lm(score ~ age + site + age:site, data = .))) %>%
  ungroup()

fit_A
summary(fit_A)
pA1 <- ggplot(df_all, aes(x = age, y = score, color = site)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ signature, scales = "free_y") +
  theme_bw() +
  labs(x = "Age (years)", y = "FB pseudo-bulk ssGSEA score",
       title = "Cross-site FB signature scores vs age (LOESS)")

pA1
ggsave('pA1.pdf',w=8,h=3.5)
fit_A_plot <- fit_A %>%
  filter(term %in% c("age", "siteeyelid", "age:siteeyelid")) %>%
  mutate(term = factor(term, levels = c("age","siteeyelid","age:siteeyelid")))

pA2 <- ggplot(fit_A_plot, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96*std.error, xmax = estimate + 1.96*std.error), height = 0.2) +
  facet_wrap(~ signature, scales = "free_x") +
  theme_bw() +
  labs(x = "Estimate (±1.96 SE)", y = NULL,
       title = "Linear model coefficients: score ~ age + site + age:site")

pA2
ggsave('pA2.pdf',w=7,h=2.5)
run_B <- function(df) {
  # df: one dataset + one signature
  d <- df %>% filter(sample_group %in% c("Y","O"))
  y <- d %>% filter(sample_group == "Y") %>% pull(score)
  o <- d %>% filter(sample_group == "O") %>% pull(score)
  
  tibble::tibble(
    nY = length(y),
    nO = length(o),
    wilcox_p = if (length(y) >= 2 && length(o) >= 2) wilcox.test(y, o)$p.value else NA_real_,
    cohens_d = if (length(y) >= 2 && length(o) >= 2) effsize::cohen.d(o, y, pooled = TRUE)$estimate else NA_real_,
    cliffs_delta = if (length(y) >= 2 && length(o) >= 2) cliff.delta(o, y)$estimate else NA_real_
  )
}

tab_B <- df_all %>%
  group_by(site, signature) %>%
  group_modify(~ run_B(.x)) %>%
  ungroup()

tab_B

pB1 <- ggplot(df_all %>% filter(sample_group %in% c("Y","O")),
              aes(x = sample_group, y = score, color = sample_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  facet_grid(signature ~ site, scales = "free_y") +
  theme_bw() +
  labs(x = "Age group", y = "FB pseudo-bulk ssGSEA score",
       title = "Young vs Old FB signature scores (external datasets)")

pB1
ggsave('pB1.pdf',w=6,h=4)

tab_B_long <- tab_B %>%
  pivot_longer(cols = c(cohens_d, cliffs_delta), names_to = "effect", values_to = "value")

pB2 <- ggplot(tab_B_long, aes(x = value, y = site)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  facet_grid(site ~ signature, scales = "free_x") +
  theme_bw() +
  labs(x = "Effect size (Old vs Young)", y = NULL,
       title = "Effect sizes for Young vs Old comparisons")

pB2
ggsave('pB2.pdf',w=5,h=3)



count_overlap <- function(pb_counts, gene_sets){
  sapply(gene_sets, function(gs) length(intersect(gs, rownames(pb_counts))))
}
count_overlap(pb_eyelid, gene_sets)
count_overlap(pb_protected, gene_sets)



summary(fit_core21)
summary(fit_senskin)
fit_core21  <- lm(core21_ssgsea_z  ~ age * site, data = df)
fit_senskin <- lm(senskin_ssgsea_z ~ age * site, data = df)
tab <- bind_rows(
  broom::tidy(fit_core21)  %>% mutate(signature = "Core21"),
  broom::tidy(fit_senskin) %>% mutate(signature = "SenSkin")
) %>% dplyr::select(signature, term, estimate, std.error, statistic, p.value)

tab

p1 <- ggplot(df, aes(age, core21_ssgsea_z, color = site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  labs(x = "Age", y = "FB pseudo-bulk Core21 ssGSEA (z)", color = "Site") +
  theme_classic()

p2 <- ggplot(df, aes(age, senskin_ssgsea_z, color = site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  labs(x = "Age", y = "FB pseudo-bulk SenSkin ssGSEA (z)", color = "Site") +
  theme_classic()

p1|p2
ggsave('test.pdf',h=3,w=10)



score_ssgsea <- function(pb_counts, gene_sets, min_gs_size = 5, kcdf = "Poisson") {
  gs2 <- lapply(gene_sets, function(gs) intersect(gs, rownames(pb_counts)))
  gs2 <- gs2[sapply(gs2, length) >= min_gs_size]
  
  gset <- GeneSetCollection(mapply(function(genes, nm) GeneSet(genes, setName = nm),
                                   gs2, names(gs2), SIMPLIFY = FALSE))
  
  scores <- GSVA::gsva(as.matrix(pb_counts), gset, method = "ssgsea", kcdf = kcdf,
                       mx.diff = TRUE, abs.ranking = FALSE, verbose = FALSE)
  # scores: signature x donor
  return(scores)
}

