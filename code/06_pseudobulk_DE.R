library(DaMiRseq)
library(limma)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(data.table)
library(ggsci) 
library(extrafont)
pdf.options(family = "Arial")
library(bioplotr)

# Function to read raw count data
read_raw_count <- function(path) {
  read.delim(path, sep = ',')
}

# Paths to raw count data
paths <- list(
  whole_skin = "whole_skin_pseudobulk_rawCount.csv",
  KC = "KC_pseudobulk_rawCount.csv",
  KC_Channel = "KC_Channel_pseudobulk_rawCount.csv",
  FB = "FB_pseudobulk_rawCount.csv",
  VEC = "VEC_pseudobulk_rawCount.csv",
  LEC = "LEC_pseudobulk_rawCount.csv",
  Pc = "Pc_pseudobulk_rawCount.csv",
  MEL = "MEL_pseudobulk_rawCount.csv",
  Schwann = "Schwann_pseudobulk_rawCount.csv",
  Lymphocyte = "Lymphocyte_pseudobulk_rawCount.csv",
  Mac_DC = "Mac_DC_pseudobulk_rawCount.csv",
  LC = "LC_pseudobulk_rawCount.csv",
  Mast = "Mast_pseudobulk_rawCount.csv",
  SGC = "SGC_pseudobulk_rawCount.csv"
)

# Read all raw count data
celltypes <- lapply(paths, read_raw_count)

sampleinfo <- read_tsv("sample_info.txt")

run_damirseq <- function(countdata, sampleinfo, group, celltype, minCounts = 10, fSample = 0.5, th.corr = 0.8, n.pred = 50) {
  if (!dir.exists(celltype)) {dir.create(celltype)}
  setwd(celltype)
  sampleinfo <- sampleinfo %>% filter(sample_group == group)
  countdata <- countdata %>% column_to_rownames("X") %>% as.data.frame()
  countdata <- countdata %>% dplyr::select(intersect(sampleinfo$sample_id, colnames(countdata)))
  keep <- rowSums(countdata) > 100
  countdata <- countdata[keep, ]

  class1 <- sampleinfo[sampleinfo$sample_id %in% colnames(countdata), ]
  class3 <- class1 %>% column_to_rownames("sample_id") %>% dplyr::rename(class = stage)
  SE <- DaMiR.makeSE(countdata, class3)

  data_norm <- DaMiR.normalization(SE, minCounts = minCounts, fSample = fSample, hyper = "no")
  if (!dir.exists('DaMiRseq_filtered')) {dir.create('DaMiRseq_filtered')}
  if (!dir.exists('DaMiRseq_adjust')) {dir.create('DaMiRseq_adjust')}
  write.table(assay(data_norm), file=paste0("DaMiRseq_filtered/", celltype, "_vst_normalized_expression.txt"), row.names = TRUE, col.names = TRUE, sep = '\t')
  data_filt <- DaMiR.sampleFilt(data_norm, th.corr = th.corr)
  sv <- DaMiR.SV(data_filt)
  data_filt$class <- as.factor(data_filt$class)
  data_adjust <- DaMiR.SVadjust(data_filt, sv)
  write.table(assay(data_adjust), file=paste0("DaMiRseq_adjust/", celltype, "_SV_adjusted_expression.txt"), row.names = TRUE, col.names = TRUE, sep = '\t', quote = FALSE)

  plot_and_save <- function(data, prefix) {
    pdf(paste0(prefix, '_all_plot.pdf'), w=10, h=10)
    DaMiR.Allplot(data, colData(data))
    dev.off()

    pdf(paste0(prefix, '_umap.pdf'), w=5, h=4)
    plot_umap(assay(data), group=data$class, size=3, label = TRUE, pal_group='d3')
    plot_umap(assay(data), group=data$class, size=3, pal_group='d3')
    dev.off()

    pdf(paste0(prefix, '_pca.pdf'), w=5, h=4)
    plot_pca(assay(data), group=data$class, size=3, label = TRUE, pal_group='d3')
    plot_pca(assay(data), group=data$class, size=3, pal_group='d3')
    plot_pca(assay(data), group=data$class, size=3, pcs = c(3L, 4L), label = TRUE, pal_group='d3')
    dev.off()
    
    pdf(paste0(prefix, '_pca_only.pdf'), w=4.5, h=3)
    plot_pca(assay(data), group=data$class, size=4, pal_group='d3')
    dev.off()
  }

  plot_and_save(data_filt, paste0('DaMiRseq_filtered/', celltype))
  plot_and_save(data_adjust, paste0('DaMiRseq_adjust/', celltype))

  data_clean <- DaMiR.transpose(assay(data_adjust))
  df <- colData(data_adjust)
  data_reduced <- DaMiR.FSelect(data_clean, df, th.corr = 0.4)
  write.table(data_reduced$data, file=paste0("DaMiRseq_adjust/", celltype, "_reduced_expression.txt"), row.names = TRUE, col.names = TRUE, sep = '\t', quote = FALSE)

  data_reduced <- DaMiR.FReduct(data_reduced$data)
  write.table(data_reduced, file=paste0("DaMiRseq_adjust/", celltype, "_reduced_expression_rmdup.txt"), row.names = TRUE, col.names = TRUE, sep = '\t', quote = FALSE)

  pdf(paste0('DaMiRseq_adjust/', celltype, '_reduced_mds.pdf'), w=10, h=8)
  DaMiR.MDSplot(data_reduced, df)
  dev.off()

  tt <- t(data_reduced)
  pdf(paste0('DaMiRseq_adjust/', celltype, '_data_reduced_umap.pdf'), w=5, h=4)
  plot_umap(tt, group=data_adjust$class, size=4, pal_group='d3')
  plot_umap(tt, group=data_adjust$class, size=4, label=TRUE, pal_group='d3')
  dev.off()

  pdf(paste0('DaMiRseq_adjust/', celltype, '_data_reduced_pca.pdf'), w=5, h=4)
  plot_pca(tt, group=data_adjust$class, size=3, pal_group='d3')
  plot_pca(tt, group=data_adjust$class, size=3, label=TRUE, pal_group='d3')
  plot_pca(tt, group=data_adjust$class, size=3, pcs=c(3L, 4L), label=TRUE, pal_group='d3')
  dev.off()

  df.importance <- DaMiR.FSort(data_reduced, df)
  write.table(df.importance, file=paste0("DaMiRseq_adjust/", celltype, "_reduced_expression_rmdup_rank.txt"), row.names=TRUE, col.names=TRUE, sep='\t', quote=FALSE)
  selected_features4 <- DaMiR.FBest(data_reduced, ranking=df.importance, autoselect='yes')
  selected_features5 <- DaMiR.FBest(data_reduced, ranking=df.importance, n.pred=length(df.importance$RReliefF))

  pdf(paste0('DaMiRseq_adjust/', celltype, '_heatmap.pdf'), h=15, w=11)
  DaMiR.Clustplot(selected_features4$data, df)
  dev.off()
}

setwd("/home/yangyu/Project/skin_age/PS/Analysis/0.all_celltype")
group = 'PS'
original_wd <- getwd()


cell_types <- c("whole_skin","KC","FB","VEC","LEC","Pc","MEL","Schwann","Lymphocyte","Mac_DC","LC","Mast"#,"KC_Channel","SGC"
               )
for (celltype in names(cell_types)) {
  print(celltype)
  run_damirseq(celltypes[[celltype]], sampleinfo, group, celltype, th.corr = 0.6)
  setwd(original_wd)
}

# limma
one_vs_one_union <- function(expr,
                             meta_data,
                             stage_col = "stage",
                             sample_col = "sample_id",
                             q_value_threshold = 0.05,
                             logFC_threshold = log2(1.5),
                             outdir = ".",
                             file_prefix = "",
                             min_n_per_group = 2) {  
  stopifnot(stage_col %in% colnames(meta_data), sample_col %in% colnames(meta_data))
  common_samples <- intersect(colnames(expr), meta_data[[sample_col]])
  if (length(common_samples) < 2) stop("可用于分析的样本不足（与表达矩阵交集 < 2）。")
  expr <- expr[, common_samples, drop = FALSE]
  meta_data_new <- meta_data[match(common_samples, meta_data[[sample_col]]), , drop = FALSE]
  
  group_labels <- droplevels(factor(meta_data_new[[stage_col]]))
  groups <- levels(group_labels)
  if (length(groups) < 2) {
    message("[skip] 该细胞类型仅有一个分组：", paste(groups, collapse = ", "))
    return(setNames(vector("list", length(groups)), groups))
  }
  
  tab_all <- table(group_labels)
  message("[sample counts] ", paste(names(tab_all), tab_all, collapse = " | "))
  
  genes_up_by_group <- setNames(vector("list", length(groups)), groups)
  for (g in groups) genes_up_by_group[[g]] <- character(0)
  
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      group1 <- groups[i]; group2 <- groups[j]
      
      subset_idx <- which(group_labels %in% c(group1, group2))
      labels_subset <- droplevels(factor(group_labels[subset_idx], levels = c(group1, group2)))
      
      if (nlevels(labels_subset) < 2) {
        message(sprintf("[skip] %s vs %s: 某组在该细胞类型中无样本", group1, group2))
        next
      }
      tab <- table(labels_subset)
      message(sprintf("[contrast] %s_vs_%s | %s", group1, group2, paste(names(tab), tab, collapse = ", ")))
      if (any(tab < min_n_per_group)) {
        message(sprintf("[skip] %s_vs_%s: 每组至少需要 >=%d 个样本", group1, group2, min_n_per_group))
        next
      }
      
      expr_subset <- expr[, subset_idx, drop = FALSE]
      
      design <- model.matrix(~ 0 + labels_subset)
      colnames(design) <- levels(labels_subset)
      contrast <- makeContrasts(contrasts = paste0(group1, "-", group2), levels = design)
      
      fit <- lmFit(expr_subset, design)
      fit2 <- tryCatch(contrasts.fit(fit, contrast),
                       error = function(e) { message("[skip contrast] ", group1, "_vs_", group2, " | ", e$message); return(NULL) })
      if (is.null(fit2)) next
      
      efit <- eBayes(fit2)
      df <- topTable(efit, number = Inf, adjust.method = "BH")
      df$gene <- rownames(df)
      
      df$change <- ifelse(
        df$adj.P.Val < q_value_threshold & df$logFC >= logFC_threshold, group1,
        ifelse(df$adj.P.Val < q_value_threshold & df$logFC <= -logFC_threshold, group2, "NOT")
      )
      df$change <- factor(df$change, levels = c("NOT", group2, group1))
      
      de_fname <- file.path(outdir, sprintf("%s%s_vs_%s_DE.txt", file_prefix, group1, group2))
      write.table(df, de_fname, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
      
      sig_df <- df %>% dplyr::filter(adj.P.Val < q_value_threshold)
      up_top10 <- sig_df %>% dplyr::filter(logFC > 0) %>% dplyr::slice_max(order_by = logFC, n = 10, with_ties = FALSE)
      down_top10 <- sig_df %>% dplyr::filter(logFC < 0) %>% dplyr::slice_min(order_by = logFC, n = 10, with_ties = FALSE)
      lab_genes <- unique(c(up_top10$gene, down_top10$gene))
      df$sign <- ifelse(df$gene %in% lab_genes, df$gene, NA_character_)
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(values = c("gray", "blue", "red"), limits = c("NOT", group2, group1)) +
        ggplot2::geom_hline(yintercept = -log10(q_value_threshold), linetype = 4) +
        ggplot2::geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = 4) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank()) +
        ggrepel::geom_label_repel(ggplot2::aes(label = sign), fontface = "bold",
                                  box.padding = grid::unit(1, "line"),
                                  point.padding = grid::unit(0.5, "line"),
                                  segment.color = "grey50", max.overlaps = 200) +
        ggplot2::labs(x = "log2FC", y = "-log10(adj.P.Val)",
                      title = sprintf("%s%s | %s vs %s", ifelse(file_prefix=="","",paste0(file_prefix," ")), "Volcano", group1, group2))
      grDevices::pdf(file.path(outdir, sprintf("%s%s_%s_volcano.pdf", file_prefix, group1, group2)), h = 5, w = 7)
      print(p); grDevices::dev.off()
      
      up_in_g1 <- df %>% dplyr::filter(adj.P.Val <= q_value_threshold, logFC >=  logFC_threshold) %>% dplyr::pull(gene)
      up_in_g2 <- df %>% dplyr::filter(adj.P.Val <= q_value_threshold, logFC <= -logFC_threshold) %>% dplyr::pull(gene)
      genes_up_by_group[[group1]] <- union(genes_up_by_group[[group1]], up_in_g1)
      genes_up_by_group[[group2]] <- union(genes_up_by_group[[group2]], up_in_g2)
    }
  }
  
  return(genes_up_by_group)
}


perform_DE_analysis <- function(file_path,
                                celltype,
                                meta_data,
                                stage_col = "stage",
                                sample_col = "sample",
                                q_value_threshold = 0.05,
                                logFC_threshold = log2(1.5),
                                outdir = ".") {
  expr <- read.table(file_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  message(sprintf("Loaded expr for %s: %d genes x %d samples", celltype, nrow(expr), ncol(expr)))
  
  cell_outdir <- file.path(outdir, celltype)
  dir.create(cell_outdir, showWarnings = FALSE, recursive = TRUE)
  
  res <- one_vs_one_union(expr,
                          meta_data = meta_data,
                          stage_col = stage_col,
                          sample_col = sample_col,
                          q_value_threshold = q_value_threshold,
                          logFC_threshold = logFC_threshold,
                          outdir = cell_outdir,
                          file_prefix = "")
  for (g in names(res)) {
    outf <- file.path(cell_outdir, sprintf("%s_genes.txt", g))
    write.table(res[[g]], outf, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  invisible(res)
}


