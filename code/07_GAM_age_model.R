library(limma)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(extrafont)
pdf.options(family = "Arial")
setwd("0.all_celltype/whole_skin/DE_continous")
# expr_mat: rows=genes, cols=samples
# meta: data.frame(sample=colnames(expr_mat), age, sex)
expr_mat <- read.table('0.all_celltype/whole_skin/DaMiRseq_adjust/whole_skin_SV_adjusted_expression.txt',sep = '\t') 
meta <- read_tsv('sample_info.txt')
meta$stage <- factor(meta$stage,levels = c('Y','M','O'))
# design matrix
design <- model.matrix(~ age, data = meta)

fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)

res <- topTable(fit, coef="age", number=Inf, sort.by="none")
res$change <- ifelse(
  res$adj.P.Val < 0.05 & res$logFC >= 0.01, 'up',
  ifelse(res$adj.P.Val < 0.05 & res$logFC <= -0.01, 'down', "NOT")
)
res$change <- ifelse(
  res$P.Value < 0.05 & res$logFC >= 0.01, 'up',
  ifelse(res$P.Value < 0.05 & res$logFC <= -0.01, 'down', "NOT")
)

res <- res%>% arrange(desc(logFC))
write.table(res,file = 'whole_skin_DEG_continous.txt',sep = '\t',quote = F)
aa <- res[rownames(res)[!grepl("\\.", rownames(res))],]
write.table(aa,file = 'whole_skin_DEG_continous.txt',sep = '\t',quote = F)

top_up   <- res %>% filter(change=="up")   %>% slice_max(logFC, n=30)
top_down <- res %>% filter(change=="down") %>% slice_min(logFC, n=30)

top_genes <- c(rownames(top_up), rownames(top_down))

expr_sub <- expr_mat[top_genes, ] 
expr_scaled <- t(scale(t(expr_sub)))
age_order <- meta$sample_id[order(meta$age)]
expr_scaled <- expr_scaled[, age_order]

meta_ordered <- meta[match(colnames(expr_scaled), meta$sample_id), ]

stage_colors <- c("Y"= "#FD7446FF",         
                  "M"= "#709AE1FF",  
                  "O"= "#FED439FF" )
age_col_fun  <- colorRamp2(c(min(meta_ordered$age), max(meta_ordered$age)),
                           c("white", "firebrick"))
ha <- HeatmapAnnotation(
  Stage = meta_ordered$stage,
  Age   = meta_ordered$age,
  col = list(
    Stage = stage_colors,
    Age   = age_col_fun
  ),
  annotation_legend_param = list(
    Stage = list(title="Stage"),
    Age   = list(title="Age")
  )
)
p <- Heatmap(expr_scaled,
        name="zscore",
        column_split = meta_ordered$stage,
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 6))
pdf('heatmap_top30.pdf',h=7.5,w=6)
print(p)
dev.off()


df_gene <- bitr(rownames(top_up) , fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db )
gene=df_gene[,2]
egoBP <- enrichGO(gene = gene,
                  OrgDb= org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.2,qvalueCutoff = 0.05,keyType = 'ENTREZID')
egoBP <- setReadable(egoBP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
go.all=data.frame(egoBP)
write.table(go.all,"up30_gene_enrich-GO_BP.txt",sep="\t",quote=F,row.names=F)
pdf('up30_gene_go_BP_barplot.pdf',h=10,w=15)
barplot(egoBP, label_format=100,showCategory=30,title='up 30') 
dev.off()

df_gene <- bitr(rownames(top_down) , fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db )
gene=df_gene[,2]
egoBP <- enrichGO(gene = gene,
                  OrgDb= org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
egoBP <- setReadable(egoBP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
go.all=data.frame(egoBP)
write.table(go.all,"down30_gene_enrich-GO_BP.txt",sep="\t",quote=F,row.names=F)
pdf('down30_gene_go_BP_barplot.pdf',h=10,w=15)
barplot(egoBP, label_format=100,showCategory=30,title='down 30') 
dev.off()

top_up   <- res %>% filter(change=="up")   %>% slice_max(logFC, n=100)
top_down <- res %>% filter(change=="down") %>% slice_min(logFC, n=100)
top_genes <- c(rownames(top_up), rownames(top_down))



