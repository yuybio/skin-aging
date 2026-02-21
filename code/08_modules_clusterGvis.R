library(dplyr)
library(ClusterGVis)
library(VennDiagram)
library(UpSetR)
setwd("DEG/FB")

O_df <- read.table('O_genes.txt')
O_genes <- O_df$V1
M_df <- read.table('M_genes.txt')
M_genes <- M_df$V1
Y_df <- read.table('Y_genes.txt')
Y_genes <- Y_df$V1

setwd('geneCluster')
data <- read.table('0.all_celltype/FB/DaMiRseq_adjust/FB_SV_adjusted_expression.txt',sep = '\t')  #135
dim(data)
sample_info <- read_tsv('sample_info.txt')

all_de_gene <- unique(c(O_genes,M_genes,Y_genes))
all_de_gene <- all_de_gene[!grepl("\\.", all_de_gene)]
de_df <- data[all_de_gene,] 
dim(de_df) #3082
order_indices <- match(colnames(de_df), sample_info$sample_id)
sample_info <- sample_info[order_indices, ]

df = as.data.frame(t(de_df)) %>% rownames_to_column('sample_id')
df <- df %>%
  left_join(sample_info,by = 'sample_id') %>% 
  column_to_rownames("sample_id")%>%
  dplyr::select(-sample_group)%>%
  dplyr::select(-age)%>%
  dplyr::select(-gender)
df_ave <- df %>%
  group_by(stage) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("stage") %>%
  as.data.frame()
df_ave <- as.data.frame(t(df_ave))
exps = df_ave
save(exps,file = 'exps_DEG.rda')
getClusters(exp = exps)

ck <- clusterData(exp = exps,
                  cluster.method = "kmeans",
                  cluster.num = 11)
save(ck,file = 'ck_DEG.rda')
pdf('cluster_DEG_11.pdf',height = 11,width = 6)
visCluster(object = ck,
           plot.type = "heatmap",
           column_names_rot = 45, add.box = T,
           sample.order = c('Y','M','O'),
           cluster.order = c(3,4,1,2,6,5,11,10,7,8,9),
           ctAnno.col = ggsci::pal_d3("category20")(11))
dev.off()


library(org.Hs.eg.db)
enrich <- enrichCluster(object = ck,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 1111)
save(enrich,file = 'enrich.rda')
pdf('termenrich.pdf',height = 15,width = 11,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,add.box = T,
           markGenes = c(FB_Y_specific_gene_aged,FB_O_specific_gene_aged),
           markGenes.side = "left",
           genes.gp = c('italic',fontsize = 5,col = "black"),
           annoTerm.data = enrich,
           line.side = "left",
           sample.order = c('Y','M','O'),
           cluster.order = c(3,4,1,2,6,5,11,10,7,8,9),
           ctAnno.col = ggsci::pal_d3("category20")(11),
           go.col = rep(ggsci::pal_d3("category20")(11),each = 5),
           go.size = "pval")
dev.off()



