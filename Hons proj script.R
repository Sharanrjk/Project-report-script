#COPD dataset DEG analysis########

setwd(" ") # Insert your working directory here
library("biomaRt")
library("dplyr")
library(ggplot2)
library(tibble)
library(tidyr)


EM = read.csv("HC_COPD_export_counts.csv")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host = "www.ensembl.org"))
genes <-  EM$ID
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), 
                  values = genes, mart= mart)

new_EM = left_join(EM, gene_IDs, by = c("ID"="ensembl_gene_id")) #there are some gene IDs that havent been converted; fix them
EM = new_EM

empty_genes = EM[EM$hgnc_symbol == '',] %>% rownames()
empty_genes = empty_genes %>% as.numeric()
EM$hgnc_symbol[empty_genes] = 'undefined_name' #find way to convert them into their respective gene IDs instead 

EM$hgnc_symbol[empty_genes] = EM$ID[empty_genes] #small problem here- somehow some gene symbols are NA so the empty_gene values also have NA. fixing them below

NA_vals = is.na(EM$hgnc_symbol) %>% which()
empty_genes = na.omit(empty_genes)
all_empty_genes = c(empty_genes,NA_vals)

#Now we can give all empty gene names their GENE IDs
EM$hgnc_symbol[all_empty_genes] = EM$ID[all_empty_genes]


#Trying the deseq now:

#trying to make rownames into genenames. below two code for simple gene IDs
# row.names(EM) = EM[,1]
# EM = EM[,-1]

row.names(EM)= make.names(EM[,"hgnc_symbol"], unique = TRUE)
EM = EM[,-c(1,17)] # getting rid of IDs and gene names
# Filtering data
any(is.na( EM )) #we dont have NAs here, but while converting gene ID we will

rm = rowMeans(EM, na.rm = TRUE)
length(rm)
length(rm[rm>0.01]) # all genes are expressed
expressed_genes = which(rm>0.01)
EM = EM[expressed_genes,]

#Covariance filter
rm = rowMeans(EM, na.rm = TRUE)
stdev = sapply(EM, sd, na.rm = TRUE)
cv = stdev/rm
length(cv[cv > 0.05])  #most genes make the cut so no need to apply it

#######
EM = readRDS("EM_COPD_processed")
#######

condition = c(rep("HC",5), rep("COPD",10)) %>% as.factor()
columns_metadata = data.frame(condition, row.names =  colnames(EM))

# creating DeSeq object
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData =  EM, #Make sure to specify correct input data & colData
                              colData = columns_metadata, 
                              design = ~  condition )



# Transform counts for data visualization & outlier detection (PCA & Heatmap)
rld <- rlog(dds, blind=TRUE)   #normalize and rlog transforms raw counts

# Plot PCA

DESeq2::plotPCA(rld, intgroup = "condition")

pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE) # these codes are if you want to modify the PCA plot w more conditions
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_text_repel(aes(label = name)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme(panel.background = element_blank(), panel.border = element_rect(fill=NA)) 

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap::pheatmap(rld_cor, annotation = columns_metadata[, c("condition"), drop=F])



# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
# Plot dispersion estimates
plotDispEsts(dds)

# Contrast
contrast_HC_COPD = c("condition","HC","COPD")
contrast_COPD_HC = c("condition","COPD","HC")  # just puts upregulated genes in positive log 

res <- results(dds, 
               contrast = contrast_COPD_HC,
               alpha = 0.05)


res <- lfcShrink(dds, 
                 contrast =  contrast_COPD_HC,
                 res=res,
                 type = "normal")


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl

# Table of results for significant genes
# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res


#Scatterplot of normalized expression of top 20 most significant genes

normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)


gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")



ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = samplename, 
                 size = 0.5), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes for selected contrast group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))

# Heatmap of all genes

# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)


library(RColorBrewer)
# Set a color palette
heat_colors <- brewer.pal(6, "GnBu")
# Run pheatmap using the metadata data frame for the annotation
a <- pheatmap::pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
                        color = heat_colors, # Find a cool color palette
                        cluster_rows = T, 
                        show_rownames = F,
                        annotation = columns_metadata[, c("condition","condition")], 
                        border_color = NA, 
                        fontsize = 10, 
                        scale = "row", 
                        fontsize_row = 10, 
                        height = 20,
                        main = "Heatmap of DEGs between HC and COPD samples")



res_table_thres <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.50)


## creating genelabels column for labelling
res_table_thres <- res_table_thres %>% arrange(padj) %>% mutate(genelabels = "")

# For specific genes, get position of them first then insert. 
interest_genes = (c("MMP9","CCR5") %>% toupper())

interest_pos = ( match(interest_genes,res_table_thres$gene) %>% na.omit() )

res_table_thres$genelabels[interest_pos] <- res_table_thres$gene[interest_pos]

library(ggrepel)
## Volcano plot
ggplot(res_table_thres, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold, alpha= 0.3, )) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +       #(this line makes plot look cleaner but labels dont work if used)
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) + theme(panel.border = element_rect(fill = NA, size = 2),
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(),
        ) + coord_cartesian(ylim = c(0,35)
        ) + scale_color_manual(values=c("black", "blue2")
        ) + geom_vline(xintercept = c(-0.5,0.5), col="red",linetype = "dashed" ) + geom_hline(yintercept = -log10(0.05), col="red", linetype = "dashed")


ggplot(res_table_thres, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold, alpha= 0.3)) +
  geom_label_repel(aes(label = genelabels), max.overlaps = 5000, box.padding = 0.7) + 
  ggtitle("Volcano plot of DEGs between RA samples and controls") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value")

# MA plot; requires 3 columns: mean, log fold change, logical table like threshold or significance

ma = res_table_thres[, c("baseMean", "log2FoldChange", "threshold")]
plotMA(ma,ylim=c(-5,5))

sig_res = sig_res_RA

# GO
library(gprofiler2)
gene_list = sig_res_RA$gene[which(sig_res_RA$log2FoldChange > 0)]
Gprofiler <- gost(a, organism = "hsapiens", ordered_query = FALSE,
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                  measure_underrepresentation = FALSE, evcodes = FALSE,
                  user_threshold = 0.05, correction_method = c("bonferroni"),#"g_SCS", "bonferroni","fdr", "false_discovery_rate", "gSCS", "analytical",
                  custom_bg = NULL,
                  numeric_ns = "", sources = NULL) #Gene set enrichment analysis

head(Gprofiler$result,10)
gostplot(Gprofiler, capped = FALSE, interactive = TRUE)
p <- gostplot(Gprofiler, capped = FALSE, interactive = FALSE) #Plot results
pp <- publish_gostplot(p, highlight_terms = c("GO:0030334","GO:0030155","GO:0010941","GO:0002694","GO:0042110","GO:0001503"), 
                       width = NA, height = NA, filename = NULL )


library(clusterProfiler)
library(enrichplot)

### Cluster profiler ###
library(BiocManager)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)

COPD_sig_genes_all_entrez = bitr(a, fromType = "SYMBOL", 
                                 toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

go_enrich_COPD = enrichGO(gene = COPD_sig_genes_all_entrez$ENTREZID, 
                          OrgDb = org.Hs.eg.db,readable = T, ont = "BP", 
                          pvalueCutoff = 0.05, qvalueCutoff = 0.10)

barplot(go_enrich_COPD, showCategory=10)
mutate(go_enrich_COPD, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

dotplot(go_enrich_COPD, showCategory = 15) 

emapplot(go_enrich_COPD)
cnetplot(go_enrich_COPD, categorySize="pvalue", foldChange = sig_res_RA$log2FoldChange)
heatplot(go_enrich_COPD)

p1 <- heatplot(go_enrich_COPD, showCategory=5)
p2 <- heatplot(go_enrich_COPD, foldChange=, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
##


### RA dataset DEG analysis
setwd(" ")
library("biomaRt")
library("dplyr")
library(ggplot2)
library(tibble)
library(tidyr)


EM = read.csv("counts_RA_exported.csv")


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host = "www.ensembl.org"))
genes <-  EM$ID
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), 
                  values = genes, mart= mart)

new_EM = left_join(EM, gene_IDs, by = c("ID"="ensembl_gene_id")) #there are some gene IDs that havent been converted; fix them
EM = new_EM

empty_genes = EM[EM$hgnc_symbol == '',] %>% rownames()
empty_genes = empty_genes %>% as.numeric()
EM$hgnc_symbol[empty_genes] = 'undefined_name' #find way to cnvert them into their respective gene IDs instead 

EM$hgnc_symbol[empty_genes] = EM$ID[empty_genes] #small problem here- somehow some gene symbols are NA so the empty_gene values also have NA. fixing them below

NA_vals = is.na(EM$hgnc_symbol) %>% which()
empty_genes = na.omit(empty_genes)
all_empty_genes = c(empty_genes,NA_vals)

#Now we can give all empty gene names their GENE IDs
EM$hgnc_symbol[all_empty_genes] = EM$ID[all_empty_genes]


#Trying the deseq now:

#trying to make rownames into genenames. below two code for simple gene IDs
# row.names(EM) = EM[,1]
# EM = EM[,-1]

row.names(EM)= make.names(EM[,"hgnc_symbol"], unique = TRUE)
EM = EM[,-c(1,16)] # getting rid of IDs and gene names
# Filtering data
any(is.na( EM )) #we dont have NAs here, but while converting gene ID we will

rm = rowMeans(EM, na.rm = TRUE)
length(rm)
length(rm[rm>0.01]) # all genes are expressed
expressed_genes = which(rm>0.01)
EM = EM[expressed_genes,]

#Covariance filter
rm = rowMeans(EM, na.rm = TRUE)
stdev = sapply(EM, sd, na.rm = TRUE)
cv = stdev/rm
length(cv[cv > 0.05])  #most genes make the cut so no need to apply it




condition = c(rep("HC",5), rep("RA",9)) %>% as.factor()
columns_metadata = data.frame(condition, row.names =  colnames(EM))

# creating DeSeq object
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData =  EM, #Make sure to specify correct input data & colData
                              colData = columns_metadata, 
                              design = ~  condition )



# Transform counts for data visualization & outlier detection (PCA & Heatmap)
rld <- rlog(dds, blind=TRUE)   #normalize and rlog transforms raw counts

# Plot PCA

DESeq2::plotPCA(rld, intgroup = "condition")

pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE) # these codes are if you want to modify the PCA plot w more conditions
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_text_repel(aes(label = name)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme(panel.background = element_blank(), panel.border = element_rect(fill=NA)) 

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap::pheatmap(rld_cor, annotation = columns_metadata[, c("condition"), drop=F])



# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
# Plot dispersion estimates
plotDispEsts(dds)

# Contrast
contrast_HC_RA = c("condition","HC","RA")
contrast_RA_HC = c("condition","RA","HC")

res <- results(dds, 
               contrast = contrast_RA_HC,
               alpha = 0.05)


res <- lfcShrink(dds, 
                 contrast =  contrast_RA_HC,
                 res=res,
                 type = "normal")


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl

# Table of results for significant genes
# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)   #No logfc cut-off here?

# Check significant genes output
sig_res


#Scatterplot of normalized expression of top 20 most significant genes

normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)


gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")



ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = samplename, 
                 size = 0.5), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes for selected contrast group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))

# Heatmap of all genes

# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)

library(RColorBrewer)
# Set a color palette
heat_colors <- brewer.pal(6, "GnBu")
# Run pheatmap using the metadata data frame for the annotation
a <- pheatmap::pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
                        color = heat_colors, # Find a cool color palette
                        cluster_rows = T, 
                        show_rownames = F,
                        annotation = columns_metadata[, c("condition","condition")], 
                        border_color = NA, 
                        fontsize = 10, 
                        scale = "row", 
                        fontsize_row = 10, 
                        height = 20,
                        main = "Heatmap of DEGs between selected contrast groups")


res_table_thres <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.50)


## creating genelabels column for labelling
res_table_thres <- res_table_thres %>% arrange(padj) %>% mutate(genelabels = "")

# For specific genes, get position of them first then insert. 
interest_genes = (c("IFI27","IL1A","CD80","CCRL2") %>% toupper())

interest_pos = ( match(interest_genes,res_table_thres$gene) %>% na.omit() )

res_table_thres$genelabels[interest_pos] <- res_table_thres$gene[interest_pos]

library(ggrepel)
## Volcano plot
ggplot(res_table_thres, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold, alpha= 0.3)) +
  geom_label_repel(aes(label = genelabels), max.overlaps = 200, box.padding = 0.5) + 
  ggtitle("Volcano plot of DEGs between RA samples and controls") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +       #(this line makes plot look cleaner but labels dont work if used)
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) + theme(panel.border = element_rect(fill = NA, size = 2),
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(),
                                                             panel.background = element_blank()
        ) + coord_cartesian(ylim = c(0,15)
        ) + scale_color_manual(values=c("black", "blue2")
        ) + geom_vline(xintercept = c(-0.5,0.5), col="red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col="red", linetype = "dashed")


# MA plot; requires 3 columns: mean, log fold change, logical table like threshold or significance

ma = res_table_thres[, c("baseMean", "log2FoldChange", "threshold")]
plotMA(ma,ylim=c(-5,5))

# GO
library(gprofiler2)
Gprofiler <- gost(sig_res$gene[which(sig_res$log2FoldChange > 0)], organism = "hsapiens", ordered_query = FALSE,
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                  measure_underrepresentation = FALSE, evcodes = FALSE,
                  user_threshold = 0.05, correction_method = c("bonferroni"),#"g_SCS", "bonferroni","fdr", "false_discovery_rate", "gSCS", "analytical",
                  custom_bg = NULL,
                  numeric_ns = "", sources = NULL) #Gene set enrichment analysis

head(Gprofiler$result,10)
gostplot(Gprofiler, capped = FALSE, interactive = TRUE)
p <- gostplot(Gprofiler, capped = FALSE, interactive = FALSE) #Plot results
p

# Enriched genes- do GO. Spatial enrichment???

#GSEA

library(clusterProfiler)
library(enrichplot)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# reading in data from deseq2
df = res_tbl

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

organism = org.Hs.eg.db
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")


require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

emapplot(gse, showCategory = 10)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)


gse_barplot <- enrichDGN(gene_list)
barplot(gse, showCategory = 10)

library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]
edo <- enrichDGN(de)
barplot(edo, showCategory=20) 

mutate(edo, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")






###### comparing RA and COPD ######################################

sig_res_RA = readRDS("sig_res_RA")
sig_res_COPD = readRDS("sig_res_COPD")
sig_res_LPS = readRDS("sig_res_LPS_1.rds") #in GSE data file

common_genes = intersect(sig_res_COPD$gene, sig_res_RA$gene)
common_genes_activated_monocytes = intersect(common_genes, sig_res_LPS$gene)
common_genes_autoimmune_sig = setdiff(common_genes, common_genes_activated_monocytes) # not a good way to do it, removes upreulation and downreg factors, look below it is done better

# Gonna do a log fold change comparison for both of them so subsetting only those first, then merging them

sig_res_RA = sig_res_RA[,c(1,3)]
sig_res_COPD = sig_res_COPD[,c(1,3)]
sig_res_LPS = sig_res_LPS[c(1,3)]

colnames( sig_res_RA ) = c('genes','Log2FC_RA')
colnames( sig_res_COPD ) = c('genes','Log2FC_COPD')
colnames(sig_res_LPS) = c('genes','Log2FC_LPS')

merged_df = merge(sig_res_RA,sig_res_COPD, by = 'genes')
merged_df_all_activatedmonocytes = merge(merged_df, sig_res_LPS, by = 'genes')

# to see highest DEGs
merged_df_thres <- merged_df %>% 
  mutate(abs(Log2FC_RA) >= 3 & abs(Log2FC_COPD) >= 3)

colnames(merged_df_thres) = c(colnames(merged_df)[1:3], "threshold")
merged_df_thres <- merged_df_thres %>% arrange(threshold) %>% mutate(genelabels = "")
merged_df_thres$genelabels[2743:2746] <- merged_df_thres$genes[2743:2746]


require(stats)
reg<- lm(merged_df$Log2FC_COPD ~ merged_df$Log2FC_RA)
reg

ggplot(merged_df_thres, aes(x = Log2FC_RA , y = Log2FC_COPD)) + geom_point() + geom_label_repel(aes(label = genelabels), max.overlaps = 200, box.padding = 0.5)+ #plotting LOG2FC plot
  theme_classic() + 
  geom_abline(intercept = 0.0002, slope = 1.4, color = 'dark red') + xlab(' Log fold change RA')  + ylab('Log fold change COPD') +
  geom_point(color = "blue", alpha = 0.3) +geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme(panel.border = element_rect(fill = NA, size = 2),
                                                                                                           panel.grid.major = element_blank(),
                                                                                                           panel.grid.minor = element_blank(),
                                                                                                           panel.background = element_blank())



plot(merged_df$Log2FC_RA, merged_df$Log2FC_COPD)  + abline(coef = c(0,0) ) + abline(coef = c(0,1)) #+ abline(lm(merged_df$Log2FC_RA ~ merged_df$Log2FC_COPD), col = 'blue')



# spearman correlation test
cor.test(merged_df$Log2FC_COPD, merged_df$Log2FC_RA, method = "spearman")

# For venn diagram of common genes

COPD_upreg = merged_df$genes[which(merged_df$Log2FC_COPD > 0)]
RA_upreg = merged_df$genes[which(merged_df$Log2FC_RA > 0)]

COPD_downreg = merged_df$genes[which(merged_df$Log2FC_COPD < 0)]
RA_downreg = merged_df$genes[which(merged_df$Log2FC_RA < 0)]

# Venn diagram for all genes
COPD_upreg1 = sig_res_COPD$genes[which(sig_res_COPD$Log2FC_COPD > 0)]
RA_upreg1 = sig_res_RA$genes[which(sig_res_RA$Log2FC_RA > 0)]
LPS_upreg1 = sig_res_LPS$genes[which(sig_res_LPS$Log2FC_LPS > 0)]

COPD_downreg1 = sig_res_COPD$genes[which(sig_res_COPD$Log2FC_COPD < 0)]
RA_downreg1 = sig_res_RA$genes[which(sig_res_RA$Log2FC_RA < 0)]
LPS_downreg1 = sig_res_LPS$genes[which(sig_res_LPS$Log2FC_LPS < 0)]

#common genes autoimmune signature
a = intersect(COPD_upreg1, RA_upreg1)
b = setdiff(a, LPS_upreg1)
c = intersect(COPD_downreg1, RA_downreg1)
d = setdiff(c, LPS_downreg1)
common_genes_autoimmune_sig = c(b,d)

# alternatively can use a package:
library(ggVennDiagram)
x <- list(COPD = brown_module, RA = midnight_blue)
y <- list(COPD = COPD_downreg1, RA = RA_downreg1, LPS = LPS_downreg1)
ggVennDiagram(x, label = "count", color = "black" ,lwd = 0.7, edge_lty = 1)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

x <- list(RA_module = brown_high, "COPD_module" = midnightblue_high, "RA-COPD signature" = common_genes_autoimmune_sig)
ggVennDiagram(x, label = "count", label_size = 6)


#### COPD coexpnet analysis
setwd(" ")
library(dplyr)
library(DESeq2)
library(tidyr)
library(WGCNA)
library(CoExpNets)
library("biomaRt")

EM = read.csv("HC_COPD_export_counts.csv")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host = "www.ensembl.org"))
genes <-  EM$ID
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), 
                  values = genes, mart= mart)

EM = left_join(EM, gene_IDs, by = c("ID"="ensembl_gene_id")) #there are some gene IDs that havent been converted; fix them

empty_genes = EM[EM$hgnc_symbol == '',] %>% rownames()
empty_genes = empty_genes %>% as.numeric()
EM$hgnc_symbol[empty_genes] = 'undefined_name' #find way to cnvert them into their respective gene IDs instead 

EM$hgnc_symbol[empty_genes] = EM$ID[empty_genes] #small problem here- somehow some gene symbols are NA so the empty_gene values also have NA. fixing them below

NA_vals = is.na(EM$hgnc_symbol) %>% which()
empty_genes = na.omit(empty_genes)
all_empty_genes = c(empty_genes,NA_vals)

#Now we can give all empty gene names their GENE IDs
EM$hgnc_symbol[all_empty_genes] = EM$ID[all_empty_genes]


row.names(EM)= make.names(EM[,"hgnc_symbol"], unique = TRUE)
EM = EM[,-c(1,17)] # getting rid of IDs and gene names

EM = t(EM)
EM = EM + 1
EM = log2(EM)
EM = EM %>% as.data.frame()


# filtering expressed genes (mean>1.5)
cm = colMeans(EM, na.rm = TRUE)
length(cm[cm>1.5])
expressed_genes = which(cm>1.5)

EM = EM[,expressed_genes]

# selecting genes with CV > 5%
cm = colMeans(EM, na.rm = TRUE)
stdev = sapply(EM, sd, na.rm = TRUE)
cv = stdev/cm
head(cv)
head(cv[cv>0.05])
length(cv[cv > 0.05])

good_gene_pos = which(cv>0.05)
EM = EM[ ,good_gene_pos]

gsg_EM = goodSamplesGenes(EM, verbose = 3);
gsg_EM$allOK


#Clustering samples and removing outliers
sampleTree = hclust(dist(EM), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#Clustering samples and removing outliers
sampleTree = hclust(dist(EM), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# remove COPD_10 here? aslo HC_5 and COPD_5 look fishy?

net = CoExpNets::getDownstreamNetwork(tissue="COPD_data",
                                      n.iterations=20,
                                      net.type = "signed",
                                      debug=F, fullAnnotation = FALSE,
                                      expr.data= EM,
                                      job.path = "results_coexpnet/")


#annotatting
net.name="netCOPD_data.29.it.20.rds"  #Note: to be able to run getdownstream we had to increase filtering (>3 mean) and decrease gene number
net = readRDS(net.name)
tissue = paste("results_coexpnet/", net.name, sep = "")

names(net$moduleColors) =  gsub("\\.[0-9]+","",names(net$moduleColors))
saveRDS(net,net.name)



#clustering of module eigengenes (which modules have the most similar expression)
which.one = "new"
egs = getNetworkEigengenes(which.one=which.one,tissue=tissue)
cat(tissue,"network was obtained from",nrow(egs),
    "samples and has",ncol(egs),"modules\n")

#And now we plot the EGs
plotEGClustering(which.one=which.one,tissue=tissue)
library(corrplot)
corrplot(cor(egs),tl.srt=45, 
         tl.col="black",
         type="upper",
         tl.cex=0.45,
         order="hclust",
         title=paste0("Eigengene correlations for ",tissue," within ",which.one),
         mar=c(0,3,3.5,2))

#GO

go = CoExpNets::getGProfilerOnNet(net.file=net.name,
                                  exclude.iea=F,
                                  out.file=paste0(net.name,"_gprof.csv"))


# Looking at which module is most enriched
e = read.csv("netCOPD_Data.29.it.20.rds_gprof.csv",stringsAsFactors=F)
sort(table(e$query.number),decreasing=T)

# looking further into module genes
midnightblue = e[e$query.number == "midnightblue",]
midnightblue[order(midnightblue$p.value)[1:10],c("term.name","domain","p.value"),]

blue = e[e$query.number == "blue" & e$domain == "keg",]
blue[order(blue$p.value)[1:10],c("term.name","domain","p.value"),]

yellow = e[e$query.number == "yellow" & e$domain == "BP",]
yellow[order(yellow$p.value)[1:10],c("term.name","domain","p.value"),]

brown = e[e$query.number == "brown" & e$domain == "BP",]
brown[order(brown$p.value)[1:10],c("term.name","domain","p.value"),]

#Look at genes in each module
(net[["moduleColors"]] == "midnightblue") %>% which() %>% length()
(net[["moduleColors"]] == "midnightblue") %>% which() %>% names()



#Get MM of all genes so you can see which module any gene is in
Allmms = CoExpNets::getMM(which.one="new",tissue=tissue,
                          genes= colnames(EM),expr.data.file = EM) 


# Correlating with categorical traits

condition = c(rep("HC",5), rep("COPD",10)) %>% as.factor()
traitData =  data.frame(condition, row.names = row.names(EM))

CoExpNets::corWithCatTraits(which.one="new",tissue= tissue,covs=traitData, covlist = colnames(traitData))
#think we have to rework the trait data file? maybe put HC in one column and COPD in the other with 1's and 0's?


#Tomplot now

library(tidyr)
expr.data = EM

mms_TOM <- function(module, tissue, Beta){
  genes = CoExpNets::getGenesFromModule(which.one="new",
                                        tissue=tissue,
                                        module=module) #get genes from module
  mms = CoExpNets::getMM(which.one="new",tissue=tissue,
                         genes= genes,expr.data.file = expr.data) #Get MM of genes
  mmsO <- mms[order(mms$mm, decreasing = TRUE),] #order by MM
  moduleExpr <- expr.data[,colnames(expr.data) %in% genes]
  createTOM = function(expr.data.file,
                       beta=Beta,
                       save.as=NULL,
                       net.type="signed",
                       debug=F){
    
    stopifnot(beta > 0 & beta < 40)
    stopifnot(net.type == "signed" | net.type == "unsigned")
    
    if(typeof(expr.data.file) == "character"){
      print(paste0("Creating matrix ",save.as," from expression data ",expr.data.file))
      expr.data <- readRDS(expr.data.file)
    }else{
      expr.data <- expr.data.file
    }
    cat("Creating TOM for",ncol(expr.data),"genes and",nrow(expr.data),"samples, beta",
        beta,"and type",net.type,"\n")
    if(debug)
      expr.data = expr.data[,1:1000]
    adjacency = adjacency(expr.data, power = beta, type = net.type )
    print("Adjacency matrix created")
    # Topological Overlap Matrix (TOM)
    # Turn adjacency into topological overlap
    print("Creating TOM")
    TOM = TOMsimilarity(adjacency)
    colnames(TOM) = colnames(expr.data)
    rownames(TOM) = colnames(TOM)
    
    if(!is.null(save.as)){
      cat("Saving TOM at",save.as,"\n")
      saveRDS(TOM,save.as)
    }
    else return(TOM)
  }
  TOM <- as.data.frame(CoExpNets::createTOM(expr.data.file = moduleExpr, beta = Beta))
  TOM$from <- rownames(TOM)
  TOM <- gather(TOM, from)
  genes <- genes[genes %in% TOM[,1]]
  TOM$to <- rep(genes)
  TOM$Zero <- rep(0)
  TOM$MOO39 <- rep("M0039")
  TOM <- TOM[,c(1, 3, 4, 5, 2)]
  MMS_tom <- list(mmsO, TOM)
  return(MMS_tom)
} #remember to change beta for different networks

Network = "midnightblue" #which module you want to investigate
Beta = 29 #the beta value chosen by co exp nets (number before ".it" in net name, must be correct to give MM)

Module <- mms_TOM(Network, paste("results_coexpnet/", net.name, sep = ""), Beta) #get MM of all genes in the module + list of TOMs'

# exporting MM and TOM text file exported
write.csv(Module[[1]], file = paste0(Network,"Module Membership", net.name,".csv")) #save MM
Module[[2]] = Module[[2]][Module[[2]]$value > 0.3,]  # filtering by adjacency

write.table(Module[[2]], file= paste0(Network,"Module TOMs 3plus", net.name, ".txt"),sep = , quote=FALSE, row.names=FALSE, col.names = FALSE) #save TOMs


#just extracting names of genes with mm > 0.9

a = Module[[1]]
midnightblue_high = a$name[which(a$mm > 0.95)]
brown_high = a$name[which(a$mm > 0.95)]

### RA coexpnets

setwd(" ")
library(dplyr)
library(DESeq2)
library(tidyr)
library(WGCNA)
library(CoExpNets)
library("biomaRt")

EM = read.csv("HC_COPD_export_counts.csv")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host = "www.ensembl.org"))
genes <-  EM$ID
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), 
                  values = genes, mart= mart)

EM = left_join(EM, gene_IDs, by = c("ID"="ensembl_gene_id")) #there are some gene IDs that havent been converted; fix them

empty_genes = EM[EM$hgnc_symbol == '',] %>% rownames()
empty_genes = empty_genes %>% as.numeric()
EM$hgnc_symbol[empty_genes] = 'undefined_name' #find way to cnvert them into their respective gene IDs instead 

EM$hgnc_symbol[empty_genes] = EM$ID[empty_genes] #small problem here- somehow some gene symbols are NA so the empty_gene values also have NA. fixing them below

NA_vals = is.na(EM$hgnc_symbol) %>% which()
empty_genes = na.omit(empty_genes)
all_empty_genes = c(empty_genes,NA_vals)

#Now we can give all empty gene names their GENE IDs
EM$hgnc_symbol[all_empty_genes] = EM$ID[all_empty_genes]


row.names(EM)= make.names(EM[,"hgnc_symbol"], unique = TRUE)
EM = EM[,-c(1,17)] # getting rid of IDs and gene names

EM = t(EM)
EM = EM + 1
EM = log2(EM)
EM = EM %>% as.data.frame()


# filtering expressed genes (mean>1.5)
cm = colMeans(EM, na.rm = TRUE)
length(cm[cm>1.5])
expressed_genes = which(cm>1.5)

EM = EM[,expressed_genes]

# selecting genes with CV > 5%
cm = colMeans(EM, na.rm = TRUE)
stdev = sapply(EM, sd, na.rm = TRUE)
cv = stdev/cm
head(cv)
head(cv[cv>0.05])
length(cv[cv > 0.05])

good_gene_pos = which(cv>0.05)
EM = EM[ ,good_gene_pos]

gsg_EM = goodSamplesGenes(EM, verbose = 3);
gsg_EM$allOK


#Clustering samples and removing outliers
sampleTree = hclust(dist(EM), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#Clustering samples and removing outliers
sampleTree = hclust(dist(EM), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# remove COPD_10 here? aslo HC_5 and COPD_5 look fishy?

net = CoExpNets::getDownstreamNetwork(tissue="COPD_data",
                                      n.iterations=20,
                                      net.type = "signed",
                                      debug=F, fullAnnotation = FALSE,
                                      expr.data= EM,
                                      job.path = "results_coexpnet/")


#annotatting
net.name="netCOPD_data.29.it.20.rds"  #Note: to be able to run getdownstream we had to increase filtering (>3 mean) and decrease gene number
net = readRDS(net.name)
tissue = paste("results_coexpnet/", net.name, sep = "")

names(net$moduleColors) =  gsub("\\.[0-9]+","",names(net$moduleColors))
saveRDS(net,net.name)



#clustering of module eigengenes (which modules have the most similar expression)
which.one = "new"
egs = getNetworkEigengenes(which.one=which.one,tissue=tissue)
cat(tissue,"network was obtained from",nrow(egs),
    "samples and has",ncol(egs),"modules\n")

#And now we plot the EGs
plotEGClustering(which.one=which.one,tissue=tissue)
library(corrplot)
corrplot(cor(egs),tl.srt=45, 
         tl.col="black",
         type="upper",
         tl.cex=0.45,
         order="hclust",
         title=paste0("Eigengene correlations for ",tissue," within ",which.one),
         mar=c(0,3,3.5,2))

#GO

go = CoExpNets::getGProfilerOnNet(net.file=net.name,
                                  exclude.iea=F,
                                  out.file=paste0(net.name,"_gprof.csv"))


# Looking at which module is most enriched
e = read.csv("netCOPD_Data.29.it.20.rds_gprof.csv",stringsAsFactors=F)
sort(table(e$query.number),decreasing=T)

# looking further into module genes
midnightblue = e[e$query.number == "midnightblue",]
midnightblue[order(midnightblue$p.value)[1:10],c("term.name","domain","p.value"),]

blue = e[e$query.number == "blue" & e$domain == "keg",]
blue[order(blue$p.value)[1:10],c("term.name","domain","p.value"),]

yellow = e[e$query.number == "yellow" & e$domain == "BP",]
yellow[order(yellow$p.value)[1:10],c("term.name","domain","p.value"),]

brown = e[e$query.number == "brown" & e$domain == "BP",]
brown[order(brown$p.value)[1:10],c("term.name","domain","p.value"),]

#Look at genes in each module
(net[["moduleColors"]] == "midnightblue") %>% which() %>% length()
(net[["moduleColors"]] == "midnightblue") %>% which() %>% names()



#Get MM of all genes so you can see which module any gene is in
Allmms = CoExpNets::getMM(which.one="new",tissue=tissue,
                          genes= colnames(EM),expr.data.file = EM) 


# Correlating with categorical traits

condition = c(rep("HC",5), rep("COPD",10)) %>% as.factor()
traitData =  data.frame(condition, row.names = row.names(EM))

CoExpNets::corWithCatTraits(which.one="new",tissue= tissue,covs=traitData, covlist = colnames(traitData))
#think we have to rework the trait data file? maybe put HC in one column and COPD in the other with 1's and 0's?


#Tomplot now

library(tidyr)
expr.data = EM

mms_TOM <- function(module, tissue, Beta){
  genes = CoExpNets::getGenesFromModule(which.one="new",
                                        tissue=tissue,
                                        module=module) #get genes from module
  mms = CoExpNets::getMM(which.one="new",tissue=tissue,
                         genes= genes,expr.data.file = expr.data) #Get MM of genes
  mmsO <- mms[order(mms$mm, decreasing = TRUE),] #order by MM
  moduleExpr <- expr.data[,colnames(expr.data) %in% genes]
  createTOM = function(expr.data.file,
                       beta=Beta,
                       save.as=NULL,
                       net.type="signed",
                       debug=F){
    
    stopifnot(beta > 0 & beta < 40)
    stopifnot(net.type == "signed" | net.type == "unsigned")
    
    if(typeof(expr.data.file) == "character"){
      print(paste0("Creating matrix ",save.as," from expression data ",expr.data.file))
      expr.data <- readRDS(expr.data.file)
    }else{
      expr.data <- expr.data.file
    }
    cat("Creating TOM for",ncol(expr.data),"genes and",nrow(expr.data),"samples, beta",
        beta,"and type",net.type,"\n")
    if(debug)
      expr.data = expr.data[,1:1000]
    adjacency = adjacency(expr.data, power = beta, type = net.type )
    print("Adjacency matrix created")
    # Topological Overlap Matrix (TOM)
    # Turn adjacency into topological overlap
    print("Creating TOM")
    TOM = TOMsimilarity(adjacency)
    colnames(TOM) = colnames(expr.data)
    rownames(TOM) = colnames(TOM)
    
    if(!is.null(save.as)){
      cat("Saving TOM at",save.as,"\n")
      saveRDS(TOM,save.as)
    }
    else return(TOM)
  }
  TOM <- as.data.frame(CoExpNets::createTOM(expr.data.file = moduleExpr, beta = Beta))
  TOM$from <- rownames(TOM)
  TOM <- gather(TOM, from)
  genes <- genes[genes %in% TOM[,1]]
  TOM$to <- rep(genes)
  TOM$Zero <- rep(0)
  TOM$MOO39 <- rep("M0039")
  TOM <- TOM[,c(1, 3, 4, 5, 2)]
  MMS_tom <- list(mmsO, TOM)
  return(MMS_tom)
} #remember to change beta for different networks

Network = "midnightblue" #which module you want to investigate
Beta = 29 #the beta value chosen by co exp nets (number before ".it" in net name, must be correct to give MM)

Module <- mms_TOM(Network, paste("results_coexpnet/", net.name, sep = ""), Beta) #get MM of all genes in the module + list of TOMs'

# exporting MM and TOM text file exported
write.csv(Module[[1]], file = paste0(Network,"Module Membership", net.name,".csv")) #save MM
Module[[2]] = Module[[2]][Module[[2]]$value > 0.3,]  # filtering by adjacency

write.table(Module[[2]], file= paste0(Network,"Module TOMs 3plus", net.name, ".txt"),sep = , quote=FALSE, row.names=FALSE, col.names = FALSE) #save TOMs


#just extracting names of genes with mm > 0.9

a = Module[[1]]
midnightblue_high = a$name[which(a$mm > 0.95)]
brown_high = a$name[which(a$mm > 0.95)]














