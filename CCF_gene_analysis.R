
### Get the globally lowly expressed genes
#####################################
###

### Arrange the global expression matrix
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
names(object[[]])
library(edgeR)
pb_obj=Seurat2PB(object, sample="age", cluster="orig.ident")
pseudobulkdata=read.delim("~/Project_PBMCage/try_for_dtw.txt", stringsAsFactors=F)
keep.genes=rownames(pseudobulkdata)
pb_obj=pb_obj[keep.genes, , keep=FALSE]
write.table(pb_obj$counts, "~/Project_PBMCage/pseudobulk_global.txt", sep="\t")

### Plot the distribution of transcriptome in PBMC from samples at various ages to get an idea of "low expression"
df=read.delim("~/Project_PBMCage/pseudobulk_global.txt")
df=t(df)
ages=as.factor(gsub("X|_cluster.*","",rownames(df)))
rownames(df)=ages

df_long=reshape2::melt(df, variable.name="Gene", value.name="Expression")
colnames(df_long)=c("Age","Gene","Expression")
df_long$Age=as.factor(df_long$Age)
p_thrdeter=
  ggplot(df_long, aes(x=Expression, color=Age)) +
  geom_density(lwd=1, fill="white", alpha=0.25) +
  scale_x_log10(labels=scales::label_log(digits=1)) +
  theme_classic() +
  theme(legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"))
# according to the plot, determine the thr of low expression
lowExpr_thr=1.5 # log10(1.5)=0.1760913

### Determine the lowly expressed genes globally or in each celltype (may due to cell specificity)
celltype.wide_lowlyExprGenes=colnames(df)[apply(df, 2, mean)<lowExpr_thr]

pseudobulkdata=read.delim("~/Project_PBMCage/try_for_dtw.txt", stringsAsFactors=F)
Celltype=gsub(".*_cluster","",colnames(pseudobulkdata))
Age=as.integer(gsub("X|_.*","",colnames(pseudobulkdata)))
pseudobulkdata_arranged=as.data.frame(t(pseudobulkdata))
pseudobulkdata_arranged=cbind(Celltype, Age, pseudobulkdata_arranged)

all_celltypes=names(table(Celltype))
lowlyExprGenes=list()
for (i in 1:length(all_celltypes)) {
  df_=subset(pseudobulkdata_arranged, Celltype==all_celltypes[i])
  ages=df_[,"Age"]
  df_=df_[,c(3:ncol(df_))]
  rownames(df_)=ages
  
  lowlyExprGenes[[i]]=colnames(df_)[apply(df_, 2, mean)<lowExpr_thr]
}
names(lowlyExprGenes)=all_celltypes
saveRDS(lowlyExprGenes, "~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_perCelltype.rds")
saveRDS(celltype.wide_lowlyExprGenes, "~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_global.rds")

### Enrich the globally lowly expressed genes
gse=enrichGO(celltype.wide_lowlyExprGenes, OrgDb="org.Hs.eg.db", keyType="SYMBOL")

p_gse=
  dotplot(gse, showCategory=8, font.size=10) +
  labs(title="Globally lowly expressed genes") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10))

### Plot the expression of globally lowly expressed genes
genes_for_plotting=c(celltype.wide_lowlyExprGenes[1:10], "HPRT1", celltype.wide_lowlyExprGenes[11:length(celltype.wide_lowlyExprGenes)])
expression_matrix_chosen=pseudobulkdata_arranged[,c("Celltype",genes_for_plotting)]
long_exp.matrix.chosen=reshape2::melt(expression_matrix_chosen, id.vars=c("Celltype"), variable.name="Gene", value.name="Expression")
p_expr=
  ggplot(long_exp.matrix.chosen, aes(x=Gene, y=Expression, color=Celltype)) +
  geom_point(shape=1) +
  scale_color_manual(values=scater:::.get_palette("tableau20")) +
  scale_x_discrete(labels=c(rep("",10), "HPRT1", rep("",length(genes_for_plotting)-10))) +
  scale_y_log10(labels=scales::label_log(digits=1)) +
  theme_void() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=NULL,
        axis.title.y=element_text(size=9, angle=90),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Globally lowly expressed genes")

pdf("~/Project_PBMCage/Plots/Globally_lowlyExprGenes.pdf", width=6, height=5)
p_thrdeter; p_gse; p_expr
dev.off()

#####################################





###
### Get the age-stable genes (稳定表达基因)
#####################################
###

### Determine globally age-stable genes
df=read.delim("~/Project_PBMCage/pseudobulk_global.txt")
celltype.wide_lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_global.rds")
df=df %>%
  subset((rownames(.) %in% celltype.wide_lowlyExprGenes)==FALSE)
df=t(df)
ages=as.factor(gsub("X|_cluster.*","",rownames(df)))
rownames(df)=ages

sd_=apply(df, 2, sd)
sd_df=data.frame(Gene=names(sd_), SD=sd_)

sd_thr=3 # with SD=3, a num is considered 99.7% confidence not by chance dropping closely around the mean of a population

p_sd=
  ggplot(sd_df, aes(x=log2(SD))) +
  geom_density(lwd=1, fill="white", alpha=0.25) +
  scale_x_continuous(labels=scales::label_math(2^.x)) +
  geom_vline(xintercept=log2(sd_thr)) +
  theme_classic()

### Determine the age-stable genes globally or in each celltype
celltype.wide_housekeepingGenes=colnames(df)[apply(df, 2, sd)<sd_thr]

pseudobulkdata=read.delim("~/Project_PBMCage/try_for_dtw.txt", stringsAsFactors=F)
Celltype=gsub(".*_cluster","",colnames(pseudobulkdata))
Age=as.integer(gsub("X|_.*","",colnames(pseudobulkdata)))
pseudobulkdata_arranged=as.data.frame(t(pseudobulkdata))
pseudobulkdata_arranged=cbind(Celltype, Age, pseudobulkdata_arranged)

lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_perCelltype.rds")

all_celltypes=names(table(Celltype))
StableGenes=list()
for (i in 1:length(all_celltypes)) {
  df_=subset(pseudobulkdata_arranged, Celltype==all_celltypes[i])
  ages=df_[,"Age"]
  df_=df_[,c(3:ncol(df_))]
  rownames(df_)=ages
  
  df_=df[, (colnames(df) %in% lowlyExprGenes[[i]])==FALSE]
  StableGenes[[i]]=colnames(df_)[apply(df_, 2, sd)<sd_thr]
}

names(StableGenes)=all_celltypes
saveRDS(StableGenes, "~/Project_PBMCage/Tempt_RDS/StableGenes_perCelltype.rds")

All_housekeepingGenes=do.call(c, StableGenes)
All_housekeepingGenes=All_housekeepingGenes[!duplicated(All_housekeepingGenes)]
expression_matrix_chosen=cbind(Celltype, Age, as.data.frame(pseudobulkdata_arranged[,All_housekeepingGenes]))
long_exp.matrix.chosen=reshape2::melt(expression_matrix_chosen, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
long_exp.matrix.chosen=long_exp.matrix.chosen %>%
  group_by(Celltype, Age) %>%
  summarise(across(Expression, mean)) %>%
  as.data.frame()
long_exp.matrix.chosen$Group="Age-stable genes"
# add aging-related genes
AgingDB=read.gmt("~/Project_PBMCage/DEMAGALHAES_AGING_UP.v2023.2.Hs.gmt")
Aging_Matrix=pseudobulkdata_arranged[, match(AgingDB$gene, colnames(pseudobulkdata_arranged))[!is.na(match(AgingDB$gene, colnames(pseudobulkdata_arranged)))]]
expression_matrix_age=cbind(Celltype, Age, as.data.frame(Aging_Matrix))
long_exp.matrix.age=reshape2::melt(expression_matrix_age, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
long_exp.matrix.age=long_exp.matrix.age %>%
  group_by(Celltype, Age) %>%
  summarise(across(Expression, mean)) %>%
  as.data.frame()
long_exp.matrix.age$Group="Age-related genes"

df_for_plot=rbind(long_exp.matrix.chosen, long_exp.matrix.age)

p_stableGeneExpr=
  ggplot(df_for_plot, aes(x=Age, y=Expression, color=Celltype, shape=Group)) +
  geom_point(alpha=0.5) +
  scale_shape_manual(values=c(1,3)) +
  scale_color_manual(values=scater:::.get_palette("tableau20")) +
  scale_y_continuous(trans="pseudo_log") +
  theme_void() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=NULL,
        axis.title.y=element_text(size=9, angle=90),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Age-stable genes") +
  geom_smooth(method='auto', show.legend=FALSE, linewidth=0.5, fill="gray")

#####################################





###
### Get the Lag-N correlation (基因时序表达的同调性及前后关系)
#####################################
###

###
### Calculate lags and correlation between each celltype per gene
pseudobulkdata=read.delim("~/Project_PBMCage/try_for_dtw.txt", stringsAsFactors=F)
Celltype=gsub(".*_cluster","",colnames(pseudobulkdata))
all_celltypes=names(table(Celltype))
compair_pair_list=c()
for (i in 1:length(all_celltypes)) {
  rest=c(1:length(all_celltypes))[-i]
  compair_pair=lapply(1:length(rest), function(x) c(all_celltypes[i], all_celltypes[rest[x]]))
  compair_pair_list=c(compair_pair_list, compair_pair)
}

lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_perCelltype.rds")
StableGenes=readRDS("~/Project_PBMCage/Tempt_RDS/StableGenes_perCelltype.rds")
Age=as.integer(gsub("X|_.*","",colnames(pseudobulkdata)))
correlation_mtx=matrix(NA, nrow=nrow(pseudobulkdata), ncol=length(compair_pair_list))
h_mtx=matrix(NA, nrow=nrow(pseudobulkdata), ncol=length(compair_pair_list))

library(TSA)
processbar=txtProgressBar(min=0, max=length(compair_pair_list), style=3, width=60, char="=")
for (i in 1:length(compair_pair_list)) {
  type1=compair_pair_list[[i]][1]
  type2=compair_pair_list[[i]][2]
  lowlyExprGenes.1=lowlyExprGenes[[which(all_celltypes==type1)]]
  lowlyExprGenes.2=lowlyExprGenes[[which(all_celltypes==type2)]]
  StableGenes.1=StableGenes[[which(all_celltypes==type1)]]
  StableGenes.2=StableGenes[[which(all_celltypes==type2)]]
  filter_out=c(lowlyExprGenes.1,lowlyExprGenes.2,StableGenes.1,StableGenes.2)[!duplicated(c(lowlyExprGenes.1,lowlyExprGenes.2,StableGenes.1,StableGenes.2))]
  keep.df=pseudobulkdata %>% subset((rownames(.) %in% filter_out)==FALSE)
  keep.df=as.data.frame(t(keep.df))
  keep.df=cbind(Celltype, Age, keep.df)
  df1=keep.df %>% subset(Celltype==type1)
  df2=keep.df %>% subset(Celltype==type2)
  for (j in 1:(ncol(df1)-2)) {
    genename=colnames(df1)[j+2]
    series1=df1[,j+2]
    series2=df2[,j+2]
    results=prewhiten(series1, series2)
    Abs.correlation=max(abs(results$ccf$acf))
    h=results$ccf$lag[which(abs(results$ccf$acf)==Abs.correlation)]
    correlation=results$ccf$acf[which(results$ccf$lag==h)]
    gene.idx=which(rownames(pseudobulkdata)==genename)
    
    correlation_mtx[gene.idx, i]=correlation
    h_mtx[gene.idx, i]=h
  }
  setTxtProgressBar(processbar, i)
}
close(processbar)

correlation_mtx=as.data.frame(correlation_mtx)
colnames(correlation_mtx)=unlist(lapply(compair_pair_list, function(x) paste0(x[1],".vs.",x[2])))
rownames(correlation_mtx)=rownames(pseudobulkdata)
h_mtx=as.data.frame(h_mtx)
colnames(h_mtx)=colnames(correlation_mtx)
rownames(h_mtx)=rownames(correlation_mtx)

saveRDS(correlation_mtx, "~/Project_PBMCage/Results/correlationMatrix.rds")
saveRDS(h_mtx, "~/Project_PBMCage/Results/hMatrix.rds")

#######################################





###
### Analyze the Corr and Lag
#####################################
###

### Arrange correlation and h matrix for plotting
correlation_mtx=readRDS("~/Project_PBMCage/Results/correlationMatrix.rds")
h_mtx=readRDS("~/Project_PBMCage/Results/hMatrix.rds")
correlation_mtx$Gene=rownames(correlation_mtx)
long_exp.matrix.cor=reshape2::melt(correlation_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="Correlation")
h_mtx$Gene=rownames(h_mtx)
long_exp.matrix.h=reshape2::melt(h_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="h")
long_exp.matrix=long_exp.matrix.cor
long_exp.matrix$h=long_exp.matrix.h$h
long_exp.matrix$abs_corr=abs(long_exp.matrix$Correlation)

### Plot distribution of correlation
pdf("~/Project_PBMCage/Plots/Correlation_Distribution_Global.pdf")
p_correlation=
  ggplot(long_exp.matrix.cor, aes(x=Correlation, color=Comparison_pair)) +
  geom_density(lwd=1, fill="transparent", alpha=0.25) +
  theme_classic() +
  theme(legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"),
        legend.position="none")
p_correlation
dev.off()

all_celltypes=gsub("\\.vs.*","",long_exp.matrix.cor$Comparison_pair)[!duplicated(gsub("\\.vs.*","",long_exp.matrix.cor$Comparison_pair))]
pdf("~/Project_PBMCage/Plots/Correlation_Distribution_PerCelltype.pdf")
for (i in 1:length(all_celltypes)) {
  df=long_exp.matrix.cor %>%
    subset(gsub("\\.vs.*","",Comparison_pair)==all_celltypes[i])
  p_correlation_eachCelltype=
    ggplot(df, aes(x=Correlation, color=Comparison_pair)) +
    geom_density(lwd=1, fill="transparent", alpha=0.25) +
    theme_classic() +
    theme(legend.key.height=unit(3,"mm"),
          legend.key.width=unit(3,"mm"))
  plot(p_correlation_eachCelltype)
}
dev.off()

### Plot scatterplot of correlation vs. h
pdf("~/Project_PBMCage/Plots/Correlation_vs_Lag_Global.pdf")
ggplot(long_exp.matrix, aes(x=Correlation, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag")
ggplot(long_exp.matrix, aes(x=abs_corr, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag",
       x="Abs(Correlation)")
dev.off()

#####################################





###
### Globally synchronized genes (defined as: a large portion of large abs positive correlation + small abs(h))
#####################################
###

### Determine the threshold for "small correlation" and "small h" based on the plot
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h<=h_thr & h>=-h_thr)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)

pdf("~/Project_PBMCage/Plots/Correlation_H_SynchronizedGenes.pdf")
ggplot(df_comparison.pair, aes(x=No.Comparison.Pair)) +
  geom_density(lwd=1, alpha=0.25) +
  theme_classic() +
  theme(legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"))
dev.off()

### Deternube the globally syn. genes
# if more than half of the comparison pairs have large abs(correlation) + small abs(h) in a gene,
# then the gene is considered globally synchronized
global_syn_genes=(df_comparison.pair %>% subset(No.Comparison.Pair>120))[["Gene"]]
saveRDS(global_syn_genes, "~/Project_PBMCage/Tempt_RDS/Global_synchronized_genes.rds")

### Enrich the globally synchronized genes
library(clusterProfiler)
library(org.Hs.eg.db)
gse=enrichGO(global_syn_genes, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
pdf("~/Project_PBMCage/Plots/Correlation_H_SynchronizedGenes_Enrichment.pdf")
gseplot0=
  dotplot(gse, showCategory=8, font.size=10) +
  labs(title="Globally synchronized genes") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10))
plot(gseplot0)
dev.off()

### Determine the genes that are expressed simultaneously in all the celltypes
range_of_h=lapply(global_syn_genes, function(x) range((df_ %>% subset(Gene==x))[["h"]]))
names(range_of_h)=global_syn_genes
determine_=lapply(range_of_h, function(x) x[1]==0 & x[2]==0)
table(determine_==TRUE) # so all these genes are expressed simultaneously in all the celltypes

### Analyze gene-gene CCF of these globally syn. genes
pseudobulkdata=read.delim("~/Project_PBMCage/try_for_dtw.txt", stringsAsFactors=F)
pseudobulk_global=read.delim("~/Project_PBMCage/pseudobulk_global.txt", stringsAsFactors=F)
pseudobulkdata_all_filtered=cbind(pseudobulkdata[global_syn_genes,], pseudobulk_global[global_syn_genes,])

Celltype=gsub(".*_cluster","",colnames(pseudobulkdata_all_filtered))
all_celltypes=names(table(Celltype))

gene_compair_pair_list=c()
for (i in 1:nrow(pseudobulkdata_all_filtered)) {
  rest=c(1:nrow(pseudobulkdata_all_filtered))[-i]
  gene_compair_pair=lapply(1:length(rest), function(x) c(global_syn_genes[i], global_syn_genes[rest[x]]))
  gene_compair_pair_list=c(gene_compair_pair_list, gene_compair_pair)
}

synGenes_correlation_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))
synGenes_h_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))

library(TSA)
processbar=txtProgressBar(min=0, max=length(gene_compair_pair_list), style=3, width=60, char="=")
for (i in 1:length(gene_compair_pair_list)) {
  gene1=gene_compair_pair_list[[i]][1]
  gene2=gene_compair_pair_list[[i]][2]
  for (j in 1:length(all_celltypes)) {
    celltype_name=all_celltypes[j]
    celltype.df=pseudobulkdata_all_filtered[, grepl(celltype_name, colnames(pseudobulkdata_all_filtered))]
    celltype.df_filtered=celltype.df %>% subset(rownames(.) %in% c(gene1,gene2))
    celltype.df_filtered=as.data.frame(t(celltype.df_filtered))
    celltype.df_filtered=cbind(Age=gsub("X|_.*","",rownames(celltype.df_filtered)), celltype.df_filtered)
    series1=celltype.df_filtered[,2]
    series2=celltype.df_filtered[,3]
    results=prewhiten(series1, series2)
    Abs.correlation=max(abs(results$ccf$acf))
    h=results$ccf$lag[which(abs(results$ccf$acf)==Abs.correlation)]
    correlation=results$ccf$acf[which(results$ccf$lag==h)]
    
    synGenes_correlation_mtx[i, j]=correlation
    synGenes_h_mtx[i, j]=h
  }
  setTxtProgressBar(processbar, i)
}
close(processbar)

synGenes_correlation_mtx=as.data.frame(synGenes_correlation_mtx)
rownames(synGenes_correlation_mtx)=unlist(lapply(gene_compair_pair_list, function(x) paste0(x[1],".vs.",x[2])))
colnames(synGenes_correlation_mtx)=all_celltypes
synGenes_h_mtx=as.data.frame(synGenes_h_mtx)
rownames(synGenes_h_mtx)=rownames(synGenes_correlation_mtx)
colnames(synGenes_h_mtx)=colnames(synGenes_correlation_mtx)

saveRDS(synGenes_correlation_mtx, "~/Project_PBMCage/Results/correlationMatrix_synGenes.rds")
saveRDS(synGenes_h_mtx, "~/Project_PBMCage/Results/hMatrix_synGenes.rds")

### Plot the expression of globally syn. genes in various celltypes
df_global.syn.gene.expr=as.data.frame(t(pseudobulkdata_all_filtered))
df_global.syn.gene.expr_age=gsub("X|_cluster.*","",rownames(df_global.syn.gene.expr))
df_global.syn.gene.expr_celltype=gsub(".*cluster","",rownames(df_global.syn.gene.expr))
df_global.syn.gene.expr=cbind(Age=as.factor(df_global.syn.gene.expr_age), Celltype=as.factor(df_global.syn.gene.expr_celltype), df_global.syn.gene.expr)
long_df_global.syn.gene=reshape2::melt(df_global.syn.gene.expr, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")

pdf("~/Project_PBMCage/Plots/Correlation_H_SynchronizedGenes_Expr.pdf")
ggplot(long_df_global.syn.gene, aes(x=Age, y=Expression, color=Celltype)) +
  geom_point(size=0.5, shape=1) +
  scale_color_manual(values=scater:::.get_palette("tableau20"),
                     labels=c(all_celltypes[1:(length(all_celltypes)-1)], "All")) +
  scale_y_log10(labels=scales::label_log(digits=1)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19))) +
  labs(title="Globally synchronized gene",
       x="Age")
dev.off()

### Analyze the gene-gene lags among globally syn. genes
synGenes_correlation_mtx=readRDS("~/Project_PBMCage/Results/correlationMatrix_synGenes.rds")
synGenes_h_mtx=readRDS("~/Project_PBMCage/Results/hMatrix_synGenes.rds")
synGenes_correlation_mtx$Comparison_pair=rownames(synGenes_correlation_mtx)
long_df_synGene_matrix=reshape2::melt(synGenes_correlation_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="Correlation")
synGenes_h_mtx$Comparison_pair=rownames(synGenes_h_mtx)
long_df_synGene_matrix_h=reshape2::melt(synGenes_h_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="h")
long_df_synGene_matrix=cbind(long_df_synGene_matrix, h=long_df_synGene_matrix_h$h)

range(long_df_synGene_matrix$h) # turns out to be all 0, so all these genes are expressed simultaneously
range(long_df_synGene_matrix$Correlation)

pdf("~/Project_PBMCage/Plots/Correlation_Distribution_SynchronizedGenes_PerCelltype.pdf")
ggplot(long_df_synGene_matrix, aes(x=Correlation, color=Celltype)) +
  geom_density(lwd=1, alpha=0) +
  scale_color_manual(values=scater:::.get_palette("tableau20"),
                     labels=c(all_celltypes[1:(length(all_celltypes)-1)], "All")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19))) +
  labs(title="Synchronicity among globally sync. genes")
dev.off()

#####################################





###
### Cell-dependent leading/lagging genes (defined as: for each celltype comparison pair, large abs positive correlation + large abs(h))
#####################################
###

### Get cell-dependent leading/lagging genes
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h>h_thr | h<(-h_thr))

gene_cluster=split(df_$Gene, paste0(df_$Comparison_pair, "__", df_$h))
gene_cluster_sort=gene_cluster[order(abs(as.numeric(gsub(".*__","",names(gene_cluster)))), decreasing=T)]

pdf("~/Project_PBMCage/Plots/Correlation_H_Leading.Lagging.Genes_Expr.pdf")
for (i in 1:length(gene_cluster_sort)) {
  cluster_=gene_cluster_sort[i]
  h_value=gsub(".*__","",names(cluster_))
  cluster=cluster_[[1]]
  type1=gsub("\\.vs.*","",names(cluster_))
  type2=gsub(".*vs\\.|__.*","",names(cluster_))
  df_byCluster.Gene=cbind(pseudobulkdata[cluster,], pseudobulk_global[cluster,])
  df_byCluster.Gene_t_age=gsub("X|_cluster.*","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t_celltype=gsub(".*cluster","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t=data.frame(Age=as.factor(df_byCluster.Gene_t_age), Celltype=as.factor(df_byCluster.Gene_t_celltype), as.data.frame(t(df_byCluster.Gene)))
  df_byCluster.Gene_t_chosen=df_byCluster.Gene_t %>%
    subset(Celltype %in% c(type1, type2))
  long_df_leading.lagging.gene.expr=reshape2::melt(df_byCluster.Gene_t_chosen, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
  
  cluster_name=ifelse(length(cluster)>5, paste0(cluster[1:5], collapse=", "), paste0(cluster, collapse=", "))
  
  p_x=
    ggplot(long_df_leading.lagging.gene.expr, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
    geom_point(shape=1, alpha=0.7) +
    scale_color_manual(values=scater:::.get_palette("tableau20")[c(1,11)]) +
    scale_y_continuous(trans="pseudo_log") +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=9),
          legend.text=element_text(size=8),
          title=element_text(size=10)) + 
    guides(color=guide_legend(override.aes=list(size=2, alpha=1))) +
    labs(title=paste0(cluster_name, ";     h=", h_value),
         x="Age") +
    geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)
  
  plot(p_x)
}
dev.off()
# turns out not to be meaningful

#####################################





###
### Negatively regulated genes (defined as: for each celltype comparison pair, large abs of negative correlation + any h)
#####################################
###

### Determine the threshold for "large abs" of negative correlation
corr_thr=0.5

df_=long_exp.matrix %>%
  subset(Correlation<=-corr_thr)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)
table(df_comparison.pair$No.Comparison.Pair)
# it turns out that for each of these neg.reg.genes, only one comparison pair of celltypes shows a strong negative corr

type1=gsub("\\.vs.*","",df_$Comparison_pair)
type2=gsub(".*vs\\.","",df_$Comparison_pair)
pdf("~/Project_PBMCage/Plots/Correlation_H_Neg.Reg.Genes_Expr.pdf")
for (i in 1:length(type1)) {
  type1_=type1[i]
  type2_=type2[i]
  gene_=df_$Gene[i]
  h_value=df_$h[i]
  corr=df_$Correlation[i]
  exprssion=pseudobulkdata_arranged %>% subset(Celltype %in% c(type1_,type2_))
  df_byCluster.Gene_t_age=exprssion$Age
  df_byCluster.Gene_t_celltype=exprssion$Celltype
  df_byCluster.Gene_t=data.frame(Age=as.factor(df_byCluster.Gene_t_age), Celltype=as.factor(df_byCluster.Gene_t_celltype), as.data.frame(exprssion[,gene_]))
  colnames(df_byCluster.Gene_t)=c("Age","Celltype","Expression")
  
  p_x=
    ggplot(df_byCluster.Gene_t, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
    geom_point(shape=1, alpha=0.7) +
    scale_color_manual(values=scater:::.get_palette("tableau20")[c(1,11)]) +
    scale_y_continuous(trans="pseudo_log") +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=9),
          legend.text=element_text(size=8),
          title=element_text(size=10)) + 
    guides(color=guide_legend(override.aes=list(size=2, alpha=1))) +
    labs(title=paste0(gene_, ";     h=", h_value, ";     correlation=", formatC(corr, digits=2)),
         x="Age") +
    geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)
  
  plot(p_x)
}
dev.off()
# turns out not to be meaningful

#####################################





###
### Globally less correlated genes (defined as: small abs(correlation) + any h)
#####################################
###

### Set the threshold based on the distribution of correlation
corr_thr=0.25

df_=long_exp.matrix %>%
  subset(abs(Correlation)<corr_thr)

### Determine the globally less correlated genes, if any
No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)

ggplot(df_comparison.pair, aes(x=No.Comparison.Pair)) +
  geom_density(lwd=1, alpha=0.25) +
  theme_classic() +
  theme(legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"))
df_comparison.pair=df_comparison.pair %>%
  subset(No.Comparison.Pair>=length(table(df_$Comparison_pair))/2)
# empty, i.e. no globally less correlated ones, which makes sense since every gene should have a shared expressed pattern between at least a couple of celltypes

### Determine the less correlated genes in various celltypes
library(clusterProfiler)
library(org.Hs.eg.db)
geneList=split(df_$Gene, df_$Comparison_pair)
GO_enrich_compare=compareCluster(geneCluster=geneList, fun=enrichGO, OrgDb="org.Hs.eg.db", pvalueCutoff=0.05)
# no term enriched, which means that the less correlated genes in each celltype comparison pair are just random genes

#####################################





###
### CCF analysis on other groups of celltypes
#####################################
###

library(edgeR)
library(reshape)
library(Seurat)
library(dplyr)

### Prepare for CCF results collection

object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
table(rownames(obj_annotation)==colnames(object)) # check
object=AddMetaData(object, metadata=obj_annotation$CD8T_NK_subsets, col.name="CD8T_NK_subsets")

names(object[[]])
table(object$CD8T_NK_subsets)
table(object$reassigned)

chose_from=names(table(object$CD8T_NK_subsets))
Celltype_groups=list(DCs=c(chose_from[grepl("DC$|DC[^CA]",chose_from)]),
                     B_cells=c(chose_from[grepl("Bnaive$|Bmem|PB$|PC$",chose_from)]),
                     Monocytes=c(chose_from[grepl("Mono$",chose_from)]),
                     CD4_T=c(chose_from[grepl("^CD4|ISAGhi T|^Treg",chose_from)]),
                     CD8_T=c(chose_from[grepl("^CD8",chose_from)]),
                     ILCs=c(chose_from[grepl("^NK\\.|ILC",chose_from)]),
                     NKTs=c(chose_from[grepl("^(NKT)",chose_from)]),
                     Other_lymph=c(chose_from[grepl("CLP|dnT|gdT|MAIT",chose_from)]),
                     Other_myelo=c(chose_from[grepl("EMP|^Eryth|Mast|Platelet",chose_from)]))

Run_my_CCF=function(group.idx) {
  cateogry_name=names(Celltype_groups[group.idx])
  celltype_group=Celltype_groups[[group.idx]]
  obj_sub=subset(object, CD8T_NK_subsets %in% celltype_group)
  pb_obj=Seurat2PB(obj_sub, sample="age", cluster="CD8T_NK_subsets")
  keep.genes=filterByExpr(pb_obj, group=pb_obj$samples$cluster, min.count=2, min.prop=2/nrow(pb_obj$samples))
  table(keep.genes)
  pb_obj=pb_obj[keep.genes, , keep=FALSE]
  write.table(pb_obj$counts, "~/Project_PBMCage/ccf_tempt_try_for_dtw.txt", sep="\t")
  
  pb_obj_all=Seurat2PB(object, sample="age", cluster="orig.ident")
  pseudobulkdata=read.delim("~/Project_PBMCage/ccf_tempt_try_for_dtw.txt", stringsAsFactors=F)
  keep.genes=rownames(pseudobulkdata)
  pb_obj_all=pb_obj_all[keep.genes, , keep=FALSE]
  write.table(pb_obj_all$counts, "~/Project_PBMCage/ccf_tempt_pseudobulk_global.txt", sep="\t")
  
  # remove lowly-expressed genes
  df=read.delim("~/Project_PBMCage/ccf_tempt_pseudobulk_global.txt")
  df=t(df)
  ages=as.factor(gsub("X|_cluster.*","",rownames(df)))
  rownames(df)=ages
  
  lowExpr_thr=1.5 # log10(1.5)=0.1760913
  celltype.wide_lowlyExprGenes=colnames(df)[apply(df, 2, mean)<lowExpr_thr]
  
  pseudobulkdata=read.delim("~/Project_PBMCage/ccf_tempt_try_for_dtw.txt", stringsAsFactors=F)
  Celltype=gsub(".*_cluster","",colnames(pseudobulkdata))
  Age=as.integer(gsub("X|_.*","",colnames(pseudobulkdata)))
  pseudobulkdata_arranged=as.data.frame(t(pseudobulkdata))
  pseudobulkdata_arranged=cbind(Celltype, Age, pseudobulkdata_arranged)
  
  all_celltypes=names(table(Celltype))
  lowlyExprGenes=list()
  for (i in 1:length(all_celltypes)) {
    df_=subset(pseudobulkdata_arranged, Celltype==all_celltypes[i])
    ages=df_[,"Age"]
    df_=df_[,c(3:ncol(df_))]
    rownames(df_)=ages
    
    lowlyExprGenes[[i]]=colnames(df_)[apply(df_, 2, mean)<lowExpr_thr]
  }
  names(lowlyExprGenes)=all_celltypes
  saveRDS(lowlyExprGenes, "~/Project_PBMCage/Tempt_RDS/ccf_tempt_Expr_lowlyExprGenes_perCelltype.rds")
  saveRDS(celltype.wide_lowlyExprGenes, "~/Project_PBMCage/Tempt_RDS/ccf_tempt_Expr_lowlyExprGenes_global.rds")
  
  # remove age-stable genes
  df=read.delim("~/Project_PBMCage/ccf_tempt_pseudobulk_global.txt")
  celltype.wide_lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/ccf_tempt_Expr_lowlyExprGenes_global.rds")
  df=df %>%
    subset((rownames(.) %in% celltype.wide_lowlyExprGenes)==FALSE)
  df=t(df)
  ages=as.factor(gsub("X|_cluster.*","",rownames(df)))
  rownames(df)=ages
  
  sd_=apply(df, 2, sd)
  sd_df=data.frame(Gene=names(sd_), SD=sd_)
  
  sd_thr=3 # with SD=3, a num is considered 99.7% confidence not by chance dropping closely around the mean of a population
  
  celltype.wide_housekeepingGenes=colnames(df)[apply(df, 2, sd)<sd_thr]
  
  pseudobulkdata=read.delim("~/Project_PBMCage/ccf_tempt_try_for_dtw.txt", stringsAsFactors=F)
  Celltype=gsub(".*_cluster","",colnames(pseudobulkdata))
  Age=as.integer(gsub("X|_.*","",colnames(pseudobulkdata)))
  pseudobulkdata_arranged=as.data.frame(t(pseudobulkdata))
  pseudobulkdata_arranged=cbind(Celltype, Age, pseudobulkdata_arranged)
  
  lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/ccf_tempt_Expr_lowlyExprGenes_perCelltype.rds")
  
  all_celltypes=names(table(Celltype))
  StableGenes=list()
  for (i in 1:length(all_celltypes)) {
    df_=subset(pseudobulkdata_arranged, Celltype==all_celltypes[i])
    ages=df_[,"Age"]
    df_=df_[,c(3:ncol(df_))]
    rownames(df_)=ages
    
    df_=df[, (colnames(df) %in% lowlyExprGenes[[i]])==FALSE]
    StableGenes[[i]]=colnames(df_)[apply(df_, 2, sd)<sd_thr]
  }
  
  names(StableGenes)=all_celltypes
  saveRDS(StableGenes, "~/Project_PBMCage/Tempt_RDS/ccf_tempt_StableGenes_perCelltype.rds")
  
  # run CCF
  pseudobulkdata=read.delim("~/Project_PBMCage/ccf_tempt_try_for_dtw.txt", stringsAsFactors=F)
  Celltype=gsub(".*_cluster","",colnames(pseudobulkdata))
  all_celltypes=names(table(Celltype))
  compair_pair_list=c()
  for (i in 1:length(all_celltypes)) {
    rest=c(1:length(all_celltypes))[-i]
    compair_pair=lapply(1:length(rest), function(x) c(all_celltypes[i], all_celltypes[rest[x]]))
    compair_pair_list=c(compair_pair_list, compair_pair)
  }
  
  lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/ccf_tempt_Expr_lowlyExprGenes_perCelltype.rds")
  StableGenes=readRDS("~/Project_PBMCage/Tempt_RDS/ccf_tempt_StableGenes_perCelltype.rds")
  Age=as.integer(gsub("X|_.*","",colnames(pseudobulkdata)))
  correlation_mtx=matrix(NA, nrow=nrow(pseudobulkdata), ncol=length(compair_pair_list))
  h_mtx=matrix(NA, nrow=nrow(pseudobulkdata), ncol=length(compair_pair_list))
  
  library(TSA)
  processbar=txtProgressBar(min=0, max=length(compair_pair_list), style=3, width=60, char="=")
  for (i in 1:length(compair_pair_list)) {
    print(i)
    type1=compair_pair_list[[i]][1]
    type2=compair_pair_list[[i]][2]
    lowlyExprGenes.1=lowlyExprGenes[[which(all_celltypes==type1)]]
    lowlyExprGenes.2=lowlyExprGenes[[which(all_celltypes==type2)]]
    StableGenes.1=StableGenes[[which(all_celltypes==type1)]]
    StableGenes.2=StableGenes[[which(all_celltypes==type2)]]
    filter_out=c(lowlyExprGenes.1,lowlyExprGenes.2,StableGenes.1,StableGenes.2)[!duplicated(c(lowlyExprGenes.1,lowlyExprGenes.2,StableGenes.1,StableGenes.2))]
    keep.df=pseudobulkdata %>% subset((rownames(.) %in% filter_out)==FALSE)
    keep.df=as.data.frame(t(keep.df))
    keep.df=cbind(Celltype, Age, keep.df)
    df1=keep.df %>% subset(Celltype==type1)
    df2=keep.df %>% subset(Celltype==type2)
    
    if (nrow(df1)==78 & nrow(df2)==78) {
      for (j in 1:(ncol(df1)-2)) {
        genename=colnames(df1)[j+2]
        series1=df1[,j+2]
        series2=df2[,j+2]
        results=prewhiten(series1, series2)
        Abs.correlation=max(abs(results$ccf$acf))
        h=results$ccf$lag[which(abs(results$ccf$acf)==Abs.correlation)]
        correlation=results$ccf$acf[which(results$ccf$lag==h)]
        gene.idx=which(rownames(pseudobulkdata)==genename)
        
        correlation_mtx[gene.idx, i]=correlation
        h_mtx[gene.idx, i]=h
      }
    } else {
      correlation_mtx[gene.idx, i]=NA
      h_mtx[gene.idx, i]=NA
    }
    
    setTxtProgressBar(processbar, i)
  }
  close(processbar)
  
  correlation_mtx=as.data.frame(correlation_mtx)
  colnames(correlation_mtx)=unlist(lapply(compair_pair_list, function(x) paste0(x[1],".vs.",x[2])))
  rownames(correlation_mtx)=rownames(pseudobulkdata)
  h_mtx=as.data.frame(h_mtx)
  colnames(h_mtx)=colnames(correlation_mtx)
  rownames(h_mtx)=rownames(correlation_mtx)
  
  saveRDS(correlation_mtx, paste0("~/Project_PBMCage/Results/correlationMatrix_", cateogry_name, ".rds"))
  saveRDS(h_mtx, paste0("~/Project_PBMCage/Results/hMatrix_", cateogry_name, ".rds"))
}

Run_my_CCF(1)
Run_my_CCF(2)
Run_my_CCF(3)
Run_my_CCF(4)
Run_my_CCF(5)
Run_my_CCF(6)
Run_my_CCF(7)
Run_my_CCF(8)
Run_my_CCF(9)

#####################################





###
### CCF analysis on DC lineage
#####################################
###

### Analyze the Corr and Lag
# Arrange correlation and h matrix for plotting
correlation_mtx=readRDS("~/Project_PBMCage/Results/correlationMatrix_DCs.rds")
h_mtx=readRDS("~/Project_PBMCage/Results/hMatrix_DCs.rds")
correlation_mtx$Gene=rownames(correlation_mtx)
long_exp.matrix.cor=reshape2::melt(correlation_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="Correlation")
h_mtx$Gene=rownames(h_mtx)
long_exp.matrix.h=reshape2::melt(h_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="h")
long_exp.matrix=long_exp.matrix.cor
long_exp.matrix$h=long_exp.matrix.h$h
long_exp.matrix$abs_corr=abs(long_exp.matrix$Correlation)
long_exp.matrix=na.omit(long_exp.matrix)
# empty, meaning that at least one of the subset of DCs are missing at some ages, so that the CCF analysis is skipped 

#####################################





###
### CCF analysis on B cell lineage
#####################################
###

### Analyze the Corr and Lag
# Arrange correlation and h matrix for plotting
correlation_mtx=readRDS("~/Project_PBMCage/Results/correlationMatrix_B_cells.rds")
h_mtx=readRDS("~/Project_PBMCage/Results/hMatrix_B_cells.rds")
compare_idx=list()
for (idx in 1:length(colnames(correlation_mtx))) {
  pair_1=colnames(correlation_mtx)[idx]
  pair_2=paste0(gsub(".*vs\\.","",pair_1), ".vs.", gsub("\\.vs.*","",pair_1))
  compare_idx1=c(1:length(colnames(correlation_mtx)))[idx]
  compare_idx2=which(colnames(correlation_mtx)==pair_2)
  if ((compare_idx2 %in% c(1:idx))==FALSE) {
    compare_idx=c(compare_idx, list(c(compare_idx1, compare_idx2)))
  }
}
idx_1_all=c()
for (i in 1:length(compare_idx)) {
  idx1=compare_idx[[i]][1]
  idx2=compare_idx[[i]][2]
  correlation_updated=sapply(1:length(correlation_mtx[,idx1]), function(x) ifelse(abs(correlation_mtx[,idx1][x])>=abs(correlation_mtx[,idx2][x]), correlation_mtx[,idx1][x], correlation_mtx[,idx2][x]))
  col_chosen=sapply(1:length(correlation_updated), function(x) ifelse(correlation_updated[x]==correlation_mtx[,idx1][x], 1, 2))
  h_updated=sapply(1:length(h_mtx[,idx1]), function(x) ifelse(col_chosen[x]==1, h_mtx[,idx1], h_mtx[,idx2]))
  
  correlation_mtx[,idx1]=correlation_updated
  h_mtx[,idx1]=h_updated
  idx_1_all=c(idx_1_all,idx1)
}
correlation_mtx=correlation_mtx[,idx_1_all]
h_mtx=h_mtx[,idx_1_all]

correlation_mtx$Gene=rownames(correlation_mtx)
long_exp.matrix.cor=reshape2::melt(correlation_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="Correlation")
h_mtx$Gene=rownames(h_mtx)
long_exp.matrix.h=reshape2::melt(h_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="h")

# remove the globally lowly expressed genes and age-stable genes
lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_perCelltype.rds"); names(lowlyExprGenes)
celltype.wide_lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_global.rds")
StableGenes=readRDS("~/Project_PBMCage/Tempt_RDS/StableGenes_perCelltype.rds"); names(StableGenes)
global_syn_genes=readRDS("~/Project_PBMCage/Tempt_RDS/Global_synchronized_genes.rds")

gene.to.remove=c(lowlyExprGenes$"Bmem",lowlyExprGenes$"Bnaive",lowlyExprGenes$"PB.PC",
                 StableGenes$"Bmem",StableGenes$"Bnaive",StableGenes$"PB.PC",
                 celltype.wide_lowlyExprGenes,
                 global_syn_genes)
gene.to.remove=gene.to.remove[!duplicated(gene.to.remove)]

long_exp.matrix.cor=long_exp.matrix.cor %>%
  subset((Gene %in% gene.to.remove)==FALSE)
long_exp.matrix.h=long_exp.matrix.h %>%
  subset((Gene %in% gene.to.remove)==FALSE)

long_exp.matrix=long_exp.matrix.cor
long_exp.matrix$h=long_exp.matrix.h$h
long_exp.matrix$abs_corr=abs(long_exp.matrix$Correlation)
long_exp.matrix=na.omit(long_exp.matrix)

# Plot distribution of correlation
p_correlation=
  ggplot(long_exp.matrix, aes(x=Correlation, color=Comparison_pair)) +
  geom_density(lwd=1, fill="transparent", alpha=0.25) +
  theme_classic() +
  theme(legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"))
p_correlation
# it seems that the correlation between Bnaive and Classical.Bmem is > that between Bnaive and DN.Bmem > that between Classical.Bmem and DN.Bmem

### Plot scatterplot of correlation vs. h
ggplot(long_exp.matrix, aes(x=Correlation, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag - B cell lineage")
ggplot(long_exp.matrix, aes(x=abs_corr, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag - B cell lineage",
       x="Abs(Correlation)")
dev.off()
# there's a portion of genes with high lags and high positive corr, meaning that there are leading/lagging genes

### B-specific synchronized genes (defined as: a large portion of large abs positive correlation + small abs(h))
# determine the threshold for "small correlation" and "small h" based on the plot
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h<=h_thr & h>=-h_thr)

df_ %>% group_by(Comparison_pair) %>% slice_max(abs_corr, n=5)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)
table(as.vector(table(df_comparison.pair$Gene)))

# determine the B-lineage common syn. genes
common_thr=length(table(df_$Comparison_pair)[unlist(lapply(table(df_$Comparison_pair), function(x) x!=0))])/2
common_syn_genes=(df_comparison.pair %>% subset(No.Comparison.Pair>common_thr))[["Gene"]]
# it turns out to be very few common synchronized gene in B cell subtypes
# "HLA-DQB1" "MARCKS"   "SRGN"     "UBALD2" 
# some endothelial cell components and plasma membrane proteins

# enrich the syn.genes in each subtypes of B cells
subtype_compair=names(table(df_$Comparison_pair)[unlist(lapply(table(df_$Comparison_pair), function(x) x!=0))])
library(clusterProfiler)
library(org.Hs.eg.db)
geneList=list()
pdf("~/Project_PBMCage/Plots/CCF_synGenes_Blineage.pdf")
for (i in 1:length(subtype_compair)) {
  df_chosen=df_ %>%
    subset(Comparison_pair==subtype_compair[i])
  geneList[[i]]=df_chosen$Gene
  tryCatch({
    gse=enrichGO(df_chosen$Gene, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p0=
      dotplot(gse, showCategory=8, font.size=10) +
      labs(title=paste0("Synchronized genes between ", paste0(strsplit(subtype_compair[i], "\\.vs\\.")[[1]], collapse=" and "))) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8),
            title=element_text(size=10))
    plot(p0)
  }, error=function(msg) {print("No enriched terms.")})
  
  geneCluster_=list(df_chosen$Gene, global_syn_genes)
  names(geneCluster_)=c(subtype_compair[i], "Global")
  tryCatch({
    GO_enrich_compare=compareCluster(geneCluster=geneCluster_, fun=enrichGO, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p1=
      dotplot(GO_enrich_compare, showCategory=8, font.size=10) +
      labs(x=NULL) +
      theme(axis.text.x=element_text(size=8, hjust=1, angle=45),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9))
    plot(p1)
  }, error=function(msg) {print("No enriched terms (with global).")})
}
geneList[[i+1]]=global_syn_genes
names(geneList)=c(subtype_compair, "Global")
tryCatch({
  GO_enrich_compare_all=compareCluster(geneCluster=geneList, fun=enrichGO, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
  p2=
    dotplot(GO_enrich_compare_all, showCategory=8, font.size=10) +
    labs(x=NULL) +
    theme(axis.text.x=element_text(size=8, hjust=1, angle=45),
          axis.text.y=element_text(size=8),
          legend.title=element_text(size=9))
  plot(p2)
}, error=function(msg) {print("No enriched terms with all.")})
dev.off()
# Bnaive vs. Classical.Bmem has some enriched terms, while Bnaive vs. DN.Bmem does not have any.

# plot the expression of the syn.genes in each subtypes of B cells
# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
# obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
# table(rownames(obj_annotation)==colnames(object)) # check
# object=AddMetaData(object, metadata=obj_annotation$CD8T_NK_subsets, col.name="CD8T_NK_subsets")
# saveRDS(object, "~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
chose_from=names(table(object$CD8T_NK_subsets))
Celltype_groups=list(DCs=c(chose_from[grepl("DC$|DC[^CA]",chose_from)]),
                     B_cells=c(chose_from[grepl("Bnaive$|Bmem|PB$|PC$",chose_from)]),
                     Monocytes=c(chose_from[grepl("Mono$",chose_from)]),
                     CD4_T=c(chose_from[grepl("^CD4|ISAGhi T|^Treg",chose_from)]),
                     CD8_T=c(chose_from[grepl("^CD8",chose_from)]),
                     ILCs=c(chose_from[grepl("^NK\\.|ILC",chose_from)]),
                     NKTs=c(chose_from[grepl("^(NKT)",chose_from)]),
                     Other_lymph=c(chose_from[grepl("CLP|dnT|gdT|MAIT",chose_from)]),
                     Other_myelo=c(chose_from[grepl("EMP|^Eryth|Mast|Platelet",chose_from)]))
cateogry_name="B_cells"
celltype_group=Celltype_groups$"B_cells"
obj_sub=subset(object, CD8T_NK_subsets %in% celltype_group)
pb_obj=edgeR::Seurat2PB(obj_sub, sample="age", cluster="CD8T_NK_subsets")
pseudobulkdata=pb_obj$counts

subtype_to_keep=as.vector(sapply(subtype_compair, function(x) strsplit(x, "\\.vs\\.")[[1]]))
subtype_to_keep=subtype_to_keep[!duplicated(subtype_to_keep)]
subtype_to_keep_for_grepl=paste0("(cluster",subtype_to_keep,")"); subtype_to_keep_for_grepl=gsub("\\."," ",subtype_to_keep_for_grepl)
subtype_to_keep_for_grepl=paste0(subtype_to_keep_for_grepl, collapse="|")
col.keep=grepl(subtype_to_keep_for_grepl, colnames(pseudobulkdata))
pseudobulkdata=pseudobulkdata[,col.keep]

gene_to_plot=df_$Gene[!duplicated(df_$Gene)]
pseudobulkdata_filtered=pseudobulkdata[gene_to_plot,]

df_syn.gene.expr=as.data.frame(t(pseudobulkdata_filtered))
df_syn.gene.expr_age=gsub("X|_cluster.*","",rownames(df_syn.gene.expr))
df_syn.gene.expr_celltype=gsub(".*cluster","",rownames(df_syn.gene.expr))
df_syn.gene.expr=cbind(Age=as.factor(df_syn.gene.expr_age), Celltype=as.factor(df_syn.gene.expr_celltype), df_syn.gene.expr)
long_df_syn.gene=reshape2::melt(df_syn.gene.expr, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")

library(ggplot2)
ggplot(long_df_syn.gene, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
  geom_point(size=0.5, shape=1, alpha=0.1) +
  scale_color_manual(values=ggpubr::get_palette("locuszoom", 3)) +
  scale_y_log10(labels=scales::label_log(digits=1)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19, alpha=1))) +
  labs(title="Synchronized genes in the B cell lineage",
       x="Age") +
  geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)

# enrich the syn.genes from all the B subtypes comparison pairs altogether
gse=enrichGO(gene_to_plot, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
p0=
  dotplot(gse, showCategory=8, font.size=10) +
  labs(title="Synchronized genes in the B cell lineage") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10))
plot(p0)

# Analyze gene-gene CCF of the syn. genes in each B-lineage comparison pair
length(gene_to_plot)>300 # too many, so set a expr cutoff:
max.expr_thr=1000
table(apply(pseudobulkdata_filtered, 1, max)>max.expr_thr) # check
gene_to_plot_chosen=rownames(pseudobulkdata_filtered)[apply(pseudobulkdata_filtered, 1, max)>max.expr_thr]

pseudobulk_global=read.delim("~/Project_PBMCage/pseudobulk_global.txt", stringsAsFactors=F)
pseudobulkdata_all_filtered=cbind(pseudobulkdata_filtered[gene_to_plot_chosen,], pseudobulk_global[gene_to_plot_chosen,])

Celltype=gsub(".*_cluster","",colnames(pseudobulkdata_all_filtered))
all_celltypes=names(table(Celltype))

gene_compair_pair_list=c()
for (i in 1:nrow(pseudobulkdata_all_filtered)) {
  rest=c(1:nrow(pseudobulkdata_all_filtered))[-i]
  gene_compair_pair=lapply(1:length(rest), function(x) c(gene_to_plot_chosen[i], gene_to_plot_chosen[rest[x]]))
  gene_compair_pair_list=c(gene_compair_pair_list, gene_compair_pair)
}

synGenes_correlation_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))
synGenes_h_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))

library(TSA)
processbar=txtProgressBar(min=0, max=length(gene_compair_pair_list), style=3, width=60, char="=")
for (i in 1:length(gene_compair_pair_list)) {
  gene1=gene_compair_pair_list[[i]][1]
  gene2=gene_compair_pair_list[[i]][2]
  for (j in 1:length(all_celltypes)) {
    celltype_name=all_celltypes[j]
    celltype.df=pseudobulkdata_all_filtered[, grepl(celltype_name, colnames(pseudobulkdata_all_filtered))]
    celltype.df_filtered=celltype.df %>% subset(rownames(.) %in% c(gene1,gene2))
    celltype.df_filtered=as.data.frame(t(celltype.df_filtered))
    celltype.df_filtered=cbind(Age=gsub("X|_.*","",rownames(celltype.df_filtered)), celltype.df_filtered)
    series1=celltype.df_filtered[,2]
    series2=celltype.df_filtered[,3]
    results=prewhiten(series1, series2)
    Abs.correlation=max(abs(results$ccf$acf))
    h=results$ccf$lag[which(abs(results$ccf$acf)==Abs.correlation)]
    correlation=results$ccf$acf[which(results$ccf$lag==h)]
    
    synGenes_correlation_mtx[i, j]=correlation
    synGenes_h_mtx[i, j]=h
  }
  setTxtProgressBar(processbar, i)
}
close(processbar)

synGenes_correlation_mtx=as.data.frame(synGenes_correlation_mtx)
rownames(synGenes_correlation_mtx)=unlist(lapply(gene_compair_pair_list, function(x) paste0(x[1],".vs.",x[2])))
colnames(synGenes_correlation_mtx)=all_celltypes
synGenes_h_mtx=as.data.frame(synGenes_h_mtx)
rownames(synGenes_h_mtx)=rownames(synGenes_correlation_mtx)
colnames(synGenes_h_mtx)=colnames(synGenes_correlation_mtx)

saveRDS(synGenes_correlation_mtx, "~/Project_PBMCage/Results/correlationMatrix_synGenes_Bsubtype.rds")
saveRDS(synGenes_h_mtx, "~/Project_PBMCage/Results/hMatrix_synGenes_Bsubtype.rds")

# Analyze the gene-gene lags among globally syn. genes
synGenes_correlation_mtx$Comparison_pair=rownames(synGenes_correlation_mtx)
long_df_synGene_matrix=reshape2::melt(synGenes_correlation_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="Correlation")
synGenes_h_mtx$Comparison_pair=rownames(synGenes_h_mtx)
long_df_synGene_matrix_h=reshape2::melt(synGenes_h_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="h")
long_df_synGene_matrix=cbind(long_df_synGene_matrix, h=long_df_synGene_matrix_h$h)

ggplot(long_df_synGene_matrix, aes(x=Correlation, color=Celltype)) +
  geom_density(lwd=1, alpha=0) +
  scale_color_manual(values=scater:::.get_palette("tableau20"),
                     labels=c(all_celltypes[1:(length(all_celltypes)-1)], "All")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19))) +
  labs(title="Synchronicity among sync. genes in the B cell lineage")
# determine which of them are leading/lagging
range(long_df_synGene_matrix$h)
range(long_df_synGene_matrix$Correlation)
long_df_synGene_matrix=long_df_synGene_matrix %>%
  mutate(Leading=gsub("\\.vs.*","",Comparison_pair)) %>%
  mutate(Lagging=gsub(".*vs\\.","",Comparison_pair)) %>%
  subset(Celltype!="SeuratProject")
leading_strength=long_df_synGene_matrix %>% count(h, Leading)
syn.Genes=(leading_strength %>% subset(h==0 & n>length(gene_compair_pair_list)/length(gene_to_plot_chosen)*3/2))$Leading
leading.Genes=(leading_strength %>% subset(h>0 & (Leading %in% syn.Genes)==FALSE))$Leading
lagging.Genes=(leading_strength %>% subset(h<0 & (Leading %in% syn.Genes)==FALSE))$Leading
# it seems that none of them are leading or lagging, but all expressed simultaneously

### Determine the B-subtype-dependent leading/lagging genes
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h>h_thr | h<(-h_thr))

gene_cluster=split(df_$Gene, paste0(df_$Comparison_pair, "__", df_$h))
gene_cluster_sort=gene_cluster[order(abs(as.numeric(gsub(".*__","",names(gene_cluster)))), decreasing=T)]

pdf("~/Project_PBMCage/Plots/Correlation_H_Leading.Lagging.Genes_Expr_Bcell.pdf")
for (i in 1:length(gene_cluster_sort)) {
  cluster_=gene_cluster_sort[i]
  h_value=gsub(".*__","",names(cluster_))
  type1=gsub("\\."," ",gsub("\\.vs.*","",names(cluster_)))
  type2=gsub("\\."," ",gsub(".*vs\\.|__.*","",names(cluster_)))
  cluster=cluster_[[1]]
  pseudobulkdata_filtered=pseudobulkdata[cluster,]
  
  df_byCluster.Gene=cbind(pseudobulkdata[cluster,], pseudobulk_global[cluster,])
  df_byCluster.Gene_t_age=gsub("X|_cluster.*","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t_celltype=gsub(".*cluster","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t=data.frame(Age=as.factor(df_byCluster.Gene_t_age), Celltype=as.factor(df_byCluster.Gene_t_celltype), as.data.frame(t(df_byCluster.Gene)))
  df_byCluster.Gene_t_chosen=df_byCluster.Gene_t %>%
    subset(Celltype %in% c(type1, type2))
  long_df_leading.lagging.gene.expr=reshape2::melt(df_byCluster.Gene_t_chosen, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
  
  cluster_name=ifelse(length(cluster)>5, paste0(paste0(cluster[1:5], collapse=", "),", ..."), paste0(cluster, collapse=", "))
  p_x=
    ggplot(long_df_leading.lagging.gene.expr, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
    geom_point(shape=1, alpha=0.3) +
    scale_color_manual(values=scater:::.get_palette("tableau20")[c(1,11)]) +
    scale_y_continuous(trans="pseudo_log") +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=9),
          legend.text=element_text(size=8),
          title=element_text(size=10)) + 
    guides(color=guide_legend(override.aes=list(size=2, alpha=1))) +
    labs(title=paste0(type1, " vs. ", type2, " (", cluster_name, "; h=", h_value, ")"),
         x="Age") +
    geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)
  
  plot(p_x)
  
  tryCatch({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    gse=enrichGO(cluster, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p0=
      dotplot(gse, showCategory=8, font.size=10) +
      labs(title=ifelse(h_value<0, paste0("Genes expressed earlier in ",type1, " than ", type2), paste0("Genes expressed earlier in ",type2, " than ", type1))) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8),
            title=element_text(size=10))
    plot(p0)
  }, error=function(msg) {print("No term enriched.")})
}
dev.off()
# we conclude that the genes related to proteolysis/neg.reg.phosphorylation/neg.reg.metabolic process/phagocytosis are expressed to the peak more quickly in DN Bmem than in Bnaive

### Analyze negatively regulated genes
corr_thr=0.5

df_=long_exp.matrix %>%
  subset(Correlation<=-corr_thr)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)
table(df_comparison.pair$No.Comparison.Pair)
split(df_$Comparison_pair, df_$Gene)
# it turns out that only one gene (CCNT1) was neg.reg in Classical.Bmem vs. DN.Bmem, which's not very meaningful

### Analyze less correlated genes
corr_thr=0.25

df_=long_exp.matrix %>%
  subset(abs(Correlation)<corr_thr)
# Determine the less correlated genes in various celltypes
library(clusterProfiler)
library(org.Hs.eg.db)
geneList=split(df_$Gene, df_$Comparison_pair)
geneList=geneList[sapply(geneList, function(x) length(x)!=0)]
GO_enrich_compare=compareCluster(geneCluster=geneList, fun=enrichGO, OrgDb="org.Hs.eg.db", pvalueCutoff=0.05)
# no term enriched, which means that the less correlated genes in each celltype comparison pair are just random genes

#####################################





###
### CCF analysis on CD8 T cell lineage
#####################################
###

### Analyze the Corr and Lag
# Arrange correlation and h matrix for plotting
correlation_mtx=readRDS("~/Project_PBMCage/Results/correlationMatrix_CD8_T.rds")
h_mtx=readRDS("~/Project_PBMCage/Results/hMatrix_CD8_T.rds")

chose_from=names(table(object$CD8T_NK_subsets))
Celltype_groups=list(DCs=c(chose_from[grepl("DC$|DC[^CA]",chose_from)]),
                     B_cells=c(chose_from[grepl("Bnaive$|Bmem|PB$|PC$",chose_from)]),
                     Monocytes=c(chose_from[grepl("Mono$",chose_from)]),
                     CD4_T=c(chose_from[grepl("^CD4|ISAGhi T|^Treg",chose_from)]),
                     CD8_T=c(chose_from[grepl("^CD8",chose_from)]),
                     ILCs=c(chose_from[grepl("^NK\\.|ILC",chose_from)]),
                     NKTs=c(chose_from[grepl("^(NKT)",chose_from)]),
                     Other_lymph=c(chose_from[grepl("CLP|dnT|gdT|MAIT",chose_from)]),
                     Other_myelo=c(chose_from[grepl("EMP|^Eryth|Mast|Platelet",chose_from)]))
all_celltypes=Celltype_groups$CD8_T
compair_pair_list=c()
for (i in 1:length(all_celltypes)) {
  rest=c(1:length(all_celltypes))[-i]
  compair_pair=lapply(1:length(rest), function(x) c(all_celltypes[i], all_celltypes[rest[x]]))
  compair_pair_list=c(compair_pair_list, compair_pair)
}

colnames(correlation_mtx)=sapply(compair_pair_list, function(x) paste0(x[1], ".vs.", x[2]))
colnames(h_mtx)=sapply(compair_pair_list, function(x) paste0(x[1], ".vs.", x[2]))
compare_idx=list()
for (idx in 1:length(colnames(correlation_mtx))) {
  pair_1=colnames(correlation_mtx)[idx]
  pair_2=paste0(gsub(".*vs\\.","",pair_1), ".vs.", gsub("\\.vs.*","",pair_1))
  compare_idx1=c(1:length(colnames(correlation_mtx)))[idx]
  compare_idx2=which(colnames(correlation_mtx)==pair_2)
  if ((compare_idx2 %in% c(1:idx))==FALSE) {
    compare_idx=c(compare_idx, list(c(compare_idx1, compare_idx2)))
  }
}
idx_1_all=c()
for (i in 1:length(compare_idx)) {
  idx1=compare_idx[[i]][1]
  idx2=compare_idx[[i]][2]
  correlation_updated=sapply(1:length(correlation_mtx[,idx1]), function(x) ifelse(abs(correlation_mtx[,idx1][x])>=abs(correlation_mtx[,idx2][x]), correlation_mtx[,idx1][x], correlation_mtx[,idx2][x]))
  col_chosen=sapply(1:length(correlation_updated), function(x) ifelse(correlation_updated[x]==correlation_mtx[,idx1][x], 1, 2))
  h_updated=sapply(1:length(h_mtx[,idx1]), function(x) ifelse(col_chosen[x]==1, h_mtx[,idx1], h_mtx[,idx2]))
  
  correlation_mtx[,idx1]=correlation_updated
  h_mtx[,idx1]=h_updated
  idx_1_all=c(idx_1_all,idx1)
}
correlation_mtx=correlation_mtx[,idx_1_all]
h_mtx=h_mtx[,idx_1_all]

correlation_mtx$Gene=rownames(correlation_mtx)
long_exp.matrix.cor=reshape2::melt(correlation_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="Correlation")
h_mtx$Gene=rownames(h_mtx)
long_exp.matrix.h=reshape2::melt(h_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="h")

# remove the globally lowly expressed genes and age-stable genes
lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_perCelltype.rds"); names(lowlyExprGenes)
celltype.wide_lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_global.rds")
StableGenes=readRDS("~/Project_PBMCage/Tempt_RDS/StableGenes_perCelltype.rds"); names(StableGenes)
global_syn_genes=readRDS("~/Project_PBMCage/Tempt_RDS/Global_synchronized_genes.rds")

gene.to.remove=c(lowlyExprGenes$"CD8.Tmem",lowlyExprGenes$"CD8.Tnaive",lowlyExprGenes$"Proliferating.NK.T",
                 StableGenes$"CD8.Tmem",StableGenes$"CD8.Tnaive",StableGenes$"Proliferating.NK.T",
                 celltype.wide_lowlyExprGenes,
                 global_syn_genes)
gene.to.remove=gene.to.remove[!duplicated(gene.to.remove)]

long_exp.matrix.cor=long_exp.matrix.cor %>%
  subset((Gene %in% gene.to.remove)==FALSE)
long_exp.matrix.h=long_exp.matrix.h %>%
  subset((Gene %in% gene.to.remove)==FALSE)

long_exp.matrix=long_exp.matrix.cor
long_exp.matrix$h=long_exp.matrix.h$h
long_exp.matrix$abs_corr=abs(long_exp.matrix$Correlation)
long_exp.matrix=na.omit(long_exp.matrix)

# Plot distribution of correlation
p_correlation=
  ggplot(long_exp.matrix, aes(x=Correlation, color=Comparison_pair)) +
  geom_density(lwd=1, fill="transparent", alpha=0.25) +
  theme_classic() +
  theme(legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm")) +
  guides(color=guide_legend(ncol=1))
p_correlation
# it seems that the correlation between CD8 Tem.GZMB+ FCGR3A+ and CD8 Tem.GZMK+ NKG7+ is > that between any of the other comparison pairs

### Plot scatterplot of correlation vs. h
ggplot(long_exp.matrix, aes(x=Correlation, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag - CD8 T cells")
ggplot(long_exp.matrix, aes(x=abs_corr, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag - CD8 T cells",
       x="Abs(Correlation)")
# there's a portion of genes with high abs(lag) and high positive corr, meaning that there are leading/lagging genes

### CD8 T-specific synchronized genes (defined as: a large portion of large abs positive correlation + small abs(h))
# determine the threshold for "small correlation" and "small h" based on the plot
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h<=h_thr & h>=-h_thr)

df_ %>% group_by(Comparison_pair) %>% slice_max(abs_corr, n=5)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)
table(as.vector(table(df_comparison.pair$Gene)))

# determine the CD8 T common syn. genes
common_thr=length(table(df_$Comparison_pair)[unlist(lapply(table(df_$Comparison_pair), function(x) x!=0))])/2
common_syn_genes=(df_comparison.pair %>% subset(No.Comparison.Pair>common_thr))[["Gene"]]
library(clusterProfiler)
library(org.Hs.eg.db)
gse=enrichGO(common_syn_genes, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
p0=
  dotplot(gse, showCategory=8, font.size=10) +
  labs(title="Common synchronized genes in CD8 T cells") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10))
plot(p0)
# it seems that the genes that related to genome replication and RNA splicing are commonly synchronized in various CD8 T cell subtypes

# enrich the subtype-specific syn.genes
library(clusterProfiler)
library(org.Hs.eg.db)
all_genes_to_enrich=split(df_$Gene,df_$Comparison_pair)[sapply(split(df_$Gene,df_$Comparison_pair), function(x) length(x)!=0)]
subtype_compair=names(all_genes_to_enrich)
geneList=list()
pdf("~/Project_PBMCage/Plots/CCF_synGenes_CD8_T.pdf")
for (i in 1:length(subtype_compair)) {
  rest_genes=unlist(all_genes_to_enrich[-i])[!duplicated(unlist(all_genes_to_enrich[-i]))]; rest_genes=as.vector(rest_genes)
  df_chosen=df_ %>%
    subset(Gene %in% all_genes_to_enrich[[i]]) %>%
    subset((Gene %in% rest_genes)==FALSE)
  geneList[[i]]=df_chosen$Gene
  tryCatch({
    gse=enrichGO(df_chosen$Gene, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p0=
      dotplot(gse, showCategory=8, font.size=10) +
      labs(title=paste0("Synchronized genes between \n", paste0(strsplit(subtype_compair[i], "\\.vs\\.")[[1]], collapse=" and "))) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8),
            title=element_text(size=10))
    plot(p0)
  }, error=function(msg) {print("No enriched terms.")})
  
  geneCluster_=list(df_chosen$Gene, global_syn_genes)
  names(geneCluster_)=c(subtype_compair[i], "Global")
  tryCatch({
    GO_enrich_compare=compareCluster(geneCluster=geneCluster_, fun=enrichGO, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p1=
      dotplot(GO_enrich_compare, showCategory=8, font.size=10) +
      labs(x=NULL) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9))
    plot(p1)
  }, error=function(msg) {print("No enriched terms (with global).")})
}
geneList[[i+1]]=common_syn_genes
geneList[[i+2]]=global_syn_genes
names(geneList)=c(subtype_compair, "Common", "Global")
geneList_deleted=geneList[names(geneList)!="CD8 Tem.GZMK+ NKG7+.vs.CD8 Tnaive.TNFAIP3+"]
tryCatch({
  GO_enrich_compare_all=compareCluster(geneCluster=geneList, fun=enrichGO, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
  p2=
    dotplot(GO_enrich_compare_all, showCategory=3, font.size=10) +
    labs(x=NULL) +
    theme(axis.text.x=element_text(size=8, hjust=1, angle=45),
          axis.text.y=element_text(size=8),
          legend.title=element_text(size=9))
  plot(p2)
}, error=function(msg) {print("No enriched terms with all.")})
dev.off()
# it seems that many comparison pairs have syn.genes that are related to biosynthesis, or response to ionizing radiation, or transport,
# or rRNA metabolism, or tRNA metabolism, apart from the common one (genome replication/mRNA splicing) among all the comparison pairs;
# CD8 Tem.GZMB+ FCGR3A+.vs.CD8 Tem.GZMK+ NKG7+ shows the most shared syn.genes, which are related to tRNA
# so we conclude that some of the CD8 T cell subtypes are expressing genes with certain functions simultaneously, but each subtype
# has its own preference except for the genome replication/mRNA splicing

# plot the expression of the syn.genes in each subtypes of CD8 T subtypes
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
chose_from=names(table(object$CD8T_NK_subsets))
Celltype_groups=list(DCs=c(chose_from[grepl("DC$|DC[^CA]",chose_from)]),
                     B_cells=c(chose_from[grepl("Bnaive$|Bmem|PB$|PC$",chose_from)]),
                     Monocytes=c(chose_from[grepl("Mono$",chose_from)]),
                     CD4_T=c(chose_from[grepl("^CD4|ISAGhi T|^Treg",chose_from)]),
                     CD8_T=c(chose_from[grepl("^CD8",chose_from)]),
                     ILCs=c(chose_from[grepl("^NK\\.|ILC",chose_from)]),
                     NKTs=c(chose_from[grepl("^(NKT)",chose_from)]),
                     Other_lymph=c(chose_from[grepl("CLP|dnT|gdT|MAIT",chose_from)]),
                     Other_myelo=c(chose_from[grepl("EMP|^Eryth|Mast|Platelet",chose_from)]))
cateogry_name="CD8_T"
celltype_group=Celltype_groups$"CD8_T"
obj_sub=subset(object, CD8T_NK_subsets %in% celltype_group)
pb_obj=edgeR::Seurat2PB(obj_sub, sample="age", cluster="CD8T_NK_subsets")
pseudobulkdata=pb_obj$counts

subtype_to_keep=as.vector(sapply(subtype_compair, function(x) strsplit(x, "\\.vs\\.")[[1]]))
subtype_to_keep=subtype_to_keep[!duplicated(subtype_to_keep)]
subtype_to_keep_for_grepl=paste0("(",gsub("\\+","\\\\+",gsub("\\.","\\\\.",subtype_to_keep)),")")
subtype_to_keep_for_grepl=paste0(subtype_to_keep_for_grepl, collapse="|")
col.keep=grepl(subtype_to_keep_for_grepl, colnames(pseudobulkdata))
pseudobulkdata=pseudobulkdata[,col.keep]

gene_to_plot=df_$Gene[!duplicated(df_$Gene)]
pseudobulkdata_filtered=pseudobulkdata[gene_to_plot,]

df_syn.gene.expr=as.data.frame(t(pseudobulkdata_filtered))
df_syn.gene.expr_age=gsub("X|_cluster.*","",rownames(df_syn.gene.expr))
df_syn.gene.expr_celltype=gsub(".*cluster","",rownames(df_syn.gene.expr))
df_syn.gene.expr=cbind(Age=as.factor(df_syn.gene.expr_age), Celltype=as.factor(df_syn.gene.expr_celltype), df_syn.gene.expr)

long_df_syn.gene=reshape2::melt(df_syn.gene.expr, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
dim(long_df_syn.gene)
# too many genes, so take the highly expressed ones for plotting
df_syn.gene.expr_thr=long_df_syn.gene %>% group_by(Celltype, Gene) %>%
  summarise(across(Expression, mean)) %>%
  group_by(Celltype) %>%
  slice_max(order_by=Expression, n=500) %>%
  as.data.frame()
genes_chosen=df_syn.gene.expr_thr$Gene[!duplicated(df_syn.gene.expr_thr$Gene)]
pseudobulkdata_filtered=pseudobulkdata[genes_chosen,]
df_syn.gene.expr=as.data.frame(t(pseudobulkdata_filtered))
df_syn.gene.expr_age=gsub("X|_cluster.*","",rownames(df_syn.gene.expr))
df_syn.gene.expr_celltype=gsub(".*cluster","",rownames(df_syn.gene.expr))
df_syn.gene.expr=cbind(Age=as.factor(df_syn.gene.expr_age), Celltype=as.factor(df_syn.gene.expr_celltype), df_syn.gene.expr)
long_df_syn.gene=reshape2::melt(df_syn.gene.expr, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
dim(long_df_syn.gene)

library(ggplot2)
p_e=
  ggplot(long_df_syn.gene, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
  geom_point(size=0.5, shape=1, alpha=0.1) +
  scale_color_manual(values=ggpubr::get_palette("locuszoom", 10)) +
  scale_y_log10(labels=scales::label_log(digits=1)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19, alpha=1))) +
  labs(title="Synchronized genes in all the subtypes of CD8 T cells",
       x="Age") +
  geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)
plot(p_e)

# enrich the syn.genes from all the CD8 T subtype comparison pairs altogether
gse=enrichGO(gene_to_plot, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
p0=
  dotplot(gse, showCategory=8, font.size=10) +
  labs(title="Synchronized genes in all the subtypes of CD8 T cells") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10))
plot(p0)
# similar to the common syn.genes among the CD8 T subtypes (RNA splicing/processing)

# Analyze gene-gene CCF of the syn. genes in each CD8 T subtypes comparison pair
length(gene_to_plot)>300 # too many, so set a expr cutoff:
max.expr_thr=100
table(apply(pseudobulkdata_filtered, 1, max)>max.expr_thr) # check
gene_to_plot_chosen=rownames(pseudobulkdata_filtered)[apply(pseudobulkdata_filtered, 1, max)>max.expr_thr]

pseudobulk_global=read.delim("~/Project_PBMCage/pseudobulk_global.txt", stringsAsFactors=F)
pseudobulkdata_all_filtered=cbind(pseudobulkdata_filtered[gene_to_plot_chosen,], pseudobulk_global[gene_to_plot_chosen,])

Celltype=gsub(".*_cluster","",colnames(pseudobulkdata_all_filtered))
all_celltypes=names(table(Celltype))

gene_compair_pair_list=c()
for (i in 1:nrow(pseudobulkdata_all_filtered)) {
  rest=c(1:nrow(pseudobulkdata_all_filtered))[-i]
  gene_compair_pair=lapply(1:length(rest), function(x) c(gene_to_plot_chosen[i], gene_to_plot_chosen[rest[x]]))
  gene_compair_pair_list=c(gene_compair_pair_list, gene_compair_pair)
}

synGenes_correlation_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))
synGenes_h_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))

library(TSA)
processbar=txtProgressBar(min=0, max=length(gene_compair_pair_list), style=3, width=60, char="=")
for (i in 1:length(gene_compair_pair_list)) {
  gene1=gene_compair_pair_list[[i]][1]
  gene2=gene_compair_pair_list[[i]][2]
  for (j in 1:length(all_celltypes)) {
    celltype_name=all_celltypes[j]
    celltype.df=pseudobulkdata_all_filtered[, grepl(gsub("\\+","\\\\+",gsub("\\.","\\\\.",celltype_name)), colnames(pseudobulkdata_all_filtered))]
    celltype.df_filtered=celltype.df %>% subset(rownames(.) %in% c(gene1,gene2))
    celltype.df_filtered=as.data.frame(t(celltype.df_filtered))
    celltype.df_filtered=cbind(Age=gsub("X|_.*","",rownames(celltype.df_filtered)), celltype.df_filtered)
    series1=celltype.df_filtered[,2]
    series2=celltype.df_filtered[,3]
    results=prewhiten(series1, series2)
    Abs.correlation=max(abs(results$ccf$acf))
    h=results$ccf$lag[which(abs(results$ccf$acf)==Abs.correlation)]
    correlation=results$ccf$acf[which(results$ccf$lag==h)]
    
    synGenes_correlation_mtx[i, j]=correlation
    synGenes_h_mtx[i, j]=h
  }
  setTxtProgressBar(processbar, i)
}
close(processbar)

synGenes_correlation_mtx=as.data.frame(synGenes_correlation_mtx)
rownames(synGenes_correlation_mtx)=unlist(lapply(gene_compair_pair_list, function(x) paste0(x[1],".vs.",x[2])))
colnames(synGenes_correlation_mtx)=all_celltypes
synGenes_h_mtx=as.data.frame(synGenes_h_mtx)
rownames(synGenes_h_mtx)=rownames(synGenes_correlation_mtx)
colnames(synGenes_h_mtx)=colnames(synGenes_correlation_mtx)

saveRDS(synGenes_correlation_mtx, "~/Project_PBMCage/Results/correlationMatrix_synGenes_CD8Tsubtype.rds")
saveRDS(synGenes_h_mtx, "~/Project_PBMCage/Results/hMatrix_synGenes_CD8Tsubtype.rds")

# Analyze the gene-gene lags among the syn. genes in all the CD8 T subtypes
synGenes_correlation_mtx$Comparison_pair=rownames(synGenes_correlation_mtx)
long_df_synGene_matrix=reshape2::melt(synGenes_correlation_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="Correlation")
synGenes_h_mtx$Comparison_pair=rownames(synGenes_h_mtx)
long_df_synGene_matrix_h=reshape2::melt(synGenes_h_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="h")
long_df_synGene_matrix=cbind(long_df_synGene_matrix, h=long_df_synGene_matrix_h$h)

ggplot(long_df_synGene_matrix, aes(x=Correlation, color=Celltype)) +
  geom_density(lwd=1, alpha=0) +
  scale_color_manual(values=scater:::.get_palette("tableau20"),
                     labels=c(all_celltypes[1:(length(all_celltypes)-1)], "All")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19))) +
  labs(title="Synchronicity among sync. genes in CD8 T cells")
# determine which of them are leading/lagging
range(long_df_synGene_matrix$h)
range(long_df_synGene_matrix$Correlation)
long_df_synGene_matrix=long_df_synGene_matrix %>%
  mutate(Leading=gsub("\\.vs.*","",Comparison_pair)) %>%
  mutate(Lagging=gsub(".*vs\\.","",Comparison_pair)) %>%
  subset(Celltype!="SeuratProject")
leading_strength=long_df_synGene_matrix %>% count(h, Leading)
syn.Genes=(leading_strength %>% subset(h==0 & n>length(gene_compair_pair_list)/length(gene_to_plot_chosen)*length(all_celltypes)/2))$Leading
leading.Genes=(leading_strength %>% subset(h<0 & (Leading %in% syn.Genes)==FALSE))$Leading
lagging.Genes=(leading_strength %>% subset(h>0 & (Leading %in% syn.Genes)==FALSE))$Leading
# it seems that none of them are leading or lagging, but all expressed simultaneously

### Determine the CD8T-subtype-dependent leading/lagging genes
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h>h_thr | h<(-h_thr))

gene_cluster=split(df_$Gene, paste0(df_$Comparison_pair, "__", df_$h))
gene_cluster_sort=gene_cluster[order(abs(as.numeric(gsub(".*__","",names(gene_cluster)))), decreasing=T)]

pdf("~/Project_PBMCage/Plots/Correlation_H_Leading.Lagging.Genes_Expr_CD8Tcell.pdf")
for (i in 1:length(gene_cluster_sort)) {
  cluster_=gene_cluster_sort[i]
  h_value=gsub(".*__","",names(cluster_))
  type1_and_2=strsplit(names(cluster_),split=".vs.")[[1]]
  type1=type1_and_2[1]
  type2=gsub("__.*","",type1_and_2[2])
  cluster=cluster_[[1]]
  pseudobulkdata_filtered=pseudobulkdata[cluster,]
  
  df_byCluster.Gene=cbind(pseudobulkdata[cluster,], pseudobulk_global[cluster,])
  df_byCluster.Gene_t_age=gsub("X|_cluster.*","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t_celltype=gsub(".*cluster","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t=data.frame(Age=as.factor(df_byCluster.Gene_t_age), Celltype=as.factor(df_byCluster.Gene_t_celltype), as.data.frame(t(df_byCluster.Gene)))
  df_byCluster.Gene_t_chosen=df_byCluster.Gene_t %>%
    subset(Celltype %in% c(type1, type2))
  long_df_leading.lagging.gene.expr=reshape2::melt(df_byCluster.Gene_t_chosen, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
  
  cluster_name=ifelse(length(cluster)>5, paste0(paste0(cluster[1:5], collapse=", "),", ..."), paste0(cluster, collapse=", "))
  p_x=
    ggplot(long_df_leading.lagging.gene.expr, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
    geom_point(shape=1, alpha=0.3) +
    scale_color_manual(values=scater:::.get_palette("tableau20")[c(1,11)]) +
    scale_y_continuous(trans="pseudo_log") +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=9),
          legend.text=element_text(size=8),
          title=element_text(size=10)) + 
    guides(color=guide_legend(override.aes=list(size=2, alpha=1))) +
    labs(title=paste0(type1, " vs. ", type2, "\n(", cluster_name, "; h=", h_value, ")"),
         x="Age") +
    geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)
  
  plot(p_x)
  
  tryCatch({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    gse=enrichGO(cluster, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p0=
      dotplot(gse, showCategory=8, font.size=10) +
      labs(title=ifelse(h_value<0, paste0("Genes expressed earlier in ",type1, "\nthan in ", type2), paste0("Genes expressed earlier in ",type2, "\nthan in ", type1))) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8),
            title=element_text(size=10))
    plot(p0)
  }, error=function(msg) {print("No term enriched.")})
}
dev.off()
# the genes with certain functions are expressed earlier in a subtype than in another, meaning that the functions are important and working properly in that celltype at young ages,
# and diminish more quickly during aging

### Analyze negatively regulated genes
corr_thr=0.5

df_=long_exp.matrix %>%
  subset(Correlation<=-corr_thr)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)
table(df_comparison.pair$No.Comparison.Pair)
split(df_$Comparison_pair, df_$Gene)
# it turns out that: APPL2, CCDC32, SLC38A6 are neg.reg in CD8 Tcm.LGALS1+.vs.CD8 Tem.GZMK+ NKG7+;
# PPP2R5D in CD8 Tcm.LGALS1+.vs.CD8 Tem.GZMB+ FCGR3A+;
# RETSAT in CD8 Tem.GZMKhi.vs.CD8 Tnaive.TNFAIP3+

### Analyze less correlated genes
corr_thr=0.25

df_=long_exp.matrix %>%
  subset(abs(Correlation)<corr_thr)
# Determine the less correlated genes in various celltypes
library(clusterProfiler)
library(org.Hs.eg.db)
geneList=split(df_$Gene, df_$Comparison_pair)
geneList=geneList[sapply(geneList, function(x) length(x)!=0)]
GO_enrich_compare=compareCluster(geneCluster=geneList, fun=enrichGO, OrgDb="org.Hs.eg.db", pvalueCutoff=0.05)
# no term enriched, which means that the less correlated genes in each celltype comparison pair are just random genes

#####################################





###
### CCF analysis on NK cell lineage
#####################################
###

### Analyze the Corr and Lag
# Arrange correlation and h matrix for plotting
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
correlation_mtx=readRDS("~/Project_PBMCage/Results/correlationMatrix_ILCs.rds")
h_mtx=readRDS("~/Project_PBMCage/Results/hMatrix_ILCs.rds")

chose_from=names(table(object$CD8T_NK_subsets))
Celltype_groups=list(DCs=c(chose_from[grepl("DC$|DC[^CA]",chose_from)]),
                     B_cells=c(chose_from[grepl("Bnaive$|Bmem|PB$|PC$",chose_from)]),
                     Monocytes=c(chose_from[grepl("Mono$",chose_from)]),
                     CD4_T=c(chose_from[grepl("^CD4|ISAGhi T|^Treg",chose_from)]),
                     CD8_T=c(chose_from[grepl("^CD8",chose_from)]),
                     ILCs=c(chose_from[grepl("^NK\\.|ILC",chose_from)]),
                     NKTs=c(chose_from[grepl("^(NKT)",chose_from)]),
                     Other_lymph=c(chose_from[grepl("CLP|dnT|gdT|MAIT",chose_from)]),
                     Other_myelo=c(chose_from[grepl("EMP|^Eryth|Mast|Platelet",chose_from)]))
all_celltypes=Celltype_groups$ILCs
compair_pair_list=c()
for (i in 1:length(all_celltypes)) {
  rest=c(1:length(all_celltypes))[-i]
  compair_pair=lapply(1:length(rest), function(x) c(all_celltypes[i], all_celltypes[rest[x]]))
  compair_pair_list=c(compair_pair_list, compair_pair)
}

colnames(correlation_mtx)=sapply(compair_pair_list, function(x) paste0(x[1], ".vs.", x[2]))
colnames(h_mtx)=sapply(compair_pair_list, function(x) paste0(x[1], ".vs.", x[2]))
compare_idx=list()
for (idx in 1:length(colnames(correlation_mtx))) {
  pair_1=colnames(correlation_mtx)[idx]
  pair_2=paste0(gsub(".*vs\\.","",pair_1), ".vs.", gsub("\\.vs.*","",pair_1))
  compare_idx1=c(1:length(colnames(correlation_mtx)))[idx]
  compare_idx2=which(colnames(correlation_mtx)==pair_2)
  if ((compare_idx2 %in% c(1:idx))==FALSE) {
    compare_idx=c(compare_idx, list(c(compare_idx1, compare_idx2)))
  }
}
idx_1_all=c()
for (i in 1:length(compare_idx)) {
  idx1=compare_idx[[i]][1]
  idx2=compare_idx[[i]][2]
  correlation_updated=sapply(1:length(correlation_mtx[,idx1]), function(x) ifelse(abs(correlation_mtx[,idx1][x])>=abs(correlation_mtx[,idx2][x]), correlation_mtx[,idx1][x], correlation_mtx[,idx2][x]))
  col_chosen=sapply(1:length(correlation_updated), function(x) ifelse(correlation_updated[x]==correlation_mtx[,idx1][x], 1, 2))
  h_updated=sapply(1:length(h_mtx[,idx1]), function(x) ifelse(col_chosen[x]==1, h_mtx[,idx1], h_mtx[,idx2]))
  
  correlation_mtx[,idx1]=correlation_updated
  h_mtx[,idx1]=h_updated
  idx_1_all=c(idx_1_all,idx1)
}
correlation_mtx=correlation_mtx[,idx_1_all]
h_mtx=h_mtx[,idx_1_all]

correlation_mtx$Gene=rownames(correlation_mtx)
long_exp.matrix.cor=reshape2::melt(correlation_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="Correlation")
h_mtx$Gene=rownames(h_mtx)
long_exp.matrix.h=reshape2::melt(h_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="h")

# remove the globally lowly expressed genes and age-stable genes
lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_perCelltype.rds"); names(lowlyExprGenes)
celltype.wide_lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_global.rds")
StableGenes=readRDS("~/Project_PBMCage/Tempt_RDS/StableGenes_perCelltype.rds"); names(StableGenes)
global_syn_genes=readRDS("~/Project_PBMCage/Tempt_RDS/Global_synchronized_genes.rds")

gene.to.remove=c(lowlyExprGenes$"CD56dim.NK",lowlyExprGenes$"CD56hi.NK",lowlyExprGenes$"Proliferating.NK.T",
                 StableGenes$"CD56dim.NK",StableGenes$"CD56hi.NK",StableGenes$"Proliferating.NK.T",
                 celltype.wide_lowlyExprGenes,
                 global_syn_genes)
gene.to.remove=gene.to.remove[!duplicated(gene.to.remove)]

library(dplyr)
long_exp.matrix.cor=long_exp.matrix.cor %>%
  subset((Gene %in% gene.to.remove)==FALSE)
long_exp.matrix.h=long_exp.matrix.h %>%
  subset((Gene %in% gene.to.remove)==FALSE)

long_exp.matrix=long_exp.matrix.cor
long_exp.matrix$h=long_exp.matrix.h$h
long_exp.matrix$abs_corr=abs(long_exp.matrix$Correlation)
long_exp.matrix=na.omit(long_exp.matrix)

# Plot distribution of correlation
library(ggplot2)
p_correlation=
  ggplot(long_exp.matrix, aes(x=Correlation, color=Comparison_pair)) +
  geom_density(lwd=1, fill="transparent", alpha=0.25) +
  theme_classic() +
  theme(legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm")) +
  guides(color=guide_legend(ncol=1))
p_correlation
# it seems that the correlation between NK.CD56dim C1orf56+ and NK.CD56dim FCER1G+ is > that between any of the other comparison pairs

### Plot scatterplot of correlation vs. h
ggplot(long_exp.matrix, aes(x=Correlation, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag - NK cells")
ggplot(long_exp.matrix, aes(x=abs_corr, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag - NK cells",
       x="Abs(Correlation)")
# there's a portion of genes with high abs(lag) and high positive corr, meaning that there are leading/lagging genes

### NK-specific synchronized genes (defined as: a large portion of large abs positive correlation + small abs(h))
# determine the threshold for "small correlation" and "small h" based on the plot
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h<=h_thr & h>=-h_thr)

df_ %>% group_by(Comparison_pair) %>% slice_max(abs_corr, n=5)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)
table(as.vector(table(df_comparison.pair$Gene)))

# determine the NK lineage common syn. genes
common_thr=length(table(df_$Comparison_pair)[unlist(lapply(table(df_$Comparison_pair), function(x) x!=0))])/2
common_syn_genes=(df_comparison.pair %>% subset(No.Comparison.Pair>common_thr))[["Gene"]]
# there seems some common syn.genes ("CFAP97"   "DONSON"   "NR4A2"    "SLC25A39" "SMCHD1"),
# but they are related to a wide range of functions - DNA replication, TF, epigenetic regulation, transporter,...

# enrich the subtype-specific syn.genes
library(clusterProfiler)
library(org.Hs.eg.db)
all_genes_to_enrich=split(df_$Gene,df_$Comparison_pair)[sapply(split(df_$Gene,df_$Comparison_pair), function(x) length(x)!=0)]
subtype_compair=names(all_genes_to_enrich)
geneList=list()
pdf("~/Project_PBMCage/Plots/CCF_synGenes_NK.pdf")
for (i in 1:length(subtype_compair)) {
  rest_genes=unlist(all_genes_to_enrich[-i])[!duplicated(unlist(all_genes_to_enrich[-i]))]; rest_genes=as.vector(rest_genes)
  df_chosen=df_ %>%
    subset(Gene %in% all_genes_to_enrich[[i]]) %>%
    subset((Gene %in% rest_genes)==FALSE)
  geneList[[i]]=df_chosen$Gene
  tryCatch({
    gse=enrichGO(df_chosen$Gene, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p0=
      dotplot(gse, showCategory=8, font.size=10) +
      labs(title=paste0("Synchronized genes between \n", paste0(strsplit(subtype_compair[i], "\\.vs\\.")[[1]], collapse=" and "))) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8),
            title=element_text(size=10))
    plot(p0)
  }, error=function(msg) {print("No enriched terms.")})
  
  geneCluster_=list(df_chosen$Gene, global_syn_genes)
  names(geneCluster_)=c(subtype_compair[i], "Global")
  tryCatch({
    GO_enrich_compare=compareCluster(geneCluster=geneCluster_, fun=enrichGO, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p1=
      dotplot(GO_enrich_compare, showCategory=8, font.size=10) +
      labs(x=NULL) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9))
    plot(p1)
  }, error=function(msg) {print("No enriched terms (with global).")})
}
geneList[[i+1]]=common_syn_genes
geneList[[i+2]]=global_syn_genes
names(geneList)=c(subtype_compair, "Common", "Global")
# geneList_deleted=geneList[names(geneList)!="CD8 Tem.GZMK+ NKG7+.vs.CD8 Tnaive.TNFAIP3+"]
tryCatch({
  GO_enrich_compare_all=compareCluster(geneCluster=geneList, fun=enrichGO, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
  p2=
    dotplot(GO_enrich_compare_all, showCategory=3, font.size=10) +
    labs(x=NULL) +
    theme(axis.text.x=element_text(size=8, hjust=1, angle=45),
          axis.text.y=element_text(size=8),
          legend.title=element_text(size=9))
  plot(p2)
}, error=function(msg) {print("No enriched terms with all.")})
dev.off()
# it seems that some comparison pairs have syn.genes that are related to mRNA processing, or RNA splicing, or rRNA processing, 
# apart from the common one (regulation in DNA metabolic process/replication) among all the comparison pairs;
# so we conclude that stimultanously expressed genes in the NK lineage are usually related to DNA/RNA processing.

# plot the expression of the syn.genes in each subtypes of NK subtypes
# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
chose_from=names(table(object$CD8T_NK_subsets))
Celltype_groups=list(DCs=c(chose_from[grepl("DC$|DC[^CA]",chose_from)]),
                     B_cells=c(chose_from[grepl("Bnaive$|Bmem|PB$|PC$",chose_from)]),
                     Monocytes=c(chose_from[grepl("Mono$",chose_from)]),
                     CD4_T=c(chose_from[grepl("^CD4|ISAGhi T|^Treg",chose_from)]),
                     CD8_T=c(chose_from[grepl("^CD8",chose_from)]),
                     ILCs=c(chose_from[grepl("^NK\\.|ILC",chose_from)]),
                     NKTs=c(chose_from[grepl("^(NKT)",chose_from)]),
                     Other_lymph=c(chose_from[grepl("CLP|dnT|gdT|MAIT",chose_from)]),
                     Other_myelo=c(chose_from[grepl("EMP|^Eryth|Mast|Platelet",chose_from)]))
cateogry_name="ILCs"
celltype_group=Celltype_groups$"ILCs"
obj_sub=subset(object, CD8T_NK_subsets %in% celltype_group)
pb_obj=edgeR::Seurat2PB(obj_sub, sample="age", cluster="CD8T_NK_subsets")
pseudobulkdata=pb_obj$counts

subtype_to_keep=as.vector(sapply(subtype_compair, function(x) strsplit(x, "\\.vs\\.")[[1]]))
subtype_to_keep=subtype_to_keep[!duplicated(subtype_to_keep)]
subtype_to_keep_for_grepl=paste0("(",gsub("\\+","\\\\+",gsub("\\.","\\\\.",subtype_to_keep)),")")
subtype_to_keep_for_grepl=paste0(subtype_to_keep_for_grepl, collapse="|")
col.keep=grepl(subtype_to_keep_for_grepl, colnames(pseudobulkdata))
pseudobulkdata=pseudobulkdata[,col.keep]

gene_to_plot=df_$Gene[!duplicated(df_$Gene)]
pseudobulkdata_filtered=pseudobulkdata[gene_to_plot,]

df_syn.gene.expr=as.data.frame(t(pseudobulkdata_filtered))
df_syn.gene.expr_age=gsub("X|_cluster.*","",rownames(df_syn.gene.expr))
df_syn.gene.expr_celltype=gsub(".*cluster","",rownames(df_syn.gene.expr))
df_syn.gene.expr=cbind(Age=as.factor(df_syn.gene.expr_age), Celltype=as.factor(df_syn.gene.expr_celltype), df_syn.gene.expr)

long_df_syn.gene=reshape2::melt(df_syn.gene.expr, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
dim(long_df_syn.gene)
# too many genes, so take the highly expressed ones for plotting
# df_syn.gene.expr_thr=long_df_syn.gene %>% group_by(Celltype, Gene) %>%
#   summarise(across(Expression, mean)) %>%
#   group_by(Celltype) %>%
#   slice_max(order_by=Expression, n=500) %>%
#   as.data.frame()
# genes_chosen=df_syn.gene.expr_thr$Gene[!duplicated(df_syn.gene.expr_thr$Gene)]
# pseudobulkdata_filtered=pseudobulkdata[genes_chosen,]
# df_syn.gene.expr=as.data.frame(t(pseudobulkdata_filtered))
# df_syn.gene.expr_age=gsub("X|_cluster.*","",rownames(df_syn.gene.expr))
# df_syn.gene.expr_celltype=gsub(".*cluster","",rownames(df_syn.gene.expr))
# df_syn.gene.expr=cbind(Age=as.factor(df_syn.gene.expr_age), Celltype=as.factor(df_syn.gene.expr_celltype), df_syn.gene.expr)
# long_df_syn.gene=reshape2::melt(df_syn.gene.expr, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
# dim(long_df_syn.gene)

library(ggplot2)
p_e=
  ggplot(long_df_syn.gene, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
  geom_point(size=0.5, shape=1, alpha=0.1) +
  scale_color_manual(values=ggpubr::get_palette("locuszoom", 10)) +
  scale_y_log10(labels=scales::label_log(digits=1)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19, alpha=1))) +
  labs(title="Synchronized genes in all the subtypes of NK cells",
       x="Age") +
  geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)
plot(p_e)

# enrich the syn.genes from all the NK subtype comparison pairs altogether
gse=enrichGO(gene_to_plot, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
p0=
  dotplot(gse, showCategory=8, font.size=10) +
  labs(title="Synchronized genes in all the subtypes of NK cells") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10))
plot(p0)
# consistent with the general enriched terms of syn.genes from many of the comparison groups (mRNA processing, or RNA splicing, or rRNA processing), 
# but not the same as the common one (regulation in DNA metabolic process/replication) among all the comparison pairs;

# Analyze gene-gene CCF of the syn. genes in each NK subtypes comparison pair
length(gene_to_plot)>300 # too many, so set a expr cutoff:
max.expr_thr=4000
table(apply(pseudobulkdata_filtered, 1, max)>max.expr_thr) # check
gene_to_plot_chosen=rownames(pseudobulkdata_filtered)[apply(pseudobulkdata_filtered, 1, max)>max.expr_thr]

pseudobulk_global=read.delim("~/Project_PBMCage/pseudobulk_global.txt", stringsAsFactors=F)
pseudobulkdata_all_filtered=cbind(pseudobulkdata_filtered[gene_to_plot_chosen,], pseudobulk_global[gene_to_plot_chosen,])

Celltype=gsub(".*_cluster","",colnames(pseudobulkdata_all_filtered))
all_celltypes=names(table(Celltype))

gene_compair_pair_list=c()
for (i in 1:nrow(pseudobulkdata_all_filtered)) {
  rest=c(1:nrow(pseudobulkdata_all_filtered))[-i]
  gene_compair_pair=lapply(1:length(rest), function(x) c(gene_to_plot_chosen[i], gene_to_plot_chosen[rest[x]]))
  gene_compair_pair_list=c(gene_compair_pair_list, gene_compair_pair)
}

synGenes_correlation_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))
synGenes_h_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))

library(TSA)
processbar=txtProgressBar(min=0, max=length(gene_compair_pair_list), style=3, width=60, char="=")
for (i in 1:length(gene_compair_pair_list)) {
  gene1=gene_compair_pair_list[[i]][1]
  gene2=gene_compair_pair_list[[i]][2]
  for (j in 1:length(all_celltypes)) {
    celltype_name=all_celltypes[j]
    celltype.df=pseudobulkdata_all_filtered[, grepl(gsub("\\+","\\\\+",gsub("\\.","\\\\.",celltype_name)), colnames(pseudobulkdata_all_filtered))]
    celltype.df_filtered=celltype.df %>% subset(rownames(.) %in% c(gene1,gene2))
    celltype.df_filtered=as.data.frame(t(celltype.df_filtered))
    celltype.df_filtered=cbind(Age=gsub("X|_.*","",rownames(celltype.df_filtered)), celltype.df_filtered)
    series1=celltype.df_filtered[,2]
    series2=celltype.df_filtered[,3]
    results=prewhiten(series1, series2)
    Abs.correlation=max(abs(results$ccf$acf))
    h=results$ccf$lag[which(abs(results$ccf$acf)==Abs.correlation)]
    correlation=results$ccf$acf[which(results$ccf$lag==h)]
    
    synGenes_correlation_mtx[i, j]=correlation
    synGenes_h_mtx[i, j]=h
  }
  setTxtProgressBar(processbar, i)
}
close(processbar)

synGenes_correlation_mtx=as.data.frame(synGenes_correlation_mtx)
rownames(synGenes_correlation_mtx)=unlist(lapply(gene_compair_pair_list, function(x) paste0(x[1],".vs.",x[2])))
colnames(synGenes_correlation_mtx)=all_celltypes
synGenes_h_mtx=as.data.frame(synGenes_h_mtx)
rownames(synGenes_h_mtx)=rownames(synGenes_correlation_mtx)
colnames(synGenes_h_mtx)=colnames(synGenes_correlation_mtx)

saveRDS(synGenes_correlation_mtx, "~/Project_PBMCage/Results/correlationMatrix_synGenes_NKsubtype.rds")
saveRDS(synGenes_h_mtx, "~/Project_PBMCage/Results/hMatrix_synGenes_NKsubtype.rds")

# Analyze the gene-gene lags among the syn. genes in all the CD8 T subtypes
synGenes_correlation_mtx$Comparison_pair=rownames(synGenes_correlation_mtx)
long_df_synGene_matrix=reshape2::melt(synGenes_correlation_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="Correlation")
synGenes_h_mtx$Comparison_pair=rownames(synGenes_h_mtx)
long_df_synGene_matrix_h=reshape2::melt(synGenes_h_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="h")
long_df_synGene_matrix=cbind(long_df_synGene_matrix, h=long_df_synGene_matrix_h$h)

ggplot(long_df_synGene_matrix, aes(x=Correlation, color=Celltype)) +
  geom_density(lwd=1, alpha=0) +
  scale_color_manual(values=scater:::.get_palette("tableau20"),
                     labels=c(all_celltypes[1:(length(all_celltypes)-1)], "All")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19))) +
  labs(title="Synchronicity among sync. genes in NK cells")
# determine which of them are leading/lagging
range(long_df_synGene_matrix$h)
range(long_df_synGene_matrix$Correlation)
long_df_synGene_matrix=long_df_synGene_matrix %>%
  mutate(Leading=gsub("\\.vs.*","",Comparison_pair)) %>%
  mutate(Lagging=gsub(".*vs\\.","",Comparison_pair)) %>%
  subset(Celltype!="SeuratProject")
leading_strength=long_df_synGene_matrix %>% count(h, Leading)
syn.Genes=(leading_strength %>% subset(h==0 & n>length(gene_compair_pair_list)/length(gene_to_plot_chosen)*length(all_celltypes)/2))$Leading
leading.Genes=(leading_strength %>% subset(h<0 & (Leading %in% syn.Genes)==FALSE))$Leading
lagging.Genes=(leading_strength %>% subset(h>0 & (Leading %in% syn.Genes)==FALSE))$Leading
# it seems that none of them are leading or lagging, but all expressed simultaneously

### Determine the NK subtype-dependent leading/lagging genes
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h>h_thr | h<(-h_thr))

gene_cluster=split(df_$Gene, paste0(df_$Comparison_pair, "__", df_$h))
gene_cluster_sort=gene_cluster[order(abs(as.numeric(gsub(".*__","",names(gene_cluster)))), decreasing=T)]

pdf("~/Project_PBMCage/Plots/Correlation_H_Leading.Lagging.Genes_Expr_NKcell.pdf")
for (i in 1:length(gene_cluster_sort)) {
  cluster_=gene_cluster_sort[i]
  h_value=gsub(".*__","",names(cluster_))
  type1_and_2=strsplit(names(cluster_),split=".vs.")[[1]]
  type1=type1_and_2[1]
  type2=gsub("__.*","",type1_and_2[2])
  cluster=cluster_[[1]]
  pseudobulkdata_filtered=pseudobulkdata[cluster,]
  
  df_byCluster.Gene=cbind(pseudobulkdata[cluster,], pseudobulk_global[cluster,])
  df_byCluster.Gene_t_age=gsub("X|_cluster.*","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t_celltype=gsub(".*cluster","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t=data.frame(Age=as.factor(df_byCluster.Gene_t_age), Celltype=as.factor(df_byCluster.Gene_t_celltype), as.data.frame(t(df_byCluster.Gene)))
  df_byCluster.Gene_t_chosen=df_byCluster.Gene_t %>%
    subset(Celltype %in% c(type1, type2))
  long_df_leading.lagging.gene.expr=reshape2::melt(df_byCluster.Gene_t_chosen, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
  
  cluster_name=ifelse(length(cluster)>5, paste0(paste0(cluster[1:5], collapse=", "),", ..."), paste0(cluster, collapse=", "))
  p_x=
    ggplot(long_df_leading.lagging.gene.expr, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
    geom_point(shape=1, alpha=0.3) +
    scale_color_manual(values=scater:::.get_palette("tableau20")[c(1,11)]) +
    scale_y_continuous(trans="pseudo_log") +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=9),
          legend.text=element_text(size=8),
          title=element_text(size=10)) + 
    guides(color=guide_legend(override.aes=list(size=2, alpha=1))) +
    labs(title=paste0(type1, " vs. ", type2, "\n(", cluster_name, "; h=", h_value, ")"),
         x="Age") +
    geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)
  
  plot(p_x)
  
  tryCatch({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    gse=enrichGO(cluster, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p0=
      dotplot(gse, showCategory=8, font.size=10) +
      labs(title=ifelse(h_value<0, paste0("Genes expressed earlier in ",type1, "\nthan in ", type2), paste0("Genes expressed earlier in ",type2, "\nthan in ", type1))) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8),
            title=element_text(size=10))
    plot(p0)
  }, error=function(msg) {print("No term enriched.")})
}
dev.off()
# the genes with certain functions are expressed earlier in a subtype than in another, meaning that the functions are important and working properly in that celltype at young ages,
# and diminish more quickly during aging

### Analyze negatively regulated genes
corr_thr=0.5

df_=long_exp.matrix %>%
  subset(Correlation<=-corr_thr)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)
table(df_comparison.pair$No.Comparison.Pair)
split(df_$Comparison_pair, df_$Gene)
# it turns out that: AGL in NK.CD56bright FCGR3A+.vs.NK.CD56dim C1orf56+
# AK5 in NK.CD56dim FCER1G+.vs.NK.CD56dim FOS+
# CAT in NK.CD56dim C1orf56+.vs.NK.CD56dim FOS+
# CNOT6, COPS2 in NK.CD56bright FCGR3A+.vs.NK.CD56dim FOS+
# GART in NK.CD56bright COTL1+.vs.NK.CD56dim FCER1G+
# NOP2, SESN2 in NK.CD56bright FCGR3A+.vs.NK.CD56dim FCER1G+
# SLC25A17 in NK.CD56bright FCGR3A+.vs.NK.CD56dim FOS+
# WAPL in NK.CD56bright FCGR3A+.vs.NK.CD56dim C1orf56+
# ZC3H8 in NK.CD56bright FCGR3A+.vs.NK.CD56dim FOS+ are neg.reg.

### Analyze less correlated genes
corr_thr=0.25

df_=long_exp.matrix %>%
  subset(abs(Correlation)<corr_thr)
# Determine the less correlated genes in various celltypes
library(clusterProfiler)
library(org.Hs.eg.db)
geneList=split(df_$Gene, df_$Comparison_pair)
geneList=geneList[sapply(geneList, function(x) length(x)!=0)]
GO_enrich_compare=compareCluster(geneCluster=geneList, fun=enrichGO, OrgDb="org.Hs.eg.db", pvalueCutoff=0.05)
# no term enriched, which means that the less correlated genes in each celltype comparison pair are just random genes

#####################################





###
### CCF analysis on CD4 T cell lineage
#####################################
###

### Analyze the Corr and Lag
# Arrange correlation and h matrix for plotting
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
correlation_mtx=readRDS("~/Project_PBMCage/Results/correlationMatrix_CD4_T.rds")
h_mtx=readRDS("~/Project_PBMCage/Results/hMatrix_CD4_T.rds")

chose_from=names(table(object$CD8T_NK_subsets))
Celltype_groups=list(DCs=c(chose_from[grepl("DC$|DC[^CA]",chose_from)]),
                     B_cells=c(chose_from[grepl("Bnaive$|Bmem|PB$|PC$",chose_from)]),
                     Monocytes=c(chose_from[grepl("Mono$",chose_from)]),
                     CD4_T=c(chose_from[grepl("^CD4|ISAGhi T|^Treg",chose_from)]),
                     CD8_T=c(chose_from[grepl("^CD8",chose_from)]),
                     ILCs=c(chose_from[grepl("^NK\\.|ILC",chose_from)]),
                     NKTs=c(chose_from[grepl("^(NKT)",chose_from)]),
                     Other_lymph=c(chose_from[grepl("CLP|dnT|gdT|MAIT",chose_from)]),
                     Other_myelo=c(chose_from[grepl("EMP|^Eryth|Mast|Platelet",chose_from)]))
all_celltypes=Celltype_groups$CD4_T
compair_pair_list=c()
for (i in 1:length(all_celltypes)) {
  rest=c(1:length(all_celltypes))[-i]
  compair_pair=lapply(1:length(rest), function(x) c(all_celltypes[i], all_celltypes[rest[x]]))
  compair_pair_list=c(compair_pair_list, compair_pair)
}

colnames(correlation_mtx)=sapply(compair_pair_list, function(x) paste0(x[1], ".vs.", x[2]))
colnames(h_mtx)=sapply(compair_pair_list, function(x) paste0(x[1], ".vs.", x[2]))
compare_idx=list()
for (idx in 1:length(colnames(correlation_mtx))) {
  pair_1=colnames(correlation_mtx)[idx]
  pair_2=paste0(gsub(".*vs\\.","",pair_1), ".vs.", gsub("\\.vs.*","",pair_1))
  compare_idx1=c(1:length(colnames(correlation_mtx)))[idx]
  compare_idx2=which(colnames(correlation_mtx)==pair_2)
  if ((compare_idx2 %in% c(1:idx))==FALSE) {
    compare_idx=c(compare_idx, list(c(compare_idx1, compare_idx2)))
  }
}
idx_1_all=c()
for (i in 1:length(compare_idx)) {
  idx1=compare_idx[[i]][1]
  idx2=compare_idx[[i]][2]
  correlation_updated=sapply(1:length(correlation_mtx[,idx1]), function(x) ifelse(abs(correlation_mtx[,idx1][x])>=abs(correlation_mtx[,idx2][x]), correlation_mtx[,idx1][x], correlation_mtx[,idx2][x]))
  col_chosen=sapply(1:length(correlation_updated), function(x) ifelse(correlation_updated[x]==correlation_mtx[,idx1][x], 1, 2))
  h_updated=sapply(1:length(h_mtx[,idx1]), function(x) ifelse(col_chosen[x]==1, h_mtx[,idx1], h_mtx[,idx2]))
  
  correlation_mtx[,idx1]=correlation_updated
  h_mtx[,idx1]=h_updated
  idx_1_all=c(idx_1_all,idx1)
}
correlation_mtx=correlation_mtx[,idx_1_all]
h_mtx=h_mtx[,idx_1_all]

correlation_mtx$Gene=rownames(correlation_mtx)
long_exp.matrix.cor=reshape2::melt(correlation_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="Correlation")
h_mtx$Gene=rownames(h_mtx)
long_exp.matrix.h=reshape2::melt(h_mtx, id.vars=c("Gene"), variable.name="Comparison_pair", value.name="h")

# remove the globally lowly expressed genes and age-stable genes
lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_perCelltype.rds"); names(lowlyExprGenes)
celltype.wide_lowlyExprGenes=readRDS("~/Project_PBMCage/Tempt_RDS/Expr_lowlyExprGenes_global.rds")
StableGenes=readRDS("~/Project_PBMCage/Tempt_RDS/StableGenes_perCelltype.rds"); names(StableGenes)
global_syn_genes=readRDS("~/Project_PBMCage/Tempt_RDS/Global_synchronized_genes.rds")

gene.to.remove=c(lowlyExprGenes$"CD4.Tmem",lowlyExprGenes$"CD4.Tnaive",lowlyExprGenes$"Proliferating.NK.T",
                 StableGenes$"CD4.Tmem",StableGenes$"CD4.Tnaive",StableGenes$"Proliferating.NK.T",
                 celltype.wide_lowlyExprGenes,
                 global_syn_genes)
gene.to.remove=gene.to.remove[!duplicated(gene.to.remove)]

library(dplyr)
long_exp.matrix.cor=long_exp.matrix.cor %>%
  subset((Gene %in% gene.to.remove)==FALSE)
long_exp.matrix.h=long_exp.matrix.h %>%
  subset((Gene %in% gene.to.remove)==FALSE)

long_exp.matrix=long_exp.matrix.cor
long_exp.matrix$h=long_exp.matrix.h$h
long_exp.matrix$abs_corr=abs(long_exp.matrix$Correlation)
long_exp.matrix=na.omit(long_exp.matrix)

# Plot distribution of correlation
library(ggplot2)
p_correlation=
  ggplot(long_exp.matrix, aes(x=Correlation, color=Comparison_pair)) +
  geom_density(lwd=1, fill="transparent", alpha=0.25) +
  theme_classic() +
  theme(legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm")) +
  guides(color=guide_legend(ncol=1))
p_correlation
# it seems that the correlation between CD4 Tem and CD4 Tnaive is > that between CD4 Tcm and CD4 Tem or between CD4 Tcm and CD4 Tnaive

### Plot scatterplot of correlation vs. h
ggplot(long_exp.matrix, aes(x=Correlation, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag - CD4 T cells")
ggplot(long_exp.matrix, aes(x=abs_corr, y=h)) +
  geom_point(shape=1, color="black", alpha=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) +
  labs(title="Correlation vs. Lag - NK cells",
       x="Abs(Correlation)")
# all the genes are with lag=0 (although some have high corr and some have low), meaning that there aren't leading/lagging genes between any comparison pair

### CD4 T-specific synchronized genes (defined as: a large portion of large abs positive correlation + small abs(h))
# determine the threshold for "small correlation" and "small h" based on the plot
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h<=h_thr & h>=-h_thr)

df_ %>% group_by(Comparison_pair) %>% slice_max(abs_corr, n=5)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)
table(as.vector(table(df_comparison.pair$Gene)))

# determine the CD4 T common syn. genes
common_thr=length(table(df_$Comparison_pair)[unlist(lapply(table(df_$Comparison_pair), function(x) x!=0))])/2
common_syn_genes=(df_comparison.pair %>% subset(No.Comparison.Pair>common_thr))[["Gene"]]
length(common_syn_genes)
library(clusterProfiler)
library(org.Hs.eg.db)
gse=enrichGO(common_syn_genes, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
p_common=
  dotplot(gse, showCategory=8, font.size=10) +
  labs(title="Common synchronized genes in CD4 T cell subtypes") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10))
plot(p_common)
# there seems to be a lot of (8794) common syn.genes in CD4 T subtypes,
# which are related to mRNA processing, RNA splicing, ribosome biogenesis, histone modification,...

# enrich the subtype-specific syn.genes
library(clusterProfiler)
library(org.Hs.eg.db)
all_genes_to_enrich=split(df_$Gene,df_$Comparison_pair)[sapply(split(df_$Gene,df_$Comparison_pair), function(x) length(x)!=0)]
subtype_compair=names(all_genes_to_enrich)
geneList=list()
pdf("~/Project_PBMCage/Plots/CCF_synGenes_CD4T.pdf")
for (i in 1:length(subtype_compair)) {
  rest_genes=unlist(all_genes_to_enrich[-i])[!duplicated(unlist(all_genes_to_enrich[-i]))]; rest_genes=as.vector(rest_genes)
  df_chosen=df_ %>%
    subset(Gene %in% all_genes_to_enrich[[i]]) %>%
    subset((Gene %in% rest_genes)==FALSE)
  geneList[[i]]=df_chosen$Gene
  tryCatch({
    gse=enrichGO(df_chosen$Gene, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p0=
      dotplot(gse, showCategory=8, font.size=10) +
      labs(title=paste0("Synchronized genes between \n", paste0(strsplit(subtype_compair[i], "\\.vs\\.")[[1]], collapse=" and "))) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8),
            title=element_text(size=10))
    plot(p0)
  }, error=function(msg) {print("No enriched terms.")})
  
  geneCluster_=list(df_chosen$Gene, global_syn_genes)
  names(geneCluster_)=c(subtype_compair[i], "Global")
  tryCatch({
    GO_enrich_compare=compareCluster(geneCluster=geneCluster_, fun=enrichGO, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p1=
      dotplot(GO_enrich_compare, showCategory=8, font.size=10) +
      labs(x=NULL) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9))
    plot(p1)
  }, error=function(msg) {print("No enriched terms (with global).")})
}
geneList[[i+1]]=common_syn_genes
geneList[[i+2]]=global_syn_genes
names(geneList)=c(subtype_compair, "Common", "Global")
# geneList_deleted=geneList[names(geneList)!="CD8 Tem.GZMK+ NKG7+.vs.CD8 Tnaive.TNFAIP3+"]
tryCatch({
  GO_enrich_compare_all=compareCluster(geneCluster=geneList, fun=enrichGO, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
  p2=
    dotplot(GO_enrich_compare_all, showCategory=3, font.size=10) +
    labs(x=NULL) +
    theme(axis.text.x=element_text(size=8, hjust=1, angle=45),
          axis.text.y=element_text(size=8),
          legend.title=element_text(size=9))
  plot(p2)
}, error=function(msg) {print("No enriched terms with all.")})
dev.off()
# it seems that CD4 Tcm and CD4 Tnaive have syn.genes that are related to DNA replication,
# while CD4 Tem does not have enriched term regarding its syn.genes with CD4 Tcm or CD4 Tnaive
# interestingly, the common syn.genes among all the three comparison groups are related to RNA processing instead.

# plot the expression of the syn.genes in each subtypes of CD4 T subtypes
# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
chose_from=names(table(object$CD8T_NK_subsets))
Celltype_groups=list(DCs=c(chose_from[grepl("DC$|DC[^CA]",chose_from)]),
                     B_cells=c(chose_from[grepl("Bnaive$|Bmem|PB$|PC$",chose_from)]),
                     Monocytes=c(chose_from[grepl("Mono$",chose_from)]),
                     CD4_T=c(chose_from[grepl("^CD4|ISAGhi T|^Treg",chose_from)]),
                     CD8_T=c(chose_from[grepl("^CD8",chose_from)]),
                     ILCs=c(chose_from[grepl("^NK\\.|ILC",chose_from)]),
                     NKTs=c(chose_from[grepl("^(NKT)",chose_from)]),
                     Other_lymph=c(chose_from[grepl("CLP|dnT|gdT|MAIT",chose_from)]),
                     Other_myelo=c(chose_from[grepl("EMP|^Eryth|Mast|Platelet",chose_from)]))
cateogry_name="CD4_T"
celltype_group=Celltype_groups$"CD4_T"
obj_sub=subset(object, CD8T_NK_subsets %in% celltype_group)
pb_obj=edgeR::Seurat2PB(obj_sub, sample="age", cluster="CD8T_NK_subsets")
pseudobulkdata=pb_obj$counts

subtype_to_keep=as.vector(sapply(subtype_compair, function(x) strsplit(x, "\\.vs\\.")[[1]]))
subtype_to_keep=subtype_to_keep[!duplicated(subtype_to_keep)]
subtype_to_keep_for_grepl=paste0("(",gsub("\\+","\\\\+",gsub("\\.","\\\\.",subtype_to_keep)),")")
subtype_to_keep_for_grepl=paste0(subtype_to_keep_for_grepl, collapse="|")
col.keep=grepl(subtype_to_keep_for_grepl, colnames(pseudobulkdata))
pseudobulkdata=pseudobulkdata[,col.keep]

gene_to_plot=df_$Gene[!duplicated(df_$Gene)]
pseudobulkdata_filtered=pseudobulkdata[gene_to_plot,]

df_syn.gene.expr=as.data.frame(t(pseudobulkdata_filtered))
df_syn.gene.expr_age=gsub("X|_cluster.*","",rownames(df_syn.gene.expr))
df_syn.gene.expr_celltype=gsub(".*cluster","",rownames(df_syn.gene.expr))
df_syn.gene.expr=cbind(Age=as.factor(df_syn.gene.expr_age), Celltype=as.factor(df_syn.gene.expr_celltype), df_syn.gene.expr)

long_df_syn.gene=reshape2::melt(df_syn.gene.expr, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
dim(long_df_syn.gene)
# too many genes, so take the highly expressed ones for plotting
df_syn.gene.expr_thr=long_df_syn.gene %>% group_by(Celltype, Gene) %>%
  summarise(across(Expression, mean)) %>%
  group_by(Celltype) %>%
  slice_max(order_by=Expression, n=1000) %>%
  as.data.frame()
genes_chosen=df_syn.gene.expr_thr$Gene[!duplicated(df_syn.gene.expr_thr$Gene)]
pseudobulkdata_filtered=pseudobulkdata[genes_chosen,]
df_syn.gene.expr=as.data.frame(t(pseudobulkdata_filtered))
df_syn.gene.expr_age=gsub("X|_cluster.*","",rownames(df_syn.gene.expr))
df_syn.gene.expr_celltype=gsub(".*cluster","",rownames(df_syn.gene.expr))
df_syn.gene.expr=cbind(Age=as.factor(df_syn.gene.expr_age), Celltype=as.factor(df_syn.gene.expr_celltype), df_syn.gene.expr)
long_df_syn.gene=reshape2::melt(df_syn.gene.expr, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
dim(long_df_syn.gene)

library(ggplot2)
p_e=
  ggplot(long_df_syn.gene, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
  geom_point(size=0.5, shape=1, alpha=0.1) +
  scale_color_manual(values=ggpubr::get_palette("locuszoom", 10)) +
  scale_y_log10(labels=scales::label_log(digits=1)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19, alpha=1))) +
  labs(title="Synchronized genes in all the subtypes of CD4 T cells",
       x="Age") +
  geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)
plot(p_e)

# enrich the syn.genes from all the CD4 T subtype comparison pairs altogether
gse=enrichGO(gene_to_plot, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
p0=
  dotplot(gse, showCategory=8, font.size=10) +
  labs(title="Synchronized genes in all the subtypes of CD4 T cells") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        title=element_text(size=10))
plot(p0)
# most of the enriched terms here are the same as those for common syn.genes in various CD4 T subtypes,
# because the common syn.genes constitute the large portion of all the syn.genes (8794 vs. 9150)

# Analyze gene-gene CCF of the syn. genes in each CD4 T subtypes comparison pair
length(gene_to_plot)>300 # too many, so set a expr cutoff:
max.expr_thr=4000
table(apply(pseudobulkdata_filtered, 1, max)>max.expr_thr) # check
gene_to_plot_chosen=rownames(pseudobulkdata_filtered)[apply(pseudobulkdata_filtered, 1, max)>max.expr_thr]

pseudobulk_global=read.delim("~/Project_PBMCage/pseudobulk_global.txt", stringsAsFactors=F)
pseudobulkdata_all_filtered=cbind(pseudobulkdata_filtered[gene_to_plot_chosen,], pseudobulk_global[gene_to_plot_chosen,])

Celltype=gsub(".*_cluster","",colnames(pseudobulkdata_all_filtered))
all_celltypes=names(table(Celltype))

gene_compair_pair_list=c()
for (i in 1:nrow(pseudobulkdata_all_filtered)) {
  rest=c(1:nrow(pseudobulkdata_all_filtered))[-i]
  gene_compair_pair=lapply(1:length(rest), function(x) c(gene_to_plot_chosen[i], gene_to_plot_chosen[rest[x]]))
  gene_compair_pair_list=c(gene_compair_pair_list, gene_compair_pair)
}

synGenes_correlation_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))
synGenes_h_mtx=matrix(NA, nrow=length(gene_compair_pair_list), ncol=length(all_celltypes))

library(TSA)
processbar=txtProgressBar(min=0, max=length(gene_compair_pair_list), style=3, width=60, char="=")
for (i in 1:length(gene_compair_pair_list)) {
  gene1=gene_compair_pair_list[[i]][1]
  gene2=gene_compair_pair_list[[i]][2]
  for (j in 1:length(all_celltypes)) {
    celltype_name=all_celltypes[j]
    celltype.df=pseudobulkdata_all_filtered[, grepl(gsub("\\+","\\\\+",gsub("\\.","\\\\.",celltype_name)), colnames(pseudobulkdata_all_filtered))]
    celltype.df_filtered=celltype.df %>% subset(rownames(.) %in% c(gene1,gene2))
    celltype.df_filtered=as.data.frame(t(celltype.df_filtered))
    celltype.df_filtered=cbind(Age=gsub("X|_.*","",rownames(celltype.df_filtered)), celltype.df_filtered)
    series1=celltype.df_filtered[,2]
    series2=celltype.df_filtered[,3]
    results=prewhiten(series1, series2)
    Abs.correlation=max(abs(results$ccf$acf))
    h=results$ccf$lag[which(abs(results$ccf$acf)==Abs.correlation)]
    correlation=results$ccf$acf[which(results$ccf$lag==h)]
    
    synGenes_correlation_mtx[i, j]=correlation
    synGenes_h_mtx[i, j]=h
  }
  setTxtProgressBar(processbar, i)
}
close(processbar)

synGenes_correlation_mtx=as.data.frame(synGenes_correlation_mtx)
rownames(synGenes_correlation_mtx)=unlist(lapply(gene_compair_pair_list, function(x) paste0(x[1],".vs.",x[2])))
colnames(synGenes_correlation_mtx)=all_celltypes
synGenes_h_mtx=as.data.frame(synGenes_h_mtx)
rownames(synGenes_h_mtx)=rownames(synGenes_correlation_mtx)
colnames(synGenes_h_mtx)=colnames(synGenes_correlation_mtx)

saveRDS(synGenes_correlation_mtx, "~/Project_PBMCage/Results/correlationMatrix_synGenes_CD4Tsubtype.rds")
saveRDS(synGenes_h_mtx, "~/Project_PBMCage/Results/hMatrix_synGenes_CD4Tsubtype.rds")

# Analyze the gene-gene lags among the syn. genes in all the CD4 T subtypes
synGenes_correlation_mtx$Comparison_pair=rownames(synGenes_correlation_mtx)
long_df_synGene_matrix=reshape2::melt(synGenes_correlation_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="Correlation")
synGenes_h_mtx$Comparison_pair=rownames(synGenes_h_mtx)
long_df_synGene_matrix_h=reshape2::melt(synGenes_h_mtx, id.vars=c("Comparison_pair"), variable.name="Celltype", value.name="h")
long_df_synGene_matrix=cbind(long_df_synGene_matrix, h=long_df_synGene_matrix_h$h)

ggplot(long_df_synGene_matrix, aes(x=Correlation, color=Celltype)) +
  geom_density(lwd=1, alpha=0) +
  scale_color_manual(values=scater:::.get_palette("tableau20"),
                     labels=c(all_celltypes[1:(length(all_celltypes)-1)], "All")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"),
        title=element_text(size=10)) + 
  guides(color=guide_legend(override.aes=list(size=1, shape=19))) +
  labs(title="Synchronicity among sync. genes in CD4 T cells")
# determine which of them are leading/lagging
range(long_df_synGene_matrix$h)
range(long_df_synGene_matrix$Correlation)
long_df_synGene_matrix=long_df_synGene_matrix %>%
  mutate(Leading=gsub("\\.vs.*","",Comparison_pair)) %>%
  mutate(Lagging=gsub(".*vs\\.","",Comparison_pair)) %>%
  subset(Celltype!="SeuratProject")
leading_strength=long_df_synGene_matrix %>% count(h, Leading)
syn.Genes=(leading_strength %>% subset(h==0 & n>length(gene_compair_pair_list)/length(gene_to_plot_chosen)*length(all_celltypes)/2))$Leading
leading.Genes=(leading_strength %>% subset(h<0 & (Leading %in% syn.Genes)==FALSE))$Leading
lagging.Genes=(leading_strength %>% subset(h>0 & (Leading %in% syn.Genes)==FALSE))$Leading
# it seems that none of them are leading or lagging, but all expressed simultaneously

### Determine the NK subtype-dependent leading/lagging genes
corr_thr=0.6
h_thr=3

df_=long_exp.matrix %>%
  subset(Correlation>=corr_thr) %>%
  subset(h>h_thr | h<(-h_thr))

gene_cluster=split(df_$Gene, paste0(df_$Comparison_pair, "__", df_$h))
gene_cluster_sort=gene_cluster[order(abs(as.numeric(gsub(".*__","",names(gene_cluster)))), decreasing=T)]

pdf("~/Project_PBMCage/Plots/Correlation_H_Leading.Lagging.Genes_Expr_NKcell.pdf")
for (i in 1:length(gene_cluster_sort)) {
  cluster_=gene_cluster_sort[i]
  h_value=gsub(".*__","",names(cluster_))
  type1_and_2=strsplit(names(cluster_),split=".vs.")[[1]]
  type1=type1_and_2[1]
  type2=gsub("__.*","",type1_and_2[2])
  cluster=cluster_[[1]]
  pseudobulkdata_filtered=pseudobulkdata[cluster,]
  
  df_byCluster.Gene=cbind(pseudobulkdata[cluster,], pseudobulk_global[cluster,])
  df_byCluster.Gene_t_age=gsub("X|_cluster.*","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t_celltype=gsub(".*cluster","",colnames(df_byCluster.Gene))
  df_byCluster.Gene_t=data.frame(Age=as.factor(df_byCluster.Gene_t_age), Celltype=as.factor(df_byCluster.Gene_t_celltype), as.data.frame(t(df_byCluster.Gene)))
  df_byCluster.Gene_t_chosen=df_byCluster.Gene_t %>%
    subset(Celltype %in% c(type1, type2))
  long_df_leading.lagging.gene.expr=reshape2::melt(df_byCluster.Gene_t_chosen, id.vars=c("Celltype","Age"), variable.name="Gene", value.name="Expression")
  
  cluster_name=ifelse(length(cluster)>5, paste0(paste0(cluster[1:5], collapse=", "),", ..."), paste0(cluster, collapse=", "))
  p_x=
    ggplot(long_df_leading.lagging.gene.expr, aes(x=Age, y=Expression, color=Celltype, group=Celltype)) +
    geom_point(shape=1, alpha=0.3) +
    scale_color_manual(values=scater:::.get_palette("tableau20")[c(1,11)]) +
    scale_y_continuous(trans="pseudo_log") +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=9),
          legend.text=element_text(size=8),
          title=element_text(size=10)) + 
    guides(color=guide_legend(override.aes=list(size=2, alpha=1))) +
    labs(title=paste0(type1, " vs. ", type2, "\n(", cluster_name, "; h=", h_value, ")"),
         x="Age") +
    geom_smooth(method='auto', show.legend=FALSE, linewidth=1, alpha=0)
  
  plot(p_x)
  
  tryCatch({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    gse=enrichGO(cluster, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
    p0=
      dotplot(gse, showCategory=8, font.size=10) +
      labs(title=ifelse(h_value<0, paste0("Genes expressed earlier in ",type1, "\nthan in ", type2), paste0("Genes expressed earlier in ",type2, "\nthan in ", type1))) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8),
            title=element_text(size=10))
    plot(p0)
  }, error=function(msg) {print("No term enriched.")})
}
dev.off()
# the genes with certain functions are expressed earlier in a subtype than in another, meaning that the functions are important and working properly in that celltype at young ages,
# and diminish more quickly during aging

### Analyze negatively regulated genes
corr_thr=0.5

df_=long_exp.matrix %>%
  subset(Correlation<=-corr_thr)

No.Comparison.Pair=unlist(lapply(split(df_$Comparison_pair, df_$Gene), length))
df_comparison.pair=data.frame(Gene=names(No.Comparison.Pair), No.Comparison.Pair=No.Comparison.Pair)
table(df_comparison.pair$No.Comparison.Pair)
split(df_$Comparison_pair, df_$Gene)
# it turns out that: AGL in NK.CD56bright FCGR3A+.vs.NK.CD56dim C1orf56+
# AK5 in NK.CD56dim FCER1G+.vs.NK.CD56dim FOS+
# CAT in NK.CD56dim C1orf56+.vs.NK.CD56dim FOS+
# CNOT6, COPS2 in NK.CD56bright FCGR3A+.vs.NK.CD56dim FOS+
# GART in NK.CD56bright COTL1+.vs.NK.CD56dim FCER1G+
# NOP2, SESN2 in NK.CD56bright FCGR3A+.vs.NK.CD56dim FCER1G+
# SLC25A17 in NK.CD56bright FCGR3A+.vs.NK.CD56dim FOS+
# WAPL in NK.CD56bright FCGR3A+.vs.NK.CD56dim C1orf56+
# ZC3H8 in NK.CD56bright FCGR3A+.vs.NK.CD56dim FOS+ are neg.reg.

### Analyze less correlated genes
corr_thr=0.25

df_=long_exp.matrix %>%
  subset(abs(Correlation)<corr_thr)
# Determine the less correlated genes in various celltypes
library(clusterProfiler)
library(org.Hs.eg.db)
geneList=split(df_$Gene, df_$Comparison_pair)
geneList=geneList[sapply(geneList, function(x) length(x)!=0)]
GO_enrich_compare=compareCluster(geneCluster=geneList, fun=enrichGO, OrgDb="org.Hs.eg.db", pvalueCutoff=0.05)
# no term enriched, which means that the less correlated genes in each celltype comparison pair are just random genes

#####################################





###
### Graph analysis on time-serial gene expression in celltypes at the detailed level
#####################################
###

library(Seurat)
### Load
THEObj=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
assigned_types=unique(THEObj$CD8T_NK_subsets2)[!grepl("lin_infidel|mix", unique(THEObj$CD8T_NK_subsets2))]
THEObj=subset(THEObj, CD8T_NK_subsets2 %in% assigned_types)

### Extract the normalized matrix from the "CD8T cells"
pseudo=AggregateExpression(THEObj, assays="RNA", return.seurat=T, group.by=c("age", "CD8T_NK_subsets2_rough"))
pseudo=NormalizeData(pseudo)
pseudo$age=gsub("_.*","",colnames(pseudo))
pseudo$CD8T_NK_subsets2_rough=gsub(".*_","",colnames(pseudo))
pseudo.CD8T=subset(pseudo, CD8T_NK_subsets2_rough=="CD8T cells")
dataMatrix=LayerData(pseudo.CD8T, assay="RNA", layer="data")
table(pseudo.CD8T$age==sort(pseudo.CD8T$age)) # check the time series

### Calculate cosine similarity between two vectors
cosine_difference=function(vector1, vector2) {
  dot_product=sum(vector1*vector2)
  magnitude1=sqrt(sum(vector1^2))
  magnitude2=sqrt(sum(vector2^2))
  cos_difference=dot_product/(magnitude1*magnitude2)
  return(cos_difference)
}

dense_matrix=as.matrix(dataMatrix)
diff_matrix=matrix(0, nrow=nrow(dense_matrix), ncol=nrow(dense_matrix))
for (i in 1:length(rownames(dense_matrix))) {
  for (j in 1:length(rownames(dense_matrix))) {
    diff_matrix[i,j]=cosine_difference(dense_matrix[i,], dense_matrix[j,])
  }
  message(paste0("====== Gene", i, " Done! ======"))
}
colnames(diff_matrix)=rownames(diff_matrix)=rownames(dense_matrix)
write.table(diff_matrix, "~/Project_PBMCage/From_and_to_Jin/TimeSerialGeneExp_CD8Tcells.txt", sep=",")

### Calculate pearson correlation between two vectors
library(qlcMatrix)
cor_result=corSparse(t(dataMatrix))
colnames(cor_result)=rownames(cor_result)=rownames(dataMatrix)
write.table(cor_result, "~/Project_PBMCage/From_and_to_Jin/TimeSerialGeneExp_CD8Tcells_spearson.txt", sep=",")

