
### Prepare the 5000 gene weights matrix
#####################################
###

### Load the weights for each donor
csv_files=list.files("~/Project_PBMCage/For_or_From_Xu/layer6_newCT/")
CSV_DATA_ALL=list()
for (i in 1:length(csv_files)) {
  file_=file.path("~/Project_PBMCage/For_or_From_Xu/layer6_newCT", csv_files[i])
  csv_data=read.csv(file_, header=F, sep=" ")
  CSV_DATA_ALL[[i]]=csv_data
}
names(CSV_DATA_ALL)=csv_files
saveRDS(CSV_DATA_ALL, "~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000_newCT.rds")

### Add up the weights from each donor
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000_newCT.rds")
# check if the genes are in the same order in each file
gene_list=Gene_weights_5000[[1]]$V2
for (i in 1:length(Gene_weights_5000)) {
  single_file=Gene_weights_5000[[i]]
  if (FALSE %in% unique(single_file$V2==gene_list)) {
    print("Unsorted genes detected.")
  } else {print("Pass.")}
}
# sum
for (i in 1:length(Gene_weights_5000)) {
  colnames(Gene_weights_5000[[i]])=c("donor_id","gene",Gene_weights_5000[[i]]$V1[1])
  rownames(Gene_weights_5000[[i]])=Gene_weights_5000[[i]]$gene
  Gene_weights_5000[[i]]=Gene_weights_5000[[i]] %>% dplyr::select(-c(donor_id, gene))
}
Gene_weights_5000_all=do.call(cbind, Gene_weights_5000)
Gene_weights_5000_all$mean=rowMeans(Gene_weights_5000_all)

### Convert ENSEMBL to SYMBOL
ConvertTable=read.delim("~/Project_PBMCage/For_or_From_Xu/GeneFeatures_ENStoSYMBOL.txt", sep="\t")
Gene_weights_5000_all$symbol=ConvertTable$feature_name[match(rownames(Gene_weights_5000_all), rownames(ConvertTable))]

### Save
saveRDS(Gene_weights_5000_all, "~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000_newCT.rds")

#####################################



### Draw Dotplots of the top50 genes
#####################################
###

library(Seurat)
library(dplyr)
library(scRNAtoolVis)

Seurat_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000_newCT.rds")
Top_genes=Gene_weights_5000 %>% 
  # subset(!grepl("^MT-|^RP[SL]", symbol)) %>%
  slice_max(mean, n=50) %>% # take the top 50 genes
  dplyr::select(symbol) %>% 
  unlist() %>% 
  unname()

plot_rough=
  jjDotPlot(object=Seurat_obj,
            gene=Top_genes,
            id="Annot.rough",
            xtree=T,
            rescale=F,
            rescale.min=0,
            rescale.max=1,
            # point.geom=F,
            tile.geom=T) +
  scale_size(range=c(4, 8))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/GeneWeight_5000gene.top50_expr_1.pdf", height=5, width=20)
plot(plot_rough)
dev.off()

plot_inter=
  jjDotPlot(object=Seurat_obj,
            gene=Top_genes,
            id="Annot.inter",
            xtree=T,
            rescale=F,
            rescale.min=0,
            rescale.max=1,
            # point.geom=F,
            tile.geom=T) +
  scale_size(range=c(4, 8))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/GeneWeight_5000gene.top50_expr_2.pdf", height=12, width=20)
plot(plot_inter)
dev.off()

plot_detailed=
  jjDotPlot(object=Seurat_obj,
            gene=Top_genes,
            id="Annot.detailed",
            xtree=T,
            rescale=T,
            rescale.min=0,
            rescale.max=1,
            # point.geom=F,
            tile.geom=T) +
  scale_size(range=c(4, 8))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/GeneWeight_5000gene.top50_expr_3.pdf", height=30, width=20)
plot(plot_detailed)
dev.off()

qpdf::pdf_combine(input=paste0("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/GeneWeight_5000gene.top50_expr_",c(1:3),".pdf"),
                  output="~/Project_PBMCage/For_or_From_Xu/Xu_Plots/GeneWeight_5000gene.top50_expr.pdf")
unlink(paste0("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/GeneWeight_5000gene.top50_expr_",c(1:3),".pdf")) # delete the tempts

#####################################



### Map the 5000 genes in the RNA.Correlation.results
#####################################
###

library(dplyr)
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000_newCT.rds")
Gene_weights_5000=Gene_weights_5000 %>% dplyr::select(symbol, mean)
colnames(Gene_weights_5000)=c("gene","weight")

### Load the Correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
Cor_analysis_DF=Cor_analysis_DF %>% subset(analysis=="all")

### Merge
merged_DF=dplyr::left_join(Cor_analysis_DF, Gene_weights_5000, by="gene")
merged_DF=merged_DF[!is.na(merged_DF$weight),]

### Plot the weight-vs.-rho dotplot
plot_=
  ggplot(merged_DF, aes(x=rho, y=weight, color=celltypes)) +
  ggrastr::geom_point_rast(size=0.5, shape=20, alpha=1) +
  facet_wrap(~celltypes, ncol=4, scales="fixed") +
  ggsci::scale_color_d3() +
  theme_classic() +
  labs(x=expression(Spearman~rho), y="Attention weight", title=NULL) +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="none",
        legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/GeneWeight_5000gene.Corr_wAges.CorPlots.pdf", height=6, width=8)
plot(plot_)
dev.off()

### Jaccard analysis
# prepare the jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}
# process the genelists
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000_newCT.rds")
XuGenes=Gene_weights_5000 %>% dplyr::select(symbol, mean) %>% arrange(desc(mean)) %>% select(symbol) %>% unlist() %>% unname()

Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
Cor_analysis_DF_sep=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05) %>% arrange(desc(abs(rho)))
Cor.Gene.list=split(Cor_analysis_DF_sep$gene, Cor_analysis_DF_sep$celltypes)
Cor_analysis_DF_sum=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05) %>% 
  group_by(gene) %>% summarize_at(colnames(.)[2:7], mean) %>% 
  arrange(desc(abs(rho))) %>% select(gene) %>% unlist() %>% unname()
Cor.Gene.list=c(Cor.Gene.list, list(All=Cor_analysis_DF_sum))
genesets_length=unlist(lapply(Cor.Gene.list, function(x) min(length(x), length(XuGenes))))
Cor.Gene.list=lapply(1:length(genesets_length), function(i) Cor.Gene.list[[i]][1:genesets_length[i]])

# run
RESULTS=list()

for (i in 1:length(genesets_length)) {
  results=list()
  cor_geneset_=Cor.Gene.list[[i]]
  for (j in 1:length(cor_geneset_)) {
    results[[j]]=jaccard(XuGenes[1:j], cor_geneset_[1:j])
  }
  RESULTS[[i]]=unlist(results)
}
JACCARD_results_list=lapply(1:length(RESULTS), function(i) data.frame(rank=1:length(RESULTS[[i]]), jaccard=RESULTS[[i]], geneset=names(genesets_length)[i]))
JACCARD_results=data.table::rbindlist(JACCARD_results_list)
write.table(JACCARD_results, "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.RNAcorResults_JaccardSimilarity.txt", sep="\t")

### Plot
# plot_=
  ggplot(JACCARD_results, aes(x=rank, y=jaccard, color=geneset)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=1) +
  # facet_wrap(~geneset, ncol=3, scales="free_x") +
  ggsci::scale_color_d3("category20c") +
  theme_classic() +
  # ylim(0,1) +
  labs(x="geneN", y="Jaccard similarity index", title=NULL) +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="right",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title="celltype"))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.RNAcorResults_JaccardSimilarity.pdf", width=5, height=3)
plot(plot_)
dev.off()

#####################################



### Plot the real correlation plots of the top50 genes to the ages
#####################################
###

library(Seurat)
bulk_rough=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")

### Plot those with clear correlation according to the RNA.correlation.DF first
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
Cor_analysis_DF_sigGenes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & cor_pval<0.05) %>%
  arrange(desc(abs(rho))) %>% select(gene) %>% unlist() %>% unname()
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
Top_genes=Gene_weights_5000 %>% 
  arrange(desc(mean)) %>%
  subset(symbol %in% Cor_analysis_DF_sigGenes) %>% select(symbol) %>% unlist() %>% unname()
Top_genes=gsub("_","-",Top_genes)

Genes_Expr=FetchData(bulk_rough, vars=c("age","sex","donor_id","Annot.rough",Top_genes), layer="data")
Genes_Expr_long=tidyr::pivot_longer(Genes_Expr, 
                                    cols=colnames(Genes_Expr)[!grepl("age|sex|donor_id|Annot.rough", colnames(Genes_Expr))], 
                                    names_to="gene", 
                                    values_to="Expression")

# for (i in 1:length(Top_genes)) { # too many genes, so plot only top30 correlated ones
for (i in 1:300) {
  plot_=
    ggplot(subset(Genes_Expr_long, gene==Top_genes[i]), aes(x=age, y=Expression, color=Annot.rough)) +
    ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
    theme_classic() +
    ggsci::scale_color_d3() +
    facet_wrap(~Annot.rough, scales="free_y") + 
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.position="none",
          legend.title=element_blank(),
          title=element_text(size=10)) +
    guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
    labs(title=gsub("_ENSG[0-9].*","",Top_genes[i])) +
    geom_smooth(method="loess", fill="gray93", alpha=0.5, show.legend=FALSE, linewidth=0.5) +
    ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5)
  
  pdf(paste0("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/5000Gene_corPlots/CorPlots_",i,"_",gsub("-","_",Top_genes[i]),".pdf"), height=4, width=6)
  plot(plot_)
  dev.off()
}

#####################################



### Map the Xu Top Weighted genes to the Sig.DMRs
#####################################
###

library(ClusterGVis)
library(GseaVis)
library(clusterProfiler)

### Load the sig.DMRs TermTOGene dataframe
TERM2GENE_df=readRDS("~/Project_PBMCage/Other_omics/Methy_results/age_methylation_results_sigDMRs.rds")
TERM2GENE_df=subset(TERM2GENE_df, orientation!="sig")
TERM2GENE_df$term=gsub("neg","hypo",TERM2GENE_df$term)
TERM2GENE_df$term=gsub("pos","hyper",TERM2GENE_df$term)

### Take the Xu top weighted genes
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
# take the top 10% genes
Top_genes=Gene_weights_5000 %>%
  slice_min(mean, n=round(nrow(Gene_weights_5000)*0.1)) %>%
  dplyr::select(symbol) %>%
  unlist() %>%
  unname()
Top_genes=gsub("_","-",Top_genes)
# or, take the whole ranking genelist
geneList=Gene_weights_5000 %>% arrange(desc(mean)) %>% mutate(mean=mean*981) %>% dplyr::select(mean) %>% unlist()
names(geneList)=Gene_weights_5000 %>% arrange(desc(mean)) %>% mutate(mean=mean*981) %>% dplyr::select(symbol) %>% unlist()

### Enrich analysis
enrich_result=enricher(Top_genes,
                       TERM2GENE=TERM2GENE_df,
                       minGSSize=-Inf,
                       maxGSSize=Inf)
# plot
show_cat=subset(enrich_result@result, p.adjust<0.05)
enrichplot=
  dotplot(enrich_result, showCategory=nrow(show_cat),
          label_format=50) +
  scale_colour_gradient(limits=quantile(show_cat$p.adjust)[c(1,5)],
                        low="pink", high="brown4",
                        labels= ~sprintf(fmt="%0.01e", .)) +
  scale_size_continuous(breaks=c(min(show_cat$Count),
                                 round((min(show_cat$Count)+max(show_cat$Count))/2),
                                 max(show_cat$Count)),
                        labels=c(min(show_cat$Count),
                                 round((min(show_cat$Count)+max(show_cat$Count))/2),
                                 max(show_cat$Count)),
                        range=c(2,6)) +
  ggtitle("DMR enrichment of Top weighted genes") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        title=element_text(size=11)) +
  scale_x_continuous(expand=expansion(mult=c(0, 0), add=c(0.01, 0.01)))

### GSEA analysis
gsea_result=GSEA(geneList,
                 TERM2GENE=TERM2GENE_df,
                 pvalueCutoff=0.1,
                 minGSSize=-Inf,
                 maxGSSize=Inf)
# plot
gsea_plot=
  gseaNb(gsea_result,
         geneSetID=gsea_result@result$ID,
         subPlot=2,
         termWidth=35,
         legend.position=c(0.8,0.8),
         addPval=T,
         pvalX=0.05,
         pvalY=0.05)

pdf("~/Project_PBMCage/For_or_From_Xu/Mapping_5000Genes.10percentTopWeights_toDMRs.pdf", height=5, width=7)
plot(enrichplot)
plot(gsea_plot)
dev.off()

#####################################



### Jaccard similarity analysis between the Xu_5000 geneset and the other previously defined age-related genesets
#####################################
###

library(dplyr)

### Prepare the jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}

### Load genesets for comparison
### !: done in "Comparison_on_donor.R"
gene_list_toCompare=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list_toCompare$Xu_5000.ranked
# cut the genesets all to the same length as Xu_5000 to compare
for (i in 1:length(gene_list_toCompare)) {
  gene_list_toCompare[[i]]=gene_list_toCompare[[i]][1:min(length(gene_list_toCompare[[i]]), length(Xu_5000.ranked))]
}
length(gene_list_toCompare) # to check

### Set the sliding window
JACCARD_results=list()
geneset_length=lapply(gene_list_toCompare, function(x) min(length(x), length(Xu_5000.ranked)))
for (i in 1:length(gene_list_toCompare)) {
  geneset_from_published=gene_list_toCompare[[i]]
  Xu_5000.ranked_Cut.To.the.Same.Length=Xu_5000.ranked[1:length(geneset_from_published)]
  results=list()
  for (j in 1:length(geneset_from_published)) {
    results[[j]]=jaccard(Xu_5000.ranked_Cut.To.the.Same.Length[1:j], geneset_from_published[1:j])
  }
  JACCARD_results[[i]]=data.frame(rank=1:length(geneset_from_published), jaccard=unlist(results), info=names(gene_list_toCompare)[i])
}
JACCARD_DF=data.table::rbindlist(JACCARD_results)
# relevel and reorder
JACCARD_DF$geneset=as.factor(JACCARD_DF$info)
levels(JACCARD_DF$geneset) # check
levels(JACCARD_DF$geneset)=c("de Magalhães et al. (2009)","GenAge Build 21","Dørum et al. (2024)","Yang et al. (2015)","Chuffa et al. (2022)",
                             "Peters et al. (2015)","Shokhirev et al. (2021); cubist-DG","Shokhirev et al. (2021); cubist-VG",
                             "Shokhirev et al. (2021); Rf-DG","Shokhirev et al. (2021); Rf-VG","Ren et al. (2020)","Zhu et al. (2023)",
                             "AgeGene5K")

write.table(JACCARD_DF, "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_JaccardSimilarity.txt", sep="\t")

### Plot
JACCARD_DF=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_JaccardSimilarity.txt", sep="\t")
JACCARD_DF=subset(JACCARD_DF, geneset!="AgeGene5K")

# plot_=
  ggplot(JACCARD_DF, aes(x=rank, y=jaccard, color=geneset)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=1) +
  facet_wrap(~geneset, ncol=3, scales="free_x") +
  ggsci::scale_color_d3("category20c") +
  theme_classic() +
  labs(x="geneN", y="Jaccard similarity index", title=NULL) +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="none",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1)))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.Published_Genesets_JaccardSimilarity.pdf", width=9, height=6)
plot(plot_)
dev.off()

#####################################



### Analyze the similarity/clusters of aging gene inflatuation
#####################################
###

library(dplyr)
library(Seurat)
library(ggrepel)
library(ClusterGVis)
if(!require("dtwclust")) install.packages("dtwclust")

### Load the data
bulk_rough=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
pseudobulkobj_all=AggregateExpression(bulk_rough, assays="RNA", return.seurat=T, group.by=c("age"))
pseudobulkobj_all$age=as.numeric(colnames(pseudobulkobj_all))

bulk_rough_females=subset(bulk_rough, sex=="female")
pseudobulkobj_females=AggregateExpression(bulk_rough_females, assays="RNA", return.seurat=T, group.by=c("age","sex"))
pseudobulkobj_females$age=as.numeric(unlist(lapply(strsplit(colnames(pseudobulkobj_females), split="_"), function(x) x[1])))
pseudobulkobj_females$sex=unlist(lapply(strsplit(colnames(pseudobulkobj_females), split="_"), function(x) x[2]))

bulk_rough_males=subset(bulk_rough, sex=="male")
pseudobulkobj_males=AggregateExpression(bulk_rough_males, assays="RNA", return.seurat=T, group.by=c("age","sex"))
pseudobulkobj_males$age=as.numeric(unlist(lapply(strsplit(colnames(pseudobulkobj_males), split="_"), function(x) x[1])))
pseudobulkobj_males$sex=unlist(lapply(strsplit(colnames(pseudobulkobj_males), split="_"), function(x) x[2]))

### Load Xu's geneset
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
Cor_analysis_DF_sigGenes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & cor_pval<0.05) %>%
  arrange(desc(abs(rho))) %>% select(gene) %>% unlist() %>% unname()
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
genes_to_analyze=Gene_weights_5000 %>% 
  arrange(desc(mean)) %>%
  subset(symbol %in% Cor_analysis_DF_sigGenes) %>% select(symbol) %>% unlist() %>% unname()

CK_list=list()
### Run ClusterGVis on both sex
pseudo_data=LayerData(pseudobulkobj_all, layer="data")
pseudo_data=pseudo_data[match(genes_to_analyze, rownames(pseudo_data)), ] # take only those genes with sig. correlation with ages
# run getClusters to determine the best cluster.num
getClusters(exp=pseudo_data) # so take cluster.num=9
ck=clusterData(exp=pseudo_data,
               cluster.method="kmeans",
               cluster.num=9)
CK_list=c(CK_list, list(ck))

### Run ClusterGVis on females only
pseudo_data=LayerData(pseudobulkobj_females, layer="data")
pseudo_data=pseudo_data[match(genes_to_analyze, rownames(pseudo_data)), ] # take only those genes with sig. correlation with ages
# run getClusters to determine the best cluster.num
getClusters(exp=pseudo_data) # so take cluster.num=9
ck=clusterData(exp=pseudo_data,
               cluster.method="kmeans",
               cluster.num=9)
CK_list=c(CK_list, list(ck))

### Run ClusterGVis on males only
pseudo_data=LayerData(pseudobulkobj_males, layer="data")
pseudo_data=pseudo_data[match(genes_to_analyze, rownames(pseudo_data)), ] # take only those genes with sig. correlation with ages
# run getClusters to determine the best cluster.num
getClusters(exp=pseudo_data) # so take cluster.num=9
ck=clusterData(exp=pseudo_data,
               cluster.method="kmeans",
               cluster.num=9)
CK_list=c(CK_list, list(ck))

### Save the ClusterGVis results on both sex, females only, or males only
names(CK_list)=paste0("Xu5000_",c("both","females","males"))
saveRDS(CK_list, "~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Xu.5000.Geneset.rds")

### Plot the clustering and map the genes in the RNA.Correlation
CK_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Xu.5000.Geneset.rds")
Analysis_list=c("all","females","males")
Subtitle_list=gsub(".*_","",names(CK_list))

plot_ckvis=plot_cluster.Cor=list()
for (i in 1:length(CK_list)) {
  ck=CK_list[[i]]

  # plot the clustering
  plot_ckvis[[i]]=
    visCluster(object=ck,
               plot.type="line", 
               mline.size=0.5) +
    labs(title=" ", subtitle=paste0("<",Subtitle_list[i],">")) +
    facet_wrap(facets="cluster_name", ncol=3) +
    theme(title=element_text(size=11), 
          plot.subtitle=element_text(size=10, face="italic"),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=9),
          strip.text.x=element_text(size=9)) +
    scale_x_discrete(breaks=seq(18, 100, 10))
  
  # extract the genes belong to each cluster
  gene_cluster_info=ck$wide.res %>% select(gene, cluster)
  
  # check those clusters in RNA.Correlation DF
  Cor_analysis_DF.sig=Cor_analysis_DF %>% subset(analysis==Analysis_list[i] & rho_pval<0.05 & cor_pval<0.05)
  merged_cluster_cor_df=dplyr::left_join(Cor_analysis_DF.sig, gene_cluster_info, by="gene") %>%
    subset(!is.na(cluster))
  merged_cluster_cor_top=merged_cluster_cor_df %>% 
    mutate(rank=rownames(.)) %>% 
    group_by(cluster) %>% 
    slice_max(abs(rho), n=5) %>% 
    ungroup() %>% 
    select(rank) %>% unlist() %>% unname()
  merged_cluster_cor_df=merged_cluster_cor_df %>% 
    mutate(mark=ifelse(rownames(.) %in% merged_cluster_cor_top, gene, NA)) %>%
    mutate(cluster=paste0("cluster ",cluster))
  
  # map the genes in each cluster to RNA.Correlation
  plot_cluster.Cor[[i]]=
    ggplot(merged_cluster_cor_df, aes(x=rho, y=-log10(rho_pval), color=celltypes)) +
    ggrastr::geom_point_rast(size=1, shape=20, alpha=0.5) +
    facet_wrap(~cluster, ncol=3) +
    theme_classic() +
    ggsci::scale_color_d3() +
    labs(y=expression(-log[10]~pval), x=expression(Spearman~rho), title=" ", subtitle=paste0("<",Subtitle_list[i],">")) +
    guides(color=guide_legend(title="celltype", override.aes=list(size=2, alpha=1))) +
    theme(title=element_text(size=11),
          plot.subtitle=element_text(size=11, face="italic"),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=9),
          axis.title.x=element_text(size=9),
          legend.position="right",
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          strip.background=element_rect(color="transparent", fill="transparent", linewidth=1.5, linetype=NULL),
          strip.text.x=element_text(size=10, color="black", face="plain")) +
    geom_smooth(aes(group=celltypes), method="loess", se=FALSE, alpha=0.1, show.legend=FALSE, linewidth=0.1) +
    geom_text_repel(aes(label=mark), size=3, show.legend=FALSE, color="black")
}

### Plot the ClusterGVis results
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/GeneWeight_5000gene.expr_ClusterGVis.pdf")
for (i in 1:length(CK_list)) {
  plot(plot_ckvis[[i]])
  plot(plot_cluster.Cor[[i]])
}
dev.off()

#####################################



### Analyze the similarity/clusters of scRNA_AgeGene gene inflatuation for comparison
#####################################
###

library(dplyr)
library(Seurat)
library(ggrepel)
library(ClusterGVis)
if(!require("dtwclust")) install.packages("dtwclust")

bulk_rough=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
pseudobulkobj_all=AggregateExpression(bulk_rough, assays="RNA", return.seurat=T, group.by=c("age"))
pseudobulkobj_all$age=as.numeric(colnames(pseudobulkobj_all))

bulk_rough_females=subset(bulk_rough, sex=="female")
pseudobulkobj_females=AggregateExpression(bulk_rough_females, assays="RNA", return.seurat=T, group.by=c("age","sex"))
pseudobulkobj_females$age=as.numeric(unlist(lapply(strsplit(colnames(pseudobulkobj_females), split="_"), function(x) x[1])))
pseudobulkobj_females$sex=unlist(lapply(strsplit(colnames(pseudobulkobj_females), split="_"), function(x) x[2]))

bulk_rough_males=subset(bulk_rough, sex=="male")
pseudobulkobj_males=AggregateExpression(bulk_rough_males, assays="RNA", return.seurat=T, group.by=c("age","sex"))
pseudobulkobj_males$age=as.numeric(unlist(lapply(strsplit(colnames(pseudobulkobj_males), split="_"), function(x) x[1])))
pseudobulkobj_males$sex=unlist(lapply(strsplit(colnames(pseudobulkobj_males), split="_"), function(x) x[2]))

### Take all the scRNA_AgeGene together, regardless which celltype they are from
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
Cor_analysis_DF_sigGenes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & cor_pval<0.05) %>%
  arrange(desc(abs(rho))) %>% select(gene) %>% unlist() %>% unname()
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000_newCT.rds")
scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
scRNA_ageGenes=scRNA_ageGenes %>% 
  mutate(Gene.name=gsub("HLA\\.","HLA-",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.AS1","-AS1",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.2HG","-2HG",Gene.name)) %>% 
  group_by(Gene.name) %>% summarize_at("Vip", mean) %>% slice_max(Vip, n=nrow(Gene_weights_5000)) %>% as.data.frame()
colnames(scRNA_ageGenes)=c("gene","weight")
genes_to_analyze=scRNA_ageGenes %>% 
  subset(gene %in% Cor_analysis_DF_sigGenes) %>% select(gene) %>% unlist() %>% unname()

CK_list=list()
CK_list_name=c()
# run ClusterGVis on both sex
pseudo_data=LayerData(pseudobulkobj_all, layer="data")
pseudo_data=pseudo_data[match(genes_to_analyze, rownames(pseudo_data)), ] # take only those genes with sig. correlation with ages
getClusters(exp=pseudo_data) # so take cluster.num=9
ck=clusterData(exp=pseudo_data,
               cluster.method="kmeans",
               cluster.num=9)
CK_list=c(CK_list, list(ck))
CK_list_name=c(CK_list_name, "alltypes||both")

# run ClusterGVis on females only
pseudo_data=LayerData(pseudobulkobj_females, layer="data")
pseudo_data=pseudo_data[match(genes_to_analyze, rownames(pseudo_data)), ] # take only those genes with sig. correlation with ages
getClusters(exp=pseudo_data) # so take cluster.num=9
ck=clusterData(exp=pseudo_data,
               cluster.method="kmeans",
               cluster.num=9)
CK_list=c(CK_list, list(ck))
CK_list_name=c(CK_list_name, "alltypes||females")

# run ClusterGVis on males only
pseudo_data=LayerData(pseudobulkobj_males, layer="data")
pseudo_data=pseudo_data[match(genes_to_analyze, rownames(pseudo_data)), ] # take only those genes with sig. correlation with ages
getClusters(exp=pseudo_data) # so take cluster.num=9
ck=clusterData(exp=pseudo_data,
               cluster.method="kmeans",
               cluster.num=9)
CK_list=c(CK_list, list(ck))
CK_list_name=c(CK_list_name, "alltypes||males")

### Take the scRNA_AgeGene in specific celltype
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
Cor_analysis_DF_sigGenes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & cor_pval<0.05) %>%
  arrange(desc(abs(rho))) %>% select(gene) %>% unlist() %>% unname()
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000_newCT.rds")
scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
scRNA_ageGenes=scRNA_ageGenes %>% 
  mutate(Gene.name=gsub("HLA\\.","HLA-",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.AS1","-AS1",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.2HG","-2HG",Gene.name)) %>% 
  group_by(Cell.type) %>% slice_max(Vip, n=nrow(Gene_weights_5000)) %>% as.data.frame()

for (i in 1:length(unique(scRNA_ageGenes$Cell.type))) {
  genes_to_analyze=scRNA_ageGenes %>% 
    subset(Cell.type==unique(scRNA_ageGenes$Cell.type)[i]) %>%
    subset(Gene.name %in% Cor_analysis_DF_sigGenes) %>% select(Gene.name) %>% unlist() %>% unname()
  
  # run ClusterGVis on both sex
  pseudo_data=LayerData(pseudobulkobj_all, layer="data")
  pseudo_data=pseudo_data[match(genes_to_analyze, rownames(pseudo_data)), ] # take only those genes with sig. correlation with ages
  getClusters(exp=pseudo_data) # so take cluster.num=9
  ck=clusterData(exp=pseudo_data,
                 cluster.method="kmeans",
                 cluster.num=9)
  CK_list=c(CK_list, list(ck))
  CK_list_name=c(CK_list_name, paste0(unique(scRNA_ageGenes$Cell.type)[i],"||both"))
  
  # run ClusterGVis on females only
  pseudo_data=LayerData(pseudobulkobj_females, layer="data")
  pseudo_data=pseudo_data[match(genes_to_analyze, rownames(pseudo_data)), ] # take only those genes with sig. correlation with ages
  getClusters(exp=pseudo_data) # so take cluster.num=9
  ck=clusterData(exp=pseudo_data,
                 cluster.method="kmeans",
                 cluster.num=9)
  CK_list=c(CK_list, list(ck))
  CK_list_name=c(CK_list_name, paste0(unique(scRNA_ageGenes$Cell.type)[i],"||females"))
  
  # run ClusterGVis on males only
  pseudo_data=LayerData(pseudobulkobj_males, layer="data")
  pseudo_data=pseudo_data[match(genes_to_analyze, rownames(pseudo_data)), ] # take only those genes with sig. correlation with ages
  getClusters(exp=pseudo_data) # so take cluster.num=9
  ck=clusterData(exp=pseudo_data,
                 cluster.method="kmeans",
                 cluster.num=9)
  CK_list=c(CK_list, list(ck))
  CK_list_name=c(CK_list_name, paste0(unique(scRNA_ageGenes$Cell.type)[i],"||males"))
}
names(CK_list)=CK_list_name

### Save all the ClusterGVis results in the scRNA_AgeGene geneset
saveRDS(CK_list, "~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_scRNA_AgeGene.Geneset.rds")

### Plot the clustering and map the genes in the RNA.Correlation
CK_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_scRNA_AgeGene.Geneset.rds")
CK_list=CK_list[1:3]
Analysis_list=c("all","females","males")
Subtitle_list=gsub(".*\\||","",names(CK_list))

plot_ckvis=plot_cluster.Cor=list()
for (i in 1:length(CK_list)) {
  ck=CK_list[[i]]
  
  # plot the clustering
  plot_ckvis[[i]]=
    visCluster(object=ck,
               plot.type="line", 
               mline.size=0.5) +
    labs(title=" ", subtitle=paste0("<",Subtitle_list[i],">")) +
    facet_wrap(facets="cluster_name", ncol=3) +
    theme(title=element_text(size=11), 
          plot.subtitle=element_text(size=10, face="italic"),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=9),
          strip.text.x=element_text(size=9)) +
    scale_x_discrete(breaks=seq(18, 100, 10))
  
  # extract the genes belong to each cluster
  gene_cluster_info=ck$wide.res %>% select(gene, cluster)
  
  # check those clusters in RNA.Correlation DF
  Cor_analysis_DF.sig=Cor_analysis_DF %>% subset(analysis==Analysis_list[i] & rho_pval<0.05 & cor_pval<0.05)
  merged_cluster_cor_df=dplyr::left_join(Cor_analysis_DF.sig, gene_cluster_info, by="gene") %>%
    subset(!is.na(cluster))
  merged_cluster_cor_top=merged_cluster_cor_df %>% 
    mutate(rank=rownames(.)) %>% 
    group_by(cluster) %>% 
    slice_max(abs(rho), n=5) %>% 
    ungroup() %>% 
    select(rank) %>% unlist() %>% unname()
  merged_cluster_cor_df=merged_cluster_cor_df %>% 
    mutate(mark=ifelse(rownames(.) %in% merged_cluster_cor_top, gene, NA)) %>%
    mutate(cluster=paste0("cluster ",cluster))
  
  # map the genes in each cluster to RNA.Correlation
  plot_cluster.Cor[[i]]=
    ggplot(merged_cluster_cor_df, aes(x=rho, y=-log10(rho_pval), color=celltypes)) +
    ggrastr::geom_point_rast(size=1, shape=20, alpha=0.5) +
    facet_wrap(~cluster, ncol=3) +
    theme_classic() +
    ggsci::scale_color_d3() +
    labs(y=expression(-log[10]~pval), x=expression(Spearman~rho), title=" ", subtitle=paste0("<",Subtitle_list[i],">")) +
    guides(color=guide_legend(title="celltype", override.aes=list(size=2, alpha=1))) +
    theme(title=element_text(size=11),
          plot.subtitle=element_text(size=11, face="italic"),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=9),
          axis.title.x=element_text(size=9),
          legend.position="right",
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          strip.background=element_rect(color="transparent", fill="transparent", linewidth=1.5, linetype=NULL),
          strip.text.x=element_text(size=10, color="black", face="plain")) +
    geom_smooth(aes(group=celltypes), method="loess", se=FALSE, alpha=0.1, show.legend=FALSE, linewidth=0.1) +
    geom_text_repel(aes(label=mark), size=3, show.legend=FALSE, color="black")
}

### Plot the ClusterGVis results
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/scRNAAgeGenes.expr.allcelltypes_ClusterGVis.pdf")
for (i in 1:length(CK_list)) {
  plot(plot_ckvis[[i]])
  plot(plot_cluster.Cor[[i]])
}
dev.off()

#####################################



### Analyze the similarity of the time-serial clusters between the scRNA_AgeGene geneset and Xu's geneset
#####################################
###

library(dplyr)

### Prepare the jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}

### Extract the genes in each time-serial clusters
scRNA_AgeGene_ckList=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_scRNA_AgeGene.Geneset.rds")
Xu_5000Gene_ckList=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Xu.5000.Geneset.rds")
JACCARD_RESULTs=list()

### Jaccard analysis on subset-both
# subset
scRNA_AgeGene_ckList_both=scRNA_AgeGene_ckList[grepl("\\|\\|both",names(scRNA_AgeGene_ckList))]
scRNA_AgeGene_ckList_Clusters=list()
for (i in 1:length(scRNA_AgeGene_ckList_both)) {
  scRNA_AgeGene_ckList_Clusters[[i]]=split(scRNA_AgeGene_ckList_both[[i]]$long.res$gene, scRNA_AgeGene_ckList_both[[i]]$long.res$cluster)
}
names(scRNA_AgeGene_ckList_Clusters)=names(scRNA_AgeGene_ckList_both)
Xu_5000Gene_ckList_both=Xu_5000Gene_ckList[grepl("both",names(Xu_5000Gene_ckList))]
Xu_5000Gene_ckList_Cluster=split(Xu_5000Gene_ckList_both$Xu5000_both$long.res$gene, Xu_5000Gene_ckList_both$Xu5000_both$long.res$cluster)

# analysis
results=results_name=set_size=c()
for (i in 1:length(Xu_5000Gene_ckList_Cluster)) {
  Xu_list=Xu_5000Gene_ckList_Cluster[[i]]
  Xu_list_name=paste0("Xu, ",names(Xu_5000Gene_ckList_Cluster)[i])
  for (j in 1:length(scRNA_AgeGene_ckList_Clusters)) {
    scRNA_LIST=scRNA_AgeGene_ckList_Clusters[[j]]
    scRNA_celltype=gsub("\\|\\|.*","",names(scRNA_AgeGene_ckList_Clusters)[j])
    for (k in 1:length(scRNA_LIST)) {
      scRNA_list=scRNA_LIST[[k]]
      scRNA_list_name=paste0(names(scRNA_AgeGene_ckList_Clusters)[j],", ",names(scRNA_LIST)[k])
      results=c(results, jaccard(Xu_list, scRNA_list))
      results_name=c(results_name, paste0(Xu_list_name, " [VS] ", scRNA_list_name))
      set_size=c(set_size, length(intersect(Xu_list, scRNA_list)[!duplicated(intersect(Xu_list, scRNA_list))]))
    }
  }
}
Jaccard_result_df=data.frame(comparison=results_name,
                             jaccard=results,
                             overlappingN=set_size)
Jaccard_result_df=Jaccard_result_df %>% 
  mutate(Xu_part=gsub(" \\[VS\\] .*","",comparison),
         scRNA_AgeGene_part=gsub(".* \\[VS\\] ","",comparison)) %>%
  mutate(Xu_cluster=gsub("Xu, ","",Xu_part),
         scRNA_AgeGene_celltype=gsub("\\|\\|.*","",scRNA_AgeGene_part),
         scRNA_AgeGene_cluster=gsub(".*, ","",scRNA_AgeGene_part)) %>%
  mutate(analysis="both") %>%
  select(analysis, comparison, Xu_part, scRNA_AgeGene_part, scRNA_AgeGene_celltype, Xu_cluster, scRNA_AgeGene_cluster, jaccard, overlappingN)

# save the subset analysis result
JACCARD_RESULTs=c(JACCARD_RESULTs, list(Jaccard_result_df))


### Jaccard analysis on subset-females
# subset
scRNA_AgeGene_ckList_both=scRNA_AgeGene_ckList[grepl("\\|\\|females",names(scRNA_AgeGene_ckList))]
scRNA_AgeGene_ckList_Clusters=list()
for (i in 1:length(scRNA_AgeGene_ckList_both)) {
  scRNA_AgeGene_ckList_Clusters[[i]]=split(scRNA_AgeGene_ckList_both[[i]]$long.res$gene, scRNA_AgeGene_ckList_both[[i]]$long.res$cluster)
}
names(scRNA_AgeGene_ckList_Clusters)=names(scRNA_AgeGene_ckList_both)
Xu_5000Gene_ckList_both=Xu_5000Gene_ckList[grepl("females",names(Xu_5000Gene_ckList))]
Xu_5000Gene_ckList_Cluster=split(Xu_5000Gene_ckList_both$Xu5000_females$long.res$gene, Xu_5000Gene_ckList_both$Xu5000_females$long.res$cluster)

# analysis
results=results_name=set_size=c()
for (i in 1:length(Xu_5000Gene_ckList_Cluster)) {
  Xu_list=Xu_5000Gene_ckList_Cluster[[i]]
  Xu_list_name=paste0("Xu, ",names(Xu_5000Gene_ckList_Cluster)[i])
  for (j in 1:length(scRNA_AgeGene_ckList_Clusters)) {
    scRNA_LIST=scRNA_AgeGene_ckList_Clusters[[j]]
    scRNA_celltype=gsub("\\|\\|.*","",names(scRNA_AgeGene_ckList_Clusters)[j])
    for (k in 1:length(scRNA_LIST)) {
      scRNA_list=scRNA_LIST[[k]]
      scRNA_list_name=paste0(names(scRNA_AgeGene_ckList_Clusters)[j],", ",names(scRNA_LIST)[k])
      results=c(results, jaccard(Xu_list, scRNA_list))
      results_name=c(results_name, paste0(Xu_list_name, " [VS] ", scRNA_list_name))
      set_size=c(set_size, length(intersect(Xu_list, scRNA_list)[!duplicated(intersect(Xu_list, scRNA_list))]))
    }
  }
}
Jaccard_result_df=data.frame(comparison=results_name,
                             jaccard=results,
                             overlappingN=set_size)
Jaccard_result_df=Jaccard_result_df %>% 
  mutate(Xu_part=gsub(" \\[VS\\] .*","",comparison),
         scRNA_AgeGene_part=gsub(".* \\[VS\\] ","",comparison)) %>%
  mutate(Xu_cluster=gsub("Xu, ","",Xu_part),
         scRNA_AgeGene_celltype=gsub("\\|\\|.*","",scRNA_AgeGene_part),
         scRNA_AgeGene_cluster=gsub(".*, ","",scRNA_AgeGene_part)) %>%
  mutate(analysis="females") %>%
  select(analysis, comparison, Xu_part, scRNA_AgeGene_part, scRNA_AgeGene_celltype, Xu_cluster, scRNA_AgeGene_cluster, jaccard, overlappingN)

# save the subset analysis result
JACCARD_RESULTs=c(JACCARD_RESULTs, list(Jaccard_result_df))


### Jaccard analysis on subset-males
# subset
scRNA_AgeGene_ckList_both=scRNA_AgeGene_ckList[grepl("\\|\\|males",names(scRNA_AgeGene_ckList))]
scRNA_AgeGene_ckList_Clusters=list()
for (i in 1:length(scRNA_AgeGene_ckList_both)) {
  scRNA_AgeGene_ckList_Clusters[[i]]=split(scRNA_AgeGene_ckList_both[[i]]$long.res$gene, scRNA_AgeGene_ckList_both[[i]]$long.res$cluster)
}
names(scRNA_AgeGene_ckList_Clusters)=names(scRNA_AgeGene_ckList_both)
Xu_5000Gene_ckList_both=Xu_5000Gene_ckList[grepl("males",names(Xu_5000Gene_ckList))]
Xu_5000Gene_ckList_Cluster=split(Xu_5000Gene_ckList_both$Xu5000_males$long.res$gene, Xu_5000Gene_ckList_both$Xu5000_males$long.res$cluster)

# analysis
results=results_name=set_size=c()
for (i in 1:length(Xu_5000Gene_ckList_Cluster)) {
  Xu_list=Xu_5000Gene_ckList_Cluster[[i]]
  Xu_list_name=paste0("Xu, ",names(Xu_5000Gene_ckList_Cluster)[i])
  for (j in 1:length(scRNA_AgeGene_ckList_Clusters)) {
    scRNA_LIST=scRNA_AgeGene_ckList_Clusters[[j]]
    scRNA_celltype=gsub("\\|\\|.*","",names(scRNA_AgeGene_ckList_Clusters)[j])
    for (k in 1:length(scRNA_LIST)) {
      scRNA_list=scRNA_LIST[[k]]
      scRNA_list_name=paste0(names(scRNA_AgeGene_ckList_Clusters)[j],", ",names(scRNA_LIST)[k])
      results=c(results, jaccard(Xu_list, scRNA_list))
      results_name=c(results_name, paste0(Xu_list_name, " [VS] ", scRNA_list_name))
      set_size=c(set_size, length(intersect(Xu_list, scRNA_list)[!duplicated(intersect(Xu_list, scRNA_list))]))
    }
  }
}
Jaccard_result_df=data.frame(comparison=results_name,
                             jaccard=results,
                             overlappingN=set_size)
Jaccard_result_df=Jaccard_result_df %>% 
  mutate(Xu_part=gsub(" \\[VS\\] .*","",comparison),
         scRNA_AgeGene_part=gsub(".* \\[VS\\] ","",comparison)) %>%
  mutate(Xu_cluster=gsub("Xu, ","",Xu_part),
         scRNA_AgeGene_celltype=gsub("\\|\\|.*","",scRNA_AgeGene_part),
         scRNA_AgeGene_cluster=gsub(".*, ","",scRNA_AgeGene_part)) %>%
  mutate(analysis="males") %>%
  select(analysis, comparison, Xu_part, scRNA_AgeGene_part, scRNA_AgeGene_celltype, Xu_cluster, scRNA_AgeGene_cluster, jaccard, overlappingN)

# save the subset analysis result
JACCARD_RESULTs=c(JACCARD_RESULTs, list(Jaccard_result_df))

### Arrange the results
jaccard_final=data.table::rbindlist(JACCARD_RESULTs)
data.table::fwrite(jaccard_final, "~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Xu.vs.scRNAgenesets_JaccardSimilarity.csv.gz", sep="\t")

### Calculate the similarity between the time-serial clusters of Xu_5000geneset vs. scRNA_AgeGeneset
jaccard_final=data.table::fread("~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Xu.vs.scRNAgenesets_JaccardSimilarity.csv.gz", sep="\t")
jaccard_final$Xu_cluster=as.factor(jaccard_final$Xu_cluster)
jaccard_final$scRNA_AgeGene_cluster=as.factor(jaccard_final$scRNA_AgeGene_cluster)
jaccard_final=jaccard_final %>% 
  subset(analysis=="both" & jaccard!=0 & scRNA_AgeGene_celltype=="alltypes")
  
ggplot(jaccard_final, aes(x=Xu_cluster, y=scRNA_AgeGene_cluster, color=jaccard, size=overlappingN)) +
  geom_point() +
  scale_color_gradient(low="white", high="brown3") +
  theme_linedraw() +
  labs(x=NULL, y=NULL, title="Time series characteristics") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        title=element_text(size=10)) +
  guides(size=guide_legend(title="overlap"))

#####################################



### Intra-gene DTW analysis of the scRNA_AgeGene geneset and Xu's geneset
#####################################
###

### Analyze DTW of scRNA_AgeGene.allcelltypes and Xu_5000 on the VerifiCombo dataset on donor_id
# Load the geneset
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
genes_Xu5000=Gene_weights_5000 %>% arrange(desc(mean)) %>% mutate(symbol=gsub("_ENSG","-ENSG",symbol)) %>% 
  select(symbol) %>% unlist() %>% unname()
scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
genes_scRNA=scRNA_ageGenes %>% 
  mutate(Gene.name=gsub("HLA\\.","HLA-",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.AS1","-AS1",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.2HG","-2HG",Gene.name)) %>% 
  group_by(Gene.name) %>%
  summarize_at("Vip", mean) %>%
  arrange(desc(Vip)) %>% 
  select(Gene.name) %>% unlist() %>% unname()
# Load the bulk counts data
pseudobulk=readRDS("~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_rough.rds")
pseudobulk=AggregateExpression(pseudobulk, assays="RNA", return.seurat=T, group.by=c("age"))
pseudobulk$age=as.numeric(colnames(pseudobulk))
pseudobulk_data=LayerData(pseudobulk, layer="counts")
pseudobulk_data=pseudobulk_data[rowSums(pseudobulk_data)!=0,]
# Set the sliding window
steps=seq(50, 4000, 50)
# Intra-gene DTW analysis
df_scRNA.Xu_list=list()
for (i in 1:length(steps)) {
  genes_Xu5000_selected=genes_Xu5000[1:steps[i]]
  genes_scRNA_selected=genes_scRNA[1:steps[i]]
  
  # scRNA_AgeGene geneset
  pseudobulk_data_scRNA.Selected=pseudobulk_data[intersect(genes_scRNA_selected, rownames(pseudobulk_data)),]
  pseudobulk_data_scRNA.Selected_scaled=apply(pseudobulk_data_scRNA.Selected, 1, scale)
  ref_series_scRNA.Selected=scale(apply(pseudobulk_data_scRNA.Selected, 2, sum))
  result_scRNA=apply(pseudobulk_data_scRNA.Selected_scaled, 2, function(x) {alignment=dtw::dtw(x, ref_series_scRNA.Selected, keep=TRUE); alignment$distance})
  df_scRNA=data.frame(geneN=steps[i], gene=names(result_scRNA), dtw_dist=result_scRNA, celltype="all", geneset="scRNA_AgeGene")
  
  # Xu_5000 geneset
  pseudobulk_data_Xu.Selected=pseudobulk_data[intersect(genes_Xu5000_selected, rownames(pseudobulk_data)),]
  pseudobulk_data_Xu.Selected_scaled=apply(pseudobulk_data_Xu.Selected, 1, scale)
  ref_series_Xu.Selected=scale(apply(pseudobulk_data_Xu.Selected, 2, sum))
  result_Xu=apply(pseudobulk_data_Xu.Selected_scaled, 2, function(x) {alignment=dtw::dtw(x, ref_series_Xu.Selected, keep=TRUE); alignment$distance})
  df_Xu=data.frame(geneN=steps[i], gene=names(result_Xu), dtw_dist=result_Xu, celltype="all", geneset="Xu_5000Gene")
  
  # merge the results
  df_scRNA.Xu_list[[i]]=data.table::rbindlist(list(df_scRNA, df_Xu))
}
DF_scRNA.Xu=data.table::rbindlist(df_scRNA.Xu_list)
data.table::fwrite(DF_scRNA.Xu, "~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Comparison_DTW.Intra.gene.csv.gz", sep="\t")

### Analyze DTW of scRNA_AgeGene.per celltype on the VerifiCombo dataset on donor_id
# Load the geneset
scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
scRNA_ageGenes=scRNA_ageGenes %>% 
  mutate(Gene.name=gsub("HLA\\.","HLA-",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.AS1","-AS1",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.2HG","-2HG",Gene.name)) %>% 
  group_by(Cell.type) %>%
  arrange(desc(Vip))
celltypes=unique(scRNA_ageGenes$Cell.type)
# Load the bulk counts data
pseudobulk=readRDS("~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_rough.rds")
pseudobulk=AggregateExpression(pseudobulk, assays="RNA", return.seurat=T, group.by=c("age"))
pseudobulk$age=as.numeric(colnames(pseudobulk))
pseudobulk_data=LayerData(pseudobulk, layer="counts")
pseudobulk_data=pseudobulk_data[rowSums(pseudobulk_data)!=0,]
# Set the sliding window
steps=seq(50, 4000, 50)
# Intra-gene DTW analysis
df_scRNA_allcelltypes_list=list()
for (j in 1:length(celltypes)) {
  
  df_scRNA_allsteps=list()
  pb=txtProgressBar(min=0, max=length(steps), width=60, style=3, char="=")
  for (i in 1:length(steps)) {
    scRNA_ageGenes_subset=subset(scRNA_ageGenes, Cell.type==celltypes[j])
    genes_scRNA_selected=scRNA_ageGenes_subset$Gene.name[1:steps[i]]
    
    pseudobulk_data_scRNA.Selected=pseudobulk_data[intersect(genes_scRNA_selected, rownames(pseudobulk_data)),]
    pseudobulk_data_scRNA.Selected_scaled=apply(pseudobulk_data_scRNA.Selected, 1, scale)
    ref_series_scRNA.Selected=scale(apply(pseudobulk_data_scRNA.Selected, 2, sum))
    result_scRNA=apply(pseudobulk_data_scRNA.Selected_scaled, 2, function(x) {alignment=dtw::dtw(x, ref_series_scRNA.Selected, keep=TRUE); alignment$distance})
    df_scRNA=data.frame(geneN=steps[i], gene=names(result_scRNA), dtw_dist=result_scRNA, celltype=celltypes[j], geneset="scRNA_AgeGene")
    
    df_scRNA_allsteps[[i]]=df_scRNA
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  message(paste0(celltypes[j], " is done!"))
  df_scRNA_allcelltypes_list[[j]]=data.table::rbindlist(df_scRNA_allsteps)
}
# combine the scRNA_AgeGene.per.celltype result with the scRNA_AgeGene.allcelltypes and Xu_5000 results
DF_scRNA=data.table::rbindlist(df_scRNA_allcelltypes_list)
DF_scRNA.Xu=data.table::fread("~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Comparison_DTW.Intra.gene.csv.gz", sep="\t")
DF_scRNA.Xu=data.table::rbindlist(list(DF_scRNA.Xu, DF_scRNA))
data.table::fwrite(DF_scRNA.Xu, "~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Comparison_DTW.Intra.gene_updated.csv.gz", sep="\t")

### Rearrange the results
DF_scRNA.Xu=data.table::fread("~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Comparison_DTW.Intra.gene_updated.csv.gz", sep="\t")
DF_scRNA.Xu$geneset=as.factor(DF_scRNA.Xu$geneset)
levels(DF_scRNA.Xu$geneset) # check
levels(DF_scRNA.Xu$geneset)=c("Zhu et al. (2023)","AgeGene5K")
DF_scRNA.Xu$detailed=as.factor(DF_scRNA.Xu$celltype)
levels(DF_scRNA.Xu$detailed) # check
levels(DF_scRNA.Xu$detailed)=c("Other.ABC","all","Mono.CD14","Mono.CD14-PPBP","Mono.CD16","CD4T.naive","CD4T.Tem","CD4T.Tm","CD4T.Treg",
                               "CD8T.CTL","CD8T.naive","CD8T.Tem","DC.cDC","Mono.inter","Other.Megakaryocytes","B.mem-CRIP1","B.mem-IL32",
                               "B.naive","B.naive-NRGN","NK.FCER1G","NK.GZMH","NK.GZMK","NK.IL7R","NK.S100A2","NK.SELL","DC.pDC","B.plasma")
detailed.celltypes=names(table(DF_scRNA.Xu$detailed))
DF_scRNA.Xu=DF_scRNA.Xu %>%
  mutate(rough=ifelse(detailed %in% c("all"), "all", detailed)) %>%
  mutate(rough=ifelse(detailed %in% detailed.celltypes[grepl("^B\\.",detailed.celltypes)], "B cells", rough)) %>%
  mutate(rough=ifelse(detailed %in% detailed.celltypes[grepl("^CD4T\\.",detailed.celltypes)], "CD4T cells", rough)) %>%
  mutate(rough=ifelse(detailed %in% detailed.celltypes[grepl("^CD8T\\.",detailed.celltypes)], "CD8T cells", rough)) %>%
  mutate(rough=ifelse(detailed %in% detailed.celltypes[grepl("^DC\\.",detailed.celltypes)], "DCs", rough)) %>%
  mutate(rough=ifelse(detailed %in% detailed.celltypes[grepl("^Mono\\.",detailed.celltypes)], "Monocytes", rough)) %>%
  mutate(rough=ifelse(detailed %in% detailed.celltypes[grepl("^NK\\.",detailed.celltypes)], "NK cells", rough)) %>%
  mutate(rough=ifelse(detailed %in% detailed.celltypes[grepl("^Other\\.",detailed.celltypes)], "Other cells", rough))
table(DF_scRNA.Xu$rough, useNA="ifany")

DF_scRNA.Xu$geneset.celltype=paste0(DF_scRNA.Xu$geneset, " | ", DF_scRNA.Xu$rough)
DF_scRNA.Xu$geneset.celltype=gsub(" \\| all","",DF_scRNA.Xu$geneset.celltype)
DF_scRNA.Xu$geneset.celltype=as.factor(DF_scRNA.Xu$geneset.celltype)
levels(DF_scRNA.Xu$geneset.celltype) # check
data.table::fwrite(DF_scRNA.Xu, "~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Comparison_DTW.Intra.gene_updated.csv.gz", sep="\t")

### Plot the results
DF_scRNA.Xu=data.table::fread("~/Project_PBMCage/For_or_From_Xu/Results/ClusterGVis_ck_results_Comparison_DTW.Intra.gene_updated.csv.gz", sep="\t")

DF_scRNA.Xu=DF_scRNA.Xu %>%
  group_by(geneN, geneset.celltype) %>%
  summarize_at("dtw_dist", mean)

plot_dtw=
  ggplot(DF_scRNA.Xu, aes(x=geneN, y=dtw_dist, color=geneset.celltype)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  ggsci::scale_color_d3("category20c") +
  theme_classic() +
  labs(x="geneN", y="DTW distance", title=NULL) +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title="geneset", nrow=3)) +
  geom_smooth(data=subset(DF_scRNA.Xu, geneset.celltype %in% c("AgeGene5K","Zhu et al. (2023)")), 
              aes(x=geneN, y=dtw_dist, group=geneset.celltype), method="loess", show.legend=FALSE, linewidth=0.75, se=FALSE, alpha=0.5) +
  geom_smooth(data=subset(DF_scRNA.Xu, !(geneset.celltype %in% c("AgeGene5K","Zhu et al. (2023)"))), 
              aes(x=geneN, y=dtw_dist, group=geneset.celltype), method="loess", show.legend=FALSE, linewidth=0.25, se=FALSE, alpha=0.5)
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.scRNAAgeGenes_DTW_lines.pdf", height=4, width=8)
plot_dtw
dev.off()

### Plot an example with top100 and 100-200 genes
# set the number of genes to choose for comparison
split_=50 # i.e., plot the mean of top split_ genes against the mean of the second split_ genes in a geneset

# Load the geneset
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
genes_Xu5000=Gene_weights_5000 %>% arrange(desc(mean)) %>% mutate(symbol=gsub("_ENSG","-ENSG",symbol)) %>% 
  select(symbol) %>% unlist() %>% unname()
genes_Xu5000_top100=genes_Xu5000[1:split_]
genes_Xu5000_200=genes_Xu5000[(split_+1):(2*split_)]
scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
genes_scRNA=scRNA_ageGenes %>% 
  mutate(Gene.name=gsub("HLA\\.","HLA-",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.AS1","-AS1",Gene.name)) %>% 
  mutate(Gene.name=gsub("\\.2HG","-2HG",Gene.name)) %>% 
  group_by(Gene.name) %>%
  summarize_at("Vip", mean) %>%
  arrange(desc(Vip)) %>% 
  select(Gene.name) %>% unlist() %>% unname()
genes_scRNA_top100=genes_scRNA[1:split_]
genes_scRNA_200=genes_scRNA[(split_+1):(2*split_)]

# Load the bulk counts data
pseudobulk=readRDS("~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_rough.rds")
pseudobulk=AggregateExpression(pseudobulk, assays="RNA", return.seurat=T, group.by=c("age"))
pseudobulk$age=as.numeric(colnames(pseudobulk))
pseudobulk_data=LayerData(pseudobulk, layer="counts")
pseudobulk_data=pseudobulk_data[rowSums(pseudobulk_data)!=0,]

# scRNA_AgeGene geneset
pseudobulk_data_scRNA.top100=pseudobulk_data[intersect(genes_scRNA_top100, rownames(pseudobulk_data)),]
pseudobulk_data_scRNA.top100=scale(apply(pseudobulk_data_scRNA.top100, 2, mean))
pseudobulk_data_scRNA.200=pseudobulk_data[intersect(genes_scRNA_200, rownames(pseudobulk_data)),]
pseudobulk_data_scRNA.200=scale(apply(pseudobulk_data_scRNA.200, 2, mean))

# Xu_5000 geneset
pseudobulk_data_Xu.top100=pseudobulk_data[intersect(genes_Xu5000_top100, rownames(pseudobulk_data)),]
pseudobulk_data_Xu.top100=scale(apply(pseudobulk_data_Xu.top100, 2, mean))
pseudobulk_data_Xu.200=pseudobulk_data[intersect(genes_Xu5000_top100, rownames(pseudobulk_data)),]
pseudobulk_data_Xu.200=scale(apply(pseudobulk_data_Xu.200, 2, mean))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.scRNAAgeGenes_DTW_example.pdf", height=5, width=8)

plot(
  dtw(pseudobulk_data_scRNA.200, pseudobulk_data_scRNA.top100, keep=TRUE,
      step=rabinerJuangStepPattern(6,"c")),
  type="twoway",
  offset=-1,
  xlab=paste0("Age Index"),
  ylab=paste0("Scaled expression"),
  col=c("black","grey50"))
title(main="Zhu et al. (2023)", adj=0, cex.main=1)
legend("topright", inset=c(0.07,0.02),
       c(paste0("top",split_," genes"),paste0("second",split_," genes")),
       col=c("grey50","black"), lty=c("dashed","solid"), seg.len=1)

plot(
  dtw(pseudobulk_data_Xu.200, pseudobulk_data_Xu.top100, keep=TRUE,
      step=rabinerJuangStepPattern(6,"c")),
  type="twoway",
  offset=-1,
  xlab=paste0("Age Index"),
  ylab=paste0("Scaled expression"),
  col=c("brown3","lightsalmon3"))
title(main="AgeGene5K", adj=0, cex.main=1)
legend("topright", inset=c(0.07,0.02),
       c(paste0("top",split_," genes"),paste0("second",split_," genes")),
       col=c("lightsalmon3","brown3"), lty=c("dashed","solid"), seg.len=1)

dev.off()

### Combine the DTW results and the examples
qpdf::pdf_combine(input=c("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.scRNAAgeGenes_DTW_lines.pdf",
                          "~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.scRNAAgeGenes_DTW_example.pdf"),
                  output="~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.scRNAAgeGenes_DTW.pdf")
unlink(c("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.scRNAAgeGenes_DTW_lines.pdf",
         "~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.scRNAAgeGenes_DTW_example.pdf"))

#####################################





