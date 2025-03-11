
#####################################

object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
library(Seurat)
query_counts=LayerData(object, assay="RNA", layer="counts", fast=TRUE)
dimnames(query_counts)=list(rownames(object@assays$RNA@counts), colnames(object@assays$RNA@counts))
query_meta=object[[]]

table(rownames(query_meta)==colnames(query_counts))

query_meta_chosen=query_meta[,c(2,3,4,5,6,8,10,11,12,45)]
head(query_meta_chosen)

counts=query_counts
metadata=query_meta_chosen

library(scrattch.io)
write_dgCMatrix_csv(counts, "~/Project_PBMCage/Matrix_20231231.csv", col1_name="gene", chunk_size=50)
write.csv(metadata, "~/Project_PBMCage/Metadata_20231231.csv", row.names=TRUE)


xxx=read.csv("~/Project_PBMCage/Metadata_20231231.csv")
View(xxx)


### Get metadata
metadata=object[[]][,c(8,ncol(object[[]]))]
metadata=metadata %>%
  mutate(Celltype_major=gsub("\\..*","",CD8T_NK_subsets2)) %>%
  mutate(Celltype_rough=gsub(" .*","",Celltype_major)) %>%
  mutate(Celltype_rough=gsub("^gdT.*","gdT",Celltype_rough)) %>%
  mutate(Celltype_rough=gsub("^B.*","B",Celltype_rough))
colnames(metadata)=c("Age","Celltype_detailed", "Celltype_major", "Celltype_rough")
names(table(metadata$Celltype_detailed))
names(table(metadata$Celltype_major))
names(table(metadata$Celltype_rough))
write.table(metadata, "~/Project_PBMCage/Results/all_cells_annoted_metadata.txt", sep="\t")


### Extract marker genes from all celltypes in the metadata
library(Seurat)

# * the object saved on my hard disk was the original one that we used for analysis
object_original=readRDS("/home/jo_79/Project_PBMCage/For_or_From_Xu/cellxgene.rds")

# # * the object was updated at https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1, so do not use this:
# object_updated=readRDS("/home/jo_79/Project_PBMCage/raw_data/raw.rds")

# Extract the ENS-SYMBOL conversion table
feature_df=object_original@assays$RNA@meta.features
write.table(feature_df, "~/Project_PBMCage/For_or_From_Xu/GeneFeatures_ENStoSYMBOL.txt", sep="\t")

metadata=read.delim("~/Project_PBMCage/Results/all_cells_annoted_metadata.txt", sep="\t")
table(rownames(metadata)==colnames(object_original)) # check
object_original=AddMetaData(object_original, metadata=metadata$Celltype_rough, col.name="Celltype_rough_")
Idents(object_original)="Celltype_rough_"
markers=FindAllMarkers(object_original, only.pos=TRUE, min.pct=0.05, logfc.threshold=0.05)
saveRDS(markers, "~/Project_PBMCage/For_or_From_Xu/Markers.rds")

# ~10000 genes
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=100000, order_by=avg_log2FC)
length(unique(top10$gene))
genes_n10000=unique(top10$gene)

# ~7000 genes
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=900, order_by=avg_log2FC)
length(unique(top10$gene))
genes_n7000=unique(top10$gene)

# 5000+ genes
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=510, order_by=avg_log2FC)
length(unique(top10$gene))
genes_n5000=unique(top10$gene)

# 3000+ genes
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=250, order_by=avg_log2FC)
length(unique(top10$gene))
genes_n3000=unique(top10$gene)

# save
df=data.frame(genes_max=c(genes_n10000),
              genes_large=c(genes_n7000, rep(NA, length(genes_n10000)-length(genes_n7000))),
              genes_medium=c(genes_n5000, rep(NA, length(genes_n10000)-length(genes_n5000))),
              genes_small=c(genes_n3000, rep(NA, length(genes_n10000)-length(genes_n3000)))
              )
write.table(df, "~/Project_PBMCage/For_or_From_Xu/MarkerGenes.txt", sep="\t")

#####################################



### Update the celltype annotation and corresponding genelist
#####################################


object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
query_meta=object[[]]
query_meta_chosen=query_meta[,c("donor_id","age","sex","Annot.detailed","Annot.inter","Annot.rough")]
write.table(query_meta_chosen, "~/Project_PBMCage/For_or_From_Xu/Metadata_20240506.txt", sep="\t")

xxx=read.delim("~/Project_PBMCage/For_or_From_Xu/Metadata_20240506.txt", sep="\t")
View(xxx)

### Extract marker genes from all celltypes in the metadata
library(Seurat)

# * the object saved on my hard disk was the original one that we used for analysis
object_original=readRDS("/home/jo_79/Project_PBMCage/For_or_From_Xu/cellxgene.rds")

# # * the object was updated at https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1, so do not use this:
# object_updated=readRDS("/home/jo_79/Project_PBMCage/raw_data/raw.rds")

# Extract the ENS-SYMBOL conversion table
feature_df=object_original@assays$RNA@meta.features
write.table(feature_df, "~/Project_PBMCage/For_or_From_Xu/GeneFeatures_ENStoSYMBOL.txt", sep="\t")

metadata=read.delim("~/Project_PBMCage/For_or_From_Xu/Metadata_20240506.txt", sep="\t")
object_original$cell_id=colnames(object_original)
object_annoted=subset(object_original, cell_id %in% rownames(metadata))
table(colnames(object_annoted)==rownames(metadata)) # check
object_annoted=AddMetaData(object_annoted, metadata=metadata$Annot.rough, col.name="Annot.rough")
Idents(object_annoted)="Annot.rough"
markers=FindAllMarkers(object_annoted, only.pos=TRUE)
saveRDS(markers, "~/Project_PBMCage/For_or_From_Xu/Markers_20240506.rds")

# ~10000 genes
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=1985, order_by=avg_log2FC)
length(unique(top10$gene))
genes_n10000=unique(top10$gene)

# ~7000 genes
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=1157, order_by=avg_log2FC)
length(unique(top10$gene))
genes_n7000=unique(top10$gene)

# 5000+ genes
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=757, order_by=avg_log2FC)
length(unique(top10$gene))
genes_n5000=unique(top10$gene)

# 3000+ genes
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=429, order_by=avg_log2FC)
length(unique(top10$gene))
genes_n3000=unique(top10$gene)

# save
df=data.frame(genes_10k=c(genes_n10000),
              genes_7k=c(genes_n7000, rep(NA, length(genes_n10000)-length(genes_n7000))),
              genes_5k=c(genes_n5000, rep(NA, length(genes_n10000)-length(genes_n5000))),
              genes_3k=c(genes_n3000, rep(NA, length(genes_n10000)-length(genes_n3000)))
)
write.table(df, "~/Project_PBMCage/For_or_From_Xu/MarkerGenes_20240506.txt", sep="\t")


##





