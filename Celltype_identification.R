# Tumor samples
immune=c("PTPRC")
epithelial_or_cancer=c("EPCAM")
stromal_fibro_or_endo=c("MME","PECAM1")



# CD4 and CD8 T cells based on functionality
naive=c("LEF1","SELL","TCF7")
effector=c("IFNG")
cytotoxicity=c("GZMB","PRF1")
early_and_general_exhaustion=c("PDCD1","CTLA4","ENTPD1")
antigen_presentation=c("CD74","HLA-DRB1","HLA-DRB5","HLA-DQA2")
ALL.ABOVE=c("LEF1","SELL","TCF7","IFNG","GZMB","PRF1","PDCD1","CTLA4","ENTPD1","CD74","HLA-DRB1","HLA-DRB5","HLA-DQA2")



# T cells based on subsets
TCR_genes=c("TRAV13-2","TRBV7-9")
naive_T=c("CCR7","LTB","SELL")
naive_CD8=c("CCR7","LEF1","CD27")
memory_CD4_ie._CD4_TEM=c("CXCR4","TNFRSF4","CCR6","RORA","IL6ST","IL17RA")
cytotoxicity_T=c("GZMH","GZMA","GZMB")
ALL.ABOVE=c("TRAV13-2","TRBV7-9","CCR7","LTB","SELL","LEF1","CD27","CXCR4","TNFRSF4","CCR6","RORA","IL6ST","IL17RA","GZMH","GZMA","GZMB")



### First-time clustering (applicable even in tumor samples)
# T Cells ("CD3D", "CD3E", "CD8A"), 
# B cells ("CD19", "CD79A", "MS4A1"), 
# Plasma cells ("IGHG1", "MZB1", "SDC1", "CD79A"), 
# Monocytes and macrophages ("CD68", "CD163", "CD14"),
# Monocytes ("S100A9", "S100A8", "MMP19"),
# mast cells ("TPSAB1", "TPSB2"),
# NK Cells ("FGFBP2", "FCG3RA", "CX3CR1", "KLRB1", "NCR1"),  
# type 3 DCs ("LAMP3", "IDO1", "IDO2"),
# type 2 DCs ("CD1E", "CD1C").

# Fibroblasts ("FGF7", "MME", "ACTA2"), 
# Endothelial cells ("PECAM1", "VWF"), 
# epi or tumor ("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24").

# immune ("PTPRC"), 
# stromal:fibro ("MME")
# stromal:endo ("PECAM1", "VWF")
# Photoreceptor cells ("RCVRN"),
# mac ("C1QA", "C1QB"),
# mouse PDAC fibo ("DCN", "LUM", "GSN"),
# Acinar_cells ("Amy1", "Amy2a2")

immune_cells=c("CD3D", "CD3E", "CD8A",
               "CD19", "CD79A", "MS4A1", 
               "IGHG1", "MZB1", "SDC1", "CD79A", 
               "CD68", "CD163", "CD14", 
               "S100A9", "S100A8", "MMP19", 
               "TPSAB1", "TPSB2", 
               "FGFBP2", "FCG3RA", "CX3CR1", "KLRB1", "NCR1", 
               "LAMP3", "IDO1", "IDO2", 
               "CD1E", "CD1C")



### In PBMC
basophil=c("TREML1","MT-ND2","PLD3","ITGA2B","SYTL3","GP1BA","FCER1A","PMM2","CCDC126","SCCPDH","SLC10A3","EDARADD","IRGM","MS4A3","NCBP2L","Z84488.2","OR6K3","HDC","RHEX","CPA3","MS4A4E","IL4","GCSAML")
eosinophil=c("CAT","PPP1R18","CHST13","GAPT","LGALS12","C10orf128","VSTM1","CCL4L2","PIK3R6","GATA1","SLC29A1","CYSLTR2","MYCT1","PYROXD2","ADORA3","CEBPE","IDO1")
neutrophil=c("IFITM2","MNDA","LITAF","VMP1","HIST2H2AA4","FPR1","HIST1H2BJ","ZFP36L1","DUSP1","QPCT","TMEM91","FCGR2A","CSF3R","USP10","RNF149","MTRNR2L11","NRBF2","MX2","LY6G6F","VNN2","ABTB1","KIAA1551","GNG10","LY96","TREM1","MPL","FPR2","R3HDM4","CEACAM4","IFIT2","PIP4P2","HIST1H2BC","SLC19A1","XPO6","SEC14L1","COQ7","HDAC7","GK","CREB5","ZNF467","NAMPT","RALBP1","FRAT2","CXCR2")

classical_mono=c("CYP1B1")
non_classical_mono=c("ICAM4","VMO1")
intermediate_mono=c("C1QC","CCL19","CCL24")

gdT=c("LIM2","KIR3DL1")
MAIT=c("CXCR6","COLQ","ELOVL4","RORC","IL23R","PRSS35")

Naive_B=c("BACH2","CTGF","AC136616.3")
Memory_B=c("RASSF6")

Myeloid_DC=c("CLEC10A","C1orf54","CD1E","CD1B","XCR1")
pDC=c("GZMB","CLIC3","IRF8","SELENOS","ATG101","PTCRA","PLD4","SMIM5","ITM2C","C12orf45","CLN8","CCDC50","IRF7","KPTN","SUSD1","AC097637.1","DNASE1L3","CLEC4C","AMIGO3","TLCD1","MYBL2","SMPD3","NLRP7","ASIP","DIAPH3","LILRA4","TLR9","ALG1L2","DPPA4","KCNK17","TNFRSF6B","ZP1")

Naive_CD4T=c("FHIT","PSMB11")
Naive_CD8T=c("SLC7A3")
Treg=c("PCLAF","PLK1","LAIR2","CARD17","UTS2","CCR4","CCR8","MELK","CTLA4","E2F8","FOXP3","BFSP2")
Memory_CD8T=c("NKX2-3","OR4A47")
# lineage
B_cells=c("IGLL5","CD79B","AC233755.1","CD79A","AC136616.2","MS4A1","FCER2","FCRLA","AC233755.2","AC141272.1","CD72","BTLA","CLECL1","BANK1","BLK","VPREB3","CXCR5","CD19","P2RX5","HLA-DOB","CD22","AIM2","CLEC17A","CR2","TLR10","FCRL1","TNFRSF13B","GRAPL","FCRL5","POU2AF1","E2F5","ELL3","TNFRSF13C","AC136428.1","FCRL2","PLEKHG1","BACH2","PAX5","CTGF","ZNF860","AC136616.3","CLEC20A","RASSF6","NLGN4X")
T_cells=c("CCR7","CCL5","GZMH","CD3D","FGFBP2","CD8B","PIK3IP1","GRAP2","CD8A","SIRPG","CD5","CD3G","KLRG1","TMEM204","CCL4","FLT3LG","S1PR5","LEF1","BATF","ZNF101","OXNAD1","THEMIS","F2R","CD8B2","TIGIT","TRABD2A","UBASH3A","CD6","ZNF683","CD40LG","FHIT","CD28","CXCR6","CAMK4","GPR171","HAPLN3","ICOS","BIRC5","AC008695.1","PCLAF","RRM2","NCR3","PLK1","PATL2","LIM2","LAIR2","NCAPG","COLQ","PDCD1","CCR5","TIMD4","RASGRP1","AURKB","TNFRSF25","AP000311.1","DLGAP5","EDAR","NOS3","BCL11B")
DCs=c("GZMB","IGFLR1","CLIC3","IRF8","SELENOS","ATG101","PTCRA","PLD4","SMIM5","ITM2C","C12orf45","CLEC10A","CLN8","CCDC50","C1orf54","IRF7","KPTN","FLT3","CD1E","SUSD1","CD1B","AC097637.1","DNASE1L3","XCR1","CLEC4C","TLCD1","AMIGO3","MYBL2","SMPD3","NLRP7","ASIP","TIFAB","DIAPH3","LILRA4","TLR9","ZP1","KCNK17","TNFRSF6B","DPPA4","ALG1L2")
Granulocytes=c("PPBP","H3F3B","GP1BB","H3F3A","TUBB1","MTRNR2L12","HIST1H3H","PF4","HIST1H2AC","IFITM2","ICAM3","MNDA","TREML1","MTRNR2L8","LITAF","MTRNR2L1","GP9","VMP1","MT-ND2","HIST2H2AA4","RIPOR2","CAT","FPR1","PLD3","RGS18","SELENOK","NINJ1","H3F3C","HIST1H2BJ","ZFP36L1","DUSP1","NFE2","PF4V1","MSRB1","RGS2","QPCT","ARHGAP9","CD244","ITGA2B","TMEM91","HIST1H4J","FCGR2A","CSF3R","NCF4","TMEM140","USP10","SMIM27","RNF149","MTRNR2L6","HIST1H2BH","STMP1","TMEM71","MTRNR2L11","CIR1","CLINT1","MTRNR2L10","NRBF2","CDA","CASP3","MPPE1","MX2","SORL1","LY6G6F","PAPSS1","RGL4","VNN2","DPEP2","TDP2","PPP1R18","ABTB1","HIST2H2BE","KIAA1551","GNG10","ADAM8","LY96","ST20-MTHFS","GPSM3","SCGB1C2","HIST1H4H","SYTL3","ST8SIA4","HIST1H1C","TREM1","MPL","FPR2","R3HDM4","GP1BA","CHST13","ARL6IP6","CEACAM4","CEP63","GAPT","IFIT2","TOPORS","WIPI1","IFIT3","FCER1A","PMM2","CCDC126","PIP4P2","LGALS12","TM6SF1","HIST1H2BC","ST20","JAK2","GFI1B","SLC19A1","C10orf128","TMEM154","VSTM1","SCCPDH","XPO6","XKR8","SEC14L1","SLC10A3","COQ7","HDAC7","B3GNT8","TCP11L2","CLC","GK","PIGX","CCL4L2","CREB5","HEATR5A","BAZ2B","PIK3R6","LPCAT2","ZNF467","CD200R1","NAMPT","DHRS9","POLM","CEP19","RALBP1","RSBN1L","CASS4","C7orf25","FRAT2","EDARADD","IKBIP","RFLNB","CXCR2","COP1")
Monocytes=c("NAAA","CYP1B1","GPBAR1","APOBEC3A","CSF1R","CD300E","MS4A7","HMOX1","SLC46A2","TNFRSF8","ICAM4","C1QA","LILRB1","C1QB","VMO1","IL10","MS4A14","C1QC","PKD2L1","TMTC1","SPIC","CCL19","CCL24","SLIT2")
NK=c("KLRF1","XCL2","AC068775.1","KLRC3","XCL1","SPTSSB","ERCC6L","DLL1","DNA2","EMID1","TNFSF11")
# Non-blood immune cells
DCs=c("LILRA4","CLEC4C","PTCRA","KRT71","PLD4","IRF7","TTC24","CXCR3","ABHD15","ENTPD7","IRF4","IL3RA","RUBCN")
Eryth=c("HBD","GYPA","GYPB","SLC4A1","HBM","AHSP","KLF1","ALAS2","TSPO2","EPB42","CLC","RHCE","GYPE","OR10Z1","RHAG","RHD","TRIM10","HBA1","SPTA1","IFIT1B","MFSD2B","AC104389.6","MPO","TRIM58","C17orf99")
NK=c("NCR1","FGFBP2","KLRF1","SH2D1B","KIR3DL1","LIM2","S1PR5","NKG7","KIR2DL1","OR2A25","GNLY","CD160","KLRD1","PRF1","FCRL6","KRT73")
Plasma_cells=c("IGLV1-50","IGKV1D-13","IGHV2-70","IGKV1D-33","IGKV2D-30","IGKV3D-7","IGHV4-39","IGKV1-16","IGKV3D-15","IGKV3OR2-268","IGHV3-74","IGHJ5","IGLV1-51","IGKV2-30","IGKV1D-39","IGHV5-10-1","IGLV1-40","IGHV3-64D","IGLV7-46","IGHV2-5","IGKV1-12","IGHV1-3","IGLV2-18","IGHV3-72","IGKV2D-40","AC233755.1","IGHV4-4","IGHV3-66","IGKV2D-28","IGLV3-19","IGHV1-69D","IGKV1-5","IGKV5-2","IGKV1D-12","IGHV4-34","IGHV4-59","IGHV4-28","IGHV4-31","IGLV6-57","IGHV3-20","IGLV2-14","IGKV2-24","IGKV2-28","IGLV3-27","IGHJ3","IGKV1-39","IGHV1-69","IGHV3-30","IGHV6-1","IGHV3-7","AC233755.2","IGHV3-15","IGKV1OR2-108","IGKV3-20","IGKV1-17","IGLV3-21","IGHV3-23","IGLV3-9","IGHV5-51","IGHV2-70D","IGKV2D-24","IGHJ4","IGHV7-4-1","IGKV1-27","IGKV3-11","IGLV2-11","IGKV3D-11","IGLV1-44","IGHV3-33","IGHV3-38","IGKV1-9","IGKV6-21","IGLV3-10","IGKV2-40","IGLV4-60","IGHV3-49","IGLV2-23","IGHV3-21","IGLV2-8","IGHV3-48","IGLV5-45","IGKV3-15","IGHV3-73","IGHV3-53","IGLV3-1","IGHV3-11","IGHJ6","IGKV1-33","IGKV1D-8","IGKV3D-20","IGLV4-69","IGLV3-25","IGKV1D-43","AC136616.2","IGLV7-43","IGHV1-2","IGHJ2","IGKV1D-42","IGHV4-61","IGHV3-43","IGKV4-1","IGLV1-36","IGKV1-6","IGHA2","IGHV3-35","IGHV1-46","IGHV3OR16-17","IGKV1-8","IGHV3-64","IGHV1-58","IGKV2D-29","IGHV1-24","IGHV3-16","IGHV1-18","IGHJ1","IGHV3-13","IGKV1D-16","IGLV4-3","IGLV10-54","AC141272.1","IGHV2-26","IGKV1-37","TNFRSF17","TXNDC5","IGHV1-45","IGHV3OR16-9","IGHV4OR15-8","IGHV3OR16-12","IGHV3OR16-13","JCHAIN","MIXL1","IGLV1-47","IGLV5-48","IGHG2","IGHA1")



### Immune system annotation in Sctype
Pro_B_cells=c("CD27","IGHD","CD24","PTPRC","PAX5","CD24","CD38","CD79A","DNTT","DEPP1","VPREB1","ARPP21","CD99","IGLL1","CD9","CD79B","TCL1A","IGLL5","HLA-DQA1","HLA-DQB1","VPREB3","IGLL5")
Pre_B_cells=c("CD19","CD27","IGHD","CD24","PTPRC","PAX5","CD24","CD38","CD79A","NSMCE1","PCDH9","ACSM3","CCDC191","TCL1A","CD79B","TCL1A","IGLL5","HLA-DQA1","HLA-DQB1","VPREB3")
Naive_B_cells=c("CD19","IGHD","CD38","CD24","MS4A1","PTPRC","PAX5","CD79A","JCHAIN","SSR4","FKBP11","SEC11C","DERL3","PRDX4","IGLL5","CD79B","TCL1A","HLA-DQA1","HLA-DQB1","SDC1","VPREB3")
Memory_B_cells=c("CD19","CD27","IGHD","CD38","CD24","MS4A1","PTPRC","PAX5","CD79A","JCHAIN","SSR4","FKBP11","SEC11C","DERL3","PRDX4","IGLL5","CD79B","TCL1A","HLA-DQA1","HLA-DQB1","SDC1","VPREB3")
Plasma_cells=c("CD27","IGHD","CD38","CD24","MS4A1","PTPRC","PAX5","CD79A","JCHAIN","SSR4","FKBP11","SEC11C","DERL3","PRDX4","IGLL5","CD79B","TCL1A","HLA-DQA1","HLA-DQB1","SDC1","VPREB3")
# Plasma cells do not express "MS4A1"

Naive_CD8_T_cells=c("CD8A","CD2","CD3D","CD3E","CD3G","CD247","SELL","CD27","IL7R","FOXP3","PTPRC","CD8B","CCR6","ITGAM","TNFRSF8","CD6","CTLA4","IL2RA","GZMB","CCR7","GNLY","TRAC","LTB","CD52","TRBC2","SHISA5","LCK","THY1","DAPL1")
# Naive CD8 T cells do not express "IL2RA","CD44","CD69","HLA-DRB1","HLA-DRB5","FAS"
Naive_CD4_T_cells=c("CD4","CD2","CD3D","CD3E","CD3G","C247","SELL","CD27","IL7R","FOXP3","PTPRC","CCR6","ITGAM","TNFRSF8","CD6","CTLA4","IL2RA","GZMB","CCR7","TRAC","LTB","CD52","TRBC2","SHISA5","LCK","THY1","DAPL1")
# Naive CD4 T cells do not express "IL2RA","CD44","CD69","HLA-DRB1","HLA-DRB5","FAS"
Memory_CD8_T_cells=c("CD8A","CD2","CD3D","CD3E","CD3G","CD247","IL2RA","PTPRC","SELL","CD27","IL7R","FOXP3","CCR7","CD8B","CCR6","ITGAM","TNFRSF8","CD6","CTLA4","GZMB","GNLY","S100A4","TRAC","LTB","CD52","TRBC2","SHISA5","LCK","THY1","DAPL1")
Memory_CD4_T_cells=c("CD4","CD2","CD3D","CD3E","CD3G","C247","IL2RA","PTPRC","SELL","CD27","IL7R","FOXP3","CCR7","CCR6","ITGAM","TNFRSF8","CD6","CTLA4","IL2RA","GZMB","GNLY","S100A4","TRAC","LTB","CD52","TRBC2","SHISA5","LCK","THY1","DAPL1")
Effector_CD8_T_cells=c("CD8A","CD2","CD3D","CD3E","CD3G","C247","IL2RA","PTPRC","SELL","CD27","IL7R","FOXP3","CCR7","CD8B","CCR6","ITGAM","TNFRSF8","CD6","CTLA4","IL2RA","GZMB","GNLY","TRAC","LTB","CD52","TRBC2","SHISA5","LCK","THY1","DAPL1")
# Effector CD8 T cells do not express "SELL","CCR7"
Effector_CD4_T_cells=c("CD4","CD2","CD3D","CD3E","CD3G","C247","IL2RA","PTPRC","SELL","CD27","IL7R","FOXP3","CCR7","CCR6","ITGAM","TNFRSF8","CD6","CTLA4","GZMB","TRAC","LTB","CD52","TRBC2","SHISA5","LCK","THY1","DAPL1")
# Effector CD4 T cells do not express "SELL","CCR7"
gdT=c("CD2","CD3D","CD3E","CD3G","CD247","IL2RA","PTPRC","SELL","CD27","IL7R","FOXP3","CCR7","CCR6","ITGAM","TNFRSF8","CD6","CTLA4","GZMB","TRDV2","TRGV9","TRGC1","TRAC","LTB","CD52","TRBC2","SHISA5","LCK","THY1","DAPL1")
# gdT cells do not express "SELL","CCR7"

Platelets=c("ITGA2B","GP1BA","ITGB3","PECAM1","PPBP","PF4","GNG11","SDPR","CLU","MPL")

CD8_NKT_like_cells=c("CD8A","NCAM1","CD2","FCGR3A","KLRD1","CD3D","CD3E","CD3G","CD247","NCR1","ITGAM","KLRB1","KLRK1","CD69","NKG7","IL2RB","KLRK1","GZMB","GZMA","GZMM","GNLY","COX6A2","ZMAT4","KIR2DL4")
CD4_NKT_like_cells=c("CD4","NCAM1","CD2","FCGR3A","KLRD1","CD3D","CD3E","CD3G","CD247","NCR1","ITGAM","KLRB1","KLRK1","CD69","NKG7","IL2RB","KLRK1","GZMB","GZMA","GZMM","GNLY","COX6A2","ZMAT4","KIR2DL4")

NK_cells=c("NCAM1","CD2","FCGR3A","KLRD1","NCR1","ITGAM","KLRB1","KLRK1","CD69","NKG7","IL2RB","KLRK1","GZMB","GZMA","GZMM","FCGR3A","GNLY","COX6A2","ZMAT4","KIR2DL4","NKG7")
# NK cells do not express "CD3D","CD3E","CD3G","CD247","CD4","CD8A","CD8B"

Eosinophils=c("ITGAM","CCR3","IL3RA","IL5RA","FUT4","SIGLEC8","CLC","GATA1","CEBPE","SEMG1","ALOX15","CCL23","PRSS41","PRSS33","THBS4","FOXI1")
Neutrophils=c("CEACAM8","ITGAM","FUT4","FCGR3A","CXCL8","FCGR3B","MNDA","CXCR2","MPO","ELANE","PRTN3","AZU1","LYZ","S100A8","S100A9","PI3","CHI3L1","ANXA3","CXCL1","TGM3","BTNL3","C4BPA","MMP9","CD24","BPI","LTF","GCA")
Basophils=c("CD63","ENPP3","IL3RA","CLC","MS4A3","TCN1","CPA3","HDC","GATA2","MS4A2","IL4","GCSAML","TPSAB1")
Mast_cells=c("KIT","ENPP3","IL2RA","SLC18A2","CD33","FCGR2A","FCER1A","TPSD1","HPGDS")

Classical_mono=c("CD14","ITGAM","CD68","HLA-DRB1","HLA-DRB5","CD33","ITGAX","IL3RA","FUT4","CD3D","CD3E","CD3G","CD247","CEACAM8","VCAN","S100A12","CXCL8","S100A8","S100A9","LYZ","CST3")
# classical monocytes do not express "NCAM1","FCGR3A"
Non_classical_mono=c("CD14","FCGR3A","ITGAM","CD68","HLA-DRB1","HLA-DRB5","CD33","ITGAX","IL3RA","FUT4","CD3D","CD3E","CD3G","CD247","CEACAM8","CDKN1C","LST1","FCER1G","MS4A7","RHOC","S100A8","S100A9","CST3","C1QC")
# nonclassical monocytes do not express "NCAM1"
Intermediate_mono=c("CD14","FCGR3A","ITGAM","CD68","HLA-DRB1","HLA-DRB5","CD33","ITGAX","IL3RA","FUT4","CD3D","CD3E","CD3G","CD247","CEACAM8","IL1B","S100A8","S100A9","CST3","C1QC")
# intermediate monocytes do not express "NCAM1"
Macrophages=c("CD68","CD163","CD14","ITGAM","MRC1","CD80","CD86","FCGR3A","FCGR1A","CCL18","CSF1R","ITGAX","FCGR2A","HLA-DRB1","HLA-DRB5","MSR1","GCA","PF4")
# macrophages do not express "NCAM1"

Megakarocytes=c("ITGB3","ITGA2B","GP1BA","ITGA2B","GP9","CXCR4","MPL")
Erythroid_cells=c("PTPRC","GYPA","RUVBL1","TFRC","FOLR1","CD36","ITGA4","HBB","GYPA","HBD","CA1")

HSC_or_MPP_cells=c("ENG","CD34","CD44","NT5E","PTPRC","ITGB1","NANOG","SOX2","PROM1","ALCAM","MCAM","PECAM1","NES","POU5F1","KIT","KDR","CXCL8","AVP","CRHBP","ALDH1A1","ITGA4","THY1","CD69","CD24","KRT19","ASPM","MME","IL3RA","ABCG2","FLT3","ITGA6","EPCAM","KRT7","ATXN1","CD14","SLAMF1","NGFR","HLA-DRB1","HLA-DRB5")
# HSC/MPP cells do not express "PTPRC","CD38"
Progenitors=c("ENG","CD34","CD44","NT5E","PTPRC","ITGB1","NANOG","SOX2","PROM1","ALCAM","MCAM","PECAM1","NES","POU5F1","KIT","KDR","AVP","CRHBP","ALDH1A1","STMN1","CD38","PTPRC","FLT3")

Myeloid_DCs=c("ITGAX","CD83","CD1C","NRP1","CLEC4C","CD86","IL3RA","CD80","CD1A","CD40","HLA-DQA1","HLA-DRB1","HLA-DRB5","HLA-DPB1","HLA-DPA1","CLEC10A","CST3","GPR31","ODF3L1","PRB2","CD207","ARSL","MRC1","EBLN1","CRIP3")
PDCs=c("ITGAX","CD83","CD1C","NRP1","CLEC4C","CD86","IL3RA","CD80","CD1A","CD40","HLA-DQA1","HLA-DRB1","HLA-DRB5","HLA-DPB1","HLA-DPA1","CLEC10A","CST3","TPM2","LRRC26","ASIP","GPM6B","KRT5","NTM","SCT","SHD","KCNA5","SCARA5","EPHA2","MYMX")

Granulocytes=c("ENPP3","FUT4","ITGAM","CD63","CEACAM8","IL3RA","FCGR3A","CD33","KIT","PTPRC","FCER1A","IL5RA","ANPEP","CD14","IL2RA","CD44","CD69","CD9","HLA-DRB1","HLA-DRB5","CCR3","CSF2RA","ITGAX","CCR3","CD24","FCGR2A","SPN","CXCL8","FCGR3B","MNDA","SIGLEC8","AZU1","MPO","CTSG","LYZ")



### Azimuth References: https://azimuth.hubmapconsortium.org/references
### Human - PBMC
### celltype.l1
Mono=c("CTSS","FCN1","NEAT1","LYZ","PSAP","S100A9","AIF1","MNDA","SERPINA1","TYROBP")
CD4_T=c("IL7R","MAL","LTB","CD4","LDHB","TPT1","TRAC","TMSB10","CD3D","CD3G")
CD8_T=c("CD8B","CD8A","CD3D","TMSB10","HCST","CD3G","LINC02446","CTSW","CD3E","TRAC")
NK_cells=c("NKG7","KLRD1","TYROBP","GNLY","FCER1G","PRF1","CD247","KLRF1","CST7","GZMB")
mature_B_cells=c("CD79A","RALGPS2","CD79B","MS4A1","BANK1","CD74","TNFRSF13C","HLA-DQA1","IGHM","MEF2C")
Other_mature_T_cells=c("CD3D","TRDC","GZMK","KLRB1","NKG7","TRGC2","CST7","LYAR","KLRG1","GZMA")
DCs=c("CD74","HLA-DPA1","HLA-DPB1","HLA-DQA1","CCDC88A","HLA-DRA","HLA-DMA","CST3","HLA-DQB1","HLA-DRB1")
### celltype.l2
Intermediate_B_cells=c("MS4A1","TNFRSF13B","IGHM","IGHD","AIM2","CD79A","LINC01857","RALGPS2","BANK1","CD79B")
Memory_B_cells=c("MS4A1","COCH","AIM2","BANK1","SSPN","CD79A","TEX9","RALGPS2","TNFRSF13C","LINC01781")
Naive_B_cells=c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3")
Plasmablasts=c("IGHA2","MZB1","TNFRSF17","DERL3","TXNDC5","TNFRSF13B","POU2AF1","CPNE5","HRASLS2","NT5DC2")
CD4_cytotoxic_T=c("GZMH","CD4","FGFBP2","ITGB1","GZMA","CST7","GNLY","B2M","IL32","NKG7")
CD4_naive_T=c("TCF7","CD4","CCR7","IL7R","FHIT","LEF1","MAL","NOSIP","LDHB","PIK3IP1")
CD4_proliferating_T=c("MKI67","TOP2A","PCLAF","CENPF","TYMS","NUSAP1","ASPM","PTTG1","TPX2","RRM2")
CD4_central_memory_T=c("IL7R","TMSB10","CD4","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL")
CD4_effector_memory_T=c("IL7R","CCL5","FYB1","GZMK","IL32","GZMA","KLRB1","TRAC","LTB","AQP3")
Treg=c("RTKN2","FOXP3","AC133644.2","CD4","IL2RA","TIGIT","CTLA4","FCRL3","LAIR2","IKZF2")
CD8_naive_T=c("CD8B","S100B","CCR7","RGS10","NOSIP","LINC02446","LEF1","CRTAM","CD8A","OXNAD1")
CD8_proliferating_T=c("MKI67","CD8B","TYMS","TRAC","PCLAF","CD3D","CLSPN","CD3G","TK1","RRM2")
CD8_central_memory_T=c("CD8B","ANXA1","CD8A","KRT1","LINC02446","YBX3","IL7R","TRAC","NELL2","LDHB")
CD8_effector_memory_T=c("CCL5","GZMH","CD8A","TRAC","KLRD1","NKG7","GZMK","CST7","CD8B","TRGC2")
AXL_pos_dendritic_cells=c("PPP1R14A","LILRA4","AXL","IL3RA","SCT","SCN9A","LGMN","DNASE1L3","CLEC4C","GAS6")
cDC1=c("CLEC9A","DNASE1L3","C1orf54","IDO1","CLNK","CADM1","FLT3","ENPP1","XCR1","NDRG2")
cDC2=c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO","PLD4","GSN","SLC38A1","NDRG2","AFF3")
pDC=c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","MZB1","SPIB","IRF4","SMPD3")
CD14_pos_monocyte=c("S100A9","CTSS","S100A8","LYZ","VCAN","S100A12","IL1B","CD14","G0S2","FCN1")
CD16_pos_monocyte=c("CDKN1C","FCGR3A","PTPRC","LST1","IER5","MS4A7","RHOC","IFITM3","AIF1","HES4")
CD56dim_NK=c("GNLY","TYROBP","NKG7","FCER1G","GZMB","TRDC","PRF1","FGFBP2","SPON2","KLRF1")
Proliferating_NK=c("MKI67","KLRF1","TYMS","TRDC","TOP2A","FCER1G","PCLAF","CD247","CLSPN","ASPM")
CD56bright_NK=c("XCL2","FCER1G","SPINK2","TRDC","KLRC1","XCL1","SPTSSB","PPP1R9A","NCAM1","TNFRSF11A")
Erythroid_cells=c("HBD","HBM","AHSP","ALAS2","CA1","SLC4A1","IFIT1B","TRIM58","SELENBP1","TMCC2")
HSPC_Hematopoietic.stem.and.progenitor.cell=c("SPINK2","PRSS57","CYTL1","EGFL7","GATA2","CD34","SMIM24","AVP","MYB","LAPTM4B")
ILC=c("KIT","TRDC","TTLL10","LINC01229","SOX4","KLRB1","TNFRSF18","TNFRSF4","IL1R1","HPGDS")
Platelet=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","HIST1H2AC","RGS18","GP9")
Double_negative_T=c("PTPN3","MIR4422HG","NUCB2","CAV1","DTHD1","GZMA","MYB","FXYD2","GZMK","AC004585.1")
gamma_delta_T=c("TRDC","TRGC1","TRGC2","KLRC1","NKG7","TRDV2","CD7","TRGV9","KLRD1","KLRG1")
MAIT_Mucosal.associated.invariant.T=c("KLRB1","NKG7","GZMK","IL7R","SLC4A10","GZMA","CXCR6","PRSS35","RBM24","NCR3")
### celltype.l3
AXL_pos_myeloid_DC=c("AXL","LILRA4","SCN9A","CLEC4C","LTK","PPP1R14A","LGMN","SCT","IL3RA","GAS6")
AXL_pos_pDC=c("LILRA4","CLEC4C","SCT","EPHB1","AXL","PROC","LRRC26","SCN9A","LTK","DNASE1L3")
Intermediate_B_with_kappa_chain=c("MS4A1","IGKC","IGHM","LINC01857","MARCKS","IGHD","TNFRSF13B","CD24","FCRL2","BANK1")
Intermediate_B_with_lambda_chain=c("MS4A1","IGLC2","IGHM","CD79A","IGLC3","IGHD","BANK1","TNFRSF13C","CD22","TNFRSF13B")
Memory_B_with_kappa_chain=c("BANK1","IGKC","LINC01781","MS4A1","SSPN","CD79A","RALGPS2","TNFRSF13C","LINC00926")
Memory_B_with_lambda_chain=c("BANK1","IGLC2","MS4A1","IGLC3","COCH","TNFRSF13C","IGHA2","BLK","TNFRSF13B","LINC01781")
Naive_B_with_kappa_chain=c("IGHM","TCL1A","IGHD","IGHG3","CD79A","IL4R","CD37","MS4A1","IGKC")
Naive_B_with_lambda_chain=c("IGHM","IGLC2","IGHD","IGLC3","CD79A","CXCR4","MS4A1","IL4R","TCL1A","CD79B")
CD14_pos_mono=c("S100A9","CTSS","LYZ","CTSD","S100A8","VCAN","CD14","FCN1","S100A12","MS4A6A")
CD16_pos_mono=c("LST1","YBX1","AIF1","FCGR3A","NAP1L1","MS4A7","FCER1G","TCF7L2","COTL1","CDKN1C")
CD4_cytotoxic_T=c("GZMH","CD4","GNLY","FGFBP2","IL7R","S100A4","GZMA","CST7","IL32","CCL5")
CD4_naive_T=c("TCF7","CD4","NUCB2","LDHB","TRAT1","SARAF","FHIT","LEF1","CCR7","IL7R")
CD4_proliferating_T=c("MKI67","TYMS","PCLAF","TOP2A","CENPF","NUSAP1","CENPM","BIRC5","ZWINT","TPX2")
CD4_TCM_1=c("LTB","CD4","FYB1","IL7R","LIMS1","MAL","TMSB4X","TSHZ2","AP3M2","TRAC")
CD4_TCM_2=c("CTLA4","MIR4435-2HG","TMSB4X","CD28","CDCA7","TMSB10","MAF","ITM2A","TRAC","CD27")
CD4_TCM_3=c("IL7R","ITGB1","LTB","S100A4","AQP3","TNFRSF4","IL32","TOB1","PDE4D","HOPX")
CD4_TEM_1=c("GZMK","IL7R","IL32","ITGB1","CCL5","GZMA","B2M","DUSP2","KLRB1","SYNE2")
CD4_TEM_2=c("GZMK","CD4","TIGIT","IFNG-AS1","CD40LG","MALAT1","CD3G","CD3D","TRAC","CD3E")
CD4_TEM_3=c("IL7R","CCL5","NOSIP","KLRB1","SERINC5","AQP3","ITGA4","IL32","TRAC","LTB")
CD4_TEM_4=c("CCR9","KDF1","DPP4","CD244","SLC4A4","KLRB1","TMIGD2","CD40LG","IL7R","ODF2L")
CD8_naive_T=c("CD8B","LDHB","LEF1","LINC02446","CD8A","S100B","ID2","TCF7","VIM","CCR7")
CD8_naive_T2=c("CD8B","CCR5","CHI3L2","SOX4","CD8A","TNFRSF9","CD38","SIRPG","LRRN3","LEF1")
CD8_proliferating_T=c("PCLAF","CD8B","TYMS","CD3D","CLSPN","CD3G","MKI67","TRAC","CHEK1","TK1")
CD8_TCM_1=c("CD8B","SELL","CD8A","LYAR","ITGB1","NELL2","DUSP2","IL7R","CCL5","LINC01871")
CD8_TCM_2=c("CD8B","C1orf162","IL7R","GATA3","YBX3","KRT1","CD8A","CTSW","INPP4B","LTB")
CD8_TCM_3=c("CD8B","KLRB1","CD8A","HOPX","IL7R","KLRD1","CCL5","SCML4","LINC02446","TRAC")
CD8_TEM_1=c("GZMK","CD8B","CD8A","CCL5","NKG7","DUSP2","CST7","IL32")
CD8_TEM_2=c("CD8A","CMC1","CD8B","CD160","GZMH","CST7","KLRD1","CCL5","TIGIT","KLRG1")
CD8_TEM_3=c("GZMK","CD8A","ITGB1","CD8B","HOPX","CCL5","KLRD1","NKG7","GNLY","YBX3")
CD8_TEM_4=c("GZMH","THEMIS","GNLY","CD8A","ITGB1","FGFBP2","CD2","GZMB","KLRD1","CD8B")
CD8_TEM_5=c("GZMH","GNLY","ZNF683","TRAC","KLRC2","TYROBP","CD8B","CD8A","TRGC2","GZMB")
CD8_TEM_6=c("TRGC2","CD8A","IFNG-AS1","CD8B","ZNF683","KLRC2","NCR3","IKZF2","DUSP2","RTKN2")
cDC1=c("WDFY4","C1orf54","CLEC9A","BATF3","CLNK","TSPAN33","FLT3","CADM1","IDO1","DNASE1L3")
cDC2_1=c("FCER1A","CD14","CLEC10A","CTSS","ENHO","CD1C","MRC1","FCGR2B","PID1","IL13RA1")
cDC2_2=c("FCER1A","BASP1","CD1C","CD74","CLEC10A","HLA-DPA1","ENHO","HLA-DPB1","PLD4","HLA-DQA1")
dnT_1=c("GZMK","NUCB2","CD8B","GPR183","TCF7","LYAR","MALAT1","C12orf57","LEF1","LDHB")
dnT_2=c("AC004585.1","GPR183","FXYD2","NUCB2","CAV1","CD27","MYB","TMSB4X","GZMK","FGFR1")
Erythoid_cells=c("HBM","ALAS2","HBD","AHSP","SLC4A1","TRIM58","SELENBP1","CA1","IFIT1B","SNCA")
gdT_1=c("TRDC","TRGC1","TRGV9","TRDV2","KLRD1","IL7R","KLRC1","DUSP2","GNLY","KLRG1")
gdT_2=c("KLRC2","CD3G","KIR3DL2","CD3D","TRDC","TRDV1","ZNF683","KLRC1","TRGC1","GZMH")
gdT_3=c("RTKN2","TRDC","TRGC2","LEF1","IKZF2","SOX4","ZNF331","ARID5B","NUCB2","CRTAM")
gdT_4=c("TRDC","TIGIT","KLRC2","TRGC2","IKZF2","GCSAM","FCRL6","TRDV1","CST7","CMC1")
HSPC=c("CDK6","SOX4","PRSS57","AC084033.3","ANKRD28","FAM30A","MYB","EGFL7","SPINK2","SMIM24")
ILC=c("KIT","TRDC","IL1R1","SOX4","TNFRSF18","TYROBP","TNFRSF4","FCER1G","IL2RA","GATA3")
MAIT=c("KLRB1","NKG7","GZMK","SLC4A10","NCR3","CTSW","IL7R","KLRG1","CEBPD","DUSP2")
Proliferating_NK_cells=c("STMN1","KLRF1","TYMS","FCER1G","PCNA","TYROBP","CLSPN","TRDC","PCLAF","SMC2")
CD56dim_NK_1=c("FGFBP2","KLRC2","GNLY","S100A4","CD3E","CST7","LGALS1","PRF1","NKG7","GZMB")
CD56dim_NK_2=c("NKG7","FCER1G","PRF1","KLRB1","SPON2","GZMB","FGFBP2","IGFBP7","CST7","B2M")
CD56dim_NK_3=c("KLRF1","CCL5","TRDC","SYNE2","KLRC1","CMC1","XCL2","KLRB1","KLRD1","IL2RB")
CD56dim_NK_4=c("XCL2","SELL","XCL1","GZMK","KLRC1","SPTSSB","KLRF1","IL2RB","TCF7","TRDC")
CD56bright_NK=c("XCL2","GPR183","SELL","IL2RB","CD44","GZMK","KLRF1","TPT1","KLRC1","XCL1")
pDCs=c("CCDC50","UGCG","TCF4","LILRA4","IRF8","IL3RA","PLD4","IRF7","SERPINF1","ITM2C")
Plasma_cells=c("MZB1","JCHAIN","TNFRSF17","ITM2C","DERL3","TXNDC5","POU2AF1","IGHA1","TXNDC11","CD79A")
Plasmablasts=c("TYMS","TNFRSF17","SHCBP1","TK1","KNL1","ASPM","TXNDC5","TPX2","RRM2","BIRC5")
Platelet=c("GNG11","PPBP","NRGN","PF4","CAVIN2","TUBB1","HIST1H2AC","PRKAR2B","CLU","F13A1")
Memory_regulatory_T=c("RTKN2","B2M","TIGIT","FCRL3","S100A4","AC133644.2","CTLA4","FOXP3","IKZF2","TRAC")
Naive_regulatory_T=c("RTKN2","LEF1","FOXP3","C12orf57","IL2RA","TOMM7","CCR7","TRAC","CD4","LDHB")



### Azimuth References: https://azimuth.hubmapconsortium.org/references
### Human - PBMC
### Plus modification

# CD14 Mono
VlnPlot(Monocytes_subset_try, features=c("S100A9","CTSS","LYZ","CTSD","S100A8","VCAN","CD14","FCN1","S100A12","MS4A6A","FCGR3A","HLA-DRB1","ITGAX","ITGAM"), pt.size=0)
# CD16 Mono
VlnPlot(Monocytes_subset_try, features=c("LST1","YBX1","AIF1","FCGR3A","NAP1L1","MS4A7","FCER1G","TCF7L2","COTL1","CDKN1C"), pt.size=0)
# Intermediate Mono
VlnPlot(Monocytes_subset_try, features=c("FCGR3A","CD14","HLA-DRB1","FCGR3A","ITGAM","ITGAX","IL3RA","FUT4","CEACAM8","CD33","GNLY"), pt.size=0)
# Granulocyte
VlnPlot(Monocytes_subset_try, features=c("CDA","CSF3R","S100A9","S100A12","S100A8","FCER1G","FCGR2A","FCGR3A","FCGR3B","GLUL"), pt.size=0)

# HSPC
VlnPlot(OtherT_subset_try, features=c("CDK6","SOX4","PRSS57","AC084033.3","ANKRD28","FAM30A","MYB","EGFL7","SPINK2","SMIM24"), pt.size=0)
# HSC: CD34+ CD38- CD90+ CD49f+
# MPP: CD34+ CD38- CD90- CD49f-
# MLP: CD34+ CD38- CD90-
# CMP: CD34+ CD38+ CD45RA-; MEP: CD34+ CD38lo; GMP: CD34+ CD38+ CD45RA+; CLP: CD34+ CD38-/lo CD45RA+ CD90-
VlnPlot(OtherT_subset_try, features=c("CD34","CD38","THY1","ITGA6","PTPRC"), pt.size=0)

# B cell
VlnPlot(B_subset_try, features=c("CD79A","RALGPS2","CD79B","MS4A1","BANK1","CD74","TNFRSF13C","HLA-DQA1","IGHM","MEF2C"), pt.size=0)
# B naive
VlnPlot(B_subset_try, features=c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3"), pt.size=0)
# lambda B naive
VlnPlot(B_subset_try, features=c("IGHM","IGHD","MS4A1","TCL1A","IL4R","IGLC6","IGLC7","IGLC3"), pt.size=0)

# B memory
VlnPlot(B_subset_try, features=c("MS4A1","COCH","AIM2","BANK1","SSPN","CD79A","TEX9","RALGPS2","TNFRSF13C","LINC01781"), pt.size=0)
# Bmem1 or 2 (1 express igd but 2 do not) - (MME & CD38 neg)
# Cmem1 or 2 - CD19-, MS4A1+, CD27-, IGHD-, IGHM-
# DN Bmem - CD19+, MS4A1+, CD27-, IGHD-, IGHMlo
VlnPlot(B_subset_try, features=c("CD19","MS4A1","CR2","CD27","IGHM","IGHD","MME","CD38", "POU2AF1","PAX5"), pt.size=0)

# Regulatory B - MS4A1neg
VlnPlot(subset_, features=c("CD1D","CD5","CD19","MS4A1","CR2","CD24","CD38","CD40"), pt.size=0)
# transitional B - CD19+IgD+CD27-CD10+
VlnPlot(B_subset_try, features=c("CD19","IGHD","CD27","IGHM","IGHD","MME"), pt.size=0)
# Transitional (CD19+IgD+CD27-CD10+); Naïve (CD19+IgD+CD27-CD10-); IgM Memory (CD19+IgD+CD27+); Classical Memory (CD19+IgD-CD27+); Double Negative (CD19+IgD-CD27-)
VlnPlot(B_subset_try, features=c("TCL1A","CD19","IGHD","CD27","MME","IGHM"), pt.size=0)
VlnPlot(B_subset_try, features=c("VPREB3","PPP1R14A","PCDH9","PLD4","IGHD","SOX4","IGLL5","TCL1A"), pt.size=0)

# B intermediate
VlnPlot(B_subset_try, features=c("MS4A1","TNFRSF13B","IGHM","IGHD","AIM2","CD79A","LINC01857","RALGPS2","BANK1","CD79B"), pt.size=0)
# non-classical Bmem
VlnPlot(CD14_Mono_try, features=c("FCRL4","CR2","CD27"), pt.size=0)
# Plasma/PBs (the first half neg, the second half pos)
VlnPlot(subset_, features=c("MZB1","CD19","IGHD","IGHM","CD24","CD27","CD38","IRF4","PRDM1","JCHAIN","SDC1","XBP1"), pt.size=0)
# Pre-plasmablast
VlnPlot(pbmc.seu_merged, features=c("PTPRC","HLA-DRA","SDC1","CXCR4","PAX5","BACH2","IRF4","PRDM1","XBP1"), pt.size=0)
# Plasmablast
VlnPlot(B_subset_try, features=c("CXCR3","HLA-DRA","KLF4","MKI67","MS4A1"), pt.size=0)
# Plasma
VlnPlot(B_subset_try, features=c("PTPRC","CXCR4","IGHA1","IGHG1","IGHA2","IGHG2","IGHG3","IGHM","TNFRSF17","SDC1"), pt.size=0)
# lambda vs. kappa
VlnPlot(B_subset_try, features=c("IGLC2","IGLC3","IGHG3","IGKC","IGHG1","IGHG2","IGHE","IGHA1","IGHA2","IGHD","IGHM"), pt.size=0)

# T vs. NK vs. NKT
VlnPlot(Monocytes_subset_try, features=c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","NKG7","KLRD1","TYROBP","GNLY","FCER1G"), pt.size=0)
# NKT
VlnPlot(CD8_subset_try, features=c("CD3D","CD3E","CD3G","NCAM1","CD8A","CD8B","FCGR3A","KLRB1","KLRD1","NKG7"), pt.size=0)
VlnPlot(CD8_subset_try, features=c("GNLY","CD7","CTSW","CCL5","CD2","ZFP36L2","CXCR4","LGALS1","CEBPB","S100A6"), pt.size=0)
# other T
VlnPlot(B_subset_try, features=c("CD3D","TRDC","GZMK","KLRB1","NKG7","TRGC2","CST7","LYAR","KLRG1","GZMA"), pt.size=0)
# dnT
VlnPlot(OtherT_subset_try, features=c("PTPN3", "MIR4422HG", "NUCB2", "CAV1", "DTHD1", "GZMA", "MYB", "FXYD2", "GZMK", "AC004585.1"), pt.size=0)
# dnT_1
VlnPlot(CD8_subset_try, features=c("GZMK","NUCB2","CD8B","GPR183","TCF7","LYAR","MALAT1","C12orf57","LEF1","LDHB"), pt.size=0)
# dnT_2
VlnPlot(CD8_subset_try, features=c("AC004585.1","GPR183","FXYD2","NUCB2","CAV1","CD27","MYB","TMSB4X","GZMK","FGFR1"), pt.size=0)

# gdT
VlnPlot(OtherT_subset_try, features=c("CD3D","CD3E","CD3G","NCAM1","CD8A","CD8B","FCGR3A","KLRB1","KLRD1","NKG7",
                                "TRDC", "TRGC1", "TRGC2", "KLRC1", "NKG7", "TRDV2", "CD7", "TRGV9", "KLRD1", "KLRG1"), pt.size=0)
# gdT_1
VlnPlot(NK_subset_try, features=c("TRDC", "TRGC1", "TRGV9", "TRDV2", "KLRD1", "IL7R", "KLRC1", "DUSP2", "GNLY", "KLRG1"), pt.size=0)
# gdT_2
VlnPlot(OtherT_subset_try, features=c("KLRC2", "CD3G", "KIR3DL2", "CD3D", "TRDC", "TRDV1", "ZNF683", "KLRC1", "TRGC1", "GZMH"), pt.size=0)
# gdT_3
VlnPlot(NK_subset_try, features=c("RTKN2", "TRDC", "TRGC2", "LEF1", "IKZF2", "SOX4", "ZNF331", "ARID5B", "NUCB2", "CRTAM"), pt.size=0)
# gdT_4
VlnPlot(OtherT_subset_try, features=c("TRDC", "TIGIT", "KLRC2", "TRGC2", "IKZF2", "GCSAM", "FCRL6", "TRDV1", "CST7", "CMC1"), pt.size=0)
VlnPlot(OtherT_subset_try, features=c("CMC1", "KLRC1", "LEF1"), pt.size=0)

# MAIT
VlnPlot(Monocytes_subset_try, features=c("KLRB1", "NKG7", "GZMK", "IL7R", "SLC4A10", "GZMA", "CXCR6", "PRSS35", "RBM24", "NCR3"), pt.size=0)
VlnPlot(Monocytes_subset_try, features=c("KLRB1", "NKG7", "GZMK", "SLC4A10", "NCR3", "CTSW", "IL7R", "KLRG1", "CEBPD", "DUSP2"), pt.size=0)

# Cycling T
VlnPlot(OtherT_subset_try, features=c("IL21", "CTLA4", "CORO1B", "MAF", "KIF15", "SMC3", "RRM2", "PDCD1", "SPOCK2", "ICOS"), pt.size=0)
VlnPlot(OtherT_subset_try, features=c("HES4","TNFRSF18","TNFRSF4","TNFRSF25","TNFRSF1B","UBXN10-AS1","STMN1","LCK","CLSPN","HEYL"), pt.size=0)

# CD4 T
VlnPlot(CD4_subset_try, features=c("IL7R","MAL","LTB","CD4", "LDHB", "TPT1", "TRAC", "TMSB10", "CD3D", "CD3G"), pt.size=0)
# CD8 T
VlnPlot(Othercells_subset_try, features=c("CD8B", "CD8A", "CD3D", "TMSB10", "HCST", "CD3G", "LINC02446", "CTSW", "CD3E", "TRAC"), pt.size=0)

# CD8 Naive
VlnPlot(CD8T_subset_try, features=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP","LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"), pt.size=0)
# CD8 Naive_1
VlnPlot(CD8T_subset_try, features=c("CD8B","LDHB","LEF1","LINC02446","CD8A","S100B","ID2","TCF7","VIM","CCR7"), pt.size=0)
# CD8 Naive_2
VlnPlot(CD8T_subset_try, features=c("CD8B","CCR5","CHI3L2","SOX4","CD8A","TNFRSF9","CD38","SIRPG","LRRN3","LEF1"), pt.size=0)
# cluster1: homeostatic proliferation (spontaneous proliferation and in partial activation)
VlnPlot(Monocytes_subset_try, features=c("LEF1","SELL","C1orf56","HNRNPH1"), pt.size=0)
# cluster3: transcriptional repressors, essential for maintaining cellular homeostasis and regulating
VlnPlot(Monocytes_subset_try, features=c("TSHZ2","FHIT","RNASET2"), pt.size=0)


# CD4 Naive
VlnPlot(CD4T_subset_try, features=c("TCF7", "CD4", "CCR7", "IL7R", "FHIT","LEF1", "MAL", "NOSIP", "LDHB", "PIK3IP1"), pt.size=0)
# cluster0: CD4 Tnaive with proliferative property (under transformation to Tcm as CD44+ SELL+)
VlnPlot(Othercells_subset_try, features=c("TCF7","SELL","LEF1","CCR7","IL7R","CD44","C1orf56","HNRNPH1","CDC42SE1"), pt.size=0)
# cluster1: CD4 Tnaive under cell regulation of events including activation and proliferation
VlnPlot(Othercells_subset_try, features=c("CD44","SELL","IL7R","COTL1","SOCS3","ANXA1","TAGLN2"), pt.size=0)
# cluster3: CD4 Tnaive under cell regulation of quiescence (repressing differentiation), reported to be downreg during aging
VlnPlot(Othercells_subset_try, features=c("TCF7","SELL","LEF1","SOX4","ID2","SATB1"), pt.size=0)

# CD4 CTL
VlnPlot(Monocytes_subset_try, features=c("GZMH", "CD4", "FGFBP2", "ITGB1", "GZMA", "CST7", "GNLY", "B2M", "IL32", "NKG7"), pt.size=0)
# CD4 sig+
VlnPlot(CD8_subset_try, features=c("ISG15", "IFI44L", "IFI44", "RSAD2", "EIF2AK2", "STAT1", "PARP14", "SAMD9L", "HERC5", "LY6E"), pt.size=0)
# Tcm or Tem or Teff (Tcm: CD44+ SELL+; Tem: CD44+ SELL- IL7R+; Teff: CD44+ SELL- IL7R-)
VlnPlot(CD4T_subset_try, features=c("CD44","SELL","IL7R"), pt.size=0)
# naive (PECAM1+CD103+), Tcm or stem cell memory T (CCR7, CD127, CD62L), Tcm (IL2RA), Tem or Teff (HLA-DR, CCR5, TBX21, GZMA), Trm (CD69,ITGAE,CTLA4)
VlnPlot(CD4T_obj_subset, features=c("PECAM1","ITGAE","CCR7","IL7R","SELL","IL2RA","HLA-DRA","HLA-DRB1","CCR5","TBX21","GZMA",
                                   "CD69","CTLA4"), pt.size=0)
# CD4 TCM
VlnPlot(CD4T_obj_subset, features=c("IL7R", "TMSB10", "CD4", "ITGB1", "LTB", "TRAC", "AQP3", "LDHB", "IL32", "MAL"), pt.size=0)
# CD4 TCM_1
VlnPlot(CD4_subset_try, features=c("LTB", "CD4", "FYB1", "IL7R", "LIMS1", "MAL", "TMSB4X", "TSHZ2", "AP3M2", "TRAC"), pt.size=0)
# CD4 TCM_2
VlnPlot(CD4_subset_try, features=c("CTLA4", "MIR4435-2HG", "TMSB4X", "CD28", "CDCA7", "TMSB10", "MAF", "ITM2A", "TRAC", "CD27"), pt.size=0)
# CD4 TCM_3
VlnPlot(CD4_subset_try, features=c("IL7R", "TMSB10", "CD4", "ITGB1", "LTB", "TRAC", "AQP3", "LDHB", "IL32", "MAL"), pt.size=0)
# cluster0: CD4 Tcm
VlnPlot(CD4_subset_try, features=c("CCR7","SELL","GPR183","PIK3IP1"), pt.size=0)
# cluster1: CD4 Tcm under cell regulation (events eg. proliferation/cytokine production/apoptosis...)
VlnPlot(CD4_subset_try, features=c("CD44","SELL","IL7R","COTL1","ITGB1","LGALS1","S100A11","ANXA1"), pt.size=0)
# cluster1: CD4 Tcm with cell proliferative property
VlnPlot(CD4_subset_try, features=c("CD44","SELL","IL7R","C1orf56","HNRNPH1","CDC42SE1"), pt.size=0)
# cluster2: CD4 Tcm TH17-like
VlnPlot(Monocytes_subset_try, features=c("CD44","SELL","IL7R","KLRB1","CTSH","AQP3","TIMP1"), pt.size=0)
# cluster1: IFN-I T at the status of Tcm
VlnPlot(Monocytes_subset_try, features=c("LIMS1","MAL","LTB","AQP3","MX1","ISG15","ISG20"), pt.size=0)
# cluster2: IFN-I T at the status of transcriptional repression
VlnPlot(CD4_subset_try, features=c("MX1","ISG15","ISG20","TSHZ2","FHIT","RNASET2"), pt.size=0)
# cluster3: CD4 Tcm with upregulated HLA-DR
VlnPlot(Monocytes_subset_try, features=c("CD44","SELL","IL7R","HLA-DRA","CD38","IL2RA"), pt.size=0)

# CD4 TEM
VlnPlot(CD4T_obj_subset, features=c("IL7R", "CCL5", "FYB1", "GZMK", "IL32", "GZMA", "KLRB1", "TRAC", "LTB", "AQP3"), pt.size=0)
# CD4 TEM_1
VlnPlot(CD4_subset_try, features=c("GZMK", "IL7R", "IL32", "ITGB1", "CCL5", "GZMA", "B2M", "DUSP2", "KLRB1", "SYNE2"), pt.size=0)
# CD4 TEM_2
VlnPlot(CD4_subset_try, features=c("GZMK","CD4","TIGIT","IFNG-AS1","CD40LG","MALAT1","CD3G","CD3D","TRAC","CD3E"), pt.size=0)
# CD4 TEM_3
VlnPlot(CD4_subset_try, features=c("IL7R", "CCL5", "NOSIP", "KLRB1", "SERINC5", "AQP3", "ITGA4", "IL32", "TRAC", "LTB"), pt.size=0)
# CD4 TEM_4
VlnPlot(CD4_subset_try, features=c("CCR9", "KDF1", "DPP4", "CD244", "SLC4A4", "KLRB1", "TMIGD2", "CD40LG", "IL7R", "ODF2L"), pt.size=0)
# cluster0: CD4 Tem with Th17 transcriptional signature (KLRB1+ GZMK+)
VlnPlot(CD4_subset_try, features=c("CD44","SELL","IL7R","GZMK","KLRB1"), pt.size=0)
# cluster2: CD4 Tem with transcriptional regulation property
VlnPlot(CD4_subset_try, features=c("JUND","MALAT1","DDX17"), pt.size=0)

# Treg
VlnPlot(CD4_subset_try, features=c("RTKN2","FOXP3","AC133644.2","CD4","IL2RA","TIGIT","CTLA4","FCRL3","LAIR2","IKZF2"), pt.size=0)

# CD8 TCM
VlnPlot(CD8T_subset_try, features=c("CD8B", "ANXA1", "CD8A", "KRT1", "LINC02446", "YBX3", "IL7R", "TRAC", "NELL2", "LDHB"), pt.size=0)
VlnPlot(CD8T_subset_try, features=c("CCR4","CCR7", "CD27", "CD3D", "CD3G", "CD3E", "CD8A", "CD8B", "SELL", "CD69", "KLRC1","KLRC2","KLRC3","IKZF2"), pt.size=0)
# CD8 TCM_1
VlnPlot(pbmc.seu_merged, features=c("CD8B","SELL","CD8A","LYAR","ITGB1","NELL2","DUSP2","CCL5","IL7R","LINC01871"), pt.size=0)
# CD8 TCM_2
VlnPlot(pbmc.seu_merged, features=c("CD8B", "C1orf162", "IL7R", "GATA3", "YBX3", "KRT1", "CD8A", "CTSW", "INPP4B", "LTB"), pt.size=0)
# CD8 TCM_3
VlnPlot(pbmc.seu_merged, features=c("CD8B", "KLRB1", "CD8A", "HOPX", "KLRD1", "IL7R", "CCL5", "SCML4", "LINC02446", "TRAC"), pt.size=0)
# CD8 TEM
VlnPlot(gdTsubset, features=c("CCL5", "GZMH", "CD8A", "TRAC", "KLRD1", "TRGC2", "NKG7", "GZMK", "CST7", "CD8B"), pt.size=0)
VlnPlot(CD8_subset_try, features=c("GZMB","FCGR3A","HNRNPH1","GZMK","NKG7","XCL1","GATA3","LYAR","MALAT1","LINC02446","TSHZ2"), pt.size=0)
# cluster0: mature cytotoxic (CD8T.Tem.GZMB_FCGR3A)
VlnPlot(CD8_subset_try, features=c("FCGR3A","PRF1","GZMB"), pt.size=0)
# cluster1: less mature, but more proliferative (CD8T.Tem.GZMB_HNRNPH1)
VlnPlot(CD8_subset_try, features=c("C1orf56","HNRNPH1","CDC42SE1"), pt.size=0)
# CD8 TEM_1
VlnPlot(CD8_subset_try, features=c("GZMK","CD8B","CD8A","CCL5","NKG7","DUSP2","CST7","IL32"), pt.size=0)
# CD8 TEM_2
VlnPlot(CD8_subset_try, features=c("CD8A","CMC1","CD8B","CD160","GZMH","CST7","KLRD1","CCL5","TIGIT","KLRG1"), pt.size=0)
# CD8 TEM_3
VlnPlot(CD8_subset_try, features=c("GZMK","CD8A","ITGB1","CD8B","HOPX","CCL5","KLRD1","NKG7","GNLY","YBX3"), pt.size=0)
# CD8 TEM_4 (CD8T.Tem.GZMB_ITGB1)
VlnPlot(CD8_subset_try, features=c("GZMH","THEMIS","GNLY","CD8A","ITGB1","FGFBP2","CD2","GZMB","KLRD1","CD8B"), pt.size=0)
# CD8 TEM_5
VlnPlot(B_subset_try, features=c("GZMH","GNLY","ZNF683","TRAC","KLRC2","TYROBP","CD8B","CD8A","TRGC2","GZMB","GZMK"), pt.size=0)
# CD8 TEM_6
VlnPlot(pbmc.seu_merged, features=c("GZMK","CD8A","ITGB1","CD8B","HOPX","CCL5","KLRD1","NKG7","GNLY","YBX3"), pt.size=0)
# CD8+ T subsets
VlnPlot(NK_subset_try, features=c("SELL","CCR7","TCF7","LEF1","GZMB","GZMK","CD69","IL7R","ITGAE","CCR9","KLRB1","CD38","HLA-DRB1","HLA-DRA","HLA-DRB5","KLRC2","NCAM1","TYROBP","ZNF683","IKZF2","CCR4","IL2RA","IL4","IL2"), pt.size=0)
# CD8 NKT vs. NKT
VlnPlot(CD8_subset_try, features=c("CD3E","KLRF1","CD8A","NCAM1","FCGR3A","CD3D","KLRB1","STMN1","TRAC","TYROBP","CD2","KLRD1","NCR3","NKG7","ITGA2"), pt.size=0)


# DC
VlnPlot(Monocytes_subset_try, features=c("CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CCDC88A", "HLA-DRA", "HLA-DMA", "CST3", "HLA-DQB1", "HLA-DRB1"), pt.size=0)
# pre-pDC
VlnPlot(DCs_subset_try, features=c("SCT", "SHD", "LILRA4", "LILRB4", "PTPRS", "TNNI2", "PLD4", "SPIB", "IRF8", "TNFRSF21"), pt.size=0)
# pre-cDC
VlnPlot(DCs_subset_try, features=c("ENHO", "CLEC10A", "RNASE2", "PLBD1", "FCER1A", "IGSF6", "MNDA", "SAMHD1", "ALDH2", "PAK1"), pt.size=0)
# ASDC_mDC
VlnPlot(DCs_subset_try, features=c("PPP1R14A","LILRA4","AXL","IL3RA","SCT","SCN9A","LGMN","CLEC4C","DNASE1L3","GAS6"), pt.size=0)
# ASDC
VlnPlot(pbmc.seu_merged, features=c("AXL", "LILRA4", "SCN9A", "CLEC4C", "LTK", "PPP1R14A", "LGMN", "SCT", "IL3RA", "GAS6"), pt.size=0)
# ASDC_pDC
VlnPlot(DCs_subset_try, features=c("LILRA4", "CLEC4C", "SCT", "EPHB1", "AXL", "PROC", "LRRC26", "SCN9A", "LTK", "DNASE1L3"), pt.size=0)
# cDC2
VlnPlot(Othercells_subset_try, features=c("FCER1A", "HLA-DQA1", "CLEC10A", "CD1C", "ENHO", "PLD4", "GSN", "SLC38A1", "NDRG2", "AFF3"), pt.size=0)
VlnPlot(DCs_subset_try, features=c("CD14", "CD163", "CLEC10A", "NOTCH2", "ITGAM", "SIRPA", "CX3CR1", "CD1C", "CD2", "ID2","IRF4","KLF4","ZBTB46"), pt.size=0)
# cDC2_1	
VlnPlot(DCs_subset_try, features=c("FCER1A", "CD14", "CLEC10A", "CTSS", "ENHO", "CD1C", "MRC1", "FCGR2B", "PID1", "IL13RA1"), pt.size=0)
# cDC2_2	
VlnPlot(DCs_subset_try, features=c("FCER1A", "BASP1", "CD1C", "CD74", "CLEC10A", "HLA-DPA1", "ENHO", "HLA-DPB1", "PLD4", "HLA-DQA1"), pt.size=0)
# cDC1
VlnPlot(Othercells_subset_try, features=c("BLTA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "BATF3","ID2","IRF8","ZBTB46"), pt.size=0)
# pDC
VlnPlot(Othercells_subset_try, features=c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","MZB1","SPIB","IRF4","SMPD3"), pt.size=0)
# AXL DC
VlnPlot(DCs_subset_try, features=c("PPP1R14A","LILRA4","AXL","IL3RA","SCT","SCN9A","LGMN","DNASE1L3","CLEC4C","GAS6"), pt.size=0)

# ILC
VlnPlot(Othercells_subset_try, features=c("KIT","TRDC","TTLL10","LINC01229","SOX4","KLRB1","TNFRSF18","TNFRSF4","IL1R1","HPGDS"), pt.size=0)
VlnPlot(pbmc.seu_merged, features=c("IL7R","KLRB1","IL1R1"), pt.size=0)
# ILC1 (TBX21/T-bet+), ILC2 (GATA3+), ILC3 (RORC+ AHR+), LTi (RORC+)
VlnPlot(Ery_dnT_ILC_try, features=c("TBX21","GATA3","RORC","AHR"), pt.size=0)
# ILC1
VlnPlot(CD4_Tem_try, features=c("TTLL10","FCER1G","SH2D1B","XCL2","XCL1","LINC00299","GNLY","TRGC1","KIT","NMUR1"), pt.size=0)
# ILC3
VlnPlot(CD4_Tem_try, features=c("TNFRSF18","C1orf141","IL23R","LINGO4","RORC","FCER1G","SH2D1B","XCL2","XCL1","LINC00298"), pt.size=0)
# NK
VlnPlot(Monocytes_subset_try, features=c("GNLY","TYROBP","NKG7","FCER1G","GZMB","TRDC","PRF1","FGFBP2","SPON2","KLRF1"), pt.size=0)
VlnPlot(B_subset_try, features=c("GNLY","TYROBP","NKG7","FCER1G","GZMB","TRDC","PRF1","FGFBP2","SPON2","KLRF1","FCER1G","KLRC2","FOS","HNRNPH1"), pt.size=0)
# CD56hi NK
VlnPlot(B_subset_try, features=c("XCL2","FCER1G","SPINK2","TRDC","KLRC1","XCL1","SPTSSB","PPP1R9A","NCAM1","TNFRSF11A","FCGR3A"), pt.size=0)
VlnPlot(pbmc.seu_merged, features=c("GNLY","KLRF1","KLRD1","NKG7","NCAM1","TYROBP","CTSW","XCL1","EOMES","AREG"), pt.size=0)
VlnPlot(NK_subset_try, features=c("GNLY","KLRF1","KLRD1","NKG7","NCAM1","TYROBP","CCL3","XCL1","COTL1","FCGR3A","S100B"), pt.size=0)
# immature NK
VlnPlot(subset_, features=c("KLRK1","NCR1","NCR2","NCR3","KLRB1"), pt.size=0)
# mature NK
VlnPlot(subset_, features=c("KLRD1","ITGB2","KIR2DL4","PRF1","IFNG","NCAM1","FCGR3A"), pt.size=0)
# proliferative NK
VlnPlot(NK_subset_try, features=c("MKI67","KLRF1","TYMS","TRDC","TOP2A","FCER1G","PCLAF","CD247","CLSPN","ASPM","MCM2"), pt.size=0)
# NK_1
VlnPlot(subset_, features=c("FGFBP2","KLRC2","GNLY","S100A4","CD3E","CST7","LGALS1","PRF1","NKG7","GZMB"), pt.size=0)
# NK_2
VlnPlot(subset_, features=c("NKG7","FCER1G","PRF1","KLRB1","SPON2","GZMB","FGFBP2","IGFBP7","CST7","B2M"), pt.size=0)
# NK_3
VlnPlot(subset_, features=c("KLRF1","CCL5","TRDC","SYNE2","KLRC1","CMC1","XCL2","KLRB1","KLRD1","IL2RB"), pt.size=0)
# NK_4
VlnPlot(B_subset_try, features=c("XCL2","SELL","XCL1","GZMK","KLRC1","SPTSSB","KLRF1","IL2RB","TCF7","TRDC"), pt.size=0)

# platelet
VlnPlot(Othercells_subset_try, features=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","HIST1H2AC","RGS18","GP9"), pt.size=0)
# Megakaryocytes
VlnPlot(Othercells_subset_try, features=c("LTBP1","ITGB3","ITGA2B","PF4","THBS1","PLEK","MED12L","PRKAR2B","RUNX1","ARHGAP6"), pt.size=0)
# Erythrocyte
VlnPlot(Othercells_subset_try, features=c("HBB","HBA1","HBA2","HBD","AHSP","ALAS2","CA1","HBM","SLC25A37","SNCA"), pt.size=0)



