library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library(RColorBrewer)


set.seed(1) 

saveext = "/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))


#Now load in other datasets
Key20307 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_SLX-20307/allQC20307.txt",sep="\t",header = T, row.names=1)
raw_counts20307 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_SLX-20307/featurecountsAll_extended_SLX-20307.csv",sep=",",header = T, row.names=1)

Key20308 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_SLX-20308/allQC20308.txt",sep="\t",header = T, row.names=1)
raw_counts20308 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_SLX-20308/featurecountsAll_extended_SLX-20308.csv",sep=",",header = T, row.names=1)

KeyCS5 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CS5/CS5Key_210122.csv",sep=",",header = T, row.names=1)
raw_countsCS5 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CS5/featurecountsCS5.csv",sep=",",header = T, row.names=1)

KeyCS6 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CS6/CS6Key_210402.csv",sep=",",header = T, row.names=1)
raw_countsCS6 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CS6/featurecountsCS6.csv",sep=",",header = T, row.names=1)

KeyCS7 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CS7/CS7Key_210217.csv",sep=",",header = T, row.names=1) 
raw_countsCS7 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CS7/featurecountsCS7.csv",sep=",",header = T, row.names=1)

KeyCRUK <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_SLX-20304/CRUKKey2.csv",sep=",",header = T, row.names=1) 
raw_countsCRUK <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_SLX-20304/featurecounts-CRUK.csv",sep=",",header = T, row.names=1)

KeyCRUK2 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_210301_A00489_0804_AH3FMGDRXY/allQC.csv",sep=",",header = T, row.names=1) 
raw_countsCRUK2 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_210301_A00489_0804_AH3FMGDRXY/featurecounts.csv",sep=",",header = T, row.names=1)

KeyCRUK3 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_210303_A00489_0808_BH3NYWDRXY/allQC.csv",sep=",",header = T, row.names=1) 
raw_countsCRUK3 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CRUK_210303_A00489_0808_BH3NYWDRXY/featurecounts.csv",sep=",",header = T, row.names=1)

KeyPre <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CS1-3/CS1-3Key.csv",sep=",",header = T, row.names=1) 
raw_countsPre <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/CS1-3/featurecountsCS1-3.csv",sep=",",header = T, row.names=1)

marmoset_data_20307 <- CreateSeuratObject(counts = raw_counts20307[,which(Key20307$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_20307) <- Key20307$Primary.Annotation.Old[which(Key20307$QC>0)]
marmoset_data_20307$Stage <- "CS6"
marmoset_data_20307$LOC <- Key20307$Location[which(Key20307$QC>0)]
marmoset_data_20307 <- subset(marmoset_data_20307, idents = c("Tb_CS6","EmDisc_Gast_CS6","ExMes_CS6","Am_CS6","SYS_CS6","PGC_CS6","VE_CS6","EmDisc_CS6","Am_CS6_EmDisc","Stalk_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
#marmoset_data_20307 <- subset(marmoset_data_20307, idents = c("Tb_CS6","EmDisc_Gast_CS6","Am_CS6","EmDisc_CS6","Am_CS6_EmDisc","VE_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_20307 <- NormalizeData(marmoset_data_20307, verbose = FALSE)
marmoset_data_20307$Dataset <- "InVivo"

marmoset_data_20308 <- CreateSeuratObject(counts = raw_counts20308[,which(Key20308$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_20308) <- Key20308$Primary.Annotation[which(Key20308$QC>0)]
marmoset_data_20308$Stage <- "CS7"
marmoset_data_20308$LOC <- Key20308$Location[which(Key20308$QC>0)]
#marmoset_data_20308 <- subset(marmoset_data_20308, idents = c("Tb_CS7","EmDisc_Gast_CS7","EmDisc_CS7","Am_CS7","Am_CS7_EmDisc","VE_CS7")) #,"ReStroma_CS7")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_20308 <- subset(marmoset_data_20308, idents = c("Tb_CS7","EmDisc_Gast_CS7","VE_CS7","EmDisc_CS7","Am_CS7","ExMes_stalk_CS7","PGC_CS7","SYS_CS7","ExMes_CS7","Am_CS7_EmDisc")) #,"ReStroma_CS7")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_20308 <- NormalizeData(marmoset_data_20308, verbose = FALSE)
marmoset_data_20308$Dataset <- "InVivo"

marmoset_data_Pre <- CreateSeuratObject(counts = raw_countsPre[,which(KeyPre$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_Pre) <- KeyPre$Primary.lineage[which(KeyPre$QC>0)]
marmoset_data_Pre$LOC <- KeyPre$Location[which(KeyPre$QC>0)]
marmoset_data_Pre <- subset(marmoset_data_Pre, idents = c("Tb_CS3","Epi_CS3","ICM_CS3","Hyp_CS3","cMor_CS3"))
#marmoset_data_Pre <- subset(marmoset_data_Pre, idents = c("Tb_CS3","Epi_CS3","ICM_CS3","Hyp_CS3","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3"))
marmoset_data_Pre <- NormalizeData(marmoset_data_Pre, verbose = FALSE)
marmoset_data_Pre$Stage <- "Pre"
marmoset_data_Pre$Dataset <- "InVivo"

marmoset_data_CRUK <- CreateSeuratObject(counts = raw_countsCRUK[,which(KeyCRUK$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CRUK) <- KeyCRUK$Primary.Annotation[which(KeyCRUK$QC>0)]
marmoset_data_CRUK$LOC <- KeyCRUK$Loc[which(KeyCRUK$QC>0)]
marmoset_data_CRUK$Stage <- KeyCRUK$Stage[which(KeyCRUK$QC>0)]
#marmoset_data_CRUK <- subset(marmoset_data_CRUK, idents = c("Tb_CS5","EmDisc_CS6","Am_CS7","EmDisc_CS7","Am_CS6","EmDisc_CS5","EmDisc_CS7_Am","Am_CS5","VE_CS6","VE_CS5")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_CRUK <- subset(marmoset_data_CRUK, idents = c("Tb_CS5","EmDisc_CS6","Am_CS7","VE_CS5","EmDisc_CS7","ExMes_CS6","ExMes_CS7","Am_CS6","SYS_CS7","EmDisc_CS5","EmDisc_CS7_Am","VE_CS6","ExMes_CS5","Am_CS5_PGC","Am_CS5","SYS_CS5")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_CRUK <- NormalizeData(marmoset_data_CRUK, verbose = FALSE)
marmoset_data_CRUK$Dataset <- "InVivo"

marmoset_data_CRUK2 <- CreateSeuratObject(counts = raw_countsCRUK2[,which(KeyCRUK2$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CRUK2) <- KeyCRUK2$Primary.Annotation[which(KeyCRUK2$QC>0)]
marmoset_data_CRUK2$LOC <- KeyCRUK2$Loc[which(KeyCRUK2$QC>0)]
marmoset_data_CRUK2$Stage <- KeyCRUK2$Stage[which(KeyCRUK2$QC>0)]
marmoset_data_CRUK2 <- subset(marmoset_data_CRUK2, idents = c("EmDisc_CS6","Am_CS6","EmDisc_CS5","Am_CS6_EmDisc","EmDisc_CS6_Gast","PGC_CS6","EmDisc_CS6_Am","ExMes_CS6","SYS_CS6","VE_CS6","Tb_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
#marmoset_data_CRUK2 <- subset(marmoset_data_CRUK2, idents = c("EmDisc_CS6","Am_CS6","EmDisc_CS5","Am_CS6_EmDisc","EmDisc_CS6_Gast","EmDisc_CS6_Am","Tb_CS6","VE_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_CRUK2 <- NormalizeData(marmoset_data_CRUK2, verbose = FALSE)
marmoset_data_CRUK2Dataset <- "InVivo"

marmoset_data_CRUK3 <- CreateSeuratObject(counts = raw_countsCRUK3[,which(KeyCRUK3$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CRUK3) <- KeyCRUK3$Primary.Annotation[which(KeyCRUK3$QC>0)]
marmoset_data_CRUK3$LOC <- KeyCRUK3$Loc[which(KeyCRUK3$QC>0)]
marmoset_data_CRUK3$Stage <- KeyCRUK3$Stage[which(KeyCRUK3$QC>0)]
#marmoset_data_CRUK3 <- subset(marmoset_data_CRUK3, idents = c("Am_CS6","Am_CS7","EmDisc_CS7","EmDisc_CS6","EmDisc_CS7_Gast","VE_CS7","VE_CS6")) #,"Stroma_CS6")) 
marmoset_data_CRUK3 <- subset(marmoset_data_CRUK3, idents = c("ExMes_CS6","Am_CS6","Am_CS7","EmDisc_CS7","EmDisc_CS6","EmDisc_CS7_Gast","PGC_CS7","VE_CS7","VE_CS6","SYS_CS7","SYS_CS6")) #,"Stroma_CS6")) 
marmoset_data_CRUK3 <- NormalizeData(marmoset_data_CRUK3, verbose = FALSE)
marmoset_data_CRUK3$Dataset <- "InVivo"

marmoset_data_CS5 <- CreateSeuratObject(counts = raw_countsCS5[,which(KeyCS5$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CS5) <- KeyCS5$Primary.Annotation[which(KeyCS5$QC>0)]
marmoset_data_CS5$LOC <- KeyCS5$Location[which(KeyCS5$QC>0)]
marmoset_data_CS5$Stage <- "CS5"
marmoset_data_CS5 <- subset(marmoset_data_CS5, idents = c("Tb_CS5","ExMes_CS5","EmDisc_CS5","SYS_CS5","Am_CS5","VE_CS5","Am_CS5_PGC","EmDisc_CS5_Am")) #","Gland_CS5","ReGland_CS5","ReStroma_CS5")) 
#marmoset_data_CS5 <- subset(marmoset_data_CS5, idents = c("Tb_CS5","EmDisc_CS5","Am_CS5","EmDisc_CS5_Am","VE_CS5")) #","Gland_CS5","ReGland_CS5","ReStroma_CS5")) 
marmoset_data_CS5 <- NormalizeData(marmoset_data_CS5, verbose = FALSE)
marmoset_data_CS5$Dataset <- "InVivo"

marmoset_data_CS6 <- CreateSeuratObject(counts = raw_countsCS6[,which(KeyCS6$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CS6) <- KeyCS6$Primary.Anotation[which(KeyCS6$QC>0)]
marmoset_data_CS6$LOC <- KeyCS6$Location[which(KeyCS6$QC>0)]
marmoset_data_CS6$Stage <- "CS6"
#marmoset_data_CS6 <- subset(marmoset_data_CS6, idents = c("Tb_CS6","EmDisc_CS6","Am_CS6","EmDisc_gast_CS6","VE_CS6")) #,"Gland_CS6","Stroma_CS6","ReStroma_CS6"))
marmoset_data_CS6 <- subset(marmoset_data_CS6, idents = c("Tb_CS6","ExMes_CS6","EmDisc_CS6","SYS_CS6","Am_CS6","VE_CS6","PGC_CS6","EmDisc_gast_CS6","Stalk_CS6")) #,"Gland_CS6","Stroma_CS6","ReStroma_CS6"))
marmoset_data_CS6 <- NormalizeData(marmoset_data_CS6, verbose = FALSE)
marmoset_data_CS6$Dataset <- "InVivo"

marmoset_data_CS7 <- CreateSeuratObject(counts = raw_countsCS7[,which(KeyCS7$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CS7) <- KeyCS7$Primary.Annotation[which(KeyCS7$QC>0)]
marmoset_data_CS7$Stage <- "CS7"
marmoset_data_CS7$LOC <- KeyCS7$Location[which(KeyCS7$QC>0)]
#marmoset_data_CS7 <- subset(marmoset_data_CS7, idents = c("Tb_CS7","EmDisc_CS7","Am_CS7")) #"Gland_CS7","ReStroma_CS7","ReGland_CS7","Myo_CS7")) 
marmoset_data_CS7 <- subset(marmoset_data_CS7, idents = c("Tb_CS7","ExMes_CS7","EmDisc_CS7","SYS_CS7","Am_CS7","Stalk_CS7_PGC","EmDisc_CS7_PGC","ExMes_stalk_CS7","Stalk_CS7")) #"Gland_CS7","ReStroma_CS7","ReGland_CS7","Myo_CS7")) 
marmoset_data_CS7 <- NormalizeData(marmoset_data_CS7, verbose = FALSE)
marmoset_data_CS7$Dataset <- "InVivo"

marmoset_dataInVivo <- merge(marmoset_data_CS5, y = c(marmoset_data_Pre,marmoset_data_CS6,marmoset_data_CS7,marmoset_data_CRUK,marmoset_data_CRUK2,marmoset_data_CRUK3,marmoset_data_20307,marmoset_data_20308), project = "merged")
marmoset_dataInVivo$Dataset <- "2) Marmoset in vivo"

#marmoset_dataInVivo2 <- subset(marmoset_dataInVivo,idents=c("Epi_CS3","Tb_CS3","Am_CS5","Am_CS6","Am_CS7","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS6","Stalk_CS7","Tb_CS5","Tb_CS6","Tb_CS7") )


newID <- as.character(Idents(marmoset_dataInVivo))
newID[which(newID=="Tb_abembryonal_CS7")] <- "Tb_CS7"
newID[which(newID=="Am_CS5_PGC")] <- "PGC_CS5"
newID[which(newID=="EmDisc_CS5_Am")] <- "EmDisc_CS5"
newID[which(newID=="EmDisc_gast_CS6 ")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_Gast_CS6")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_stalk_CS6")] <- "Stalk_CS6"
newID[which(newID=="EmDisc_stalk_CS7")] <- "Stalk_CS7"
newID[which(newID=="EmDisc_CS6_PGC")] <- "PGC_CS6"
newID[which(newID=="Am_CS6_PGC")] <- "PGC_CS6"
newID[which(newID=="Stalk_CS7_PGC")] <- "PGC_CS7"
newID[which(newID=="EmDisc_CS7_PGC")] <- "PGC_CS7"
newID[which(newID=="ExMes_stalk_CS7")] <- "Stalk_CS7"
newID[which(newID=="EmDisc_CS7_Am")] <- "EmDisc_CS7"
newID[which(newID=="EmDisc_CS7_Am")] <- "EmDisc_CS7"
newID[which(newID=="Am_CS6_EmDisc")] <- "Am_CS6"
newID[which(newID=="EmDisc_CS6_Gast")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_gast_CS6")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_CS6_Am")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_Stalk_CS6")] <- "Stalk_CS6"
newID[which(newID=="EmDisc_Gast_CS7")] <- "EmDisc_CS7"
newID[which(newID=="EmDisc_CS7_Gast")] <- "EmDisc_CS7"
newID[which(newID=="EmDisc_CS6_Gast")] <- "EmDisc_CS6"
newID[which(newID=="Am_CS7_EmDisc")] <- "Am_CS7"
newID[which(newID=="cMor_CS3")] <- "cMor_CS2"


marmoset_dataInVivo$Lab <- newID
Idents(marmoset_dataInVivo) <- newID

KeyInVitro1 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/InVitro/InVitroKey_merged.csv",sep=",",header = T, row.names=1) 
raw_countsInVitro1 <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/InVitro/featurecountsInVitro_merged.csv",sep=",",header = T, row.names=1)

KeyClara <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/ClaraMetaUpdated.csv",sep=",",header = T, row.names=1)
raw_countsClara <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Clara_samples_3.csv",sep=",",header = T, row.names=1)

KeyPAVS<- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/InVitro/newTbKey.csv",sep=",",header = T, row.names=1)
raw_countsPAVS <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/Old Data/InVitro/newTb_fixed.csv",sep=",",header = T, row.names=1)

KeyAmnioid4F<- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/QCClaraAll2_v2.csv",sep=",",header = T, row.names=1)
raw_countsAmnioid4F <- read.table("/Users/Clara Munger/OneDrive - University of Cambridge/Documents/PhD/Bioinformatics/Single cell expression/ECM paper/Data/featurecountsClaraAll2_v2.txt",header = T, row.names=1)




marmoset_data_InVitro1 <- CreateSeuratObject(counts = raw_countsInVitro1[,which(KeyInVitro1$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_InVitro1) <- KeyInVitro1$Primary.lineage[which(KeyInVitro1$QC>0)]
marmoset_data_InVitro1$LOC <- "Invitro"
marmoset_data_InVitro1 <- subset(marmoset_data_InVitro1, idents = c("PLAXA","newTSP3","Conventional","Am","BMP_noMEF","Amnoid_bead","EmDisc","EmD","ActA_MEF","CM_TSC","CM_TSF2","OKAEP5esc","BMP_noMEF","ActA_noMEF","BMP_MEF","CHIR_MEF","SB43_MEF"))
marmoset_data_InVitro1 <- NormalizeData(marmoset_data_InVitro1, verbose = FALSE)
marmoset_data_InVitro1$Stage <- "Invitro"
marmoset_data_InVitro1$Dataset <- "InVitro"

marmoset_data_PAVS <- CreateSeuratObject(counts = raw_countsPAVS[,which(KeyPAVS$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_PAVS) <- KeyPAVS$Primary.lineage[which(KeyPAVS$QC>0)]
marmoset_data_PAVS <- NormalizeData(marmoset_data_PAVS, verbose = FALSE)
marmoset_data_PAVS$Stage <- "Invitro"
marmoset_data_PAVS$Dataset <- "InVitro"


marmoset_data_Clara <- CreateSeuratObject(counts = raw_countsClara[,which(KeyClara$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_Clara) <- KeyClara$Type2[which(KeyClara$QC>0)]
marmoset_data_Clara <- NormalizeData(marmoset_data_Clara, verbose = FALSE)
marmoset_data_Clara$Stage <- "Invitro"
marmoset_data_Clara$Dataset <- "InVitro"

marmoset_data_Amnioid4F <- CreateSeuratObject(counts = raw_countsAmnioid4F[,which(KeyAmnioid4F$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_Amnioid4F) <- KeyAmnioid4F$Type2[which(KeyAmnioid4F$QC>0)]
marmoset_data_Amnioid4F <- NormalizeData(marmoset_data_Amnioid4F, verbose = FALSE)
marmoset_data_Amnioid4F$Stage <- "Invitro"
marmoset_data_Amnioid4F$Dataset <- "InVitro"


marmoset_data_InVitroAlt1 <- merge(marmoset_data_InVitro1, y = c(marmoset_data_Clara, marmoset_data_PAVS, marmoset_data_Amnioid4F), project = "merged")
marmoset_data_InVitroAlt1$Dataset <- "InVitro"
marmoset_data_InVitroAlt1$Lab <- Idents(marmoset_data_InVitroAlt1)


all_marm <- merge(marmoset_data_InVitroAlt1, y = marmoset_dataInVivo, project = "merged")
all_marm$newID <- Idents(all_marm)

#Idents(all_marm,cells=WhichCells(all_marm,idents=c("EmDisc","EmD",'Ctlr MEF'))) <- "Epispheroids" #relabel indent
#Idents(all_marm,cells=WhichCells(all_marm,idents=c("Am","Amnoid_bead","Amnioid"))) <- "Amspheroids" #relabel indent
Idents(all_marm,cells=WhichCells(all_marm,idents=c("Tb_CS5","Tb_CS6",'Tb_CS7'))) <- "Tb" #relabel indent
Idents(all_marm,cells=WhichCells(all_marm,idents=c("VE_CS5","VE_CS6",'VE_CS7'))) <- "VE" #relabel indent
Idents(all_marm,cells=WhichCells(all_marm,idents=c("SYS_CS5","SYS_CS6",'SYS_CS7'))) <- "SYS" #relabel indent
Idents(all_marm,cells=WhichCells(all_marm,idents=c("ExMes_CS5","ExMes_CS6",'ExMes_CS7'))) <- "ExMes" #relabel indent
#Idents(all_marm,cells=WhichCells(all_marm,idents=c("PGC_CS5","PGC_CS6",'PGC_CS7'))) <- "PGC" #relabel indent
#Idents(all_marm,cells=WhichCells(all_marm,idents=c("Am_CS5","Am_CS6",'Am_CS7'))) <- "Am" #relabel indent
#Idents(all_marm,cells=WhichCells(all_marm,idents=c("EmDisc_CS5","EmDisc_CS6",'EmDisc_CS7'))) <- "EmDisc" #relabel indent


#mammal.anchors <- FindIntegrationAnchors(object.list = list(Marmoset.InVivo.Test,Marmoset.InVitro.Test), 
#                                         dims = 1:20, anchor.features = 2000, k.filter = 50)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

# to remove the one VTCN1 high ExMes cell
Idents(all_marm,cells="X3536STDY8612558Aligned.sortedByCoord.out.bam") <- "Mixed"
all_marm <- subset(all_marm,idents="Mixed",invert=TRUE)



#-------------------------------------------------------------
# Integration 

Marmoset.InVivo.Test =  subset(all_marm, idents =
                                 c('Epi_CS3','Hyp_CS3','Tb_CS3',
                                   'EmDisc_CS5','EmDisc_CS6', 'EmDisc_CS7',
                                   'Am_CS5','Am_CS6','Am_CS7',
                                   'ExMes',
                                   'Tb',
                                   'SYS',
                                   'VE'))
Marmoset.InVitro.Test =  subset(all_marm, idents =
                                  c('Ctlr MEF','cmTSCPAVS','PLAXA','newTSP3','Conventional', 'Amnioid', '24h post','3D post','6D post','3D post BMP4','Epi BMP4'))

mammal.anchors <- FindIntegrationAnchors(object.list = list(Marmoset.InVivo.Test,Marmoset.InVitro.Test), 
                                         dims = 1:20, anchor.features = 2000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)

avg_expr <- AverageExpression(mammal.combined, assays = "integrated", slot = 'scale.data', return.seurat = FALSE)
avg_expr_matrix <- avg_expr$integrated

#-------------------------------------------------------------
#FIGURE 2: Correlation analysis on subset of interest: Epispheroids and Conventional

InVitro <- avg_expr_matrix[, c('Conventional','Ctlr MEF','Amnioid')]
InVivo <- avg_expr_matrix[, c('Epi_CS3','EmDisc_CS5','EmDisc_CS6', 'EmDisc_CS7',
                              'Am_CS5','Am_CS6','Am_CS7',
                              'VE','SYS','ExMes','Tb')]

correlation_matrix <- cor(InVitro, InVivo ,method = 'pearson' )
pheatmap(correlation_matrix, 
         cluster_rows = FALSE,  # Disable automatic clustering
         cluster_cols = FALSE,  # Disable automatic clustering
         color = colorRampPalette(c("#00B0F0", "#FFFFFF", "#FF0B07"))(50),
         main = "Pearson Correlation - scaled")


#FIGURE 2: Correlation analysis on subset of interest: Epi_BMP4

InVitro <- avg_expr_matrix[, c('Ctlr MEF','Epi BMP4')]
InVivo <- avg_expr_matrix[, c('EmDisc_CS5','EmDisc_CS6', 'EmDisc_CS7',
                              'Am_CS5','Am_CS6','Am_CS7',
                              'VE','SYS','ExMes','Tb')]

correlation_matrix <- cor(InVitro, InVivo ,method = 'pearson' )
pheatmap(correlation_matrix, 
         cluster_rows = FALSE,  # Disable automatic clustering
         cluster_cols = FALSE,  # Disable automatic clustering
         color = colorRampPalette(c("#00B0F0", "#FFFFFF", "#FF0B07"))(50),
         main = "Pearson Correlation - scaled")


#-------------------------------------------------------------
# FIGURE3 Correlation analysis on subset of interest: Posteriorised structures 

InVitro <- avg_expr_matrix[, c('Ctlr MEF','24h post','3D post','6D post')]
InVivo <- avg_expr_matrix[, c('EmDisc_CS5','EmDisc_CS6', 'EmDisc_CS7',
                              'Am_CS5','Am_CS6','Am_CS7',
                              'VE','SYS','ExMes','Tb')]

correlation_matrix <- cor(InVitro, InVivo ,method = 'pearson' )
pheatmap(correlation_matrix, 
         cluster_rows = FALSE,  # Disable automatic clustering
         cluster_cols = FALSE,  # Disable automatic clustering
         color = colorRampPalette(c("#00B0F0", "#FFFFFF", "#FF0B07"))(50),
         main = "Pearson Correlation - scaled")




InVitro <- avg_expr_matrix[, c('3D post','3D post BMP4')]
InVivo <- avg_expr_matrix[, c('EmDisc_CS5','EmDisc_CS6', 'EmDisc_CS7',
                              'Am_CS5','Am_CS6','Am_CS7',
                              'VE','SYS','ExMes','Tb')]

correlation_matrix <- cor(InVitro, InVivo ,method = 'pearson' )
pheatmap(correlation_matrix, 
         cluster_rows = FALSE,  # Disable automatic clustering
         cluster_cols = FALSE,  # Disable automatic clustering
         color = colorRampPalette(c("#00B0F0", "#FFFFFF", "#FF0B07"))(50),
         main = "Pearson Correlation - scaled")



#-------------------------------------------------------------
#Correlation PGCLC without subclustering

# subclustering PGCLC

Marmoset.PGCLC =  subset(all_marm, idents =c('PGC','PGC_CS5','PGC_CS6','PGC_CS7'))
Marmoset.PGCLC$newID <- Idents(Marmoset.PGCLC)
Marmoset.PGCLC <- ScaleData(Marmoset.PGCLC, verbose = FALSE)
Marmoset.PGCLC <- FindVariableFeatures(Marmoset.PGCLC, selection.method = "vst", nfeatures = 20000)
Marmoset.PGCLC <- RunPCA(Marmoset.PGCLC, npcs = 20, verbose = FALSE)
Marmoset.PGCLC <- RunUMAP(Marmoset.PGCLC, reduction = "pca", dims = 1:20)
Marmoset.PGCLC <- FindNeighbors(Marmoset.PGCLC, dims = 1:20)
#Marmoset.PGCLC <- RunTSNE(Marmoset.PGCLC, reduction = "pca", dims = 1:20)
Marmoset.PGCLC <- FindClusters(Marmoset.PGCLC, resolution = 2)

DimPlot(Marmoset.PGCLC, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
FeaturePlot(Marmoset.PGCLC, reduction = "umap", features = "PRDM1", cols =  c("lightgrey", "black"), pt.size = 4)

# based on marker expression:
# cluster 1, 5: PGCLC
# cluster 2: Am
# cluster 0,4,3: ExMes 


PGCL <- WhichCells(Marmoset.PGCLC, expression = newID == "PGC" & seurat_clusters %in% c(1,5))
PGCLC_ExMes <- WhichCells(Marmoset.PGCLC, expression = newID == "PGC" & seurat_clusters %in% c(0,3,4))
PGCLC_Am <- WhichCells(Marmoset.PGCLC, expression = newID == "PGC" & seurat_clusters %in% c(2))
Marmoset.PGCLC$newID <- factor(Marmoset.PGCLC$newID, levels = c(levels(Marmoset.PGCLC$newID), 'PGCL',"PGCLC_ExMes",'PGCLC_Am'))
Marmoset.PGCLC$newID[PGCL] <- "PGCL"
Marmoset.PGCLC$newID[PGCLC_ExMes] <- "PGCLC_ExMes"
Marmoset.PGCLC$newID[PGCLC_Am] <- "PGCLC_Am"
Idents(Marmoset.PGCLC) <- Marmoset.PGCLC$newID


# Integration 

Marmoset.InVivo.Test =  subset(all_marm, idents =
                                 c('Epi_CS3','Hyp_CS3','Tb_CS3',
                                   'EmDisc_CS5','EmDisc_CS6', 'EmDisc_CS7',
                                   'Am_CS5','Am_CS6','Am_CS7',
                                   'ExMes',
                                   'Tb',
                                   'SYS',
                                   'VE','PGC_CS5','PGC_CS6','PGC_CS7'))

Marmoset.InVitro.Test =  subset(all_marm, idents =
                                  c('Ctlr MEF','cmTSCPAVS','PLAXA','newTSP3','Conventional', 'Amnioid', '24h post','3D post','6D post','3D post BMP4','LN MEF','LN no MEF','Epi BMP4'))
Marmoset.PGCLC =  subset(Marmoset.PGCLC, idents =
                                  c('PGCL','PGCLC_ExMes','PGCLC_Am'))
Marmoset.InVitro.Test <- merge(Marmoset.InVitro.Test, y = Marmoset.PGCLC, project = "merged")



mammal.anchors <- FindIntegrationAnchors(object.list = list(Marmoset.InVivo.Test,Marmoset.InVitro.Test), 
                                         dims = 1:20, anchor.features = 2000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)

avg_expr <- AverageExpression(mammal.combined, assays = "integrated", slot = 'scale.data', return.seurat = FALSE)
avg_expr_matrix <- avg_expr$integrated


InVitro <- avg_expr_matrix[, c('PGCL','PGCLC_ExMes','PGCLC_Am')]
InVivo <- avg_expr_matrix[, c('EmDisc_CS5','EmDisc_CS6', 'EmDisc_CS7',
                              'Am_CS5','Am_CS6','Am_CS7',
                              'VE','SYS','ExMes','Tb','PGC_CS5','PGC_CS6','PGC_CS7')]

correlation_matrix <- cor(InVitro, InVivo ,method = 'pearson' )
pheatmap(correlation_matrix, 
         cluster_rows = FALSE,  # Disable automatic clustering
         cluster_cols = FALSE,  # Disable automatic clustering
         color = colorRampPalette(c("#00B0F0", "#FFFFFF", "#FF0B07"))(50),
         main = "Pearson Correlation - scaled")

