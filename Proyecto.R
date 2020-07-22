library(dplyr)
library(tidyr)
library(vroom)
library(tidyverse)


metadata <- vroom::vroom(file = "Metadata-TCGA-All-18116-Samples.csv")

kraken_matched <- vroom::vroom(file = "Kraken-TCGA-Matched2SHOGUN-Raw-Data.csv")

# implemented a pipeline that converted discrete taxonomic counts into log-counts
# per million (log-cpm) per sample using Voom.
# Data available in https://drive.google.com/drive/folders/18V2ON-Go5AeEtZLe1f9EeJToWOhg81ab

kraken_TCGA_matched2SHOGUN_Voomsnm <- vroom::vroom(file = "Kraken-TCGA-Matched2SHOGUN-Voom-SNM-Quantile-Data.csv")

COAD  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma") 
# Individual with COAD
# 1017 samples 42 variables

Stage_IandII  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                pathologic_stage_label %in% c("Stage I","Stage IA","Stage II",
                                                              "Stage IIA","Stage IIB","Stage IIC")) # 576 Samples Stage I and II 


COAD_PT_StIandII <- Stage_IandII %>% filter(sample_type == "Primary Tumor")  # 473 samples primary tumor

COAD_BDN_StIandII <- Stage_IandII %>% filter(sample_type == "Blood Derived Normal") # 64 samples blood derived normal

COAD_STN_StIandII <- Stage_IandII %>% filter(sample_type == "Solid Tissue Normal") # 38 samples solid tissue normal


Stage_IIIandIV <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                        pathologic_stage_label %in% c ("Stage III","Stage IIIA",
                                                                       "Stage IIIB","Stage IIIC",
                                                                       "Stage IV", "Stage IVA",
                                                                       "Stage IVB"))
# 423 observaciones Stage III y IV

COAD_PT_StIIIandIV <- Stage_IIIandIV %>% filter(sample_type == "Primary Tumor") #348 samples primary tumor

COAD_BDN_StIIIandIV <- Stage_IIIandIV %>% filter(sample_type == "Blood Derived Normal") # 43 samples Blood derived normal

COAD_STN_StIIIandIV <- Stage_IIIandIV %>% filter(sample_type == "Solid Tissue Normal") #32 samples solid tissue normal

#15 Phatologic stage label: not available
#3  Phatologic stage label: NA
#1  Recurrent Tumor
#1  Metastatic

kraken_TCGA_Voom_Novirus <- kraken_TCGA_matched2SHOGUN_Voomsnm %>% select(- contains("Viroid"))
namecolums <- colnames(kraken_TCGA_Voom_Novirus)
grep(pattern = "Viroid",x = namecolums) # 0 viroids


COAD_StIandII_biom_PT <- dplyr::semi_join(x = kraken_TCGA_Voom_Novirus,
                                          y = COAD_PT_StIandII,
                                          by = c("...1" = "...1")) #190 samples

COAD_StIandII_biom_BDN <- dplyr::semi_join(x = kraken_TCGA_Voom_Novirus,
                                           y = COAD_BDN_StIandII,
                                           by = c("...1" = "...1")) # 29 samples

COAD_StIandII_biom_STN <- dplyr::semi_join(x = kraken_TCGA_Voom_Novirus,
                                           y = COAD_STN_StIandII,
                                           by = c("...1" = "...1")) #15 samples

COAD_StIIIandIV_biom_PT <- dplyr::semi_join(x = kraken_TCGA_Voom_Novirus,
                                          y = COAD_PT_StIIIandIV,
                                          by = c("...1" = "...1")) #150 samples

COAD_StIIIandIV_biom_BDN <- dplyr::semi_join(x = kraken_TCGA_Voom_Novirus,
                                           y = COAD_BDN_StIIIandIV,
                                           by = c("...1" = "...1")) # 17 samples

COAD_StIIIandIV_biom_STN <- dplyr::semi_join(x = kraken_TCGA_Voom_Novirus,
                                           y = COAD_STN_StIIIandIV,
                                           by = c("...1" = "...1")) #9 samples



# Matrices de conteo de abundancia

COAD_StIandII_biom_PT_1 <- COAD_StIandII_biom_PT %>% select(-...1)
Samples1 <- COAD_StIandII_biom_PT$...1
AbMat_COAD_StIandII_biom_PT <- as.matrix(t(COAD_StIandII_biom_PT_1))  # StageIandII Primary Tumor
colnames(AbMat_COAD_StIandII_biom_PT) <- Samples1

COAD_StIandII_biom_BDN_1 <- COAD_StIandII_biom_BDN %>% select(-...1)
Samples2 <- COAD_StIandII_biom_BDN$...1
AbMat_COAD_StIandII_biom_BDN <- as.matrix(t(COAD_StIandII_biom_BDN_1)) # StageIandII Blood Derived Normal
colnames(AbMat_COAD_StIandII_biom_BDN) <- Samples2

COAD_StIandII_biom_STN_1 <- COAD_StIandII_biom_STN %>% select(-...1)
Samples3 <- COAD_StIandII_biom_STN$...1
AbMat_COAD_StIandII_biom_STN <- as.matrix(t(COAD_StIandII_biom_STN_1)) # StageIandII Solid Tissue Normal
colnames(AbMat_COAD_StIandII_biom_STN) <- Samples3

COAD_StIIIandIV_biom_PT_1 <- COAD_StIIIandIV_biom_PT %>% select(-...1)
Samples4 <- COAD_StIIIandIV_biom_PT$...1
AbMat_COAD_StIIIandIV_biom_PT <- as.matrix(t(COAD_StIIIandIV_biom_PT_1)) # StageIIIandIV Primary Tumor
colnames(AbMat_COAD_StIIIandIV_biom_PT) <- Samples4

COAD_StIIIandIV_biom_BDN_1 <- COAD_StIIIandIV_biom_BDN %>% select(-...1)
Samples5 <- COAD_StIIIandIV_biom_BDN$...1
AbMat_COAD_StIIIandIV_biom_BDN <- as.matrix(t(COAD_StIIIandIV_biom_BDN_1)) # StageIIIandIV Blood Derived Normal
colnames(COAD_StIIIandIV_biom_BDN_1) <- Samples5

COAD_StIIIandIV_biom_STN_1 <- COAD_StIIIandIV_biom_STN %>% select(-...1)
Samples6 <- COAD_StIIIandIV_biom_STN$...1
AbMat_COAD_StIIIandIV_biom_STN <- as.matrix(t(COAD_StIIIandIV_biom_STN_1)) # StageIIIandIV Solid Tissue Normal
colnames(COAD_StIIIandIV_biom_STN_1) <- Samples6
