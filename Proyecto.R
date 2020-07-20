metadata <- vroom::vroom(file = "Metadata-TCGA-All-18116-Samples.csv")

kraken_matched <- vroom::vroom(file = "Kraken-TCGA-Matched2SHOGUN-Raw-Data.csv")

COAD  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma") 
# Pacientes con cancer COAD
# 1017 observaciones 42 variables

Stage_IandII  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                pathologic_stage_label %in% c("Stage I","Stage IA","Stage II",
                                                              "Stage IIA","Stage IIB","Stage IIC")) 
# 576 observaciones Stage I and II

COAD_PT_StIandII <- Stage_IandII %>% filter(sample_type == "Primary Tumor")  # 473 muestras primary tumor

COAD_BDN_StIandII <- Stage_IandII %>% filter(sample_type == "Blood Derived Normal") # 64 muestras blood derived normal

COAD_STN_StIandII <- Stage_IandII %>% filter(sample_type == "Solid Tissue Normal") # 38 muestras solid tissue normal


Stage_IIIandIV <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                        pathologic_stage_label %in% c ("Stage III","Stage IIIA",
                                                                       "Stage IIIB","Stage IIIC",
                                                                       "Stage IV", "Stage IVA",
                                                                       "Stage IVB"))
# 423 observaciones Stage III y IV

COAD_PT_StIIIandIV <- Stage_IIIandIV %>% filter(sample_type == "Primary Tumor") #348 muestras primary tumor

COAD_BDN_StIIIandIV <- Stage_IIIandIV %>% filter(sample_type == "Blood Derived Normal") # 43 muestras Blood derived normal

COAD_STN_StIIIandIV <- Stage_IIIandIV %>% filter(sample_type == "Solid Tissue Normal") #32 muestras solid tissue normal

#15 Phatologic stage label: not available
#3  Phatologic stage label: NA
#1 Recurrent Tumor
#1 Metastatic

COAD_StIandII_biom_PT <- dplyr::semi_join(x = kraken_matched,
                                          y = COAD_PT_StIandII,
                                          by = c("...1" = "...1")) %>%
                                select(contains("Bacteria"))
#190 muestras

COAD_StIandII_biom_BDN <- dplyr::semi_join(x = kraken_matched,
                                           y = COAD_BDN_StIandII,
                                           by = c("...1" = "...1")) %>%
                                select(contains("Bacteria"))
# 29 muestras

COAD_StIandII_biom_STN <- dplyr::semi_join(x = kraken_matched,
                                           y = COAD_STN_StIandII,
                                           by = c("...1" = "...1")) %>%
                                select(contains("Bacteria"))
#15 muestras
