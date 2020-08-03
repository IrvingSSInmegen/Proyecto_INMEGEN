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

# StageIandII Primary Tumor
COAD_StIandII_biom_PT_trans <- t(COAD_StIandII_biom_PT) %>% row_to_names(row_number = 1)

# StageIandII Blood Derived Normal
COAD_StIandII_biom_BDN_trans <- t(COAD_StIandII_biom_BDN) %>% row_to_names(row_number = 1)

# StageIandII Solid Tissue Normal
COAD_StIandII_biom_STN_trans <- t(COAD_StIandII_biom_STN) %>% row_to_names(row_number = 1)

# StageIIIandIV Primary Tumor
COAD_StIIIandIV_biom_PT_trans <- t(COAD_StIIIandIV_biom_PT) %>% row_to_names(row_number = 1)

# StageIIIandIV Blood Derived Normal
COAD_StIIIandIV_biom_BDN_trans <- t(COAD_StIIIandIV_biom_BDN) %>% row_to_names(row_number = 1)

# StageIIIandIV Solid Tissue Normal
COAD_StIIIandIV_biom_STN_trans <- t(COAD_StIIIandIV_biom_STN) %>% row_to_names(row_number = 1)

##########################
#######            #######
######   GRÁFICAS   ################################
#######            #######
##########################

# Idea sacar la media de las columnas nos dira cual es mas abundante??
# Si es así entonces acomodar en max- min
# presentar un gráfica de pastel para ver cierto porcentaje
# Vero como hacer los geom_violin()
#  Encontrar una manera de solo gráficas los mas importantes
#  o en su defecto gráficar el del interes propio.

library(ggplot2)

#  Como las muestras de tejido sólido tienen mas de 100 muestras
#   son las que se analizarán

# Data
# Matrices de conteo de abundancia 

# COAD_StIandII_biom_PT <- Matriz.primary tumor

#COAD_StIIIandIV_biom_PT <- Matriz.primary tumor 


Media_StIandII_PT <- COAD_StIandII_biom_PT %>% select(-...1) %>%
                      summarise_all(mean) %>%
                      t() %>%
                      as_tibble(rownames = "Bac_Arch")%>%
                      arrange(desc(V1))

Media_StIIIandIV_PT <- COAD_StIIIandIV_biom_PT %>% select(-...1) %>%
                        summarise_all(mean) %>%
                        t() %>%
                        as_tibble(rownames = "Bac_Arch") %>%
                        arrange(desc(V1))

StIandII_PT_vplot <- COAD_StIandII_biom_PT %>% select(...1,k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales.f__Halomonadaceae.g__Cobetia)

StIandII_PT_vplot_pivot %>% StIandII_PT_vplot %>% tidyr::pivot_longer(cols = -"...1",
                                                                     names_to = "Bacteria",
                                                                     values_to = "Valores") %>%
                                                 select(-...1) 



StIandII_PT_vplot %>% tidyr::pivot_longer(cols = -"...1",
                                          names_to = "Bacteria",
                                          values_to = "Valores")%>%
                              ggplot(mapping = aes(x = "Valores",
                                                   y = "Bacteria",
                                                   fill = as.factor(Bacteria))) +
                              geom_boxplot() + 
                              facet_wrap(~as.factor(Bacteria),scales = "free")

