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



Difference <- dplyr::left_join(x = Media_StIandII_PT,
                                y = Media_StIIIandIV_PT,
                                suffix = c ("_StIandII","_StIIIandIV"),
                                by = c("Bac_Arch" = "Bac_Arch"))
diferencia <- Difference %>% mutate(diferencias = abs(V1_StIandII - V1_StIIIandIV),
                                    Estado = case_when(diferencias == 0 ~ "Nula",
                                                       diferencias < .5 ~ "Baja",
                                                       diferencias > .5 ~ "Alta"))
Diferencia_alta <- diferencia %>% filter(Estado == "Alta")

nombres <- Diferencia_alta %>% select(Bac_Arch)
bacarchea <- as.vector(nombres)                    



# Violin plots  Stage I and II

StIandII_PT_vplot <- COAD_StIandII_biom_PT %>% select(...1,
                                                      k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Gallibacterium,
                                                      k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Aggregatibacter,
                                                      k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Leptotrichiaceae.g__Leptotrichia,
                                                      k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Buchnera,
                                                      k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Limnohabitans,
                                                      k__Archaea.p__Euryarchaeota.c__Methanococci.o__Methanococcales.f__Methanococcaceae.g__Methanococcus,
                                                      k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium,
                                                      k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Bradyrhizobiaceae.g__Rhodopseudomonas,
                                                      k__Bacteria.p__Cyanobacteria.o__Nostocales.f__Nostocaceae.g__Anabaena,
                                                      k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Listeriaceae.g__Brochothrix,
                                                      k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Pelomonas)
                                                    

StIandII_PT_vplot_pivot<-  StIandII_PT_vplot %>% tidyr::pivot_longer(cols = -"...1",
                                                                     names_to = "Bacteria",
                                                                     values_to = "Valores") %>%
                                                 select(-...1) 



StIandII_PT_vplot_pivot$Bacteria <- factor(StIandII_PT_vplot_pivot$Bacteria,
                                           labels = c("k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Gallibacterium",
                                                      "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Aggregatibacter",
                                                      "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Leptotrichiaceae.g__Leptotrichia",
                                                      "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Buchnera",
                                                      "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Limnohabitans",
                                                      "k__Archaea.p__Euryarchaeota.c__Methanococci.o__Methanococcales.f__Methanococcaceae.g__Methanococcus",
                                                      "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium",
                                                      "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Bradyrhizobiaceae.g__Rhodopseudomonas",
                                                      "k__Bacteria.p__Cyanobacteria.o__Nostocales.f__Nostocaceae.g__Anabaena",
                                                      "k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Listeriaceae.g__Brochothrix",
                                                      "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Pelomonas"))
Ploteo_StIandII <- ggplot(StIandII_PT_vplot_pivot,aes(x = Bacteria,
                                                      y = Valores,
                                                      color = Bacteria)) + 
                            geom_violin()+ 
                            theme(legend.position = "none")
Ploteo_StIandII

# Violin plots  Stage III and IV

StIIIandIV_PT_vplot <- COAD_StIIIandIV_biom_PT %>% select(...1,
                                                      k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Gallibacterium,
                                                      k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Aggregatibacter,
                                                      k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Leptotrichiaceae.g__Leptotrichia,
                                                      k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Buchnera,
                                                      k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Limnohabitans,
                                                      k__Archaea.p__Euryarchaeota.c__Methanococci.o__Methanococcales.f__Methanococcaceae.g__Methanococcus,
                                                      k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium,
                                                      k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Bradyrhizobiaceae.g__Rhodopseudomonas,
                                                      k__Bacteria.p__Cyanobacteria.o__Nostocales.f__Nostocaceae.g__Anabaena,
                                                      k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Listeriaceae.g__Brochothrix,
                                                      k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Pelomonas)


StIIIandIV_PT_vplot_pivot<-  StIIIandIV_PT_vplot %>% tidyr::pivot_longer(cols = -"...1",
                                                                     names_to = "Bacteria",
                                                                     values_to = "Valores") %>%
                                                     select(-...1) 



StIIIandIV_PT_vplot_pivot$Bacteria <- factor(StIIIandIV_PT_vplot_pivot$Bacteria,
                                           labels = c("k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Gallibacterium",
                                                      "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Aggregatibacter",
                                                      "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Leptotrichiaceae.g__Leptotrichia",
                                                      "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae.g__Buchnera",
                                                      "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Limnohabitans",
                                                      "k__Archaea.p__Euryarchaeota.c__Methanococci.o__Methanococcales.f__Methanococcaceae.g__Methanococcus",
                                                      "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium",
                                                      "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Bradyrhizobiaceae.g__Rhodopseudomonas",
                                                      "k__Bacteria.p__Cyanobacteria.o__Nostocales.f__Nostocaceae.g__Anabaena",
                                                      "k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Listeriaceae.g__Brochothrix",
                                                      "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Pelomonas"))
Ploteo_StIIIandIV <- ggplot(StIIIandIV_PT_vplot_pivot,aes(x = Bacteria,
                                                          y = Valores,
                                                          color = Bacteria)) + 
                             geom_violin()+
                             theme(legend.position = "none")
Ploteo_StIIIandIV


