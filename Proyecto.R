###############################
#######                 #######
#######     CODIGO      #################
#######                 #######
###############################
library(vroom)
# Data

metadata <- vroom::vroom(file = "Metadata-TCGA-All-18116-Samples.csv")

kraken_TCGA_matched2SHOGUN <- vroom::vroom(file = "Kraken-TCGA-Matched2SHOGUN-Raw-Data.csv")

# implemented a pipeline that converted discrete taxonomic counts into log-counts
# per million (log-cpm) per sample using Voom, and performed supervised normalization (SNM)
# Data available in https://drive.google.com/drive/folders/18V2ON-Go5AeEtZLe1f9EeJToWOhg81ab

kraken_TCGA_matched2SHOGUN_VoomSNM <- vroom::vroom(file = "Kraken-TCGA-Matched2SHOGUN-Voom-SNM-Quantile-Data.csv")

# Filter
COAD_metadata  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma") # Individual with COAD
                                                                              # 1017 samples 42 variables

# Stage Early
COAD_metadata_StageIandII  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                                  pathologic_stage_label %in% c("Stage I","Stage IA",
                                                                                "Stage II","Stage IIA",
                                                                                "Stage IIB","Stage IIC")) # 576 Samples Stage I and II

COAD_metadata_StageI  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                             pathologic_stage_label %in% c("Stage I","Stage IA")) # 180 samples only stage I ( 1 sample Stage IA)
COAD_metadata_StageII  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                              pathologic_stage_label %in% c("Stage II","Stage IIA",
                                                                            "Stage IIB","Stage IIC")) # 396  samples only stage II ( 1 sample Stage IIC)

COAD_metadata_StIandII_PT <- COAD_metadata_StageIandII %>% filter(sample_type == "Primary Tumor") # 473 samples primary tumor
COAD_metadata_StIandII_BDN <- COAD_metadata_StageIandII %>% filter(sample_type == "Blood Derived Normal") # 64 samples blood derived normal
COAD_metadata_StIandII_STN <- COAD_metadata_StageIandII %>% filter(sample_type == "Solid Tissue Normal") # 38 samples solid tissue normal


# Stage Lated

COAD_metadata_StageIIIandIV <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                                   pathologic_stage_label %in% c ("Stage III","Stage IIIA",
                                                                                  "Stage IIIB","Stage IIIC",
                                                                                  "Stage IV", "Stage IVA","Stage IVB")) # 423 observaciones Stage III y IV

COAD_metadata_StageIII <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                              pathologic_stage_label %in% c ("Stage III","Stage IIIA",
                                                                             "Stage IIIB","Stage IIIC")) # 284 samples only Stage III
COAD_metadata_StageIV <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                             pathologic_stage_label %in% c ("Stage IV", "Stage IVA",
                                                                                        "Stage IVB")) # 139 sample only Stage IV
                                                                                                      # (2 samples Stage IVB)

COAD_metadata_StIIIandIV_PT <- COAD_metadata_StageIIIandIV %>% filter(sample_type == "Primary Tumor") #348 samples primary tumor
COAD_metadata_StIIIandIV_BDN <- COAD_metadata_StageIIIandIV %>% filter(sample_type == "Blood Derived Normal") # 43 samples Blood derived normal
COAD_metadata_StIIIandIV_STN <- COAD_metadata_StageIIIandIV %>% filter(sample_type == "Solid Tissue Normal") #32 samples solid tissue normal

#15 Phatologic stage label: not available
#3  Phatologic stage label: NA
#1  Recurrent Tumor
#1  Metastatic

#Filter virus in dataframe

kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus <- kraken_TCGA_matched2SHOGUN_VoomSNM %>% select(- contains("Viroid"))
namecolums <- colnames(kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus)
grep(pattern = "Viroid",x = namecolums) # 0 viroids

# Bacteria/archea Stage I and II
COAD_StIandII_biom_PT <- dplyr::semi_join(x = kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus,
                                          y = COAD_metadata_StIandII_PT,
                                          by = c("...1" = "...1"))  #190 samples
COAD_StIandII_biom_BDN <- dplyr::semi_join(x = kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus,
                                           y = COAD_metadata_StIandII_BDN,
                                           by = c("...1" = "...1"))  # 29 samples
COAD_StIandII_biom_STN <- dplyr::semi_join(x = kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus,
                                           y = COAD_metadata_StIandII_STN,
                                           by = c("...1" = "...1"))  #15 samples

# Bacteria/archea Stage III and IV
COAD_StIIIandIV_biom_PT <- dplyr::semi_join(x = kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus,
                                            y = COAD_metadata_StIIIandIV_PT,
                                            by = c("...1" = "...1")) #150 samples
COAD_StIIIandIV_biom_BDN <- dplyr::semi_join(x = kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus,
                                             y = COAD_metadata_StIIIandIV_BDN,
                                             by = c("...1" = "...1")) # 17 samples
COAD_StIIIandIV_biom_STN <- dplyr::semi_join(x = kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus,
                                             y = COAD_metadata_StIIIandIV_STN,
                                             by = c("...1" = "...1")) #9 samples
# Matrices de conteo de abundancia

# StageIandII Primary Tumor , Blood derived normal , Solid tissue normal.
COAD_StIandII_biom_PT_trans <- t(COAD_StIandII_biom_PT) %>% row_to_names(row_number = 1)
COAD_StIandII_biom_BDN_trans <- t(COAD_StIandII_biom_BDN) %>% row_to_names(row_number = 1)
COAD_StIandII_biom_STN_trans <- t(COAD_StIandII_biom_STN) %>% row_to_names(row_number = 1)

# StageIIIandIV Primary Tumor , Blood derived normal , Solid tissue normal
COAD_StIIIandIV_biom_PT_trans <- t(COAD_StIIIandIV_biom_PT) %>% row_to_names(row_number = 1)
COAD_StIIIandIV_biom_BDN_trans <- t(COAD_StIIIandIV_biom_BDN) %>% row_to_names(row_number = 1)
COAD_StIIIandIV_biom_STN_trans <- t(COAD_StIIIandIV_biom_STN) %>% row_to_names(row_number = 1)


# Metadata COAD Primary tumor Stage I and II match2biom
COAD_metadata_StgIandII_PT_match2biom <- dplyr::semi_join(x = COAD_metadata_StIandII_PT,
                                                          y = COAD_StIandII_biom_PT,
                                                          by = c("...1" = "...1"))

COAD_metadata_StgIandII_PT_match2biom_filter <- COAD_metadata_StgIandII_PT_match2biom %>% 
                                                       filter(!...1 %in% c("s13650","s13242","s13786",
                                                                           "s13575","s13843","s13719",
                                                                           "s13574","s13695","s13548",
                                                                           "s13679","s13542","s13831",
                                                                           "s13543","s13699","s13720",
                                                                           "s13685","s13696","s13683",
                                                                           "s13098","s13674","s13550",
                                                                           "s13072","s13756","s13727",
                                                                           "s13768","s13547","s13670",
                                                                           "s13515","s13672","s13689",
                                                                           "s13667","s13511","s13532",
                                                                           "s13563","s13781","s13663",
                                                                           "s13537","s13666","s13038",
                                                                           "s13557","s13544","s13726",
                                                                           "s13753","s13516","s13519",
                                                                           "s13716","s13040","s13509",
                                                                           "s13759","s13538","s13517",
                                                                           "s13677","s13513"))  # Result: 137 samples 

Diccionario_STIandII_PT <- COAD_metadata_StgIandII_PT_match2biom_filter %>% select(...1,
                                                                                  case_uuid,
                                                                                  aliquot_uuid,
                                                                                  sample_uuid,
                                                                                  filename,
                                                                                  gdc_file_uuid) %>% arrange(case_uuid)

#  Matriz de genes Stage I and II
# setwd("C:/) Necesito cambiar de directorio para cargar solo los
# archivos de etapa temprana

archivos = list.files(pattern = "*.gz")
files <- as_tibble(archivos)

mi_lista <- lapply(1:nrow(files), function(i){
  
              mi_archivo = files[["value"]][i]
              mi_caseid  = Diccionario_STIandII_PT[["case_uuid"]][i]
              mi_tabla = vroom::vroom(file = mi_archivo, col_names = F)
              mis_colnames = c("gene", mi_caseid) 
              colnames(mi_tabla) <- mis_colnames
  
            return(mi_tabla)
            }) 

StageIandII_genes <- mi_lista %>% reduce(full_join, by = "gene")

# Metadata COAD Primary tumor Stage III and Iv match2biom
COAD_metadata_StgIIIandIV_PT_match2biom <- dplyr::semi_join(x = COAD_metadata_StIIIandIV_PT,
                                                            y = COAD_StIIIandIV_biom_PT,
                                                            by = c("...1" = "...1"))

COAD_metadata_StgIIIandIV_PT_match2biom_filter <- COAD_metadata_StgIIIandIV_PT_match2biom %>% 
                                                       filter(!...1 %in% c("s13682","s13701","s13842",
                                                                           "s13556","s13724","s13527",
                                                                           "s13567","s13562","s13822",
                                                                           "s13680","s13669","s12880",
                                                                           "s13840","s13090","s13566",
                                                                           "s13540","s13505","s13652",
                                                                           "s13649","s13572","s13013",
                                                                           "s13565","s13523","s13723",
                                                                           "s13353","s13518","s13651",
                                                                           "s13728","s13194","s13580",
                                                                           "s13573","s13764","s13096",
                                                                           "s13678","s13065","s13551",
                                                                           "s13721","s13788","s12914",
                                                                           "s13657","s13525","s13237",
                                                                           "s13681","s13839","s13700",
                                                                           "s12852")) # 104 samples

Diccionario_STIIIandIV_PT <- COAD_metadata_StgIIIandIV_PT_match2biom_filter %>% select(...1,
                                                                                       case_uuid,
                                                                                       aliquot_uuid,
                                                                                       sample_uuid,
                                                                                       filename,
                                                                                       gdc_file_uuid) %>% arrange(case_uuid)

#  Matriz de genes Stage III and IV
# setwd("C:/) Necesito cambiar de directorio para cargar solo los
# archivos de etapa tardía

archivos.1 = list.files(pattern = "*.gz")
files.1 <- as_tibble(archivos.1)

mi_lista.1 <- lapply(1:nrow(files.1), function(i){
  
                mi_archivo.1 = files.1[["value"]][i]
                mi_caseid.1  = Diccionario_STIIIandIV_PT[["case_uuid"]][i]
                mi_tabla.1 = vroom::vroom(file = mi_archivo.1, col_names = F)
                mis_colnames.1 = c("gene", mi_caseid.1) 
                colnames(mi_tabla.1) <- mis_colnames.1
  
                return(mi_tabla.1)
               }) 

StageIIIandIV_genes <- mi_lista.1 %>% reduce(full_join, by = "gene")

#  Matrices de microbioma filtradas
#   Stage I and II
nombres_StIandII <- Diccionario_STIandII_PT %>% select(...1,case_uuid) %>% arrange(...1)

StageIandII_biom <- dplyr::full_join(x = nombres_StIandII,
                                     y = kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus,
                                     by = c("...1" = "...1")) %>%
                                      select(-...1)
StageIandII_biom_filter <- dplyr:: right_join(x = StageIandII_biom,
                                              y = nombres_StIandII,
                                              by = c("case_uuid" = "case_uuid")) %>%
                                              t() %>% row_to_names(row_number = 1)

#  Stage III and IV
nombres_StIIIandIV <- Diccionario_STIIIandIV_PT %>% select(...1,case_uuid) %>% arrange(...1)

StageIIIandIV_biom <- dplyr::full_join(x = nombres_StIIIandIV,
                                       y = kraken_TCGA_matched2SHOGUN_VoomSNM_Novirus,
                                       by = c("...1" = "...1")) %>% 
                                        select(-...1)
StageIIIandIV_biom_filter <- dplyr::right_join(x = StageIIIandIV_biom,
                                               y = nombres_StIIIandIV,
                                               by = c("case_uuid" = "case_uuid")) %>%
                                              t() %>% row_to_names(row_number = 1)

#  Guardar las matrices en archivo 
write.csv(StageIIIandIV_biom_filter, file="StageIIIandIV_bac_arch.csv")
write.csv(StageIandII_biom_filter, file="StageIandII_bac_arch.csv")
write.csv(StageIandII_genes, file="StageIandII_gene.csv")
write.csv(StageIIIandIV_genes, file="StageIIIandIV_gene.csv")

##########################################
#######                            #######
######    Differetial expression    ##################
#######                            #######
##########################################
# Preparando los datos

metadata_StIandII <- COAD_metadata_StgIandII_PT_match2biom_filter %>% 
                                    select(...1,case_uuid,pathologic_stage_label)%>%
                                         mutate(condition = case_when(pathologic_stage_label == "Stage I" ~ "early",
                                                                      pathologic_stage_label == "Stage IA" ~ "early",
                                                                      pathologic_stage_label == "Stage II" ~ "early",
                                                                      pathologic_stage_label == "Stage IIA" ~ "early",
                                                                      pathologic_stage_label == "Stage IIB" ~ "early",
                                                                      pathologic_stage_label == "Stage IIC" ~ "early"))%>%
                                                        select(...1,case_uuid,condition)

metadata_StIIIandIV <- COAD_metadata_StgIIIandIV_PT_match2biom_filter %>%
                                    select(...1,case_uuid,pathologic_stage_label) %>%
                                         mutate(condition = case_when(pathologic_stage_label == "Stage III" ~ "late",
                                                                      pathologic_stage_label == "Stage IIIA" ~ "late",
                                                                      pathologic_stage_label == "Stage IIIB" ~ "late",
                                                                      pathologic_stage_label == "Stage IIIC" ~ "late",
                                                                      pathologic_stage_label == "Stage IV" ~ "late",
                                                                      pathologic_stage_label == "Stage IVA" ~ "late",
                                                                      pathologic_stage_label == "Stage IVB" ~ "late"))%>%
                                                        select(...1,case_uuid,condition)

metadata_COAD <- rbind(metadata_StIandII,metadata_StIIIandIV)
Metadada_COAD <- metadata_COAD %>% select(-...1)

bac_early <- metadata_COAD %>% select(...1,case_uuid)

# Bacterias sin normalizar
kraken_TCGA_matched2SHOGUN_Novirus <- kraken_TCGA_matched2SHOGUN %>% select(- contains("Viroid"))
namecolums <- colnames(kraken_TCGA_matched2SHOGUN_Novirus)
grep(pattern = "Viroid",x = namecolums) # 0 viroids   

Bacteria_countdata <- dplyr::left_join(x = bac_early,
                                       y = kraken_TCGA_matched2SHOGUN_Novirus,
                                       by = c("...1" = "...1"))%>% select(-...1)%>% 
                                   t()%>% row_to_names(row_number = 1)

Bacteria_countdata.1 <- t(sapply(Bacteria_countdata, as.numeric))

# La matriz de Bacteria_countdata de hace de caracteres necesito hacerla numérica
# Solución

write.csv(Metadada_COAD,file = "metadata_uuid.csv",row.names = FALSE)
write.csv(Bacteria_countdata,file = "Countdata_Bacteria.csv")

Countdata_uuid_COAD <- read.csv(file = "Countdata_Bacteria.csv",row.names = 1)%>%as.matrix()
Metadata_uuid_COAD <- read.csv(file = "metadata_uuid.csv", row.names = 1)

#Ahora si aplicamos el paquete DESeq2
library(DESeq2)

Countdata_uuid_COAD = Countdata_uuid_COAD[rowSums(Countdata_uuid_COAD)>1,] 

dds = DESeqDataSetFromMatrix(countData = Countdata_uuid_COAD,
                             colData = Metadata_uuid_COAD,
                             design =~condition)
dds = DESeq(dds)

res = results(dds, contrast = c("condition","early","late"))
res = res[order(res$pvalue),]
summary(res)

head(res,10)

##########################
#######            #######
######    GRÁFICAS    ##################
#######            #######
##########################
library(igraph)

# Gráficas analizadas con python
# Early
G_early = read.graph(file = "early_analyzed.graphml",format = c("graphml"))
vertex_attr_names(G_early)
V(G_early)
list_attr_early <- get.data.frame(G_early,what = c("vertices"))
# Quitar nodos de grado cero
aislados_early = which(degree(G_early) == 0)
G2_early = delete.vertices(G_early,aislados_early)
nodos_nonulos_early <- get.data.frame(G2_early,what = c("vertices"))

write.graph(G2_early,file = "early_analyzed_positive.gml",format = "gml")

#Late
G_late = read.graph(file = "late_analyzed.graphml",format = c("graphml"))
vertex_attr_names(G_late)
V(G_late)
list_attr_late <- get.data.frame(G_late,what = c("vertices"))
# Quitar nodos de grado cero
aislados_late = which(degree(G_late) == 0)
G2_late = delete.vertices(G_late,aislados_late)
nodos_nonulos_late <- get.data.frame(G2_late,what = c("vertices"))

write.graph(G2_late,file = "late_analized_positive.gml",format = "gml")



G_early = read_graph(file = "early.graphml",format = c("graphml"))
vertex_attr_names(G_early)
V(G_early)
nodos <- get.data.frame(G_early,what = c("vertices"))


##########################
#######                #######
######   VIOLIN PLOTS   ##################
#######                #######
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
                             theme(legend.position = "none") + 
                             coord_flip()
                             
                             
Ploteo_StIIIandIV



### SEPARADOS


StIandII_PT_vplot.1 <- COAD_StIandII_biom_PT %>% select(...1,k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Gallibacterium)

StIandII_PT_vplot_pivot.1 <-  StIandII_PT_vplot.1 %>% tidyr::pivot_longer(cols = -"...1",
                                                                          names_to = "Bacteria",
                                                                          values_to = "Valores") %>% select(-...1)
StIandII_PT_vplot_pivot.1$Bacteria <- factor(StIandII_PT_vplot_pivot.1$Bacteria,
                                               labels = c("k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Gallibacterium"))

Ploteo_StIandII.1 <- ggplot(StIandII_PT_vplot_pivot.1,aes(x = Bacteria,
                                                          y = Valores,
                                                          )) + 
                                     geom_violin(fill = "darkgreen")+
                                     theme(legend.position = "none") + 
                                     labs(title = "Pacientes con COAD en estados I y II.",
                                          subtitle = "tipo de muestra: tumor primario.",
                                          x = "Bacteria",
                                          y = "Normalized Abundance (log2-cpm)")
                                     
                                                       
                                           
Ploteo_StIandII.1

StIIIandIV_PT_vplot.1 <- COAD_StIIIandIV_biom_PT %>% select(...1,k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Gallibacterium)


StIIIandIV_PT_vplot_pivot.1<-  StIIIandIV_PT_vplot.1 %>% tidyr::pivot_longer(cols = -"...1",
                                                                             names_to = "Bacteria",
                                                                             values_to = "Valores") %>% select(-...1) 




StIIIandIV_PT_vplot_pivot.1$Bacteria <- factor(StIIIandIV_PT_vplot_pivot.1$Bacteria,
                                             labels = c("k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Gallibacterium"))

Ploteo_StIIIandIV.1 <- ggplot(StIIIandIV_PT_vplot_pivot.1,aes(x = Bacteria,
                                                              y = Valores)) + 
                                                              
                                       geom_violin(fill = "firebrick1")+
                                       theme(legend.position = "none") + 
                                       labs(title = "Pacientes con COAD en estados III y IV.",
                                            subtitle = "tipo de muestra: tumor primario.",
                                            x = "Bacteria",
                                            y = "Normalized Abundance (log2-cpm)")
                                


Ploteo_StIIIandIV.1


def.par <- par(no.readonly = TRUE)
Conf2x2 <- matrix(c(1:4), 2, 2, byrow = TRUE)
Conf2x2
layout.show(1)
