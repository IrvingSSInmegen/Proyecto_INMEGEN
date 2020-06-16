# Filtrar las observaciones que sean correspondientes a cancer de colon
colon  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma") 
# Total de 1017 observaciones que cumplen "Colon Adenocarcinoma"

# Posibles candidatos para filtrar: investigation == "TCGA-COAD"
#                                   histological_diagnosis_label == "Colon Adenocarcinoma"
#                                   primary_site == "Colorectal"

Stage_IandII  <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                pathologic_stage_label %in% c("Stage I","Stage IA","Stage II",
                                                              "Stage IIA","Stage IIB","Stage IIC")) 

Stage_IIIandIV <- metadata %>% filter(disease_type == "Colon Adenocarcinoma" & 
                                        pathologic_stage_label %in% c ("Stage III","Stage IIIA",
                                                                       "Stage IIIB","Stage IIIC",
                                                                       "Stage IV", "Stage IVA",
                                                                       "Stage IVB"))
# Total de 1017 observaciones que cumple "Colon Adenocarcinoma"
# Total de 576 que cumplen Stage I and II
# Total de 423 que cumplen Stage III and IV

# Conclusión Hay 18 observaciones que no contienen la informacion del Stage
# ¿es necesario cambiar nuestro filtro?
# No podemos filtrar por T o N ya que solo marca el crecimiento del tumon
# Y o si ya se paso a diferentes linfocitos cercanos
# o ya hizo metastasis

# Una vez separado los datos en dos tablas queremos ver
# que observaciones coinciden con las que SI tenemos
# conocimiento de viruses archeas o bacterias
# tenemos un total de 18116 observaciones
# de las cuales solo de 13517 se tiene conocimiento del microbioma

matchcolIyII <- dplyr::semi_join(x = kraken_matched,
                                 y = Stage_IandII,
                                 by = c("...1" = "...1")) %>%
                select(contains("Bacteria"))

matchcolIIIyIV <- dplyr::semi_join(x = kraken_matched,
                                   y = Stage_IIIandIV,
                                   by = c("...1" = "...1")) %>%
                         select(contains("Bacteria"))

# tenemos lectura de 234 observaciones para Stage I and  II
# tenemos lectura de 177 observaciones para Stage III and IV

# Necesitamos saber que especies microbiana son mas abundantes en cada caso
# Primero hay que normalizar?

# recuento logaritmico por millon?

# Normalización cuantil?
#
#  Voom?

# Diagrama de cajas?

# Tengo que revisar con cuidado ??
#https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
#https://www.niehs.nih.gov/research/resources/software/biostatistics/pvca/index.cfm
