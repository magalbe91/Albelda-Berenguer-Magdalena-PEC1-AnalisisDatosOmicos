# PEC1 - Analisis de Datos Omicos ----

# Eleccion del Dataset ----

library(readr) 

# Dataset de: KANG-AUTISM-2018 (Kang, Dae-Wook, et al. "Differences in fecal microbial metabolites and microbiota of children with autism spectrum disorders." Anaerobe 49 (2018): 121-131)
# https://github.com/borenstein-lab/microbiome-metabolome-curated-data

# Cargamos los datos:
library(readr)

# Metadata
metadata <- read_tsv("~/Desktop/UOC_Master/Asignaturas/Semestre 2/1- Analisis de datos omicos/PEC1/Data microbiome/data/processed_data/KANG_AUTISM_2017/metadata.tsv")
head(metadata)

# Metabolitos
metabolites <- read_tsv("~/Desktop/UOC_Master/Asignaturas/Semestre 2/1- Analisis de datos omicos/PEC1/Data microbiome/data/processed_data/KANG_AUTISM_2017/mtb.tsv")
head(metabolites)

# Extra: Abundancia relativa de generos bacterianos
genera <- read_tsv("~/Desktop/UOC_Master/Asignaturas/Semestre 2/1- Analisis de datos omicos/PEC1/Data microbiome/data/processed_data/KANG_AUTISM_2017/genera.tsv")
head(genera)

# Creacion del SummarizedExperiment ----
library(SummarizedExperiment)

# 1- Assays: ----
# Convertimos los datos en matrices
metabolites_matrix <- as.matrix(metabolites[, -1]) # quitamos la primera columna con los nombres de las muestras
rownames(metabolites_matrix) <- metabolites$Sample # nombramos las filas de la matriz

# SummarizedExperiment quiere las muestras en las columnas
# Transponemos las matrices
if (nrow(metabolites_matrix) == nrow(metadata)) {
  metabolites_matrix <- t(metabolites_matrix)
}

# Creamos lista de assays
assay_list <- list(metabolites = metabolites_matrix)

# 2- colData:----
library(dplyr)
library(tibble)

# Asignamos el identificador de las muestras como rownames
colData_df <- metadata %>%
  column_to_rownames("Sample")

# Convertimos a DataFrame
class(colData_df)
colData_df <- DataFrame(colData_df)

# 3- rowData: ----
# DataFrame con nombres de los metabolitos
rowData_df <- DataFrame(Metabolite = rownames(metabolites_matrix))

# 4- SummarizedExperiment: ----
# Creamos el SummarizedExperiment
se <- SummarizedExperiment(assays = assay_list,
                           colData = colData_df,
                           rowData = rowData_df)

# 5- Añadir informacion
metadata(se) <- list(Experimenter = "Kang, Dae-Wook, et al.",
                     Title = "Differences in fecal microbial metabolites and microbiota of children with autism spectrum disorders",
                     Journal = "Anaerobe",
                     Year = 2018,
                     DOI = "10.1016/j.anaerobe.2017.12.007",
                     PMID = 29274915)

# Visualizamos que todo este correcto
se # todo el SummarozedExperiment
head(assay(se)) # datos (amount metabolitesXsamples)
head(colData(se)) # metadata de cada muestra
head(rowData(se)) # metabolitos
metadata(se) # informacion del estudio

# Creacion archivo .rda ----
save(se, file = "summarized_experiment_kang_data.rda")



