# PEC1 - Analisis de Datos Omicos ----

# Eleccion del Dataset ----

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
# Transponemos la matriz
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
colData_df <- DataFrame(colData_df)

# 3- rowData: ----
# DataFrame con nombres de los metabolitos
rowData_df <- DataFrame(Metabolite = rownames(metabolites_matrix))

# 4- SummarizedExperiment: ----
# Creamos el SummarizedExperiment
se <- SummarizedExperiment(assays = assay_list,
                           colData = colData_df,
                           rowData = rowData_df)

# 5- Añadir informacion del estudio 
metadata(se) <- list(Experimenter = "Kang, Dae-Wook, et al.",
                     Title = "Differences in fecal microbial metabolites and microbiota of children with autism spectrum disorders",
                     Journal = "Anaerobe",
                     Year = 2018,
                     DOI = "10.1016/j.anaerobe.2017.12.007",
                     PMID = 29274915,
                     Subjects_Information = "21 neurotypical children; 23 children with Autism Spectrum Disorders (ASD), age between 4 to 17 years old",
                     Samples_Information = "Fecal samples from the 44 subjects",
                     Metabolomics_Approach = "Targeted, 1H-NMR metabolite data",
                     Metabolites_Unit = "μmole per g-dry stool")

# Visualizamos que todo este correcto
se # todo el SummarizedExperiment
head(assay(se)) # datos (metabolitesXsamples)
head(colData(se)) # metadata de cada muestra
head(rowData(se)) # metabolitos
metadata(se) # informacion del estudio

# Creacion archivo .rda ----
save(se, file = "summarized_experiment_kang_data.rda")

# Analisis exploratorio ----
# 1- Resumen estadístico de metabolitos
summary(t(assay(se))) 

# 2- Visualización de todos los metabolitos por grupo (Neurotypical vs ASD)
library(reshape2)
library(ggplot2)

# Transformamos la estructura de la matriz
metabolites_long <- melt(assay(se)) 
colnames(metabolites_long) <- c("Metabolite", "Sample", "Concentration") # nombramos las columnas de la nueva tabla
metabolites_long <- merge(metabolites_long, colData(se), by.x = "Sample", by.y = "row.names") # añadimos el metadata
metabolites_long_filtered <- as.data.frame(metabolites_long) %>%
  filter(!Metabolite %in% c("pH", "DSS-d6 (Chemical Shape Indicator)")) # eliminamos columnas que no son metabolitos a estudiar

# Realizamos el grafico con ggplot
ggplot(metabolites_long_filtered, 
       aes(x = Study.Group, y = Concentration, fill = Study.Group)) +
  geom_boxplot() +
  facet_wrap(~ Metabolite, scales = "free") +
  labs(title = "Metabolites per Study Group",
       y = "μmole metabolite per g-dry stool") +
  theme_minimal() +
  theme(legend.position = "none") 

# 3- PCA
# Extraemos la matriz de metabolitos del SummarizedExperiment
metabolites_matrix <- assay(se)

# Filtramos los datos de la matriz
# Eliminamos datos que no corresponden a metabolitos
metabolites_matrix_filtered <- metabolites_matrix[!rownames(metabolites_matrix) %in% c("pH", "DSS-d6 (Chemical Shape Indicator)"), ]
# Buscamos si hay valores NA (y filtramos si hay)
sum(is.na(metabolites_matrix_filtered))
# Eliminamos metabolitos no detectados (valor = 0, para todas las muestras)
zero_variance_columns <- apply(metabolites_matrix_filtered, 1, var) == 0 # buscamos los metabolitos que tengan una varianza 0
metabolites_matrix_filtered_no_zero_variance <- metabolites_matrix_filtered[!zero_variance_columns, ] # eliminamos los metabolitos de la matriz

# Transponemos la matriz de metabolitos (para que coincida con la organizacion de los metadatos)
t_metabolites_matrix <- t(metabolites_matrix_filtered_no_zero_variance)

# Normalizamos los datos
metabolites_centered <- scale(t_metabolites_matrix, 
                              center = TRUE, # datos centrados (matriz de covarianzas)
                              scale = FALSE) # datos en la misma unidad (no escalamos)

# Realizamos la PCA
pca_result <- prcomp(metabolites_centered)
summary(pca_result)

# Visualizamos los resultados
# Extraemos los porcentajes de varianza explicada por cada PC
explained_variance <- summary(pca_result)$importance[2, ] * 100  

# Extraemos los metadatos
colData <- colData(se)

# Creamos un data frame con las coordenadas de los primeros dos PCs y porcentajes
pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], 
                     StudyGroup = colData$Study.Group,
                     PC1_var = explained_variance[1],
                     PC2_var = explained_variance[2])

# Hacemos el grafico 
ggplot(pca_df, 
       aes(x = PC1, y = PC2, 
           color = StudyGroup)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(pca_df)), 
            vjust = -0.5, size = 3) +
  stat_ellipse(aes(group = StudyGroup), 
               level = 0.95, linetype = 2, size = 1) +
  labs(title = "PCA Metabolites per Study Group",
       x = paste("PC1 (", round(explained_variance[1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(explained_variance[2], 2), "%)", sep = "")) +
  theme_minimal()

# Que metabolitos afectan más a los dos primeros PCs?
# Obtenemos los vectores de carga (loadings) para cada PC
loadings_df <- data.frame(
  PC1_Loading = pca_result$rotation[, 1],
  PC2_Loading = pca_result$rotation[, 2]
)

# Ordenamos por la magnitud absoluta en PC1 y PC2
top_pc1 <- loadings_df[order(abs(loadings_df$PC1_Loading), decreasing = TRUE), ]
top_pc2 <- loadings_df[order(abs(loadings_df$PC2_Loading), decreasing = TRUE), ]

# Mostramos los 10 metabolitos que más afectan a PC1 y PC2
head(top_pc1[, "PC1_Loading", drop = FALSE], 10)
head(top_pc2[, "PC2_Loading", drop = FALSE], 10)














