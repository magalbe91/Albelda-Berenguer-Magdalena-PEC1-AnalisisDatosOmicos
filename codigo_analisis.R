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

# SummarizedExperiment ----
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
colData_DF <- DataFrame(colData_df)

# 3- rowData: ----
# DataFrame con nombres de los metabolitos
rowData_DF <- DataFrame(Metabolite = rownames(metabolites_matrix))

# 4- SummarizedExperiment: ----
# Creamos el SummarizedExperiment
se <- SummarizedExperiment(assays = assay_list,
                           colData = colData_DF,
                           rowData = rowData_DF)

# 5- Añadir informacion del estudio 
metadata(se) <- list(Researcher = "Kang, Dae-Wook, et al.",
                     Title = "Differences in fecal microbial metabolites and microbiota of children with autism spectrum disorders",
                     Journal = "Anaerobe",
                     Year = 2018,
                     DOI = "10.1016/j.anaerobe.2017.12.007",
                     PMID = 29274915,
                     Subjects_Information = "21 neurotypical children; 23 children with Autism Spectrum Disorders (ASD), age between 4 to 17 years old",
                     Samples_Information = "Fecal samples from the 44 subjects",
                     Metabolomics_Approach = "Targeted, 1H-NMR metabolite data",
                     Metabolites_Unit = "μmol per g-dry stool",
                     Data_Download = "https://github.com/borenstein-lab/microbiome-metabolome-curated-data?tab=readme-ov-file")

# Visualizamos que todo este correcto
se # todo el SummarizedExperiment
head(assay(se)) # datos (metabolitesXsamples)
head(colData(se)) # metadata de cada muestra
head(rowData(se)) # metabolitos
metadata(se) # informacion del estudio


# Archivo .rda ----
save(se, file = "summarized_experiment_kang_data.rda")


# Analisis exploratorio ----
# 1- Resumen estadístico de metabolitos ----
stats <- summary(t(assay(se))) 
stats
write.csv(stats, file = "stats_result.csv")

# 2- Visualización de todos los metabolitos por grupo (Neurotypical vs ASD) ----
library(reshape2)
library(ggplot2)

# Extraemos los datos del SummarizedExperiment y transformamos la estructura de la matriz
metabolites_long <- melt(assay(se)) 
colnames(metabolites_long) <- c("Metabolite", "Sample", "Concentration") # nombramos las columnas del nuevo dataframe
metabolites_long <- merge(metabolites_long, colData(se), by.x = "Sample", by.y = "row.names") # añadimos el metadata
metabolites_long_filtered <- as.data.frame(metabolites_long) %>%
  filter(!Metabolite %in% c("pH", "DSS-d6 (Chemical Shape Indicator)")) # eliminamos columnas que no son metabolitos a estudiar

# Realizamos el grafico con ggplot
ggplot(metabolites_long_filtered, 
       aes(x = Study.Group, y = Concentration, fill = Study.Group)) +
  geom_boxplot() +
  facet_wrap(~ Metabolite, scales = "free") +
  labs(title = "Metabolites per Study Group",
       y = "μmol / g-dry stool") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x=element_blank()) 

ggsave("metabolites_boxplot.png", bg = "white", width = 13, height = 9, dpi = 300)

# Wilcoxon test ----
# Calculamos diferencias significativas de cada metabolito entre los grupos de estudio
# Realizamos el test de Shapiro-Wilk para evaluar la normalidad antes de aplicar el test de Wilcoxon 
shapiro_results <- metabolites_long_filtered %>%
  group_by(Metabolite) %>%
  summarise(
    # Verificamos si la concentración tiene más de un valor único antes de aplicar el Shapiro-Wilk
    p_value_shapiro = ifelse(length(unique(Concentration)) > 1,
                             shapiro.test(Concentration)$p.value,
                             NA)  # Asignamos NA si todos los valores son iguales
  )

# Mostramos los p-value del Shapiro-Wilk
head(shapiro_results)

# Los datos no siguen una distribución normal: Realizamos el test de Wilcoxon sin exactitud
wilcoxon_results <- metabolites_long_filtered %>%
  group_by(Metabolite) %>%
  summarise(p_value = wilcox.test(Concentration ~ Study.Group, exact = FALSE)$p.value)

# Aplicamos el ajuste de Benjamini-Hochberg a los valores p
wilcoxon_results$adjusted_p_value <- p.adjust(wilcoxon_results$p_value, method = "BH")

# Mostramos los p-value de Wilcoxon y ajustados
head(wilcoxon_results)

write.csv(wilcoxon_results, file = "wilcoxon_results.csv")


# 3- PCA ----
# Extraemos de nuevo la matriz de metabolitos del SummarizedExperiment
metabolites_matrix <- assay(se)

# Filtramos los datos de la matriz
sum(is.na(metabolites_matrix)) # buscamos si hay valores NA (y filtramos si hay)
metabolites_matrix_filtered <- metabolites_matrix[!rownames(metabolites_matrix) %in% 
                                                    c("pH", "DSS-d6 (Chemical Shape Indicator)"), ] # eliminamos datos que no corresponden a metabolitos
metabolites_matrix_filtered <- metabolites_matrix_filtered[rowSums(metabolites_matrix_filtered > 0) > 3, ] # eliminamos metabolitos no detectados (valor = 0, para todas las muestras) y metabolitos detectados solo en hasta 3 muestras

# Transponemos la matriz de metabolitos (Obtenemos: muestras en filas vs metabolitos en columnas)
t_metabolites <- t(metabolites_matrix_filtered)

# Normalizamos los datos
t_metabolites_centered <- scale(t_metabolites, 
                              center = TRUE, # datos centrados (matriz de covarianzas)
                              scale = FALSE) # datos en la misma unidad (no escalamos)

# Realizamos la PCA
pca_result <- prcomp(t_metabolites_centered)
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
       aes(x = PC1, y = PC2, color = StudyGroup)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(pca_df)), vjust = -0.8, size = 3) +
  stat_ellipse(aes(group = StudyGroup), level = 0.95, linetype = 2, linewidth = 0.5) +
  labs(title = "PCA",
       x = paste("PC1 (", round(explained_variance[1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(explained_variance[2], 2), "%)", sep = "")) +
  theme_minimal()

ggsave("pca_plot.png", bg = "white", width = 8, height = 5, dpi = 300)

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


# 4- Clustering ----
library(pheatmap)

# A) Agrupamos los metabolitos segun los grupos de estudio (Autistic vs Neurotypical)
# Escalamos la matriz filtrada (esta vez centramos y escalamos)
# Trasponemos para escalar y luego transponemos otra vez para mantener los metabolitos en las filas
metabolites_scaled <- t(scale(t(metabolites_matrix_filtered), center = TRUE, scale = TRUE))

# Ordenamos las muestras según el grupo de estudio
ordered_samples <- rownames(colData(se))[order(colData(se)$Study.Group)]
metabolites_scaled_ordered <- metabolites_scaled[, match(ordered_samples, colnames(metabolites_scaled))]

# Cambiamos el nombre de la anotación
# Creamos la anotación con el mismo orden que ordered_samples
annotation_col <- data.frame(Study_Group = colData(se)$Study.Group[match(ordered_samples, colnames(se))])
rownames(annotation_col) <- ordered_samples

# Visualizamos el mediante un heatmap como se agrupan los metabolitos segun el grupo de estudio
heatmap_plot <- pheatmap(metabolites_scaled_ordered, 
                         annotation_col = annotation_col, 
                         show_rownames = TRUE,  
                         show_colnames = TRUE,
                         fontsize_row = 6,  
                         cluster_cols = FALSE,  
                         clustering_distance_rows = "euclidean",
                         clustering_method = "ward.D2",
                         main = "Heatmap - Metabolites per Study Group")

# Usamos ggsave() para guardarlo
library(gridExtra)
ggsave("metabolites_heatmap.png", plot = heatmap_plot[[4]], width = 8, height = 6, dpi = 300)

# B) Agrupamos los sujetos del estudio 
# Calculamos la distancia entre sujetos (sobre la matriz escalada)
distancia_muestras <- dist(t(metabolites_scaled))

# Realizamos el clustering jerárquico usando el método de Ward
clustering <- hclust(distancia_muestras, method="ward.D2")

# Visualizamos el dendrograma (y guardamos la imagen)
png("dendrograma_clustering.png", width = 800, height = 600, res = 100, bg = "white")

plot(clustering, main = "Hierarchical Clustering Dendrogram")

dev.off()


# 5- Correlacion ----
# Calculamos la correlación entre los metabolitos (utilizamos la matriz filtrada)
cor_matrix <- cor(t_metabolites)
head(cor_matrix)

# Visualizar la matriz de correlación como un mapa de calor
correlation_plot <- pheatmap(cor_matrix, 
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             clustering_distance_rows = "euclidean",
                             clustering_method = "ward.D2",
                             main = "Heatmap - Correlation Metabolites")

# Usamos ggsave() para guardarlo
ggsave("correlation_heatmap.png", plot = correlation_plot[[4]], width = 9, height = 9, dpi = 300)


# Archivos de texto ----
# Guardamos la matriz de datos metabolicos en formato texto
write.table(assay(se), file = "metabolites_data.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Guardamos el metadata de las muestras
write.table(cbind(Sample = rownames(colData(se)), as.data.frame(colData(se))), 
            file = "metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Marckdown ----
# Desrrollamos el archivo Markdown explicativo
markdown_text <- "
### **PEC1 - Analisis de Datos Omicos**

#### **Introduction**
In this repository, you will find the following files:

  - **`metabolites_data.tsv`**: Metabolomic data from Kang et al. (2018), which investigated fecal microbial metabolites in children with Autism Spectrum Disorder (ASD) and neurotypical children.  
  The dataset was obtained from [The Curated Gut Microbiome Metabolome Data Resource](https://github.com/borenstein-lab/microbiome-metabolome-curated-data).
  - **`metadata.tsv`**: Sample metadata, also obtained from [The Curated Gut Microbiome Metabolome Data Resource](https://github.com/borenstein-lab/microbiome-metabolome-curated-data).
  - **`codigo_analisis.R`**: Contains the code used to process the data, create the SummarizedExperiment object, and perform data analysis.
  - **`Informe.pdf`**: A detailed report describing data selection, processing, analysis, and interpretation of results.
  - **`summarized_experiment_kang_data.rda`**: The SummarizedExperiment object created from the dataset.

---
  
#### **Dataset Description**

###### **Metabolites Data (`metabolites_data.tsv`)**
- **Format**: Tab-separated values (.tsv).
- **Rows**: 61 metabolites and stool parameters (e.g., pH).
- **Columns**: 44 samples (subjects).
- **Values**: Metabolite concentrations (μmol per gram of dry stool).

This dataset contains metabolomic profiles derived from targeted 1H-NMR analysis of fecal samples.
Each row represents a metabolite or stool parameter, and each column corresponds to a sample.

###### **Sample Metadata (`metadata.tsv`)**
- **Format**: Tab-separated values (.tsv).
- **Rows**: 44 samples (subjects).
- **Columns**: Metadata variables.

The metadata includes information about the 44 study participants, such as:

  - Study group (Neurotypical or Autistic),
  - Age, gender, and gastrointestinal symptoms,
  - Clinical assessments (Autism Treatment Evaluation Checklist (ATEC), Pervasive Developmental Disorder Behavior Inventory (PDD-BI)).

---
  
References: 

  1. The Curated Gut Microbiome Metabolome Data Resource (https://github.com/borenstein-lab/microbiome-metabolome-curated-data).
  Muller E, Algavi YM, Borenstein E. The gut microbiome-metabolome dataset collection: a curated resource for integrative meta-analysis. NPJ Biofilms Microbiomes. 2022 Oct 15;8(1):79. doi: 10.1038/s41522-022-00345-5.
  2. Kang DW, Ilhan ZE, Isern NG, Hoyt DW, Howsmon DP, Shaffer M, Lozupone CA, Hahn J, Adams JB, Krajmalnik-Brown R. Differences in fecal microbial metabolites and microbiota of children with autism spectrum disorders. Anaerobe. 2018 Feb;49:121-131. doi: 10.1016/j.anaerobe.2017.12.007. 
"

# Guardamos el documento Markdown
writeLines(markdown_text, "README.md")

# Guradamos tambien en html el archivo Markdown
rmarkdown::render("README.md", output_format = "html_document")







