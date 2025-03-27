
### **PEC1 - Analisis de Datos Omicos**

#### **Introduction**
In this repository, you will find the following files:

  - **`metabolites_data.tsv`**: Metabolomic data from Kang et al. (2018), which investigated fecal microbial metabolites in children with Autism Spectrum Disorder (ASD) and neurotypical children.  
  The dataset was obtained from [The Curated Gut Microbiome Metabolome Data Resource](https://github.com/borenstein-lab/microbiome-metabolome-curated-data).
  - **`metadata_samples.tsv`**: Sample metadata, also obtained from [The Curated Gut Microbiome Metabolome Data Resource](https://github.com/borenstein-lab/microbiome-metabolome-curated-data).
  - **`codigo_analisis.R`**: Contains the code used to process the data, create the SummarizedExperiment object, and perform data analysis.
  - **`Informe.Rmd`**: A detailed report describing data selection, processing, analysis, and interpretation of results.
  - **`summarized_experiment_kang_data.rda`**: The SummarizedExperiment object created from the dataset.

---
  
#### **Study Information**
**Original Study**: 

  - **Authors**: Dae-Wook Kang, Zehra Esra Ilhan, Nancy G. Isern, David W. Hoyt, Daniel P. Howsmon, Michael Shaffer, Catherine A. Lozupone, Juergen Hahn, James B. Adams, Rosa Krajmalnik-Brown  
  - **Title**: Differences in fecal microbial metabolites and microbiota of children with autism spectrum disorders 
  - **Journal**: Anaerobe  
  - **Year**: 2018  
  - **DOI**: 10.1016/j.anaerobe.2017.12.007

**Curated Data Source**: 

  - **Authors**: Efrat Muller, Yadid M. Algavi, and Elhanan Borenstein  
  - **Title**: The gut microbiome-metabolome dataset collection: a curated resource for integrative meta-analysis
  - **Journal**: NPJ Biofilms and Microbiomes  
  - **Year**: 2022  
  - **DOI**: 10.1038/s41522-022-00345-5  

---
  
#### **Dataset Description**

###### **Metabolites Data (`metabolites_data.tsv`)**
- **Format**: Tab-separated values (.tsv).
- **Rows**: 61 metabolites and stool parameters (e.g., pH).
- **Columns**: 44 samples (subjects).
- **Values**: Metabolite concentrations (μmole per gram of dry stool).

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

