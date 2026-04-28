# TrioGCTA_EDU
# Scripts for  
"Parental Influence on Children's Educational Achievement: Analysing Direct and Indirect Genetic Effects through Trio-GCTA"

## Overview
This repository contains the analysis scripts used in the above study. The code is organized to reflect the full workflow, from data preparation to genetic analysis and final model estimation.

Due to data protection regulations related to MoBa and Statistics Norway (SSB), detailed comments and some data-specific annotations have been removed. However, all scripts remain functionally complete and correspond exactly to the original analyses, enabling full replication for authorized users.

---

## Data Access
The analyses rely on restricted-access data from:

- The Norwegian Mother, Father and Child Cohort Study (MoBa)
- Statistics Norway (SSB)

Researchers with approved access to these data sources should be able to reproduce the analyses using the scripts provided.  
For questions regarding data preparation or variable construction, please feel free to get in touch.

---

## Repository Structure

### 01_GPA_National_Assessments_Grades_Preparation
Data cleaning and preprocessing of:
- GPA (Grunnskolepoeng)
- National assessments (e.g., reading, math, English)

To comply with data protection requirements, detailed comments have been removed.    
The scripts also include extensive code for inspecting and understanding the data structure, which may make them appear relatively verbose. 
In addition, many processing steps are closely tailored to the specific structure of the underlying datasets. 
Users are therefore encouraged to adapt the data cleaning procedures based on their own data structure and analytical logic.

Nevertheless, the scripts are complete and reflect the original data cleaning pipeline.

---

### 02_Match_Grades_with_Genes_Get_Trios
Scripts for merging:
- Phenotypic data (education outcomes)
- Genetic data
- Trio

This step constructs the trio dataset linking children with their parents.

---

### 03_GRM_Julia
Julia scripts for computing the Genetic Relatedness Matrix (GRM).

⚠️ These computations are resource-intensive and typically require:
- High-performance computing (HPC) environments or clusters

---

### 04_TrioGCTA_Julia
Julia scripts for:
- Sample selection compatible with the Trio-GCTA model (threshold of GRM)
- Running the Trio-GCTA analyses
- Generating results used in the manuscript (tables and estimates)

⚠️ For large-scale data, execution on HPC/cluster systems is recommended.

---

## Reproducibility Note
All scripts are consistent with the original analysis pipeline used in the study.  
Although comments are minimized due to data protection constraints, the code is sufficient for full replication given appropriate data access.

---

## Contact
For questions, collaboration, or clarification regarding the scripts or data preparation, please contact:

Qiyuan Peng  
📧 qiyuanp@uio.no

---

## Disclaimer
This repository does not include any raw or processed data.  
All analyses must be conducted in compliance with the data access agreements of MoBa and Statistics Norway
