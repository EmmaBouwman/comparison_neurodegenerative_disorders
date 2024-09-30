## Comparing Neurodegenerative Disorders: Discovering patterns through protein-protein interaction network analysis

## Abstract
 Alzheimer’s disease (AD), Huntington’s disease (HD), and spinocerebellar ataxia (SCA),
 especially SCA1 and SCA3, are characterized by neuronal degeneration. This thesis aimed to
 identify overlapping patterns and common mechanisms among these disorders using protein
protein interaction networks (PPINs). PPINs for AD, HD, and SCA were constructed
 using data from KEGG, STRING, and WikiPathways. Network analysis, clustering, and
 functional enrichment analysis revealed shared biological processes, including responses to
 misfolded proteins and the endoplasmic reticulum-associated degradation (ERAD) pathway,
 across all three disorders. A consistent cluster was found in three different PPINs. So,
 despite their distinct genetic origins and clinical differences, these diseases share common
 mechanisms. This thesis demonstrates that PPINs can uncover patterns across different
 diseases that may help identify shared drug targets.


## Introduction
This repository contains all the code, external scripts, cytoscape files, results and images for this bachelors thesis. 

 ## External Automated Workflow
 For this thesis an automated workflow was used to create the different C-PPIN. This automated workflow was originally created by Aster de Boer (A.d.B). For this thesis the automated workflow is modified for each individual C-PPIN. The Jupyter notebook files with the name '[disease]_ppin_automated_workflow.ipynb' are tailord to with the specific data that was needed for each disease. 
 
 ## External Python Scripts
 This repository also contains scripts developed by other authors. These scripts can be found in the 'external_scripts' folder.
 - "external_scripts\kegg_data_converter_modified.py" - original author: Aster de Boer (A.d.B). This python script was originally created by A.d.B but is slightly modified for this study. 
 - "external_scripts\wikipathways_convert_uniprotIDs_NoriaYousufi.py" - original author: Noria Yousufi (N.Y.).  
