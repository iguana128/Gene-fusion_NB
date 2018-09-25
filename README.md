# Gene-fusion_NB
This project aims to unravel the gene fusion profiles in 498 neuroblastoma patients for regrouping patients into different risk group to look for better survival and potential drug repositioning opportunities.
# Prerequisites:
- The sequencing analysis for gene fusion detection were carried out under HPC environmental.<br/>
The details requirement and scripts for identifying the fusions were listed in <br/>
**Gene_fusion_detection_command.txt** 

- The downstream statistical analysis was conducted under R (version 3.3.2) program. 
The required packages including: <br/>
**ggplot2, reshape2, RColorBrewer, survminer, survival, venneuler, VennDiagram, and NMF**

# file/fold desciption: 
- **Data**: the gene fusion calling results from three calling algorithms incluing ChimeraScan, SOAPfuse and TopHat-Fusion
- **DEG**: the deffiential expressed genes (DEGs) for the refined high-risk group patients for each gen fusion algorithm
- **Processed_results**: the processed gene fusion results based on 
