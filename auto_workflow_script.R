##Run standard scRNAseq workflow on HTC
##TC Feb 6 2020

library(rmarkdown)
library(knitr)
render(input="seurat_general_workflow_htc.rmd",output_file="seurat_workflow_htc_results.html")
