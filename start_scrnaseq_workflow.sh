##Run standard scRNAseq workflow locally
##TC March 6 2020

R --vanilla -f auto_workflow_min_r_script.R --args integration=YES metadata=metadata_file.txt
