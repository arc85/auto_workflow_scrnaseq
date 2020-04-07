##Run standard scRNAseq workflow locally
##TC March 6 2020

R --vanilla -f auto_workflow_script.R --args integration=NO metadata=metadata_file.txt unhashing=cell_unhashing_identities.csv
