##Run standard scRNAseq workflow locally
##TC March 6 2020

R --vanilla -f auto_workflow_script.R --args integration=NO \
  integrate_by=NONE \
  metadata=./ancillary_files/metadata_file.csv \
  unhashing=YES \
  unhashing_file=./ancillary_files/cell_unhashing_identities.csv \
  tcr_bcr_addition=NO
