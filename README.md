# CSI-Flagship-Manuscript-Scripts

 wes_pipeline - directory containing all files required for joint genotyping CSI Flagship cohort
 wes_pipeline/scripts/hgsc_batch_processing_hg19.snakemake - snakefile detailing processing steps
 wes_pipeline/scripts/hgsc_batch_processing_hg19.sh - shell script for running and submitting master job for WES pipeline
 wes_pipeline/resources/processing_references_hg19_locus.json - json containing reference files for running WES pipeline
 wes_pipeline/resources/processing_cluster_locus.json - json containing cluster parameters for WES pipeline steps
 pc_plot.py - makes the ancestry PCA panel necessary for the ancestry plot (Fig. 1B). To run this script, update the directory of your ancestry PC spreadsheet (generated using Somalier) in line 16 and update output directory at line 52 and 66. Than call pc_plot.py in your terminal.
