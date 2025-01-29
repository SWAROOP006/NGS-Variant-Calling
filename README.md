# NGS Variant Calling Workflow
## Description
Analysis of whole-genome sequencing (WGS) data from BioProject [PRJEB62494](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB62494) (tongue cancer samples and cell lines) using:
- **Tools**: GATK, BWA, Samtools, SRA Toolkit
- **Reference Genome**: hg38 (chromosomes 6 & 7)

## Workflow Steps
1. Download raw data using SRA Toolkit.
2. Download Refrence genome and index with BWA.
3. Quality control with FastQC/fastp.
4. Alignment with BWA-MEM.
5. Variant calling with GATK HaplotypeCaller/Mutect2.
6. Annotation using ENSEMBL VEP.

[View Full Workflow Code →](workflow/scripts/workflow-code.sh)
