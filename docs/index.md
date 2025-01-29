# NGS Workflow Documentation
## Overview
- **Goal**: This project focuses on identifying somatic variants in whole-genome sequencing (WGS) data from tongue cancer samples and cell lines. The workflow leverages industry-standard tools for alignment, quality control, and variant calling, specifically targeting **chromosomes 6 and 7** of the human reference genome (hg38).
- **Data Source**: [PRJEB62494](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB62494).
- **Tools Used**:  
  ![GATK](https://img.shields.io/badge/GATK-4.0-blue) ![BWA](https://img.shields.io/badge/BWA-0.7.17-green)
  - **BWA**: Burrows-Wheeler Aligner for read alignment.
  - **Samtools**: Processing and indexing BAM files.
  - **GATK (v4)**: Variant calling with HaplotypeCaller.
  - **SRA Toolkit**: Downloading raw sequencing data.
  - **FastQC**: Quality control reports.
- **Data**:
  - **Samples**: Three accessions (`ERR11468775`, `ERR11468776`, `ERR11468777`).
  - **Reference Genome**: hg38 (chromosomes 6 and 7).



### Key Steps
1. **Data Acquisition**:
   - Download raw FASTQ files from SRA (500,000 reads per sample).
   - Example command:
     ```bash
     fastq-dump --split-files --gzip -X 500000 ERR11468775
     ```
     

2. **Reference Genome Setup**:
   - Chromosomes 6 and 7 from UCSC hg38.
   - Merged and indexed with BWA.
   

3. **Quality Control**:
   - FastQC reports for raw and processed data.
   - Example QC output:  
     ![FastQC Report](workflow/results/Picture1.png)

4. **Alignment**:
   - BWA-MEM for paired-end alignment.
   - Samtools for sorting and indexing BAM files.
   - Example command:
     ```bash
     bwa mem ref_genome/hg38_chr6_7.fa sample_1.fastq.gz sample_2.fastq.gz | samtools sort -o sample.bam
     ```

5. **Variant Calling**:
   - GATK HaplotypeCaller for germline variants, Mutect2 for somatic variants.
   - Example output:  
     ![IGV Screenshot](workflow/results/Picture1.png)

6. **Annotation**: ENSEMBL VEP.

[Explore the Code →](https://github.com/SWAROOP006/NGS-Variant-Calling/tree/main/workflow/scripts)
