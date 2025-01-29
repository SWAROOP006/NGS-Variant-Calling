#!/bin/bash
set -e  # Exit script on any error

# =========================================
# NGS Variant Calling Workflow
# BioProject: PRJEB62494 (Tongue Cancer)
# Reference Genome: hg38 (chr6 & chr7)
# =========================================

# --------------------------
# 1. Setup Environment
# --------------------------
# Create directory structure
mkdir -p \
  ref_genome \
  raw_data_mapping \
  qc_reports \
  aligned_data \
  variants \
  logs

# Install system dependencies
sudo apt-get update
sudo apt-get install -y \
  docker.io \
  wget \
  git \
  fastqc \
  samtools \
  bcftools

# Initialize log file
exec > logs/workflow.log 2>&1

# --------------------------
# 2. Reference Genome Setup
# --------------------------
echo "=== SETTING UP REFERENCE GENOME ==="
cd ref_genome

# Download chr6 and chr7
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz
gunzip *.fa.gz

# Merge chromosomes
cat chr6.fa chr7.fa > hg38_chr6_7.fa

# BWA Indexing
git clone https://github.com/lh3/bwa.git
cd bwa && make
cd ..
bwa index hg38_chr6_7.fa
cd ..

# --------------------------
# 3. Download Raw Data
# --------------------------
echo "=== DOWNLOADING RAW DATA ==="
cd raw_data_mapping

# Configure SRA Toolkit
if [ ! -d "../sra_toolkit" ]; then
  wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
  tar -xzvf sratoolkit.current-ubuntu64.tar.gz
  mv sratoolkit.* ../sra_toolkit
  ../sra_toolkit/bin/vdb-config -i  # Interactively set permissions
fi

# Download samples (500,000 reads each)
for ACC in ERR11468775 ERR11468776 ERR11468777; do
  ../sra_toolkit/bin/fastq-dump \
    --split-files \
    --gzip \
    -X 500000 \
    $ACC
done
cd ..

# --------------------------
# 4. Quality Control
# --------------------------
echo "=== RUNNING QUALITY CONTROL ==="
fastqc raw_data_mapping/*.fastq.gz -o qc_reports/

# Adapter trimming with fastp
wget http://opengene.org/fastp/fastp
chmod a+x fastp

for R1 in raw_data_mapping/*_1.fastq.gz; do
  R2=${R1/_1.fastq.gz/_2.fastq.gz}
  BASE=$(basename $R1 _1.fastq.gz)
  
  ./fastp \
    -i $R1 -I $R2 \
    -o raw_data_mapping/${BASE}_trimmed_1.fastq.gz \
    -O raw_data_mapping/${BASE}_trimmed_2.fastq.gz \
    --detect_adapter_for_pe \
    -h qc_reports/${BASE}_fastp.html
done

# --------------------------
# 5. Alignment with BWA
# --------------------------
echo "=== ALIGNING READS ==="
for R1 in raw_data_mapping/*_trimmed_1.fastq.gz; do
  R2=${R1/_trimmed_1.fastq.gz/_trimmed_2.fastq.gz}
  BASE=$(basename $R1 _trimmed_1.fastq.gz)

  bwa mem -t 4 \
    ref_genome/hg38_chr6_7.fa \
    $R1 $R2 \
    | samtools sort -@ 4 -o aligned_data/${BASE}.bam

  samtools index aligned_data/${BASE}.bam
done

# --------------------------
# 6. GATK Preprocessing
# --------------------------
echo "=== GATK PREPROCESSING ==="
# Download known sites
wget -P ref_genome/ \
  https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

# Run in Docker for version control
docker run -v $(pwd):/data broadinstitute/gatk:latest bash -c "
  cd /data/aligned_data

  for BAM in *.bam; do
    # Add read groups
    gatk AddOrReplaceReadGroups \
      -I \$BAM \
      -O \${BAM%.bam}_withRG.bam \
      -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM \${BAM%.bam}

    # Mark duplicates
    gatk MarkDuplicates \
      -I \${BAM%.bam}_withRG.bam \
      -O \${BAM%.bam}_markedDups.bam \
      -M \${BAM%.bam}_metrics.txt

    # Base recalibration
    gatk BaseRecalibrator \
      -I \${BAM%.bam}_markedDups.bam \
      -R /data/ref_genome/hg38_chr6_7.fa \
      --known-sites /data/ref_genome/Homo_sapiens_assembly38.dbsnp138.vcf \
      -O \${BAM%.bam}_recal_data.table

    gatk ApplyBQSR \
      -R /data/ref_genome/hg38_chr6_7.fa \
      -I \${BAM%.bam}_markedDups.bam \
      --bqsr-recal-file \${BAM%.bam}_recal_data.table \
      -O \${BAM%.bam}_recal.bam
  done
"

# --------------------------
# 7. Variant Calling
# --------------------------
echo "=== VARIANT CALLING ==="
# Somatic variants (Mutect2)
for BAM in aligned_data/*_recal.bam; do
  BASE=$(basename $BAM _recal.bam)
  
  gatk Mutect2 \
    -R ref_genome/hg38_chr6_7.fa \
    -I $BAM \
    -O variants/${BASE}_somatic.vcf.gz
done

# Germline variants (Trios)
gatk HaplotypeCaller \
  -R ref_genome/hg38_chr6_7.fa \
  -I aligned_data/sample1_recal.bam \
  -O variants/sample1_germline.g.vcf.gz \
  -ERC GVCF

# Repeat for sample2 and sample3...

# --------------------------
# 8. Variant Filtering
# --------------------------
echo "=== FILTERING VARIANTS ==="
for VCF in variants/*_somatic.vcf.gz; do
  BASE=$(basename $VCF _somatic.vcf.gz)
  
  gatk FilterMutectCalls \
    -R ref_genome/hg38_chr6_7.fa \
    -V $VCF \
    -O variants/${BASE}_filtered.vcf.gz

  bcftools filter \
    -i 'QUAL>20' \
    variants/${BASE}_filtered.vcf.gz \
    > variants/${BASE}_highconf.vcf
done

echo "=== WORKFLOW COMPLETED ==="
