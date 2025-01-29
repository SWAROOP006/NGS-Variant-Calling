# NGS Workflow Documentation
## Overview
- **Goal**: Identify somatic/germline variants in tongue cancer WGS data.
- **Data Source**: [PRJEB62494](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB62494).
- **Tools Used**:  
  ![GATK](https://img.shields.io/badge/GATK-4.0-blue) ![BWA](https://img.shields.io/badge/BWA-0.7.17-green)

## Workflow Diagram
![Workflow](https://mermaid.ink/svg/pako:eNpVjz0LwjAQhP_KkVMoKf5cuoigoJ0c3Eo4k5NcLBeDSP53E0EEL3u7M7M3S1hYwoZ7rHjEhhU9ZtzwxAM3vPDKjR1B6Iq9w8XkAqUoKkFjVJ7qNk9jX5VZ6lq3HtH5g1F8G2H4Zg8Lw-6TfYyBqWcXeO5Gm9mNn9hRj9Yj8aQlHrRjRY8FG3o8cY0b4M6YkF8m4wPvjC7m0P3U)

### Key Steps
1. **Data Download**: 500,000 reads/sample from SRA.
2. **Alignment**: BWA-MEM for chr6/chr7.
3. **Variant Calling**: GATK HaplotypeCaller for germline variants, Mutect2 for somatic variants.
4. **Annotation**: ENSEMBL VEP.

[Explore the Code →](https://github.com/yourusername/ngs-variant-calling/tree/main/workflow/scripts)
