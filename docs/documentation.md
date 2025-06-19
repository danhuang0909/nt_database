# BAM File Processing Tool Documentation

This tool analyzes sequencing data in BAM format, supporting functions such as filtering high-quality reads (best40), identifying RNA splice junctions, and processing mate junctions. It is designed for transcriptome data analysis scenarios.


## 1. Environment Requirements
- **Python 3.7+**
- **Dependencies**: `pysam`, `numpy`, `re`, `collections` (install via `pip install pysam numpy`)
- **External tools**: 
  - `samtools` (for BAM file operations, must be system-accessible)
  - `fastqc` (for quality control, optional)


## 2. Quick Start

### 2.1 Script Parameters
| Parameter | Type | Description | Required |
|-----------|------|-------------|----------|
| `-h`      | Flag | Display help message and exit | No |
| `-b`      | String (file path) | Input BAM file path (e.g., `input.bam`) | Yes |
| `-o`      | String (directory path) | Root output directory (results stored in `step1_filter_current_junction_read` subdirectory) | Yes |
| `-r`      | String (identifier) | Identifier for the target splice junction (optional but required for junction-specific filtering) | No |


### 2.2 Example Command
```bash
python script.py -b ./data/sample.bam -o ./results -r target_junction

### Explanation  
This command processes `sample.bam`, saves results to `./results/step1_filter_current_junction_read`, and filters reads related to `target_junction`.  


### 3. Core Function Details  

#### 3.1 High-Quality Read Filtering (best40)  
Uses `samtools` to filter reads with mapping quality (MAPQ) ≥ 40, generating `best40.bam`. If the file exists and its size is ≥ 10GB, it skips regeneration; otherwise, it creates the file and runs quality control (FastQC).  


#### 3.2 Splice Junction Identification (`filter_current_junction`)  
Analyzes RNA splicing events using CIGAR strings from BAM files, outputting:  
- `_all_new_junction_read.bed`: Reads supporting novel splice junctions.  
- `_all_current_read_count.bed`: Read counts for known splice junctions.  
- `_statistical_num.txt`: Statistical report (total reads, junction reads, etc.).  


#### 3.3 Mate Junction Processing (`get_mate_junction`)  
Analyzes insert sizes of paired-end reads, generating:  
- `_good_reads_isize_num.txt`: Insert sizes of high-quality paired reads.  
- `_good_reads_quantile.txt`: Quantile statistics of insert sizes (10%, 50%, 90%).  


### 4. Notes  
- Input BAM files must be sorted and indexed (`samtools index`).  
- The output directory will be created automatically if it does not exist. Use absolute paths to avoid parsing issues.  
- The `-r` parameter must match identifiers in known splice junction files (e.g., BED format); otherwise, junctions will be marked as novel.
