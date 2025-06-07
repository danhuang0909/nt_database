# ntRNA Junction Database of Non-productive Transcript RNA (nt-RNA)

## Overview
This repository contains a comprehensive database of non-productive transcript RNA (nt-RNA) signature junctions identified from The Cancer Genome Atlas (TCGA) pan-cancer analysis. These data represent the first systematic characterization of unproductive splicing events (nt-RNA) across multiple cancer types, providing a valuable resource for understanding alternative splicing dysregulation in cancer. By focusing only on junctional reads, our method circumvents the shortcomings of non-directional mRNA sequencing in the TCGA database, enabling a unique whole-genome survey of nt-RNA.

## Background
Non-productive transcripts (nt-RNAs) are isoform transcripts produced from protein-coding genes that contain frame-shifts or premature termination codons due to alternative splicing. These transcripts can be produced by various alternative splicing mechanisms (e.g., exon-skipping, exon extension) and cannot be translated into functional proteins. They may undergo nonsense-mediated decay (NMD). Our analysis of a large number of cancer transcriptomes (RNAseq) across 13 cancer types from TCGA revealed that nt-RNA expression is widespread, with most protein-coding genes producing nt-RNAs at levels typically representing up to 10% or more of the steady-state amount of gene transcripts.

## Data Content
### Known nt-RNA Signature Junctions
- Curated catalog of nt-RNA signature junctions from NCBI, Ensembl, and NONCODE databases.
- Expression profiles across 13 cancer types and normal tissues.
- Quality-controlled junction reads with stringent filtering criteria.

### Novel Signature Junctions
- Previously unannotated junctions supported by high-quality reads.
- Frame-shift analysis relative to known protein-coding isoforms.
- Cancer-specific alternative splicing events.

### Expression Quantification
Expression values are provided as percentiles (50th, 75th, 90th, 95th, and maximum) of raw junction read counts across discovery datasets for each cancer type, enabling comparative analysis of junction expression patterns.

## Data Visualization
### UCSC Genome Browser Integration
Interactive visualization of all nt-RNA signature junctions is available through custom tracks on the UCSC Genome Browser:
- **UCSC hg38 version**: [Access hg38 tracks](https://genome.ucsc.edu/s/dandan_0909/hg38_5_26)
- **UCSC hg37 version**: [Access hg19 tracks](https://genome.ucsc.edu/s/dandan_0909/hg19_version)

### GitHub Bulk Data Download
To facilitate your own analysis, you can download the data directly via the following link:  
[Data Download](https://github.com/username/repository-name/archive/main.zip)

Each junction displays five expression values (median, 75th, 90th, 95th, and max), allowing researchers to study expression patterns of nt-RNA signature junctions across different cancer types.

## Data Files
### 1. BedGraph Files
**Format**: UCSC bedGraph (4 columns)  
**Column Description**:  
- `chrom`: Chromosome (e.g., chr1, chrX)  
- `start`: Start position (0-based)  
- `end`: End position (0-based, interval: [start, end))  
- `value`: Expression percentile value  

**Position Encoding**: Each junction is represented by 5 consecutive base pairs centered around the junction midpoint, with each position representing a different expression percentile:  
- 50th percentile: `med_pos - 3` to `med_pos - 2`  
- 75th percentile: `med_pos - 2` to `med_pos - 1`  
- 90th percentile: `med_pos - 1` to `med_pos`  
- 95th percentile: `med_pos` to `med_pos + 1`  
- Max value: `med_pos + 1` to `med_pos + 2`  

**Example**:
chr1 10000 10001 15.2 # 50th percentile of junction chr1:9997-10003
chr1 10001 10002 22.3 # 75th percentile
chr1 10002 10003 28.1 # 90th percentile
chr1 10003 10004 31.5 # 95th percentile
chr1 10004 10005 45.7 # Maximum value

### 2. Annotation Files
**Format**: Tab-delimited (6 columns)  
**Column Description**:  
- `chrom`: Chromosome  
- `start`: BedGraph interval start position  
- `end`: BedGraph interval end position  
- `junction_tag`: Junction metadata: `{cancer_type};{chr};{start};{end}_{percentile}`  
- `percentile`: Expression percentile (50th, 75th, 90th, 95th, max)  
- `value`: Numeric expression value  

**Example**:
chr2 88592137 88592138 LUSC;2;88890526;88892790_50th 50th 15.2
chr2 88592138 88592139 LUSC;2;88890526;88892790_75th 75th 22.3

## Usage Instructions
### UCSC Genome Browser Visualization
1. Navigate to the [UCSC Genome Browser](https://genome.ucsc.edu/).
2. Select the appropriate genome assembly (hg19 or hg38).
3. Access our custom tracks through the provided links:
   - [hg38 tracks](https://genome.ucsc.edu/s/dandan_0909/hg38_5_26)
   - [hg19 tracks](https://genome.ucsc.edu/s/dandan_0909/hg19_version)
4. Navigate to your region of interest to view nt-RNA junction expression.

### Local Data Analysis
1. Download bedGraph and annotation files from this repository.
2. Use the annotation file to map bedGraph coordinates to original junction metadata.
3. Cross-reference expression values across different cancer types and percentiles.

### Data Integration
The junction expression data can be integrated with:
- Gene expression profiles
- Clinical data
- Mutation data
- Survival analysis workflows

## Applications
This database enables researchers to:
- Identify cancer-specific nt-RNA expression patterns.
- Discover potential biomarkers based on alternative splicing.
- Study the relationship between unproductive splicing and cancer progression.
- Investigate tissue-specific alternative splicing events.
- Characterize novel transcript isoforms.

## Acknowledgment
The results shown here are in whole or part based upon data generated by the TCGA Research Network: [https://www.cancer.gov/tcga](https://www.cancer.gov/tcga).
