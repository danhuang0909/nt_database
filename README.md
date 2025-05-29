# [ntRNA Database] Research Content and Data Resources


## UCSC Data Visualization

The data related to this research has been visualized on the UCSC Genome Browser for intuitive viewing and analysis. Click the link below to access the data:

[UCSC hg38 version](https://genome.ucsc.edu/s/dandan_0909/hg38_5_26)

[UCSC hg37 version](https://genome.ucsc.edu/s/dandan_0909/hg19_version)
## GitHub Data Download

To facilitate academic communication and data reuse, the original data and processed results of this research have been uploaded to the GitHub repository. You can download the data directly via the following link for further analysis and research:

[Data Download](https://github.com/danhuang0909/nt_database/tree/main/data)



### 1. BedGraph File (`*.bedgraph`)  
**Format**: UCSC bedGraph (4 columns)  
| Column       | Description                                                                 |  
|--------------|-----------------------------------------------------------------------------|  
| `chrom`      | Chromosome (e.g., `chr1`, `chrX`)                                         |  
| `start`      | Start position (0-based, exclusive)                                        |  
| `end`        | End position (0-based, exclusive; interval: `[start, end)` )               |  
| `value`      | Expression value (numeric, e.g., percentile or max value)                   |  

**Position Calculation**:  
- Each interval is **1 base pair wide** and centered around the junction's midpoint (`med_pos`):  
  - **50th percentile**: `med_pos - 3` to `med_pos - 2`  
  - **75th percentile**: `med_pos - 2` to `med_pos - 1`  
  - **90th percentile**: `med_pos - 1` to `med_pos`  
  - **95th percentile**: `med_pos` to `med_pos + 1`  
  - **Max value**: `med_pos + 1` to `med_pos + 2`  

**Example**:  
```bedgraph
chr1    10000    10001    15.2    # 50th percentile of junction chr1:9997-10003 (med_pos=10000)
chr1    10001    10002    22.3    # 75th percentile

## Contact Information

If you have any questions or wish to discuss potential collaborations, please feel free to contact us at [danhuang2018dana@gmail.com]. We look forward to hearing from you!

Thank you for your interest in our research. We hope these data resources will contribute to your academic pursuits.
    
