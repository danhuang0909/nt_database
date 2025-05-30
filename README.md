# [ntRNA Database] Research Content and Data Resources


## UCSC Data Visualization

The data related to this research has been visualized on the UCSC Genome Browser for intuitive viewing and analysis. Click the link below to access the data:

[UCSC hg38 version](https://genome.ucsc.edu/s/dandan_0909/hg38_5_26)

[UCSC hg37 version](https://genome.ucsc.edu/s/dandan_0909/hg19_version)
## GitHub Data Download

To facilitate academic communication and data reuse, the original data and processed results of this research have been uploaded to the GitHub repository. You can download the data directly via the following link for further analysis and research:

[Data Download](https://github.com/danhuang0909/nt_database/tree/main/data)



## 1. BedGraph File   
**Format**: UCSC bedGraph (4 columns)  
| Column       | Description                                                                 |  
|--------------|-----------------------------------------------------------------------------|  
| `chrom`      | Chromosome (e.g., `chr1`, `chrX`)                                         |  
| `start`      | Start position (0-based, exclusive)                                        |  
| `end`        | End position (0-based, exclusive; interval: `[start, end)` )               |  
| `value`      | Expression value (numeric, e.g., percentile or max value)                   |  

**Position Calculation**:  
- Each interval is **1 base pair wide** and centered around the junction's midpoint (`med_pos`).med_pos = (junction_start + junction_end) // 2:  
  - **50th percentile**: `med_pos - 3` to `med_pos - 2`  
  - **75th percentile**: `med_pos - 2` to `med_pos - 1`  
  - **90th percentile**: `med_pos - 1` to `med_pos`  
  - **95th percentile**: `med_pos` to `med_pos + 1`  
  - **Max value**: `med_pos + 1` to `med_pos + 2`  

**Example**:  
```bedgraph
chr1    10000    10001    15.2    # 50th percentile of junction chr1:9997-10003 (med_pos=10000)
chr1    10001    10002    22.3    # 75th percentile

```

## 2. Corresponding annotation file

## File Format Overview
This annotation file maps **bedGraph positions** to **original junctions**, enabling cross-referencing between genomic coordinates and junction metadata. It contains **6 columns**:

| Column              | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| `chrom`             | Chromosome (e.g., `chr2`)                                                   |
| `start`             | Start position of the bedGraph interval (0-based, exclusive)               |
| `end`               | End position of the bedGraph interval (0-based, exclusive)                 |
| `junction_tag`      | Original junction metadata formatted as `{cancer_type};{chr};{start};{end}_{percentile}` |
| `percentile`        | Expression percentile (e.g., `50th`, `max`)                                |
| `value`             | Numeric expression value                                                    |


## Mapping BedGraph Positions to Junctions

### Step 1: Identify BedGraph Coordinates
In your bedGraph file, each line contains:
- **Chromosome**: `chr2`
- **Start Position**: `88592137`
- **End Position**: `88592138`


### Step 2: Find Matching Entries in Annotation File
Use the bedGraph's **first three columns** (`chrom`, `start`, `end`) to query the annotation file.  
**Example Query**:
chr2 88592137 88592138 [look up in annotation file]
### Step 3: Extract Junction Metadata
The matching entry in the annotation file will contain:

chr2 88592137 88592138 LUSC;2;88890526;88892790_50th 50th 0.0



### How to Use in UCSC Genome Browser

1. **Upload the BedGraph File**  
   - Go to the [UCSC Genome Browser](https://genome.ucsc.edu/).  
   - Select the appropriate genome assembly (e.g., hg38).  
   - Click **Add Custom Track** → **Paste Data** and upload your bedGraph file.  

2. **Configure Track Settings**  
   - **Track Type**: `bedGraph`  
   - **Visibility**: `full` (to show all percentiles as separate bars)  
   - **Color**: Customize (e.g., `0,99,99` for teal).  

3. **Submit to Visualize**  
   - Click **Submit** to visualize the data.  

## Contact Information

If you have any questions or wish to discuss potential collaborations, please feel free to contact us at [danhuang2018dana@gmail.com]. We look forward to hearing from you!

Thank you for your interest in our research. We hope these data resources will contribute to your academic pursuits.

## Citation
Huang, Dan. "Profiling of Alternative Splicing in Cancer." Doctoral dissertation, The Chinese University of Hong Kong, 2018.
