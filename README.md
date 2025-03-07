# GoldMiner: Gene Origin, Duplication, and Loss Analysis Tool

GoldMiner is a tool for analyzing and comparing homology gene clusters (HOC) across different genomes. It provides a complete analysis pipeline from gene cluster identification and multi-genome comparison to evolutionary event inference.

## Features

- **TDGFinder**: Identify clusters in each genome (N = 1)
- **PairLink**: Connect clusters between pairwise genomes (N = 2)
- **MultiCluster**: Build clusters colinear network in all genomes (N ≥ 3)
- **OdlRecon**: Infer cluster origin, loss, and duplication events

## Requirements

- R ≥ 4.0
- R packages: `igraph`, `reshape2`, `dplyr`, `parallel`, `data.table`
- Python 3.x
- MCL tool (for clustering analysis)

## Usage

### Basic Command Format

```bash
goldminer <command> [options]
```

### Available Subcommands

```
[ pipeline ]
     TDGFinder      Step2: Identify clusters in each genome (N = 1)
     PairLink       Step3: Connect clusters between pairwise genomes (N = 2)
     MultiCluster   Step4: Build clusters colinear network in all genomes (N ≥ 3)
     OdlRecon       Step5&6: Clusters origin, loss and duplication inference

[ tool ]
     HocAliPlot     Dotplot of HoC alignment between two genomes
```

### Complete Analysis Pipeline

#### Step 1: Genome Comparison Analysis
This tool relies on **genetribe** for core gene comparison between genomes. genetribe serves as a foundational framework for gene collinearity and orthology analysis. 

For more details about genetribe and its functionalities, visit [genetribe Website](https://chenym1.github.io/genetribe/) and [TGT database](http://wheat.cau.edu.cn/TGT/).
```bash
genetribe core -l rice -f rice
genetribe core -l aet -f aet
genetribe core -l aet -f rice
```

#### Step 2: Identify Gene Clusters
Use TDGFinder to identify gene clusters in each genome:
```bash
goldminer TDGFinder -f data_path -b genome_bed_path -o output_path -s subgenome_info_file -t num_threads -d max_distance
```

Parameters:
- `-f/--one2many`: one2many file path
- `-b/--bed`: bed file path
- `-o/--output`: output file path
- `-t/--threads`: threads to use in parallel
- `-d/--distance`: the max distance of tandem genes
- `-s/--subgenome`: the file of subgenome

#### Step 3: Connect Clusters Between Pairwise Genomes
Use PairLink to connect clusters between pairwise genomes:
```bash
goldminer PairLink -l first_genome_prefix -f second_genome_prefix -o output_path -c cluster_path -i iblocks_path -b blocks_path -t num_threads
```

Parameters:
- `-l/--first`: Prefix name of first file
- `-f/--second`: Prefix name of second file
- `-t/--threads`: threads to use in parallel
- `-o/--out`: link file output path
- `-c/--clu`: cluster file path
- `-i/--iblocks`: iblocks file path
- `-b/--blocks`: blocks file path

#### Step 4: Build Clusters Colinear Network
Use MultiCluster to build clusters colinear network in all genomes:
```bash
goldminer MultiCluster -d output_directory -i input_directory -t genome_list_file -p output_prefix -c cluster_files_directory
```

Parameters:
- `-d/--dir`: directory of output files
- `-i/--input`: directory of input files
- `-c/--cluster`: directory of cluster files
- `-t/--table`: genome list file
- `-p/--prefix`: prefix of output file

#### Step 5&6: Infer Evolutionary Events
Use OdlRecon to infer cluster origin, loss, and duplication events:
```bash
goldminer OdlRecon data_path file_prefix selected_groups
```

Parameters:
- `path`: Path to the data files
- `prefix`: Prefix for the data files
- `selected_groups`: Comma-separated string of selected groups (e.g., 'D,A,B,Thinopyrum,Secale,Hordeum,Avena,Brachypodium,Oryza,Zea')

### Output Files

- `.cludb` files: Gene cluster database files
- `.clu` files: Gene cluster information files
- `.link` files: Mapping files between genome clusters
- `.matrix` files: Wide-format matrices for gene cluster comparison

## Citation

If you use GoldMiner in your research, please cite:

```

```

## License

MIT License

## Contact

Author: Xiaoming Xie  
Email: xiexm@cau.edu.cn
