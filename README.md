# GoldMiner: Gene Origin, Duplication, and Loss Analysis Tool

GoldMiner is a tool for analyzing and comparing homology gene clusters (HOC) across different genomes. It provides a complete analysis pipeline from gene cluster identification and multi-genome comparison to evolutionary event inference.

## 📚 Documentation Overview

- **Quick Checklist**: [INSTALL_SUMMARY.md](INSTALL_SUMMARY.md) - Step-by-step installation checklist
- **Getting Started**: [Installation Guide](#installation-guide) below
- **Windows Users**: See [WINDOWS_INSTALL.md](WINDOWS_INSTALL.md) for detailed Windows installation with WSL2, Cygwin, or Docker
- **Quick Start**: Run the [Installation Check Scripts](#quick-installation-check)
- **R Dependencies**: Use `Rscript install_dependencies.R` to install all R packages

## Features

- **TDGFinder**: Identify clusters in each genome (N = 1)
- **PairLink**: Connect clusters between pairwise genomes (N = 2)
- **MultiCluster**: Build clusters colinear network in all genomes (N ≥ 3)
- **OdlRecon**: Infer cluster origin, loss, and duplication events

## Requirements

- R ≥ 4.0
- R packages: `igraph`, `reshape2`, `dplyr`, `parallel`, `data.table`, `optparse`, `ggplot2`, `Cairo`
- Python 3.x
- MCL tool (for clustering analysis)
- genetribe (for genome comparison analysis)

## Installation Guide

### 📋 Quick Installation Check

We've provided helper scripts to check your system for required dependencies:

- **Windows**: Run `powershell -ExecutionPolicy Bypass -File install_windows.ps1`
- **Linux/macOS**: Run `bash install_check.sh`

These scripts will check if all required software is installed and provide specific instructions for missing dependencies.

### 📖 Detailed Installation Instructions

Follow our comprehensive installation guide:

- **Windows Users**: See [WINDOWS_INSTALL.md](WINDOWS_INSTALL.md) for Windows-specific instructions with WSL2, Cygwin, and Docker options
- **Linux/macOS Users**: Follow the detailed steps below

### 🚀 Quick Install R Dependencies

To quickly install all R packages required for GoldMiner:

```bash
Rscript install_dependencies.R
```

This will automatically detect and install any missing R packages.

### Step 1: Install R (≥ 4.0)

#### Windows Users

1. **Download R for Windows**
   - Visit the official CRAN website: https://cloud.r-project.org/
   - Click "Download R-4.x.x for Windows" (where 4.x.x is the latest version)
   - Save the `.exe` installer file

2. **Install R**
   - Double-click the downloaded `.exe` file
   - Follow the installation wizard:
     - Select language (English)
     - Click "Next" through the welcome screen
     - Choose installation directory (default: `C:\Program Files\R\R-4.x.x`)
     - Select components to install (keep all default selections)
     - Choose startup options (accept defaults)
     - Select additional tasks:
       - ☑ Create a desktop shortcut
       - ☑ Save version number in registry
       - ☑ Associate R with .RData files
   - Click "Next" and wait for installation to complete

3. **Verify R Installation**
   - Open Command Prompt or PowerShell
   - Type: `R --version`
   - You should see R version 4.0 or higher

#### macOS Users

1. **Download R for macOS**
   - Visit: https://cloud.r-project.org/bin/macosx/
   - Download the latest `.pkg` file (R-4.x.x.pkg)

2. **Install R**
   - Double-click the `.pkg` file
   - Follow the installation wizard
   - Enter your password when prompted

3. **Verify Installation**
   - Open Terminal
   - Type: `R --version`

#### Linux Users (Ubuntu/Debian)

```bash
# Add CRAN repository
sudo apt update
sudo apt install software-properties-common
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Install R
sudo apt update
sudo apt install r-base r-base-dev

# Verify installation
R --version
```

### Step 2: Install R Packages

#### Method 1: Automatic Installation (Recommended)

Open R console and run:

```r
# Install all required packages in one command
install.packages(c("igraph", "reshape2", "dplyr", "data.table", "optparse", "ggplot2", "Cairo"))
```

#### Method 2: Manual Installation

Open R console and install each package individually:

```r
# Install packages one by one
install.packages("igraph")
install.packages("reshape2")
install.packages("dplyr")
install.packages("data.table")
install.packages("optparse")
install.packages("ggplot2")
install.packages("Cairo")
```

#### Method 3: Using RScript (Command Line)

Create a file named `install_r_packages.R`:

```r
#!/usr/bin/env Rscript
packages <- c("igraph", "reshape2", "dplyr", "data.table", "optparse", "ggplot2", "Cairo")
for (p in packages) {
  if (!require(p, character.only = TRUE)) install.packages(p, repos = "https://cloud.r-project.org/")
}
```

Then run: `Rscript install_r_packages.R`

### Step 3: Install Python 3.x

#### Windows Users

1. **Download Python**
   - Visit: https://www.python.org/downloads/windows/
   - Download the latest Python 3.x installer (Windows x86-64 executable installer)

2. **Install Python**
   - Run the installer
   - **IMPORTANT**: Check "Add Python 3.x to PATH" at the bottom of the installer window
   - Click "Install Now"
   - Wait for installation to complete

3. **Verify Python Installation**
   - Open Command Prompt or PowerShell
   - Type: `python --version`
   - Should show: Python 3.x.x

#### macOS Users

Python 3 comes pre-installed on macOS. To verify:

```bash
python3 --version
```

If not installed, install via Homebrew:

```bash
brew install python
```

#### Linux Users (Ubuntu/Debian)

```bash
sudo apt update
sudo apt install python3 python3-pip
python3 --version
```

### Step 4: Install MCL (Markov Cluster Algorithm)

#### Windows Users (WSL Recommended)

The easiest way to use MCL on Windows is through Windows Subsystem for Linux (WSL):

1. **Install WSL2**
   - Open PowerShell as Administrator
   - Run: `wsl --install`
   - Restart your computer when prompted
   - Set up your Linux username and password

2. **Install MCL in WSL**
   ```bash
   sudo apt update
   sudo apt install mcl
   mcl --version  # Verify installation
   ```

#### Alternative for Windows (Native)

MCL doesn't have official Windows support, but you can:

1. Install Cygwin (https://www.cygwin.com/)
2. During installation, search for and install "mcl" package
3. Or use Docker container with MCL installed

#### macOS Users

```bash
brew install mcl
mcl --version
```

#### Linux Users (Ubuntu/Debian)

```bash
sudo apt update
sudo apt install mcl
mcl --version
```

### Step 5: Install genetribe

genetribe is a foundational framework for gene collinearity and orthology analysis.

#### Installation

Visit the official genetribe website for installation instructions:

- Website: https://chenym1.github.io/genetribe/
- GitHub: https://github.com/chenym1/genetribe

Basic installation:

```bash
# Clone the repository
git clone https://github.com/chenym1/genetribe.git
cd genetribe

# Follow the installation instructions in the README
# (genetribe typically requires additional setup)
```

#### genetribe Usage for GoldMiner

Before using GoldMiner, you need to run genetribe core analysis:

```bash
# Example commands for pairwise genome comparison
genetribe core -l genome1 -f genome1
genetribe core -l genome2 -f genome2
genetribe core -l genome1 -f genome2
```

### Step 6: Install GoldMiner

#### Option 1: Using install.sh Script (Unix-like systems)

```bash
# Clone the repository
git clone https://github.com/yourusername/goldminer.git
cd goldminer

# Run the installation script
bash install.sh

# Add to PATH
export PATH=$(pwd):$PATH
```

#### Option 2: Manual Installation (All systems)

```bash
# Clone the repository
git clone https://github.com/yourusername/goldminer.git
cd goldminer

# Create bin directory and make scripts executable
mkdir -p bin

# Copy and make all scripts executable
for script in src/*R src/*py src/*sh; do
    name=$(basename "$script" | cut -d'.' -f1)
    cp "$script" "bin/$name"
    chmod +x "bin/$name"
done

# Add to PATH (Unix/Linux/macOS)
export PATH="$(pwd)/bin:$PATH"

# For Windows PowerShell
$env:PATH = "$(pwd)\bin;" + $env:PATH
```

#### Option 3: Set Permanent PATH (Windows)

1. Open System Properties:
   - Right-click "This PC" → Properties
   - Click "Advanced system settings"
   - Click "Environment Variables"

2. Edit PATH:
   - Under "System variables", find and select "Path"
   - Click "Edit"
   - Click "New" and add the GoldMiner directory path (e.g., `F:\\xiexm\\Downloads\\goldminer`)
   - Click OK on all windows

3. Restart PowerShell or Command Prompt

### Step 7: Verify Installation

Run the following commands to verify all dependencies are installed:

```bash
# Check R
R --version
# Expected: R version 4.x.x

# Check Python
python --version
# Expected: Python 3.x.x

# Check MCL (if using Unix-like system or WSL)
mcl --version
# Expected: mcl 12.x.x

# Check GoldMiner
goldminer --help
# Should show GoldMiner usage information
```

### Troubleshooting

#### Common Issues

1. **"R not found" error**
   - Make sure R is added to PATH
   - Windows: Reinstall R and check "Add R to PATH"

2. **"python not found" error**
   - Make sure Python is added to PATH
   - Windows: Run Python installer again and select "Add Python to PATH"

3. **"mcl not found" error (Windows)**
   - Install WSL and use Linux environment
   - Or use Docker container with MCL

4. **R package installation fails**
   - Update R to the latest version
   - Install Rtools for Windows: https://cran.r-project.org/bin/windows/Rtools/
   - Run R as Administrator when installing packages

5. **Permission denied errors**
   - On Windows: Run PowerShell as Administrator
   - On Unix/Linux: Use `sudo` for system-wide installations

#### Getting Help

If you encounter issues during installation:

1. Check the [GoldMiner GitHub Issues page](https://github.com/xiexm1/goldminer/issues)
2. Contact: Xiaoming Xie (xiexm@cau.edu.cn)
3. Include the following information in your help request:
   - Operating system and version
   - R version (`R --version`)
   - Python version (`python --version`)
   - Exact error message
   - Steps to reproduce the issue

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
