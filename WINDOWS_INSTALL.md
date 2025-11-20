# GoldMiner Windows Installation Guide

This guide provides detailed instructions for installing GoldMiner on Windows 10/11 systems.

## Overview

Windows users face unique challenges when installing GoldMiner, primarily because:
1. **MCL** (Markov Cluster Algorithm) is a Linux/Unix tool without native Windows support
2. Some scripts expect Unix-like shell environments

We provide three installation options for Windows users.

---

## Option 1: WSL2 (Recommended) ⭐

Windows Subsystem for Linux 2 provides a full Linux environment on Windows, making all tools work seamlessly.

### Prerequisites

- Windows 10 version 2004 or higher, or Windows 11
- (For older Windows) Enable WSL2 manually

### Installation Steps

#### Step 1: Install WSL2

1. **Open PowerShell as Administrator**
   - Press `Win + X`
   - Select "Windows PowerShell (Admin)" or "Windows Terminal (Admin)"

2. **Run the WSL installation command**
   ```powershell
   wsl --install
   ```

3. **Restart your computer** when prompted

4. **Set up your Linux distribution**
   - After restart, an Ubuntu terminal will open automatically
   - Create a username and password
   - Remember this password - you'll need it for `sudo` commands

#### Step 2: Install Dependencies in WSL

1. **Open WSL terminal** (Ubuntu)
   - You can type "Ubuntu" in the Start menu
   - Or run: `wsl` in PowerShell

2. **Update package lists**
   ```bash
   sudo apt update
   ```

3. **Install dependencies**
   ```bash
   sudo apt install -y r-base python3 python3-pip mcl git
   ```

4. **Verify installations**
   ```bash
   R --version          # Should show R 4.x.x
   python3 --version    # Should show Python 3.x.x
   mcl --version        # Should show mcl version
   git --version        # Should show git version
   ```

#### Step 3: Install R Packages

In the WSL terminal, run:

```bash
sudo R -e "install.packages(c('igraph','reshape2','dplyr','data.table','optparse','ggplot2','Cairo'), repos='https://cloud.r-project.org/')"
```

#### Step 4: Install GoldMiner in WSL

1. **Clone GoldMiner repository**
   ```bash
   cd ~
   git clone https://github.com/yourusername/goldminer.git
   cd goldminer
   ```

2. **Run GoldMiner install script**
   ```bash
   bash install.sh
   ```

3. **Add to PATH (optional)**
   Add to your `~/.bashrc` file:
   ```bash
   echo 'export PATH="$HOME/goldminer/bin:$PATH"' >> ~/.bashrc
   source ~/.bashrc
   ```

#### Step 5: Setting Up File Access

WSL can access your Windows files at `/mnt/`

- **C: drive** → `/mnt/c/`
- **D: drive** → `/mnt/d/`
- Any other drive → `/mnt/[drive letter]/`

**Example**: Accessing `D:\\mydata` from WSL:
```bash
cd /mnt/d/mydata
```

**Recommended workflow**:
1. Store your data in Windows (e.g., `D:\\goldminer-data`)
2. Access it from WSL at `/mnt/d/goldminer-data`
3. Run GoldMiner commands in WSL

#### How to Run GoldMiner

**From PowerShell or Command Prompt**:
```powershell
wsl goldminer TDGFinder -f /mnt/d/your-data/one2many -b /mnt/d/your-data/bed -o /mnt/d/your-data/output -s /mnt/d/your-data/subgenome
```

**From WSL terminal**:
```bash
goldminer TDGFinder -f /mnt/d/your-data/one2many -b /mnt/d/your-data/bed -o /mnt/d/your-data/output -s /mnt/d/your-data/subgenome
```

---

## Option 2: Cygwin (Alternative)

Cygwin provides Unix-like environment on Windows without full virtualization.

### Installation Steps

#### Step 1: Install Cygwin

1. **Download Cygwin**
   - Visit: https://www.cygwin.com/
   - Download and run `setup-x86_64.exe`

2. **Install Cygwin**
   - Run the installer
   - Choose "Install from Internet"
   - Select Root Directory (default: `C:\cygwin64`)
   - Select Package Directory
   - Choose Internet connection settings
   - Select download site (choose nearby mirror)

3. **Select Packages to Install**
   - In the package selection screen
   - Search for and install these packages:
     - `R` (Math category)
     - `python3` (Python category)
     - `mcl` (Math category)
     - `git` (Devel category)
   - Click "Next" and complete installation

#### Step 2: Install Additional R Packages

Open Cygwin terminal and run:

```bash
R -e "install.packages(c('igraph','reshape2','dplyr','data.table','optparse','ggplot2','Cairo'), repos='https://cloud.r-project.org/')"
```

#### Step 3: Install GoldMiner in Cygwin

1. **Clone and install GoldMiner**
   ```bash
   cd /cygdrive/c/Users/YourUsername  # Or any directory
   git clone https://github.com/yourusername/goldminer.git
   cd goldminer
   bash install.sh
   ```

2. **Access Windows files**
   - Windows paths are accessible under `/cygdrive/` in Cygwin
   - C:\Users → `/cygdrive/c/Users`
   - D:\data → `/cygdrive/d/data`

#### Step 4: Running GoldMiner

**From Cygwin terminal**:
```bash
goldminer TDGFinder -f /cygdrive/d/your-data/one2many -b /cygdrive/d/your-data/bed -o /cygdrive/d/your-data/output -s /cygdrive/d/your-data/subgenome
```

---

## Option 3: Docker (Advanced)

Create a portable, containerized environment with all dependencies.

### Prerequisites

- Docker Desktop for Windows: https://www.docker.com/products/docker-desktop
- WSL2 backend enabled in Docker

### Installation Steps

#### Step 1: Create Dockerfile

Create a `Dockerfile` in your project directory:

```dockerfile
FROM ubuntu:22.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    r-base \
    python3 \
    python3-pip \
    mcl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('igraph','reshape2','dplyr','data.table','optparse','ggplot2','Cairo'), repos='https://cloud.r-project.org/')"

# Install GoldMiner
WORKDIR /opt
RUN git clone https://github.com/yourusername/goldminer.git
WORKDIR /opt/goldminer
RUN bash install.sh

# Set environment variables
ENV PATH="/opt/goldminer/bin:${PATH}"

# Set working directory
WORKDIR /data

CMD ["bash"]
```

#### Step 2: Build Docker Image

From PowerShell in the directory with your Dockerfile:

```powershell
docker build -t goldminer .
```

This may take 5-10 minutes on first build.

#### Step 3: Run GoldMiner with Docker

**Before running**, prepare your data in a folder (e.g., `D:\goldminer-data`)

```powershell
# Run GoldMiner in Docker
docker run -it --rm -v D:\goldminer-data:/data goldminer goldminer TDGFinder -f /data/one2many -b /data/bed -o /data/output -s /data/subgenome

# Or get an interactive shell
docker run -it --rm -v D:\goldminer-data:/data goldminer
```

**Explanation**:
- `-v D:\goldminer-data:/data` → Mounts your Windows folder to `/data` in the container
- `--rm` → Removes container after it exits
- `-it` → Interactive terminal

---

## Troubleshooting Windows Installation

### Problem: "mcl not found" error

**Solution**: You need WSL2 or Cygwin. MCL is Linux-only.

### Problem: File path issues

**In WSL**:
- Use Linux paths: `/mnt/d/mydata/file.txt`
- NOT Windows paths: `D:\\mydata\\file.txt`

**In PowerShell calling WSL**:
```powershell
# Wrong:
wsl goldminer TDGFinder -f "D:\\mydata\\file.txt"

# Correct:
wsl goldminer TDGFinder -f /mnt/d/mydata/file.txt
```

### Problem: Permission denied

**Solution**: Run PowerShell/Cygwin as Administrator, or ensure your user has permission to access the files.

### Problem: R packages fail to install

**Solutions**:
1. Run R as Administrator (Windows)
2. Or install packages from WSL/Cygwin with sudo
3. Update R to latest version
4. Install Rtools from: https://cran.r-project.org/bin/windows/Rtools/

### Problem: Large file processing is slow

**Solution**: Store data on SSD if possible. Consider splitting large datasets.

---

## Performance Comparison

| Method | Pros | Cons | Recommendation |
|--------|------|------|----------------|
| **WSL2** | ✅ Native Linux performance<br>✅ Full tool compatibility<br>✅ Easy file access | ⚠️ Requires Windows 10 2004+<br>⚠️ Slight setup time | ⭐⭐⭐⭐⭐ **Best overall** |
| **Cygwin** | ✅ Works on older Windows<br>✅ Lighter than WSL | ⚠️ Slower performance<br>⚠️ Some compatibility issues | ⭐⭐⭐ Good alternative |
| **Docker** | ✅ Portable<br>✅ Isolated environment<br>✅ Reproducible | ⚠️ Steep learning curve<br>⚠️ More disk space | ⭐⭐⭐⭐ For advanced users |

---

## Quick Reference: Common Commands

### WSL2
```powershell
# Open WSL
wsl

# Install missing packages in WSL
sudo apt update && sudo apt install -y r-base python3 mcl

# Run goldminer on Windows files
goldminer TDGFinder -f /mnt/d/data/input -b /mnt/d/data/bed -o /mnt/d/data/output -s /mnt/d/data/subgenome
```

### Cygwin
```bash
# Open Cygwin terminal
# Install packages via Cygwin installer
# Access Windows files
cd /cygdrive/c/Users/YourUsername/Documents
goldminer TDGFinder -f ./input -b ./bed -o ./output -s ./subgenome
```

### Docker
```powershell
# Build image
docker build -t goldminer .

# Run with data folder
docker run -it --rm -v D:\data:/data goldminer goldminer TDGFinder -f /data/one2many -b /data/bed -o /data/output -s /data/subgenome
```

---

## Getting Help

If you encounter issues:

1. **Run the installation check script**
   ```powershell
   # Windows
   powershell -ExecutionPolicy Bypass -File install_windows.ps1
   ```

2. **Check the main README troubleshooting section**

3. **Common Windows-specific issues**:
   - Path format problems
   - WSL not enabled
   - Permission errors
   - R packages compilation errors

4. **Contact**: Xiaoming Xie (xiexm@cau.edu.cn)
   - Include: Windows version, WSL/Cygwin/Docker choice, error messages

---

## Next Steps

Once installation is complete:

1. **Visit the main [README.md](README.md)** for detailed usage instructions
2. **Explore the [examples/](examples/) directory** for sample data and workflows
3. **Run test analysis** with provided example data to verify installation

Happy analyzing!
