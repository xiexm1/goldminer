# GoldMiner Installation Summary

This document provides a quick reference checklist for installing GoldMiner.

## Installation Checklist

### ¡ Step 0: Check Requirements

Run the appropriate installation check script:

**Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy Bypass -File install_windows.ps1
```

**Linux/macOS:**
```bash
bash install_check.sh
```

**Quick Fix for R packages:**
```bash
Rscript install_dependencies.R
```

---

### ¡ Step 1: Install R (e 4.0)

- **Windows**: Download from https://cloud.r-project.org/
- **macOS**: Download from https://cloud.r-project.org/bin/macosx/
- **Linux**: `sudo apt install r-base r-base-dev`

Verify: `R --version`

---

### ¡ Step 2: Install R Packages

Open R console and run:
```r
install.packages(c("igraph", "reshape2", "dplyr", "data.table",
                   "optparse", "ggplot2", "Cairo"))
```

Or use the automated script:
```bash
Rscript install_dependencies.R
```

Verify all packages:
```bash
Rscript -e "print(installed.packages()[, 'Package'])" | grep -E "igraph|reshape2|dplyr|data.table|optparse|ggplot2|Cairo"
```

---

### ¡ Step 3: Install Python 3.x

- **Windows**: Download from https://www.python.org/downloads/
  - CHECK "Add Python to PATH" during installation
- **macOS**: `brew install python`
- **Linux**: `sudo apt install python3 python3-pip`

Verify: `python --version` or `python3 --version`

---

### ¡ Step 4: Install MCL

#### Windows Users (Choose one):
- **Option A (Recommended)**: Install WSL2
  ```powershell
  wsl --install
  # Then in WSL: sudo apt install mcl
  ```
- **Option B**: Use Cygwin with mcl package
- **Option C**: Use Docker (see WINDOWS_INSTALL.md)

#### macOS/Linux:
- **macOS**: `brew install mcl`
- **Linux**: `sudo apt update && sudo apt install mcl`

Verify: `mcl --version`

---

### ¡ Step 5: Install genetribe

```bash
git clone https://github.com/chenym1/genetribe.git
cd genetribe
# Follow installation instructions in README
```

Visit: https://chenym1.github.io/genetribe/

---

### ¡ Step 6: Install GoldMiner

#### Download/Clone:
```bash
git clone https://github.com/yourusername/goldminer.git
cd goldminer
```

#### Install:
```bash
bash install.sh
```

#### Add to PATH (Optional):
```bash
export PATH="$(pwd):$PATH"
# For permanent setup, add to ~/.bashrc or ~/.zshrc
```

Verify: `goldminer --help`

---

## Platform-Specific Instructions

### =Ì Windows Users

**Key Points:**
- MCL requires WSL2, Cygwin, or Docker
- Best option: WSL2 (Windows 10 2004+ or Windows 11)
- See [WINDOWS_INSTALL.md](WINDOWS_INSTALL.md) for detailed WSL2/Cygwin/Docker setup

**Quick WSL2 Setup:**
```powershell
# PowerShell as Administrator
wsl --install
# Restart, then in WSL: sudo apt install r-base python3 mcl
```

### =Ì macOS Users

**Key Points:**
- Homebrew recommended for MCL and Git
- R packages may need Xcode command line tools

```bash
# Install Homebrew if not installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install mcl git
brew install r
```

### =Ì Linux Users (Ubuntu/Debian)

**Key Points:**
- All tools available via package manager
- May need to add CRAN repository for latest R

```bash
# Add CRAN repository for R
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Install all dependencies
sudo apt update
sudo apt install -y r-base python3 python3-pip mcl git
```

---

## Troubleshooting Quick Guide

### =4 "R not found"
- Add R to PATH during installation
- Windows: Re-run R installer, check "Add R to PATH"

### =4 "python not found"
- Add Python to PATH during installation
- Check both `python` and `python3` commands

### =4 "mcl not found" (Windows)
- Install WSL2 and run: `sudo apt install mcl`
- Alternative: Use Cygwin or Docker

### =4 "R package installation fails"
- Update R to latest version
- Install Rtools (Windows)
- Run R as Administrator

### =4 "Permission denied"
- Use sudo for system installations
- Check file/folder permissions

---

## Verification Commands

Run these to verify installation:

```bash
# Check R
R --version

# Check Python
python --version

# Check MCL (Linux/macOS/WSL only)
mcl --version

# Check GoldMiner
goldminer --help

# Check all at once (run from goldminer directory)
Rscript install_dependencies.R
bash install_check.sh  # or install_windows.ps1 on Windows
```

---

## Need Help?

1. Run installation check script first
2. Check [README.md](README.md) Troubleshooting section
3. Check [WINDOWS_INSTALL.md](WINDOWS_INSTALL.md) for Windows-specific issues
4. Contact: Xiaoming Xie <xiexm@cau.edu.cn>

Include in your message:
- Operating system and version
- Output from installation check script
- Specific error message
- Steps to reproduce

---

**Next Step:** Proceed to [Usage](#usage) section in README.md