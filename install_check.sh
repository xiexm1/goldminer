#!/bin/bash
# GoldMiner Installation Check Script for Unix/Linux/macOS
# This script checks if all dependencies are installed

echo "====================================="
echo "GoldMiner Installation Check"
echo "====================================="
echo ""

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to check command existence
check_command() {
    if command -v "$1" &> /dev/null; then
        echo -e "${GREEN}✓${NC} $1 is installed"
        return 0
    else
        echo -e "${RED}✗${NC} $1 is NOT installed"
        return 1
    fi
}

# Function to check R packages
check_r_package() {
    if Rscript -e "if (!require('$1', quietly=TRUE)) quit(status=1)" &> /dev/null; then
        echo -e "${GREEN}✓${NC} R package '$1' is installed"
        return 0
    else
        echo -e "${RED}✗${NC} R package '$1' is NOT installed"
        return 1
    fi
}

echo "[*] Checking R installation..."
if check_command "R"; then
    R --version | head -n 1
fi
echo ""

echo "[*] Checking R packages..."
r_packages=("igraph" "reshape2" "dplyr" "data.table" "optparse" "ggplot2" "Cairo")
missing_r_packages=()

for pkg in "${r_packages[@]}"; do
    if ! check_r_package "$pkg"; then
        missing_r_packages+=("$pkg")
    fi
done
echo ""

# Generate install command for missing R packages
if [ ${#missing_r_packages[@]} -gt 0 ]; then
    echo -e "${YELLOW}To install missing R packages, run:${NC}"
    echo "Rscript -e \"install.packages(c('${missing_r_packages[@]}'), repos='https://cloud.r-project.org/')\""
    echo ""
fi

echo "[*] Checking Python installation..."
if check_command "python3"; then
    python3 --version
elif check_command "python"; then
    python --version
else
    echo -e "${RED}✗${NC} Python is NOT installed"
fi
echo ""

echo "[*] Checking MCL installation..."
if check_command "mcl"; then
    mcl --version 2>&1 | head -n 1
fi
echo ""

echo "[*] Checking Git installation..."
check_command "git"
echo ""

echo "====================================="
echo "Installation Summary"
echo "====================================="
echo ""
echo "If any dependencies are missing:"
echo "1. For R and Python: Install from their official websites"
echo "2. For MCL on Linux/macOS: Run 'sudo apt install mcl' or 'brew install mcl'"
echo "3. For R packages: Use the install command shown above"
echo "4. For genetribe: git clone https://github.com/chenym1/genetribe.git"
echo ""
echo "To install MCL on Ubuntu/Debian:"
echo "  sudo apt update && sudo apt install mcl"
echo ""
echo "To install MCL on macOS with Homebrew:"
echo "  brew install mcl"
echo ""
