# GoldMiner Windows Installation Helper Script
# This script helps Windows users verify and install dependencies

Write-Host "=====================================" -ForegroundColor Green
Write-Host "GoldMiner Installation Helper (Windows)" -ForegroundColor Green
Write-Host "=====================================" -ForegroundColor Green
Write-Host ""

# Function to check if command exists
function Test-CommandExists {
    param($Command)
    $null = Get-Command $Command -ErrorAction SilentlyContinue
    return $?
}

# Function to check R packages
function Test-RPackage {
    param($Package)
    $test = Rscript -e "print(require($Package, quietly=T))" 2>$null
    return $test -match "TRUE"
}

# Check R installation
Write-Host "[*] Checking R installation..." -ForegroundColor Cyan
if (Test-CommandExists "R") {
    $rVersion = R --version 2>&1 | Select-String -Pattern "R version" | Select-Object -First 1
    Write-Host "  ✓ R is installed: $rVersion" -ForegroundColor Green
} else {
    Write-Host "  ✗ R is NOT installed" -ForegroundColor Red
    Write-Host "    Please download and install R from: https://cloud.r-project.org/" -ForegroundColor Yellow
    Write-Host "    Then re-run this script." -ForegroundColor Yellow
    Write-Host ""
}

# Check R packages
Write-Host "[*] Checking R packages..." -ForegroundColor Cyan
$rPackages = @("igraph", "reshape2", "dplyr", "data.table", "optparse", "ggplot2", "Cairo")
$rScript = @"
packages <- c("$($rPackages -join '","')")
missing <- packages[!sapply(packages, requireNamespace, quietly=TRUE)]
if(length(missing) > 0) {
  cat("Missing packages:", paste(missing, collapse=", "), "\n")
} else {
  cat("All packages installed\n")
}
"@

$packageStatus = $rScript | Rscript - 2>&1
Write-Host "  $packageStatus" -ForegroundColor Yellow

# Check Python installation
Write-Host "[*] Checking Python installation..." -ForegroundColor Cyan
if (Test-CommandExists "python") {
    $pyVersion = python --version 2>&1
    Write-Host "  ✓ Python is installed: $pyVersion" -ForegroundColor Green
} elseif (Test-CommandExists "python3") {
    $pyVersion = python3 --version 2>&1
    Write-Host "  ✓ Python is installed: $pyVersion" -ForegroundColor Green
} else {
    Write-Host "  ✗ Python is NOT installed" -ForegroundColor Red
    Write-Host "    Please download and install Python from: https://www.python.org/downloads/" -ForegroundColor Yellow
    Write-Host "    Make sure to check 'Add Python to PATH' during installation." -ForegroundColor Yellow
    Write-Host ""
}

# Check WSL for MCL
Write-Host "[*] Checking MCL (via WSL)..." -ForegroundColor Cyan
if (Test-CommandExists "wsl") {
    Write-Host "  ✓ WSL is available" -ForegroundColor Green
    $mclCheck = wsl which mcl 2>&1
    if ($mclCheck -match "mcl") {
        Write-Host "  ✓ MCL is installed in WSL" -ForegroundColor Green
    } else {
        Write-Host "  ! MCL is NOT installed in WSL" -ForegroundColor Yellow
        Write-Host "    To install MCL, run in WSL: sudo apt install mcl" -ForegroundColor Yellow
    }
} else {
    Write-Host "  ! WSL is NOT installed" -ForegroundColor Yellow
    Write-Host "    To use MCL on Windows, install WSL2:" -ForegroundColor Yellow
    Write-Host "    Open PowerShell as Administrator and run: wsl --install" -ForegroundColor Yellow
    Write-Host "    Then restart your computer and install MCL in WSL." -ForegroundColor Yellow
}

# Check Git (for genetribe installation)
Write-Host "[*] Checking Git installation..." -ForegroundColor Cyan
if (Test-CommandExists "git") {
    $gitVersion = git --version 2>&1
    Write-Host "  ✓ Git is installed: $gitVersion" -ForegroundColor Green
} else {
    Write-Host "  ! Git is NOT installed" -ForegroundColor Yellow
    Write-Host "    Git is optional but recommended for installing genetribe." -ForegroundColor Yellow
    Write-Host "    Download from: https://git-scm.com/download/win" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "=====================================" -ForegroundColor Green
Write-Host "Installation Summary" -ForegroundColor Green
Write-Host "=====================================" -ForegroundColor Green
Write-Host ""
Write-Host "To complete the installation:" -ForegroundColor Cyan
Write-Host "1. Install any missing dependencies shown above" -ForegroundColor White
Write-Host "2. For R packages, open R console and run:" -ForegroundColor White
Write-Host "   install.packages(c('igraph','reshape2','dplyr','data.table','optparse','ggplot2','Cairo'))" -ForegroundColor Magenta
Write-Host "3. For MCL on Windows: Install WSL2 and run 'sudo apt install mcl'" -ForegroundColor White
Write-Host "4. To install genetribe: git clone https://github.com/chenym1/genetribe.git" -ForegroundColor White
Write-Host ""
Write-Host "Press any key to exit..." -ForegroundColor Gray
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
