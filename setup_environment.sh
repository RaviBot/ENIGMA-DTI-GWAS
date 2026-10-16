#!/bin/bash
###############################################################################
# setup_environment.sh
#
# Creates a conda environment with ALL dependencies needed for the
# ENIGMA-DTI GWAS pipeline (R, Python, PLINK2, bcftools, etc.)
#
# Usage:
#   bash setup_environment.sh
#
# After setup, activate with:
#   conda activate enigma_dti_gwas
###############################################################################

set -uo pipefail

ENV_NAME="enigma_dti_gwas"

#-- Where to install the conda environment.
#   Option 1 (recommended): Use conda's default env location.
#     Set ENV_PREFIX="" and the env will be created by name ("enigma_dti_gwas").
#     You can then activate it with:  conda activate enigma_dti_gwas
#   Option 2: Install in a specific directory (e.g., on a scratch filesystem).
#     Set ENV_PREFIX to the full path, e.g.:
#     ENV_PREFIX="/scratch/mylab/envs/enigma_dti_gwas"
ENV_PREFIX=""

# If ENV_PREFIX is empty, use conda's default location
if [[ -z "${ENV_PREFIX}" ]]; then
    ENV_PREFIX="$(conda info --base 2>/dev/null)/envs/${ENV_NAME}"
fi
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=================================================================="
echo " ENIGMA-DTI GWAS Environment Setup"
echo " Environment: ${ENV_NAME}"
echo " Location:    ${ENV_PREFIX}"
echo " Started:     $(date)"
echo "=================================================================="
echo ""

# Check conda
if ! command -v conda &>/dev/null; then
    echo "[ERROR] conda not found. Please install Anaconda/Miniconda first."
    exit 1
fi

# Initialize conda for this shell
eval "$(conda shell.bash hook)" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh"

# Check if environment already exists
if [[ -d "${ENV_PREFIX}" ]]; then
    echo "[INFO] Environment already exists at ${ENV_PREFIX}."
    read -p "  Recreate it? (y/n): " answer
    if [[ "${answer}" == "y" ]]; then
        echo "  Removing old environment..."
        conda deactivate 2>/dev/null || true
        conda env remove --prefix "${ENV_PREFIX}" -y 2>&1 | tail -3
    else
        echo "  Keeping existing environment. Activating..."
        conda activate "${ENV_PREFIX}"
        echo "  Done. Run your scripts with: conda activate ${ENV_PREFIX}"
        exit 0
    fi
fi

echo ""
echo "--- Step 1: Creating conda environment (core packages) ---"
echo "  This may take 5-10 minutes..."
echo ""

# Create environment with packages that resolve cleanly together.
# numpy<2.0:  eDTI_outliers uses numpy int types with str.join(), which broke
#             in numpy 2.0 (numpy.int64 no longer inherits from Python int).
# plink2:     installed as a direct binary (Step 3) -- bioconda's plink2
#             package uses different version numbering and is often outdated.
# vcftools:   installed separately (Step 4) to avoid perl dependency conflicts.
# xlsxwriter: goes through pip (Step 5) to avoid channel priority issues.
conda create \
    --prefix "${ENV_PREFIX}" \
    -c conda-forge -c bioconda \
    -y \
    python=3.10 \
    r-base=4.3 \
    pandas \
    "numpy<2.0" \
    openpyxl \
    bcftools \
    samtools \
    htslib \
    tabix \
    git \
    wget \
    unzip \
    2>&1 | tail -20

echo ""
echo "--- Step 2: Activating environment ---"
conda activate "${ENV_PREFIX}"

echo ""
echo "--- Step 3: Installing PLINK2 (direct binary, v2.0 Oct 2026) ---"
# Download the exact binary rather than relying on bioconda's package,
# which may lag behind or use incompatible version numbering.
PLINK2_DIR="${ENV_PREFIX}/bin"
if [[ -f "${PLINK2_DIR}/plink2" ]]; then
    echo "  plink2 already in env bin, checking version..."
    cur_ver=$("${PLINK2_DIR}/plink2" --version 2>/dev/null | head -1 || echo "unknown")
    echo "  Current: ${cur_ver}"
fi

echo "  Downloading plink2_linux_x86_64_20261007.zip..."
TMPDIR_PL=$(mktemp -d)
wget -q "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20261007.zip" \
    -O "${TMPDIR_PL}/plink2.zip" 2>&1 || \
wget -q "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip" \
    -O "${TMPDIR_PL}/plink2.zip" 2>&1 || \
curl -sL "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip" \
    -o "${TMPDIR_PL}/plink2.zip" 2>&1

if [[ -f "${TMPDIR_PL}/plink2.zip" ]]; then
    unzip -o -q "${TMPDIR_PL}/plink2.zip" -d "${TMPDIR_PL}/" 2>&1
    if [[ -f "${TMPDIR_PL}/plink2" ]]; then
        cp "${TMPDIR_PL}/plink2" "${PLINK2_DIR}/plink2"
        chmod +x "${PLINK2_DIR}/plink2"
        echo "  [OK] Installed: $("${PLINK2_DIR}/plink2" --version 2>/dev/null | head -1)"
    else
        echo "  [FAIL] plink2 binary not found in zip"
    fi
    rm -rf "${TMPDIR_PL}"
else
    echo "  [FAIL] Could not download plink2. Install manually from:"
    echo "         https://www.cog-genomics.org/plink/2.0/"
fi

echo ""
echo "--- Step 4: Installing vcftools ---"
# vcftools has perl dependencies that can conflict with the main env solve.
# Install separately so a failure here doesn't block everything else.
conda install --prefix "${ENV_PREFIX}" -c bioconda -c conda-forge -y vcftools \
    2>&1 | tail -10
if command -v vcftools &>/dev/null; then
    echo "  [OK] vcftools installed"
else
    echo "  [WARN] vcftools failed via conda."
    echo "         It may already be available as a system module (try: module load vcftools)"
    echo "         Or install from source: https://github.com/vcftools/vcftools"
    echo "         The pipeline only uses vcftools for optional checks; it is not critical."
fi

echo ""
echo "--- Step 5: Installing pip packages (xlsxwriter, eDTI_outliers) ---"
pip install --no-cache-dir xlsxwriter 2>&1 | tail -3
pip install --no-cache-dir git+https://github.com/lizhaddad/eDTI_additional_QC.git 2>&1 | tail -5

echo ""
echo "--- Step 6: Verifying installation ---"

check_tool() {
    local tool=$1
    if command -v "${tool}" &>/dev/null; then
        local loc=$(command -v "${tool}")
        echo "  [OK]   ${tool} -- ${loc}"
        return 0
    else
        echo "  [FAIL] ${tool} -- not found"
        return 1
    fi
}

FAIL=0
echo ""
echo "  Tools:"
check_tool plink2       || FAIL=1
check_tool R            || FAIL=1
check_tool Rscript      || FAIL=1
check_tool python3      || FAIL=1
check_tool bcftools     || FAIL=1
check_tool vcftools     || echo "  (vcftools is optional -- see Step 4 above)"
check_tool tabix        || FAIL=1
check_tool git          || FAIL=1
check_tool eDTI_outliers || FAIL=1

echo ""
echo "  R version:      $(R --version | head -1)"
echo "  Python version: $(python3 --version)"
echo "  PLINK2 version: $(plink2 --version 2>&1 | head -1)"

# Check Python packages
echo ""
echo "  Python packages:"
python3 -c "import pandas; print(f'  [OK]   pandas {pandas.__version__}')" 2>/dev/null || echo "  [FAIL] pandas"
python3 -c "import numpy; print(f'  [OK]   numpy {numpy.__version__}')" 2>/dev/null || echo "  [FAIL] numpy"
python3 -c "import openpyxl; print(f'  [OK]   openpyxl {openpyxl.__version__}')" 2>/dev/null || echo "  [FAIL] openpyxl"
python3 -c "import xlsxwriter; print(f'  [OK]   xlsxwriter {xlsxwriter.__version__}')" 2>/dev/null || echo "  [FAIL] xlsxwriter"

if [[ ${FAIL} -eq 1 ]]; then
    echo ""
    echo "  [WARNING] Some tools failed to install."
    echo "  You may need to install them manually."
fi

echo ""
echo "=================================================================="
echo " Setup complete!"
echo "=================================================================="
echo ""
echo " To use this environment:"
echo ""
echo "   conda activate ${ENV_NAME}"
echo ""
echo " Then update your run0/run1/run2 scripts:"
echo ""
echo "   Rbin=\$(which R)"
echo ""
echo " The environment includes:"
echo "   - R 4.3 with base packages"
echo "   - Python 3.10 with pandas, numpy, openpyxl, xlsxwriter"
echo "   - PLINK2 v2.0 (Oct 2026 binary from cog-genomics.org)"
echo "   - bcftools, tabix, samtools"
echo "   - vcftools (if install succeeded)"
echo "   - eDTI_outliers (ENIGMA-DTI statistical QC)"
echo ""
echo " All scripts should use tools FROM this environment."
echo " No more module conflicts or GLIBC errors."
echo "=================================================================="
