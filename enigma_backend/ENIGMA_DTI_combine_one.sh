#!/bin/bash
###############################################################################
# ENIGMA_DTI_combine_one.sh
#
# Backend: Combine per-chromosome results for ONE phenotype×model.
# Called as a qsub array task — do not run directly.
#
# Arguments:
#   $1 = group list file (one group name per line)
#   $2 = sumstats directory
#   $3 = combined output directory
#
# Uses SGE_TASK_ID to pick which group to process.
###############################################################################

set -uo pipefail

group_list="$1"
sumstats_dir="$2"
combined_dir="$3"

mkdir -p "${combined_dir}"

# Get this task's group name
grp=$(sed -n "${SGE_TASK_ID}p" "${group_list}")

if [[ -z "${grp}" ]]; then
    echo "ERROR: No group for task ${SGE_TASK_ID}"
    exit 1
fi

combined_file="${combined_dir}/${grp}_allchr.txt.gz"

# Skip if already done
if [[ -f "${combined_file}" ]]; then
    echo "[DONE] ${grp} already combined"
    exit 0
fi

# Find chromosome files for this group
chr_files=$(ls "${sumstats_dir}/${grp}_chr"*.glm.linear.gz 2>/dev/null | sort -V)
n_chr=$(echo "${chr_files}" | grep -c "\.gz$" || echo 0)

if [[ ${n_chr} -eq 0 ]]; then
    echo "[SKIP] No files for ${grp}"
    exit 0
fi

echo "Combining: ${grp} (${n_chr} chr files)"

# Combine: cleaned header from first file, data from all
first=1
for f in ${chr_files}; do
    if [[ ${first} -eq 1 ]]; then
        # Clean header: remove # from #CHROM, ? from PROVISIONAL_REF?
        zcat "${f}" | head -1 | sed 's/^#//; s/PROVISIONAL_REF?/PROVISIONAL_REF/g'
        zcat "${f}" | tail -n +2
        first=0
    else
        zcat "${f}" | tail -n +2
    fi
done | gzip > "${combined_file}"

n_lines=$(($(zcat "${combined_file}" | wc -l) - 1))
echo "  -> ${n_lines} variants in ${combined_file}"
