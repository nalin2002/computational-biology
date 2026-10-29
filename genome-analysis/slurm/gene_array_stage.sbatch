#!/usr/bin/env bash
# AI_CODEX: Generated the dependent Slurm array stages.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Expected: STAGE FASTA OUTPUT_DIR TASK_COUNT" >&2
    exit 2
fi

stage="$1"
fasta_path="$2"
output_dir="$3"
tasks="$4"
repo_dir="${SLURM_SUBMIT_DIR:?Submit from the repository root}"
export PYTHONPATH="$repo_dir/src${PYTHONPATH:+:$PYTHONPATH}"
cd "$repo_dir"

case "$stage" in
    index)
        python3 -m module_2_assignment.parallel_gene_finder \
            index "$fasta_path" "$output_dir" --tasks "$tasks"
        ;;
    scan)
        python3 -m module_2_assignment.parallel_gene_finder \
            scan "$output_dir" "${SLURM_ARRAY_TASK_ID:?Missing array task ID}"
        ;;
    prepare)
        python3 -m module_2_assignment.parallel_gene_finder \
            prepare "$output_dir"
        ;;
    filter)
        python3 -m module_2_assignment.parallel_gene_finder \
            filter "$output_dir" "${SLURM_ARRAY_TASK_ID:?Missing array task ID}"
        ;;
    finalize)
        python3 -m module_2_assignment.parallel_gene_finder \
            finalize "$output_dir"
        ;;
    *)
        echo "Unknown stage: $stage" >&2
        exit 2
        ;;
esac
