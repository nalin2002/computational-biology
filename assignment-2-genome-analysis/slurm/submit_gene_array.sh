#!/usr/bin/env bash
# AI_CODEX: Generated the five-job array submission workflow.
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: bash slurm/submit_gene_array.sh FASTA [OUTPUT_DIR]" >&2
    exit 2
fi

if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch is unavailable; run this script on the Slurm login node" >&2
    exit 1
fi

repo_dir="$(pwd -P)"
if [[ ! -f "$repo_dir/src/module_2_assignment/parallel_gene_finder.py" ]]; then
    echo "Run from the repository root" >&2
    exit 2
fi

fasta_path="$1"
if [[ "$fasta_path" != /* ]]; then
    fasta_path="$repo_dir/$fasta_path"
fi
if [[ ! -f "$fasta_path" ]]; then
    echo "FASTA not found: $fasta_path" >&2
    exit 2
fi

if [[ $# -eq 2 ]]; then
    output_dir="$2"
    if [[ "$output_dir" != /* ]]; then
        output_dir="$repo_dir/$output_dir"
    fi
    mkdir -p "$output_dir"
else
    output_dir="$(mktemp -d "$repo_dir/data/NC_000021_parallel_run_XXXXXX")"
fi
if [[ -e "$output_dir/index.json" ]]; then
    echo "Output directory already contains a run: $output_dir" >&2
    exit 2
fi
mkdir -p "$output_dir/logs"

tasks=50
array="0-$((tasks - 1))%10"
stage_script="$repo_dir/slurm/gene_array_stage.sbatch"
log_pattern="$output_dir/logs/%x-%A_%a.out"

index_job="$(sbatch --parsable --chdir="$repo_dir" \
    --job-name=gene-index --output="$log_pattern" \
    "$stage_script" index "$fasta_path" "$output_dir" "$tasks")"
scan_job="$(sbatch --parsable --chdir="$repo_dir" \
    --job-name=gene-scan --output="$log_pattern" \
    --dependency="afterok:$index_job" --array="$array" \
    "$stage_script" scan "$fasta_path" "$output_dir" "$tasks")"
prepare_job="$(sbatch --parsable --chdir="$repo_dir" \
    --job-name=gene-prepare --output="$log_pattern" \
    --dependency="afterok:$scan_job" \
    "$stage_script" prepare "$fasta_path" "$output_dir" "$tasks")"
filter_job="$(sbatch --parsable --chdir="$repo_dir" \
    --job-name=gene-filter --output="$log_pattern" \
    --dependency="afterok:$prepare_job" --array="$array" \
    "$stage_script" filter "$fasta_path" "$output_dir" "$tasks")"
finalize_job="$(sbatch --parsable --chdir="$repo_dir" \
    --job-name=gene-finalize --output="$log_pattern" \
    --dependency="afterok:$filter_job" \
    "$stage_script" finalize "$fasta_path" "$output_dir" "$tasks")"

echo "Index job: $index_job"
echo "Scan array: $scan_job ($array)"
echo "Global CpG merge and stop index: $prepare_job"
echo "Filter array: $filter_job ($array)"
echo "Final output job: $finalize_job"
echo "Output directory: $output_dir"
echo "After completion, total elapsed time: $output_dir/summary.json"
