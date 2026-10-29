# Parallel chromosome 21 gene finder using Slurm arrays

Submit from the repository root on a Slurm login node. The uncompressed FASTA
and output directory must be on a shared filesystem visible to compute nodes.

```bash
bash slurm/submit_gene_array.sh data/NC_000021.fasta
```

To choose a new output directory:

```bash
bash slurm/submit_gene_array.sh data/NC_000021.fasta data/NC_000021_my_run
```

The submitter creates five dependent Slurm jobs:

1. `gene-index`: index the single FASTA record and establish 50 base ranges.
2. `gene-scan`: a 50-element job array (`0-49%10`) scans codons and CpG windows.
3. `gene-prepare`: merge CpG regions globally and index stops by frame.
4. `gene-filter`: a 50-element job array (`0-49%10`) pairs starts with the
   first downstream stop and checks overlap with the merged CpG regions.
5. `gene-finalize`: combine candidates and write results and timing.

The `%10` array concurrency cap means no more than 10 array elements run at a
time. Each element requests one node and one core, so this stays within the
10-node and 48-core-per-node limits. It does not require all 10 nodes to be
available at once. Add your cluster's required account or partition to the
`sbatch` calls in the submitter if needed.

Outputs:

- `gene_candidates.tsv`: start, stop, length in base pairs, and strand.
- `gene_candidates.txt`: the notebook's human-readable candidate format.
- `summary.json`: candidate and merged-island counts, total elapsed time from
  the start of indexing through final output, and each stage's time.
- `logs/`: one Slurm log per array element and reduction job.

The wall time in `summary.json` includes waiting between dependent jobs, but
not the queue wait before the index job begins. For full time from submission,
use Slurm's `sacct` on the printed job IDs.

The workflow requires a cluster `python3` with access to this repository's
`src` directory; the batch file sets `PYTHONPATH` itself. It uses only Python's
standard library and the project's existing CpG and merge functions.

This workflow treats ambiguous bases as `N` when computing the reverse
complement, preserving genomic coordinates. The older `ProteinSyn` method
deletes the three `M`/`R` bases in this chromosome's FASTA, which shifts many
reverse-strand coordinates in the saved notebook output by three bases. The
working-tree `gene_finder.py` also currently has its merge call commented out;
this workflow always merges CpG windows globally as requested.
