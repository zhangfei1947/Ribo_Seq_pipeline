# Grace deployment

The workflow uses Snakemake's `cluster-generic` executor and Singularity. Scientific
software is not loaded from Grace modules inside rules.

Before running:

1. Build containers outside Grace or pull/convert them in a Grace compute job.
2. Store immutable SIF files on shared scratch and record SHA256 checksums.
3. Create a pinned Snakemake controller environment with the cluster-generic plugin.
4. Edit `profiles/grace/config.yaml` and set the correct `slurm_account` and
   partitions for the allocation.
5. Keep raw data, references, containers, and final results outside cleanup scope.

Grace jobs are single-node, single-task jobs. Rule `threads` map directly to
`--cpus-per-task`. Jobs expected to exceed two hours must use `medium` rather than
the default `short` partition.

The `cluster-generic` executor is intentional: it exercises only Grace's stable
`sbatch`, `sacct`, and `squeue` interface and does not require the Slurm executor
plugin. The status helper prevents a submitted job from being treated as successful
before SLURM has recorded its terminal state.

Suggested controller setup and run sequence:

```bash
module load GCC/13.3.0 Singularity
conda env create -f environment-controller.yml
conda activate riboseq-snakemake

cp config/config.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
cp config/contrasts.example.tsv config/contrasts.tsv
cp containers/lock.example.tsv containers/lock.tsv
# Replace both checksum placeholders with the transferred SIF SHA-256 values.

snakemake -s workflow/Snakefile --profile profiles/grace qc
cp results/qc/qc_decisions.template.tsv config/qc_decisions.tsv
# Review every row and replace all REVIEW_REQUIRED/pending fields.
snakemake -s workflow/Snakefile --profile profiles/grace analysis
snakemake -s workflow/Snakefile --profile profiles/grace cleanup_manifest
```

Container building is deliberately separate. Run `bin/build_containers` on a host
that supports `singularity build --fakeroot`, transfer the SIFs to Grace, then copy
`containers/lock.example.tsv` to `containers/lock.tsv` and record their checksums.
