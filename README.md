# Drosophila Ribo-seq pipeline

Reproducible, two-stage Ribo-seq/RNA-seq workflow for *Drosophila melanogaster*
FlyBase r6.68 on the TAMU HPRC Grace cluster.

The workflow deliberately separates quality review from inference:

1. `qc` validates inputs, preprocesses and aligns libraries, and produces Ribo-seq
   quality metrics plus a reviewable QC decision template.
2. `analysis` requires reviewed QC decisions, produces gene-level count matrices,
   RNA differential expression, differential ribosome occupancy, paired differential
   translation efficiency (DTE), reports, and provenance.

## Scientific defaults

- FBgn `gene_id` is the machine key; `gene_symbol` is included in every exported
  gene-level matrix/result and is used for plot labels.
- Ribo-seq v1 supports single-end 20–35 nt monosome footprints.
- Only uniquely mapped genome alignments are used in primary counting.
- Multi-mappers are reported but are not fractionally or probabilistically assigned.
- Coordinate-only duplicate removal is disabled. UMI-aware deduplication is used
  only when an explicit UMI pattern is supplied.
- Ribo-seq is quantified as whole-CDS, length-specific P-site counts.
- RNA-seq produces whole-CDS counts for DTE and exon/gene counts for RNA DE.
- DTE uses matched Ribo/RNA libraries and an assay-by-condition interaction.
- Default normalization is relative. Global translation changes require valid
  experimental spike-ins.
- At least three biological replicates per condition are considered formal,
  two are exploratory, and unreplicated comparisons are descriptive only.

See [docs/design.md](docs/design.md) for the frozen v1 behavior and
[docs/metadata.md](docs/metadata.md) for project metadata, and
[docs/grace.md](docs/grace.md) for Grace deployment.

## Repository layout

- `workflow/`: Snakemake workflow, scripts, and reusable Python code.
- `config/`: example project configuration and metadata templates.
- `schemas/`: versioned machine-readable input/output contracts.
- `profiles/grace/`: cluster-generic SLURM profile for Grace.
- `containers/`: image definitions/lock template.
- `tests/`: dependency-free unit and smoke tests.

## Development checks

The core validation and reference utilities use only the Python standard library:

```bash
python3 -m unittest discover -s tests -v
Rscript tests/smoke_differential.R
```

Snakemake and scientific tools are expected to run from the locked controller
environment and Singularity images on Grace.

## Start a project

```bash
cp config/config.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
cp config/contrasts.example.tsv config/contrasts.tsv
```

Edit the three copied files, verify the reference and SIF paths, then run:

```bash
snakemake -s workflow/Snakefile --profile profiles/grace qc
```

Review `results/qc/multiqc/multiqc_report.html`, every per-library Ribo QC table,
and `results/qc/replicate_correlations.tsv`. Copy
`results/qc/qc_decisions.template.tsv` to `config/qc_decisions.tsv`; replace every
`REVIEW_REQUIRED`, fill reviewer/date/reason fields, and freeze the accepted RPF
lengths and `length:offset` pairs. Then run:

```bash
snakemake -s workflow/Snakefile --profile profiles/grace analysis
```

If a library has too few reads for a stable sample-specific offset, the template
falls back to the pooled recommendation from the same `study_id` and
`ribo_protocol_group`; it is still marked for human review.

If an inferred adapter is rejected, put the corrected sequence in `samples.tsv`
and rerun `qc`. Do not approve a decision table generated from the wrong trim.

## Outputs and cleanup

Each study is written under `results/STUDY/ANALYSIS_ID/`. Count matrices and all
differential tables start with `gene_id` and `gene_symbol`. Plots use gene symbols;
duplicated symbols are displayed as `symbol [FBgn...]`.

After all studies have a `FINALIZED` marker, checksums, and provenance:

```bash
snakemake -s workflow/Snakefile --profile profiles/grace cleanup_manifest
bin/cleanup_intermediates \
  --manifest results/cleanup/analysis_001.tsv \
  --analysis-id analysis_001
```

The cleanup command is a dry run unless both `--execute` and
`--confirm-analysis-id analysis_001` are supplied. It can delete only exact files
listed under `work/preprocessed`, `work/libraries`, `work/depleted`, and
`work/alignment`; counts, raw FASTQs, references, SIFs, logs, and results are never
eligible.
