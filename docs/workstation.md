# Ubuntu workstation deployment

This deployment leaves the Grace/SLURM version unchanged.  It reuses the same
Snakefile, scientific defaults, SIF images, and QC review gate, while replacing
SLURM submission with Snakemake's local executor.

## Chosen execution model

- Snakemake is the process scheduler.  No SLURM service or executor plugin is
  installed.
- The profile exposes 12 of 16 logical CPUs and 28 GiB of 32 GiB RAM.  Linux
  schedules work across the i5-13400 P- and E-cores.
- An `io_slots=4` global resource limits aggregate FASTQ/BAM I/O on the single
  NVMe hot disk.  STAR consumes three slots; lighter readers may overlap when
  CPU and memory also permit it.
- STAR uses 10 threads and a 14 GiB reservation.  Drosophila r6.68 indexes are
  much smaller than the 48 GiB conservative Grace allocation.
- The RTX 3060 is not requested.  STAR, Bowtie2, featureCounts, DESeq2, and the
  P-site code in this workflow are CPU applications.

These are scheduler reservations, not Linux hard memory limits.  Observe the
first full dataset with `htop`; if peak resident memory approaches 30 GiB, lower
`cores` or increase the affected rule's `mem_mb` reservation so Snakemake runs
fewer jobs concurrently.

## Why the default still uses SIF containers

Root access on the workstation makes native package installation possible, but
it does not make it more reproducible.  Apptainer is the maintained
Singularity-compatible runtime and runs the same immutable SIF files used on
Grace.  This keeps STAR, Bowtie2, R, Bioconductor, and system libraries identical
between machines and avoids a second dependency solution.

Miniforge3 manages only the small Snakemake controller environment in
`environment-controller-workstation.yml`.  The local executor is built into
Snakemake, so no SLURM executor plugin is included.

If containers become impossible, the fallback is Miniforge3 with Conda's
libmamba solver plus per-rule environment files and linux-64 lockfiles.  A single
large unlocked Conda environment is not an equivalent production replacement;
the Snakefile would also need `conda:` directives before enabling that mode.

## Storage layout

Use the NVMe as the only analysis hot tier: repository, staged raw FASTQs,
references, SIF images, `work/`, `results/`, `.snakemake/`, and temporary sorting
files all live there while a project is active.  The HDD and NAS are cold tiers
for source copies, completed projects, and checksummed archives.

Do not read active inputs directly from CIFS and do not place the project root or
`work/` there.  CIFS metadata latency, transient reconnects, and file-lock
semantics can make a local Snakemake run slow or appear incomplete.  Stage the
current project's raw inputs from HDD/NAS to SSD before starting Snakemake.

Because the 1 TB SSD also contains Ubuntu, keep at least 200 GiB free before a
new run and avoid exceeding roughly 80% sustained usage.  Stage one project (or
one manageable batch of libraries) at a time.  Archive finalized outputs and run
the existing safe cleanup workflow before staging the next large project.

## One-time setup

Install the non-setuid Apptainer package using the current official Ubuntu PPA
instructions.  Install Miniforge3, then create the controller environment:

```bash
conda env create -f environment-controller-workstation.yml
conda activate riboseq-snakemake-workstation
```

Copy or build the exact SIF images.  Building locally requires a working
fakeroot configuration:

```bash
bin/build_containers_workstation
```

Copy `containers/lock.example.tsv` to `containers/lock.tsv`, replace both digest
placeholders with the printed SHA-256 values, and keep the SIFs immutable for a
production analysis.

## Project setup and execution

```bash
cp config/config.workstation.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
cp config/contrasts.example.tsv config/contrasts.tsv
```

Edit `/CHANGE_ME` reference paths and metadata paths.  If raw data or references
are mounted outside the home/project tree, expose their mount roots to Apptainer,
for example:

```bash
export APPTAINER_BINDPATH=/data,/reference
```

Run preflight and a dry run:

```bash
bin/preflight_workstation
bin/run_workstation -n qc
```

Start long runs inside `tmux` so an SSH or Tailscale interruption does not stop
the controlling Snakemake process:

```bash
tmux new -s riboseq
conda activate riboseq-snakemake-workstation
bin/run_workstation qc
```

After reviewing and freezing `config/qc_decisions.tsv` as described in the main
README, run:

```bash
bin/run_workstation analysis
```

The existing cleanup command remains unchanged.  Generate its manifest with the
workstation wrapper, inspect it, and only then execute the guarded deletion:

```bash
bin/run_workstation cleanup_manifest
bin/cleanup_intermediates --manifest results/cleanup/analysis_001.tsv --analysis-id analysis_001
```
