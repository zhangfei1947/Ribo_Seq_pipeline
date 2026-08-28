# Container policy

Production runs use prebuilt, immutable SIF files. Do not build or pull images from
Grace login nodes. Build OCI images in CI or on a workstation, transfer the SIF
files to the project, record SHA-256 digests in `containers/lock.tsv`, and make the
image directory read-only for production.

The Grace module system remains responsible for the Snakemake controller,
Singularity/Apptainer, and SLURM commands. Scientific executables run inside SIFs.

