# Frozen v1 analysis behavior

## Scope

The workflow analyzes one published experiment (`study_id`) at a time. Technical
runs are grouped into logical libraries; matched Ribo-seq and RNA-seq libraries
share a `bio_sample_id`. Cross-study joint normalization, isoform-specific TE,
multi-mapper allocation, novel ORFs, uORFs, readthrough, disomes, and codon-level
pausing are outside v1. The supplied Salmon decoy-aware index is retained for a
future optional transcript-level branch; v1 TE deliberately uses genome-aligned,
gene-level RNA CDS counts so its feature definition matches Ribo CDS counts.

## Project isolation

All mutable state is project-owned. `directories.project_root` is the Snakemake
working directory and therefore owns `.snakemake/` and temporary files;
`directories.generated` owns reference/container validation records and derived
gene/CDS/exon annotations. Pipeline source, profiles, and immutable containers may
be shared, but no mutable symlink or generated output is stored in the pipeline
checkout. This permits sequential or concurrent projects without cross-project
state reuse.

## Annotation

FlyBase r6.68 FBgn identifiers are primary keys. Gene symbols are display-only
annotations and never keys. Gene span length is calculated from GTF coordinates as
`end - start + 1`. Exon and CDS lengths are union lengths across all transcripts.
Gene-level counting uses the union of annotated CDS intervals. A deterministic
representative transcript is used only for frame, offset, and metagene QC.

## Mapping and counting

Ribo-seq is adapter-trimmed, screened against the contaminant index, and aligned
end-to-end to the short-read STAR index. RNA-seq is aligned to the 150-nt STAR
index. Primary counting retains `NH:i:1` alignments only. Cross-gene ambiguous
features are excluded.

## Quality review

Automated metrics do not silently remove biological libraries. The QC stage writes
metrics and a decision template. Formal analysis requires decisions with reviewer,
date, selected RPF lengths, and offset provenance. A changed QC fingerprint
invalidates old decisions.

If an automatically inferred adapter is rejected, edit `samples.tsv` to provide the
explicit adapter and rerun `qc`; do not mark a replacement against counts generated
with the old adapter. The resulting changed QC fingerprint invalidates the old
decision table.

## Statistics

RNA exon/gene counts and Ribo CDS P-site counts are modeled separately with
DESeq2. Matched DTE uses combined RNA CDS and Ribo CDS counts with a paired design:

`~ bio_sample_id + assay + assay:condition`

Because condition is nested in `bio_sample_id`, the implementation uses the
full-rank equivalent `~ bio_sample_id + assay + dte_interaction`, where
`dte_interaction` is 1 only for Ribo libraries in the test condition. Its
coefficient is the assay-by-condition interaction (DTE). FDR uses
Benjamini-Hochberg correction. Direction/effect
summaries use shrunken log2 fold changes when available. One replicate per condition
produces descriptive results without p-values or significance classes.
