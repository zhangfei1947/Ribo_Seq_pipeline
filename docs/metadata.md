# Metadata contract

`samples.tsv` has one row per sequencing run. Multiple rows may share a
`library_id` only when they are technical runs of the same logical library. A
matched Ribo/RNA pair shares `study_id`, `bio_sample_id`, and `condition`, but has
different `library_id` values.

## Required operational fields

| Field | Accepted values / meaning |
|---|---|
| `study_id` | Independent experiment; never jointly normalized across studies. |
| `run_id` | Unique FASTQ/run identifier. |
| `library_id` | Technical runs with this value are concatenated after preprocessing. |
| `bio_sample_id` | Biological replicate and Ribo/RNA pairing key. |
| `assay` | `ribo` or `rna`. |
| `condition` | Experimental group used by contrasts. |
| `fastq1`, `fastq2` | Existing paths; `fastq2` is blank for SE. |
| `layout` | Ribo v1: `SE`; RNA: `SE` or `PE`. |
| `adapter_3p` | Explicit sequence, `none`, or `infer`. Inference requires later review. |
| `strandedness` | `forward`, `reverse`, or RNA-only `unstranded`; never silently inferred. |
| `umi_pattern` | `none` or a valid `umi_tools extract --bc-pattern`. |
| `ribo_protocol_group` | Required for Ribo. Pooled offset evidence is restricted to the same study and this group. |

`batch`, accession/publication, nuclease, cycloheximide, size selection, and notes
are retained as provenance. A batch term is included in single-assay DE only when
it is non-confounded and the design matrix remains full rank.

If `adapter_3p=infer`, Fastp performs the preliminary trim. The reviewer either
accepts the inferred adapter in `qc_decisions.tsv`, or replaces `adapter_3p` in
`samples.tsv` and reruns QC. UMI extraction and UMI-aware BAM deduplication happen
only for an explicit non-`none` pattern; genomic-coordinate-only deduplication is
never performed.

## Contrasts

`contrasts.tsv` declares test and reference conditions and independently enables
Ribo differential occupancy, RNA differential expression, and DTE. Conditions must
exist within the same study. Replication is assessed separately for Ribo, RNA, and
matched Ribo/RNA pairs.

## Reviewed QC decisions

The QC stage creates the template. Analysis accepts it only when every library has
a decision, the fingerprint matches current QC outputs, reviewer/date are present,
and all inclusion flags are `yes` or `no`. Ribo inclusion additionally requires
explicit selected lengths and `length:offset` pairs. Exclusion, caution, or failure
requires a reason. `include_te=yes` must be symmetric for both assays in a matched
pair.

