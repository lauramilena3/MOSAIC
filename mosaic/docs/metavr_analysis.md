# MetaVR environment, viral-taxonomy and host-metadata report

Enable the existing `metavr_blast: True` setting, or use the wrapper's
`--metavr-blast` option. The annotation targets (`runAnnotation`, `runWorkflow`,
and the default `all` target) then request the MetaVR BLAST search and
the report. With the flag off, the report is not automatically requested.
An explicit report-file target can also run it and any missing prerequisites.

The rule is `METAVR_analysis` in `rules/08_resultsParsing.smk`. It uses the same
`input`, `output`, `params`, `log: notebook=...`, and `notebook:` conventions as
the existing reporting rules. The source is `notebooks/08_METAVR_analysis.py.ipynb`;
it runs in the workflow's notebook environment, using pandas, NumPy, matplotlib,
Biopython and IPython. No new host-prediction environment is installed.

## Inputs and scope

- The shared `filtered_<REPRESENTATIVE_CONTIGS_BASE>.tot.fasta` catalogue.
- Its 11-column MetaVR BLAST output; headerless or exact-header TSV is accepted.
- The existing merged CheckV quality summary.
- `METAVR_main_table.parquet` and `IMG_full_metadata.tsv.gz` inside `METAVR_db`.

The report does not change the catalogue, clustering or filtering. The BLAST search uses the prebuilt MetaVR nucleotide database.
It includes all filtered vOTUs and a separate subset whose **query** CheckV quality
is `Complete` or `High-quality`. This subset is not defined by MetaVR subject quality.
It does not require a separate high-quality FASTA or project-specific sample IDs.

## Download and metadata preparation

With `metavr_blast: True`, `downloadMetaVR` fetches `METAVR_blastdb.tar.zst`,
`METAVR_main_table.parquet` and `IMG_full_metadata.tsv.gz` from
[the NERSC MetaVR directory](https://portal.nersc.gov/cfs/m342/METAVR/),
then checks the extracted BLAST database. Set `METAVR_db` to the desired local
storage path, and `metavr_download_url` to a mirror if needed. The BLAST archive
alone is about 65 GB compressed; plan for substantially more space after extraction.
The `selectMetaVRMetadata` rule reads the Parquet file in row groups and writes
only metadata for subjects in the BLAST hits to a small TSV consumed by the
notebook. It uses the existing `env7.yaml` environment, while the download uses
`env5.yaml`. All workflow flags, rules, report paths and exported MetaVR field names now use
`METAVR`. Isolate-relative searches use separate nucleotide and protein FASTAs
(`METAVR_reference_fasta` and `METAVR_protein_db`), downloaded only when those
isolate targets are requested. These are much larger than the metadata files;
allow substantial disk space and download time.

## Outputs under results_dir

- `NOTEBOOKS/08_METAVR_analysis.tot.ipynb`: executed notebook with tables and plots.
- `FIGURES_AND_TABLES/08_METAVR_analysis.tot.html`: report with embedded figures.
- `FIGURES_AND_TABLES/08_METAVR_analysis.tot/`: detailed hit pairs, per-vOTU and
  sample summaries, environment/viral-taxonomy/host-taxonomy/method counts,
  PNG/SVG figures, selected metadata cache and provenance JSON.

The per-vOTU summary retains no-hit queries and explicitly labels strongest-hit
metadata with `top_hit_` prefixes. The metadata cache is validated against subject
IDs, source paths, source sizes/mtimes and schema before reuse. Metadata are read
in chunks; source databases are never edited. Missing metadata remain visible,
with warnings and `metavr_metadata_found` / `source_metadata_found` indicators.
The report raises clear errors for incompatible schemas or duplicate selected IDs.

## Interpretation

BLAST alignment fragments (HSPs) are aggregated once per query–subject pair.
Query coverage is stratified into `<10%`, `10–<50%`, `50–<85%`, and `85–100%`.
HSP-weighted identity is descriptive, not an ANI estimate or a species assignment.
Statistics describe the matches saved by the existing BLAST search, not every
possible MetaVR relative. Pair counts are not independent samples or read abundance.
Taxonomy tables also report counts of unique MetaVR UViGs within each query subset.

Recorded `host_taxonomy` and `host_taxonomy_method` are preserved together, with
parsed ranks and explicit unassigned categories. `\N` database missing-value
markers are treated as missing. Host taxonomy belongs to the matched MetaVR genome:
it is not a demonstrated host of the query. Physical-source environments are kept
separate from host assignments. No iPHoP, TaxMyPhage or VIRIDIC jobs are launched
by this report; it does not import the additional Bangladesh iPHoP analyses.

The sample panels use **representative assembly origin**, inferred from configured
sample prefixes, not read-mapping presence or abundance. Unmatched identifiers
(including external references) remain under `Unassigned origin`. Every configured
short-read sample appears in both views, with one of these statuses:

- `No vOTUs in this subset`
- `No MetaVR hits`
- `No hits ≥50% coverage`
- `Hits ≥50% coverage`

Zero-hit samples get explanatory labels rather than fabricated percentage bars.
Having no representatives originating from a sample does not prove absence of
viruses: shared vOTUs can be represented by another sample's contig. Rare viral
families/host genera are pooled only in figures; exported tables retain every taxon.

## Checks

Run the synthetic notebook-cell tests in an environment containing the dependencies
above. Set `MOSAIC_SNAKEMAKE` to additionally validate the workflow dry-run:

```bash
MOSAIC_SNAKEMAKE=/path/to/snakemake MPLBACKEND=Agg \
    python -m unittest discover -s mosaic/tests -p test_metavr_analysis.py -v
```

Tests cover headered/headerless input, repeated HSPs across chunks, metadata-cache
reuse, host assignment provenance, unknown/missing metadata, quality subsets,
sample-prefix attribution, all missing-sample statuses and empty inputs.
