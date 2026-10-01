# Cenote-Taker3 on the most abundant assembled contigs

Run from the `mosaic/` directory with your usual configuration, adding:

```bash
snakemake --use-conda -p runWorkflow --config \
  input_dir=/path/to/00_RAW_DATA \
  map_to_all_assembled=True run_cenote=True all_assembled_top_n=100 \
  -j 32 --rerun-incomplete
```

The automatic Cenote branch requires both flags. `run_cenote=False` is the
default and leaves the existing top-contig analysis unchanged. RNA enrichment
is independent: when enabled, RNA assemblies can contribute to the same
abundance-ranked selection. These are selected contigs, not necessarily complete
or viral genomes.

## Rules and environment

`downloadCenoteDB` downloads the core HMM, hallmark-taxonomy, RefSeq-taxonomy,
MMseqs CDD and viral-domain databases into `cenote_db`, by default
`db/cenote-taker3`. This directory is ignored by Git. `download_date.txt` records
the successful installation day in UTC. The download is approximately 3 GB
decompressed; optional HHsuite databases are not downloaded.

`annotate_cenote` runs discovery plus annotation on
`03_CONTIGS/ALL_ASSEMBLED/AllAssembled_top_contigs_tot.fasta`. It uses virion and
RdRP hallmark evidence, not forced annotation-only mode. The generic rule also
uses the existing `annotation_fastas` mechanism, but the workflow automatically
requests it only for the top selection.

Both rules use `envs/cenote.yaml`, pinning Cenote-Taker3 3.4.4 and Python 3.11.
The current Pharokka environment cannot share this version without changing its
existing Mash/GSL pins. The separate environment leaves those pins untouched.

Configuration:

```yaml
run_cenote: False
cenote_min_contig_length: 1000
cenote_threads: 8
cenote_mem_mb: 16000
cenote_prune_prophage: False
cenote_db: "db/cenote-taker3"
```

The memory value is a Snakemake scheduling reservation, not a hard memory limit
inside Cenote. Its circularity check requires at least 1,000 nt; the linear
length cutoff follows `cenote_min_contig_length`.

## Results

Results are in `07_ANNOTATION/Cenote_AllAssembled_top_contigs.tot/`, including:

- `mosaic_ct3_virus_summary.tsv`: reported viruses, hallmark counts, taxonomy and
  annotation information, including the original identifier in `input_name`.
- `final_genes_to_contigs_annotation_summary.tsv`: gene annotations using Cenote's
  internal identifiers.
- `contig_name_map.tsv`: internal identifiers mapped to input FASTA headers.
- `sequin_and_genome_maps/`: GenBank files and related annotation outputs, when
  produced by Cenote.
- `mosaic_ct3_prune_summary.tsv`: region coordinates, when pruning produces regions.
- Run arguments, filtered sequences and the tool's own log, when produced.

The outer log is `07_ANNOTATION/Cenote_AllAssembled_top_contigs.tot.log`.
Both rules also have their own benchmarks.

`06_MAPPING/ALL_ASSEMBLED/AllAssembled_top_contigs_metadata_tot.tsv` gains
`Cenote_` columns when this branch is enabled. Results are joined using
`input_name`, not Cenote's renamed contig IDs. Every selected contig remains in
the table. Short contigs are marked `below length cutoff` in `Cenote_assessed`;
unreported classifications and annotations say `not reported`, not “non-viral”.
Multiple reported regions are retained as ` | `-separated values, not silently
reduced to the first hit.

The assembly catalogue, top FASTA, raw RPKM ranking, mapping and normalization
remain unchanged. Rotation is disabled, but Cenote can still trim terminal
repeats in its own output sequences. Metadata distinguishes whole-contig,
processed-contig and region evidence. Do not transfer gene coordinates from a
processed sequence directly to the original assembly without accounting for
that change.

Cenote's sequence-submission exports use its default `DNA` molecule-type label.
This is a supplied export label, not an inferred genome type. In mixed RNA/DNA
analyses, inspect that field before using the GenBank files for submission or
genome-type summaries; RNA assembly origin alone does not establish an RNA
virus genome type. The hallmark/taxonomy results are collected independently.

Temporary working files are cleaned on success and failure. The identifier map
is retained, while `ct_processing` is removed. An empty selection, a selection
below the length cutoff, or a completed search without qualifying hallmark hits
produces header-only annotation tables. Missing outputs from an unsuccessful
tool invocation are not silently treated as a negative result.

Upstream installation, modes and output documentation:
[Cenote-Taker3](https://github.com/mtisza1/Cenote-Taker3).
