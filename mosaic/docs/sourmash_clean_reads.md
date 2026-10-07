# Clean-read Sourmash profiles for phage isolates

Enable the optional screen with:

```bash
snakemake --use-conda -p phage_isolates \
  --config input_dir=/path/to/00_RAW_DATA isolates=True metagenome=False sourmash_clean_reads=True \
  -j 32 --rerun-incomplete
```

The corresponding defaults in `config.yaml` are:

```yaml
sourmash_clean_reads: False
sourmash_clean_min_shared_bp: 50000
sourmash_clean_min_reference_fraction: 0.10
```

This flag is independent of the existing trimmed-read `sourmash` option and
the cluster-representative `sourmash_contig_catalogue` screen. It adds inputs
and outputs to `phage_isolates_summary`; other workflow modes do not request
the new branch.

## Reads and pooling

For each sample, `sourmash_sketch_clean_reads` uses all three final-clean files:

```
02_CLEAN_DATA/{sample}_forward_paired_clean.tot.fastq.gz
02_CLEAN_DATA/{sample}_reverse_paired_clean.tot.fastq.gz
02_CLEAN_DATA/{sample}_unpaired_clean.tot.fastq.gz
```

It does not use the 2M subset or normalized reads. DNA k-mer sketches use
`k=31`, `scaled=1000`, with abundance tracking. `sourmash_pool_clean_reads`
merges every sample sketch, including any negative control. No additional
pooled FASTQ is created, and the negative control has no special treatment.
The pooled profile is weighted by each library's actual read contribution,
not an equal-weight average of sample profiles.

The existing `sourmash_gather` and `sourmash_tax` rules accept the new sketch
locations. The existing environment and configured GTDB RS226 RocksDB and
taxonomy download rules are reused. The new locations use standard `sourmash
gather`, which preserves abundance information; the older QC branch keeps
`fastmultigather`. The older QC sketch, gather and kreport
filenames remain unchanged. Each profiling rule has its own benchmark.

## Evidence and database scope

Gather searches return reference genomes with at least the configured shared
base-pair overlap. A returned match is labelled `supported` when both apply:

- `unique_intersect_bp >= sourmash_clean_min_shared_bp`;
- `f_match_orig >= sourmash_clean_min_reference_fraction`.

Other returned matches are labelled `weak`, not discarded. Matches below the
gather threshold are not searched/reported. This is a support screen rather
than a claim of confirmed contamination or species identification.

The report includes accession/name, lineage, family/genus/species, estimated
shared bases and hashes, reference k-mer fraction/percent, abundance-weighted
query fraction, and mean/median matching k-mer abundance. Shared base pairs
are estimated from sampled hashes. Reference k-mer containment is not measured
read-mapping breadth, and k-mer abundance is not sequencing depth.

The configured GTDB database screens bacteria and archaea, not all possible
contaminants. Unmatched reads may come from phages, eukaryotes, PhiX, or taxa
absent from the database. In isolate mode these are all fastp-passed reads:
biological read removal is bypassed. The screen does not remove reads or change
retained/excluded contig decisions. No six-category biological classification
is imposed.

## Outputs

Under `07_ANNOTATION/SOURMASH_CLEAN/`:

```
SAMPLES/{sample}_sourmash.sig.zip
SAMPLES/{sample}_gather_sourmash.csv
SAMPLES/{sample}_gather_sourmash.with-lineages.csv
SAMPLES/{sample}_sourmash.summarized.csv
SAMPLES/{sample}_sourmash.kreport.txt
POOLED/pooled_sourmash.sig.zip
POOLED/pooled_gather_sourmash.csv
POOLED/pooled_gather_sourmash.with-lineages.csv
POOLED/pooled_sourmash.summarized.csv
POOLED/pooled_sourmash.kreport.txt
```

Final report tables are under `FIGURES_AND_TABLES/`, not the raw Sourmash directory:

```text
08_phage_isolates_sourmash_clean_reads.tot.tsv
08_phage_isolates_sourmash_clean_taxonomy.tot.tsv
```

The first table contains reference-level evidence, support status and profile
scope (`sample` or `pooled`); the second contains the Sourmash taxonomic profiles.
There is no requirement that a read-profile match is present in the assembled
contig catalogue.

The existing executed `NOTEBOOKS/08_phage_isolates_summary.tot.ipynb` displays
profile status, reference matches, taxonomic summaries, a pooled genus
barplot, and a per-sample genus heatmap. Plotting code is explicit in the
notebook. Plots show up to 30 genera; all returned references remain in the TSV.
A heatmap star means at least one reference in that genus passed both support
filters. Blank stars do not imply absence, and no special NC comparison is made.

Figures are saved as PNG/SVG under `FIGURES_AND_TABLES` using prefixes
`08_phage_isolates_sourmash_clean_pooled.tot` and
`08_phage_isolates_sourmash_clean_genus.tot`. They are also included in the
existing isolate HTML report. Empty/no-match profiles produce header-only
tables and explanatory placeholder plots.

See the [isolate workflow guide](isolate_purity.md) for inputs, read accounting,
the output index and the Snakemake dependency graphs.
