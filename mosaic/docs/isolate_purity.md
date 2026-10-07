# Phage isolate recovery and purity

The isolate workflow assesses genome recovery, purity and relatedness. It does not require a geNomad viral-positive call.

Run from `mosaic/`:

```bash
snakemake --use-conda -p phage_isolates \
  --config \
    input_dir=/home/lmf/EMG/NOVOGENE_260921/00_RAW_DATA \
    isolates=True metagenome=False remove_euk=False RNA_enriched=False \
    map_to_RefSeq=True metavr_blast=False \
    sourmash_clean_reads=True sourmash_contig_catalogue=True \
    host_identification_test=False \
  -j 128 -k --rerun-incomplete -n
```

Remove `-n` to execute. RefSeq read mapping and the Sourmash flags remain optional; METAVR BLAST can be enabled later.

## Reads and contigs

With `isolates=True`, biological read removal is bypassed, including configured-contaminant and eukaryotic removal. Kraken classifies the QC-passed reads once for reporting; no post-decontamination Kraken job is required. All fastp-passed paired reads and orphans enter mapping, never the 2M subset. PCR duplicate removal/normalisation is for assembly, not purity mapping.

Existing saved SPAdes assemblies retain their filenames and >=1 kb filtering. Full-read own-assembly mapping supplies contig depth. Retention requires:

```text
depth >= isolate_min_depth
AND (length >= isolate_min_length_bp OR depth >= isolate_short_min_depth)
AND no qualifying host chromosome/host viral-region match
```

Defaults are 5x, 4,000 bp and 10x. Host BLAST requires >=90% identity AND >=90% query coverage, controlled by `isolate_host_min_identity` and `isolate_host_min_query_coverage`. Only matches to the sample's assigned host chromosome or its predicted viral regions drive exclusion. Supported host-prophage-associated genomes are excluded from the retained set but preserved and flagged. Other-host matches remain available in metadata and the all-host BLAST table/heatmap; they never automatically exclude or relabel a contig.

Isolate host and RefSeq/METAVR BLAST outputs append `sstart send bitscore btop` to the existing eleven columns. Query/reference breadth uses merged inclusive coordinate intervals, without double-counting overlapping alignments. Host identity is calculated from the same selected query portions, using the highest-scoring alignment where HSPs overlap and counting mismatches/gaps from BTOP. Each query/reference pair is assessed separately. Historical eleven-column files remain readable with identity explicitly labelled `legacy HSP identity estimate` and unavailable reference breadth left missing; regenerate BLAST for exact traceback-based identity. Other annotation datasets retain their existing BLAST output format.

geNomad, CheckV, terminal repeats and Sourmash remain evidence, not mandatory viral selection gates. Six biological categories are not imposed. RNA assembly, virome-positive filtering, VIBRANT, VirSorter, Pharokka, Phynteny and legacy AAI analysis are not required by the isolate target.

## Priority read accounting

```text
own retained -> masked host chromosome -> host viral regions
             -> own excluded -> other samples' retained
             -> other samples' excluded -> unexplained
```

Only reads without a passing alignment advance. Every stage uses >=95% identity and >=85% aligned read length with deterministic Bowtie2 mapping. The first stage reuses the own-assembly filtered BAM.

Reads, not alignment records or pairs, are counted. Broken pairs become orphans with original mate identity recorded. Single-mate, discordant and non-unique evidence remains explicit. Passing alternative references are retained, but one deterministic placement per read contributes to stage coverage. Both full-QC and entering-stage percentages are reported.

Host FASTAs are discovered under `HOST/`. When host FASTAs exist and `isolates=True`, `host_mapping_file.tsv` is mandatory at the project root or under `HOST/`, with `sample` and `host` columns. Every analysed sample must have one unambiguous assignment to an available host FASTA. Missing assignments, unavailable assigned references and conflicting rows/files stop DAG construction before jobs run. If both metadata locations exist, their assignments must agree. Sample names are not used to guess the host.

Residual host stages always use the assigned host; there is no all-host fallback. Runs without any host FASTAs remain possible and record host origin as not assessed, not as evidence of purity. Only reported viral coordinates are masked in isolate mode, without extra flanks. Whole-contig viral candidates are moved out of the chromosome reference into the host-viral reference.

Cross-sample references reuse pooled exact representatives and membership tables. Self-only representatives are ineligible; representatives shared with other samples are eligible. A cross-sample match is sequence sharing, not proof of its source.

`host_identification_test=True` separately enables full-read mappings against every masked/unmasked host. Breadth/depth/count comparisons are not summed into the purity fractions and never automatically relabel samples. Minimum evidence defaults to 100 reads and 1% masked-host breadth; equal best-host measurements remain ambiguous.

## Shared clustering and outputs

All saved contigs retain the hierarchy `original_contig -> exact_rep -> mosaic_cluster -> cluster_rep`. Existing concatenation, exact MMseqs dereplication, provenance, per-assembly geNomad/CheckV and CheckV ANI calculation are reused, including DNA-only sharing links.

Every clustering command uses the original installed CheckV tool:

```bash
aniclust --fna ... --ani ... --out ... \
  --min_ani 95 --min_tcov 85 --min_qcov 0
```

ANI and target coverage pass independently. The local modified product-based `scripts/aniclust_checkv.py` is not invoked or edited.

| Location | Outputs |
| --- | --- |
| `03_CONTIGS/ISOLATES/` | Retained/excluded FASTAs and shared reference pools |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/` | `all_contig_metadata.tsv`, `cluster_metadata.tsv`, `sample_host_assignments.tsv` |
| `06_MAPPING/ISOLATES/{sample}/` | Stage summaries, read-level evidence, counts, coverage/BAMs and unexplained FASTQs |
| `NOTEBOOKS/` | Executed mapping-statistics, isolate-summary, host-summary and catalogue notebooks |
| `FIGURES_AND_TABLES/` | Tables/HTML and purity, contig-count, relatedness, host-match and depth plots |

Plots are explicit notebook cells using existing MOSAIC conventions. Cluster RPKM plots show maximum member own-assembly raw RPKM, not remapped or summed cluster abundance. Figure limits do not truncate tables.

## Host prophage activity

`host_prophage_activity=True` (default) adds `NOTEBOOKS/08_host_prophage_activity.tot.ipynb` to isolate runs with host references. Set it to `False` to skip this optional report. It reuses full-read `map_to_host` for each sample's assigned host only; it does not enable the all-host identification test. The existing rule produces both masked and unmasked host mappings in the same job; only the unmasked BAM contributes to activity.

The activity depth rule excludes unmapped, secondary and supplementary alignment records, then writes compressed per-base coverage under `06_MAPPING/HOST/ACTIVITY/`. All QC-passed paired and orphan reads are used. The notebook restores uncovered positions to zero, masks scaffold ends and compares each embedded prophage with its own scaffold outside the union of all embedded prophages. Masked positions do not contribute to means or breadth denominators; assessed lengths are included. Display binning does not affect statistics.

Named activity thresholds use [PropagAtE v1.1 defaults](https://github.com/AnantharamanLab/PropagAtE#flag-descriptions): ratio >=2, Cohen's d >=0.70, prophage mean depth >=1x, breadth >=50% at >=1x, scaffold-end mask 150 bp and minimum prophage length 1,000 bp. The effect uses the unweighted pooled variance in the original implementation. MOSAIC retains its 95% identity/85% aligned-read-length policy, rather than the tool's 97% identity default: this is PropagAtE-style analysis, not an identical stock run.

The TSV `08_host_prophage_activity.tot.tsv` contains all assigned-host sample/region combinations, depth, breadth, background support, ratios, effects and individual threshold checks. Core calls are `active`, `dormant`, `ambiguous` and `not_present`. Zero host background, short regions, fully masked regions and whole-contig viral candidates receive explicit `not_assessed_*` calls. `active` means coverage enrichment compatible with replication, not confirmed induction; `dormant` does not prove inactivity. The main heatmap and per-region coverage profiles are saved and displayed inside the notebook. This evidence never changes contig retention or purity accounting. No new environment or external tool installation is required.

The host report includes CheckM, embedded prophages versus whole-contig viral candidates, CheckV and terminal-repeat evidence. The terminal-repeat candidate flag requires explicit geNomad `DTR`/`ITR` topology or the circularity rule's threshold-qualified `terminal_repeat_type` (using `circularity_min_repeat_bp`). `No terminal repeats` and sub-threshold raw overlaps do not qualify. Whole-host-contig repeat evidence is not assigned to an embedded prophage. BACPHLIP is assessed only for CheckV-complete Caudoviricetes genomes; other sequences remain not assessed. Terminal repeats alone do not establish physical circularity, completeness or biological activity.

Review warnings default to retained-contig count !=1, <70% own-retained mapping, >10% host-associated reads, >10% unexplained reads and supported host-prophage exclusions. These warnings do not erase a recovered genome.

Existing project reports are not edited in place: rerun the workflow to regenerate them. Previously product-clustered outputs must be regenerated using `--forcerun vOUTclustering` when migrating. The general `07_Normalise.py.ipynb` and virome-positive filtering remain unchanged; independent original 95/85 clustering is the intentional repository-wide change.
