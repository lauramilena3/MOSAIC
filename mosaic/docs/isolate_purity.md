# Phage isolate recovery and purity

The isolate workflow assesses genome recovery, purity and relatedness across a collection of sequenced phage isolates. The expected result is one dominant viral genome per sample, but fragmentation, host DNA and possible mixtures remain visible for review. A warning does not automatically invalidate a recovered genome.

Unlike the virome workflow, it does not require a geNomad viral-positive call. All saved SPAdes contigs enter the catalogue; retained/excluded is a coverage-and-assigned-host decision. Novel or unclassified sequences are not removed simply because a viral classifier did not recognise them.

## Quick start

Run from `mosaic/`:

```bash
snakemake --use-conda -p phage_isolates \
  --config \
    input_dir=/path/to/project/00_RAW_DATA \
    isolates=True metagenome=False \
    sourmash_clean_reads=True sourmash_contig_catalogue=True \
  -j 128 -k --rerun-incomplete -n
```

Remove `-n` to execute. This example enables both optional Sourmash screens, with the default 5,000 bp cutoff for cluster representatives. Omit those two flags for the core isolate analysis. `isolates=True` enforces DNA assembly mode and bypasses biological read removal, regardless of the virome defaults for `metagenome`, `RNA_enriched` or `remove_euk`.

Add `map_to_RefSeq=True` for independent RefSeq read mapping/detection, `metavr_blast=True` for METAVR similarity searches, or `host_identification_test=True` for the all-host comparison. These are not prerequisites for the core recovery/purity assessment. RefSeq BLAST of cluster representatives is already part of the core report; it is separate from the optional read-mapping branch.

## Inputs and project layout

```text
project/
├── 00_RAW_DATA/
│   ├── IsolateA_R1.fastq.gz
│   ├── IsolateA_R2.fastq.gz
│   ├── IsolateB_R1.fastq.gz
│   └── IsolateB_R2.fastq.gz
├── HOST/
│   ├── HostA.fasta
│   └── HostB.fasta
└── host_mapping_file.tsv
```

`host_mapping_file.tsv` is tab-separated:

```tsv
sample	host
IsolateA	HostA
IsolateB	HostB
```

| Input | Requirement |
| --- | --- |
| Paired Illumina FASTQs | Required. Default names are `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`; change `forward_tag`/`reverse_tag` if needed. Sample IDs are discovered from the reverse-read filenames. |
| `HOST/{host}.fasta` | Supply the available bacterial host assemblies directly under `HOST/`. The host ID is the filename without `.fasta`. |
| `host_mapping_file.tsv` | Required whenever host FASTAs exist. Every analysed sample needs one assignment to an available host. Several samples can use the same host. The file can instead be under `HOST/`; if both copies exist they must agree. |
| `config.yaml` and Conda | Use the repository configuration and `--use-conda`. The command-line options override the YAML defaults. |
| Databases and external tools | Missing default databases/tools are provided by their existing download rules. First execution needs network access and sufficient disk space; existing installations are reused. Custom database paths must refer to usable installations or outputs supported by a download rule. |

By default, results are written beside `00_RAW_DATA`, in its parent project directory. `results_dir=/path/to/project` can override that location; place `HOST/` and the host assignment table under the selected results directory. Runs without host FASTAs are allowed, but host origin is explicitly unassessed and the isolate decision requires review. There is no inferred host assignment from a sample name.

## Workflow overview

![Phage isolate workflow: reads, assembled contigs and assigned hosts](figures/phage_isolates_overview.png)

[Zoomable overview (SVG)](figures/phage_isolates_overview.svg)

The diagrams show the biological steps, not individual Snakemake jobs. Blue is reads, teal is contigs, purple is hosts, amber is exclusions/review, and grey is unexplained reads/reports. Each major step has its own diagram below. Thresholds shown are the repository defaults; a run with config overrides can use different values. The saved SPAdes filter is **1,000 bp (1 kb)**, not 1,000 kb.

| Step | What happens |
| --- | --- |
| [1. Read QC](#1-read-qc) | fastp trims adapters/poly-G and low-quality sequence. FastQC/MultiQC, read counts, SuperDeduper statistics and Kraken provide QC/contamination evidence. Biological read removal is bypassed. |
| [2. Assembly and contig evidence](#2-assembly-and-contig-evidence) | BBnorm-normalised QC-passed reads enter SPAdes. Save contigs >=1 kb, assign persistent IDs, then collect geNomad, CheckV and terminal-repeat evidence. Assembly reads are not the purity-mapping denominator. |
| [3. Host characterisation](#3-host-characterisation) | CheckM, geNomad, CheckV and eligible BACPHLIP predictions describe assigned hosts and viral regions. Optional activity testing measures coverage enrichment. |
| [4. Retained and excluded contigs](#4-retained-and-excluded-contigs) | Full-read own-assembly mapping supplies depth. Apply the configured length/depth rule and assigned-host chromosome/viral-region BLAST screen; preserve both FASTAs. |
| [5. Dereplication and clustering](#5-dereplication-and-clustering) | Concatenate all saved contigs, exactly dereplicate with MMseqs, then cluster at independent >=95% ANI and >=85% target coverage. Extract representatives and search RefSeq; optional Sourmash/METAVR screens add evidence. |
| [6. Priority read accounting](#6-priority-read-accounting) | Assign each QC-passed read once through six stages. Preserve alternative/single-mate evidence and the final unexplained reads. |
| [7. Reports and review](#7-reports-and-review) | Sample decisions, contig/cluster metadata, host summaries and figures. Optional full-read Sourmash, RefSeq mapping and host identification remain separate comparisons. |

The actual Snakemake graphs and database-download rules are in the [technical appendix](#technical-appendix-snakemake-dependency-graphs). To regenerate the compact diagrams from the repository root:

```bash
conda activate Mosaic_workflow
python mosaic/docs/generate_isolate_diagrams.py
```

This documentation-only generator uses matplotlib and PyYAML from the workflow environment. It reads `config.yaml` for configurable thresholds and saves PNG plus editable-text SVG files under `mosaic/docs/figures/`. It does not call Snakemake, execute analysis, download databases or change project results. Layouts and shared alignment/clustering policies are explicit in the generator, so update those labels if the workflow policy changes.

See [clean-read Sourmash](sourmash_clean_reads.md), [cluster-representative Sourmash](sourmash_contig_catalogue.md), [contig identifiers](contig_identifiers.md) and [optional RefSeq mapping](refseq_mapping.md) for the branch-specific details.

## 1. Read QC

![Read QC: fastp, raw-read duplicate reporting and all QC-passed reads](figures/phage_isolates_01_read_qc.png)

[Zoomable read-QC diagram (SVG)](figures/phage_isolates_01_read_qc.svg)

With `isolates=True`, biological read removal is bypassed, including configured-contaminant and eukaryotic removal. Kraken classifies the QC-passed reads once for reporting; no post-decontamination Kraken job is required. All fastp-passed paired reads and orphans enter mapping, never the 2M subset. BBnorm normalisation is for assembly, not purity mapping. SuperDeduper measures PCR duplicates on the raw pairs for QC reporting; the current assembly path does not consume its deduplicated output, and neither fastp deduplication nor a separate PCR-duplicate filter is enabled.

## 2. Assembly and contig evidence

![Assembly: BBnorm, SPAdes, saved contigs and annotation evidence](figures/phage_isolates_02_assembly_evidence.png)

[Zoomable assembly diagram (SVG)](figures/phage_isolates_02_assembly_evidence.svg)

Existing saved SPAdes assemblies retain their filenames and >=1 kb filtering. Persistent renamed IDs identify contigs throughout the analysis; neighbouring `.ids.tsv` sidecars preserve old assembler names for provenance only. geNomad, CheckV and terminal-repeat checks run on the saved assembly. Their results are evidence, not mandatory viral selection gates. Six biological categories are not imposed. RNA assembly, virome-positive filtering, VIBRANT, VirSorter, Pharokka, Phynteny and legacy AAI analysis are not required by the isolate target.

## 3. Host characterisation

![Assigned hosts: classification, viral regions and optional activity testing](figures/phage_isolates_03_hosts.png)

[Zoomable host diagram (SVG)](figures/phage_isolates_03_hosts.svg)

Host FASTAs require explicit sample-to-host assignments. CheckM describes host assembly quality; geNomad supplies predicted embedded prophages and whole-contig viral candidates, kept as separate evidence scopes. CheckV and terminal-repeat evidence describe those viral regions. BACPHLIP is limited to CheckV-complete Caudoviricetes genomes.

The assigned host supplies the masked chromosome and viral-region references used below. The optional all-host comparison is report-only. [Host activity testing](#host-prophage-activity) uses the assigned unmasked host mapping, compares each embedded prophage with its whole parent-scaffold background, and does not change retention. It is enabled by default for isolate runs with hosts.

## 4. Retained and excluded contigs

![Selection: coverage support and assigned-host matches, with both FASTAs preserved](figures/phage_isolates_04_selection.png)

[Zoomable selection diagram (SVG)](figures/phage_isolates_04_selection.svg)

Full-read own-assembly mapping supplies contig depth. Retention requires:

```text
depth >= isolate_min_depth
AND (length >= isolate_min_length_bp OR depth >= isolate_short_min_depth)
AND no qualifying host chromosome/host viral-region match
```

Defaults are 5x, 4,000 bp and 10x. Host BLAST requires >=90% identity AND >=90% query coverage, controlled by `isolate_host_min_identity` and `isolate_host_min_query_coverage`. Only matches to the sample's assigned host chromosome or its predicted viral regions drive exclusion. Supported host-prophage-associated genomes are excluded from the retained set but preserved and flagged. Other-host matches remain available in metadata and the all-host BLAST table/heatmap; they never automatically exclude or relabel a contig.

Isolate host and RefSeq/METAVR BLAST outputs append `sstart send bitscore btop` to the existing eleven columns. Query/reference breadth uses merged inclusive coordinate intervals, without double-counting overlapping alignments. Host identity is calculated from the same selected query portions, using the highest-scoring alignment where HSPs overlap and counting mismatches/gaps from BTOP. Each query/reference pair is assessed separately. Historical eleven-column files remain readable with identity explicitly labelled `legacy HSP identity estimate` and unavailable reference breadth left missing; regenerate BLAST for exact traceback-based identity. Other annotation datasets retain their existing BLAST output format.

## 5. Dereplication and clustering

![All-contig relatedness: exact representatives, independent 95/85 clustering and provenance](figures/phage_isolates_05_clustering.png)

[Zoomable clustering diagram (SVG)](figures/phage_isolates_05_clustering.svg)

All saved contigs retain the hierarchy `original_contig -> exact_rep -> mosaic_cluster -> cluster_rep`. Existing concatenation, exact MMseqs dereplication, provenance, per-assembly geNomad/CheckV and CheckV ANI calculation are reused, including DNA-only sharing links. Both retained and excluded sequences enter this global catalogue, without a top-N restriction. MMseqs uses 100% identity and 100% target coverage.

Every clustering command uses the original installed CheckV tool:

```bash
aniclust --fna ... --ani ... --out ... \
  --min_ani 95 --min_tcov 85 --min_qcov 0
```

ANI and target coverage pass independently. The local modified product-based `scripts/aniclust_checkv.py` is not invoked or edited. Optional Sourmash screens representatives >=5 kb but never removes shorter clusters from the catalogue.

## 6. Priority read accounting

![Read accounting: six priority mapping stages and final unexplained reads](figures/phage_isolates_06_read_accounting.png)

[Zoomable read-accounting diagram (SVG)](figures/phage_isolates_06_read_accounting.svg)

```text
own retained -> masked host chromosome -> host viral regions
             -> own excluded -> other samples' retained
             -> other samples' excluded -> unexplained
```

Only reads without a passing alignment advance. Every stage uses >=95% identity and >=85% aligned read length with deterministic Bowtie2 mapping. The first stage reuses the own-assembly filtered BAM, which is retained in isolate mode so a revised host/selection does not require repeating that alignment.

Reads, not alignment records or pairs, are counted. Broken pairs become orphans with original mate identity recorded. Single-mate, discordant and non-unique evidence remains explicit. Passing alternative references are retained, but one deterministic placement per read contributes to stage coverage. Both full-QC and entering-stage percentages are reported.

Host FASTAs are discovered under `HOST/`. When host FASTAs exist and `isolates=True`, `host_mapping_file.tsv` is mandatory at the project root or under `HOST/`, with `sample` and `host` columns. Every analysed sample must have one unambiguous assignment to an available host FASTA. Missing assignments, unavailable assigned references and conflicting rows/files stop DAG construction before jobs run. If both metadata locations exist, their assignments must agree. Sample names are not used to guess the host.

Residual host stages always use the assigned host; there is no all-host fallback. Runs without any host FASTAs remain possible and record host origin as not assessed, not as evidence of purity. Only reported viral coordinates are masked in isolate mode, without extra flanks. Whole-contig viral candidates are moved out of the chromosome reference into the host-viral reference.

The startup sample/directory summary is printed only by the main Snakemake process. Contig extraction, reference pooling, residual-read extraction and BACPHLIP input selection use explicit inline Python in `shell:` blocks, avoiding additional Snakemake workers for these lightweight steps. They reuse `env5.yaml` for pandas/Biopython and `env1_mapping.yaml` for residual-read extraction; no new environment is required. Other existing `run:` or shadow rules elsewhere can still start normal Snakemake workers, without repeating the startup summary or restarting the full workflow.

Cross-sample references reuse pooled exact representatives and membership tables. Self-only representatives are ineligible; representatives shared with other samples are eligible. A cross-sample match is sequence sharing, not proof of its source.

`host_identification_test=True` separately enables full-read mappings against every masked/unmasked host. Breadth/depth/count comparisons are not summed into the purity fractions and never automatically relabel samples. Minimum evidence defaults to 100 reads and 1% masked-host breadth; equal best-host measurements remain ambiguous.

## 7. Reports and review

![Reports: PASS, REVIEW, FAIL and separate QC warnings](figures/phage_isolates_07_reports.png)

[Zoomable reports diagram (SVG)](figures/phage_isolates_07_reports.svg)

Start with the sample summary, then use contig/cluster tables and host reports to review the evidence. `REVIEW` flags a concern without automatically rejecting a recovered genome. Read-QC warnings remain separate. The [full decision definitions](#recoverypurity-decisions-and-separate-qc-warnings) distinguish review triggers from unusable data or invalid accounting.

| Location | Outputs |
| --- | --- |
| `03_CONTIGS/ISOLATES/` | Retained/excluded FASTAs and shared reference pools |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/` | `all_contig_metadata.tsv`, `cluster_metadata.tsv`, `sample_host_assignments.tsv` |
| `06_MAPPING/ISOLATES/{sample}/` | Stage summaries, read-level evidence, counts, coverage/BAMs and unexplained FASTQs |
| `NOTEBOOKS/` | Executed mapping-statistics, isolate-summary, host-summary and catalogue notebooks |
| `FIGURES_AND_TABLES/` | Tables/HTML and purity, contig-count, relatedness, host-match and depth plots |

### Which outputs to open first

All paths below are relative to the project results directory. `tot` identifies the full-data analysis, not a 2M mapping subset.

| File or directory | What to use it for |
| --- | --- |
| `FIGURES_AND_TABLES/08_phage_isolates_summary.tot.html` | Start here for the combined figures, sample decisions, evidence and review notes. |
| `FIGURES_AND_TABLES/08_phage_isolates_summary.tot.csv` | Main machine-readable table: one row per sample, with recovery/purity decision, QC, mapping fractions and dominant-contig evidence. |
| `NOTEBOOKS/08_phage_isolates_summary.tot.ipynb` | Executed report with explicit analysis/plotting code and displayed results. |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/all_contig_metadata.tsv` | One row per saved original contig, including retained/excluded status, host evidence, annotation and cluster membership. |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/cluster_metadata.tsv` | One row per global cluster: member/sample counts, annotations, retention counts and optional Sourmash evidence. |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates_cluster_representatives.tot.fasta` | Representative sequence from every global all-contig cluster, not only viral-positive or retained clusters. |
| `03_CONTIGS/ISOLATES/{sample}_{retained,excluded}.tot.fasta` | Review the sequences kept for own-retained mapping versus those excluded by the length/depth/assigned-host policy. Excluded does not mean deleted. |
| `03_CONTIGS/ISOLATES/SELECTION/{sample}.tot.tsv` | Per-sample selection decision and exclusion reason before global reporting. |
| `FIGURES_AND_TABLES/08_phage_isolates_read_accounting.tot.tsv` | Sample/stage read totals, fractions and mapped-reference counts; includes the terminal unexplained group. |
| `06_MAPPING/ISOLATES/{sample}/` | Detailed stage summaries, reference-level counts, read evidence, BAM/coverage files and final unexplained FASTQs. Some intermediates are temporary and may already be cleaned. |
| `FIGURES_AND_TABLES/08_phage_isolates_sample_clusters.tot.tsv` | Sample-by-cluster membership, retained/excluded counts, representative breadth and annotation scope. |
| `FIGURES_AND_TABLES/08_phage_isolates_clusters.tot/` | Paged cluster heatmaps and individual Complete-representative plots; `complete_cluster_figures.tsv` indexes the latter. |
| `FIGURES_AND_TABLES/08_hosts_summary.tot.html` and `.tsv` | Host CheckM quality and numbers/types of predicted viral regions. The executed notebook is `NOTEBOOKS/08_hosts_summary.tot.ipynb`. |
| `FIGURES_AND_TABLES/08_host_viral_regions.tot.tsv` | Host-region coordinates, embedded versus whole-contig scope, CheckV, terminal-repeat and eligible BACPHLIP evidence. |
| `FIGURES_AND_TABLES/08_host_prophage_activity.tot.tsv` | Optional assigned-host sample/region coverage-enrichment calls. The matching notebook and coverage profiles explain the measurements. |
| `FIGURES_AND_TABLES/08_phage_isolates_sourmash_clean_reads.tot.tsv` and `08_phage_isolates_sourmash_clean_taxonomy.tot.tsv` | Optional per-sample and pooled read-profile evidence and taxonomy. Raw Sourmash outputs are under `07_ANNOTATION/SOURMASH_CLEAN/`. |
| `06_MAPPING/REFSEQ_VIRAL/` | Independent RefSeq count/RPKM/breadth/depth matrices and detection tables, only with `map_to_RefSeq=True`. |

### Main sample-table columns

| Column(s) | Meaning |
| --- | --- |
| `sample`, `expected_host`, `host_assessment` | Library ID, explicit host assignment and whether an assigned host reference was assessed. |
| `decision`, `warnings`, `failure_reasons` | Recovery/purity verdict and measured reasons. `FAIL` takes precedence over `REVIEW`. |
| `recovery_state` | `single_retained_contig`, `fragmented_or_mixed`, `no_retained_genome`, or `supported_host_associated_genome_excluded`. These describe recovery, not a definitive biological diagnosis. |
| `qc_status`, `qc_warnings` | Independent read-QC assessment. A PCR-duplicate warning does not automatically change the recovery/purity decision. |
| `original_reads` | All QC-passed paired and orphan reads entering accounting. Each mate is one read, not one pair. |
| `assembled_contigs`, `retained_contigs`, `excluded_contigs` | Saved >=1 kb contig counts before and after the coverage/assigned-host decision. |
| `all_contig_clusters`, `retained_contig_clusters` | Distinct global cluster representatives represented by all saved contigs or retained contigs in this sample. |
| `{stage}_reads`, `{stage}_percent` | Reads assigned at one of stages `01_own_retained` through `06_other_excluded`, and percent of `original_reads`. The six fractions plus `unexplained_percent` sum to 100%. |
| `{stage}_mapped_references` | Distinct reference sequences receiving assigned reads at that stage; the heatmap cell's `n`. Not an additive genome count. |
| `retained_mapping_percent` | Alias of `01_own_retained_percent`, used for the recovery warning. |
| `unexplained_reads`, `unexplained_percent` | Reads with no accepted match after all six stages. These are not automatically contaminants. |
| `dominant_contig`, `dominant_contig_reads`, `dominant_contig_length_bp` | Most read-supported retained contig, its stage-01 assigned reads and sequence length. |
| `dominant_percent_full_qc`, `dominant_percent_retained` | Dominant-contig fraction of all QC-passed reads versus only own-retained assigned reads. The denominators differ. |
| `dominant_contig_genomad_*`, `dominant_contig_checkv_*` | Evidence for that original retained contig, not annotations borrowed from its cluster representative. |
| `retained_genomad_virus_calls`, `retained_genomad_plasmid_calls` | Non-exclusive counts of retained original contigs with each geNomad call. |
| `supported_host_prophage_contigs` | Coverage-supported original contigs excluded because they match their assigned host's predicted viral regions. |
| `host_prophage_activity_assessment`, `active_host_prophage_count`, `active_host_prophage_regions`, `active_host_prophage_details` | Assigned-host activity evidence, including region IDs/depth/breadth. Unassessed counts are missing, not zero. |
| `host_identification_status`, `host_identification_best_host` | Optional all-host comparison; never an automatic reassignment. The best-host column may be absent when no comparison supplies a best hit. |
| `read_accounting_conserved`, `read_accounting_valid` | Total-count conservation versus the stronger check of valid counts, consistent denominators and stage transitions. |

In the per-contig table, `original_id` is the saved, renamed assembly ID before dereplication, not the original SPAdes header. `exact_rep`, `mosaic_cluster` and `cluster_rep` record the clustering hierarchy. `coverage_supported`, `host_match`, `host_prophage_match`, `retained`, `contig_set` and `exclusion_reason` explain selection. `other_host_*` evidence does not drive exclusion. `genomad_*`, `CheckV_*`, terminal-repeat evidence and optional `sourmash_*` fields describe their recorded evidence scopes. The old assembler header can be looked up in the neighbouring `.ids.tsv` sidecar; downstream tools do not need that table.

In the cluster table, `number_exact_representatives`, `number_original_contigs` and `number_samples` use different units. `number_retained`/`number_excluded` count original contigs. `number_genomad_virus`, `number_genomad_plasmid` and their fractions are non-exclusive original-member evidence. `best_CheckV_quality`/`maximum_CheckV_completeness` are the best member measurements, not proof that every member or the representative has that quality. `number_host_matching` and `number_phage_plasmid_candidate` retain the underlying evidence counts. Use the per-contig table to identify the contributing members.

Text annotations without evidence use `not reported`; missing numeric values remain missing. Neither means that the sequence was proven non-viral or host-free.

### Figures and reporting dependencies

Plots are explicit notebook cells using existing MOSAIC conventions. Cluster RPKM plots show maximum member own-assembly raw RPKM, not remapped or summed cluster abundance. Figure limits do not truncate tables.

Catalogue generation uses `notebooks/08_isolate_contig_catalogue.py.ipynb`; the final report uses `notebooks/08_phage_isolates_summary.py.ipynb`. The existing executed notebook paths and all table names are unchanged. Report/plot edits therefore invalidate only the report, not catalogue selection, retained/excluded FASTAs or the six-stage read accounting. Changes to the shared catalogue/selection notebook source or actual selection inputs/parameters still invalidate the relevant biological results. Global annotation tables are no longer mapping dependencies.

Bowtie indexes, non-filtered assembly-mapping intermediates, residual FASTQs and Sourmash sketches/gather intermediates remain temporary. The isolate own-assembly filtered BAM and stage read-evidence tables remain available for reuse. Missing temporary files alone do not require rerunning an up-to-date consumer. Residual FASTQs are reconstructed from the full QC-passed reads and cumulative assigned-read evidence by `extract_isolate_remaining_reads`; their deletion does not require repeating earlier alignments. Report-only refreshes must not force the catalogue, disable provenance checking globally or touch biological outputs.

## Incremental assembly and host updates

`select_isolate_contigs` runs the catalogue notebook in sample-selection mode, using only that sample's assembly, own-assembly coverage, assigned host and assigned-host viral-region BLAST. It writes `03_CONTIGS/ISOLATES/SELECTION/{sample}.tot.tsv`. The global catalogue collects these decisions alongside clustering and annotation results; it no longer determines the per-sample mapping references. The same selection and BLAST parsers are reused in both notebook modes.

Host-stage references and indexes are shared only among samples with the same assigned host, under `03_CONTIGS/ISOLATES/REFERENCES/HOST/{host}/`. Runs without hosts use empty `REFERENCES/UNASSIGNED/` references. Host assignments are tracked as per-sample rule parameters, not as an input dependency on a globally regenerated catalogue table. Per-sample host BLAST reuses `run_BLASTn_host` with `-subject`, avoiding concurrent writes to the same host BLAST database; existing all-host BLAST output paths remain available for reporting.

An improved assembly reruns its own coverage/selection and early accounting stages. A changed host reruns host-specific analysis and early accounting for its assigned samples. Changes to retained/excluded sets or exact dereplication can still require cross-sample stages 05/06 for every sample, because those shared references really change. Their pools depend on the selected FASTAs and existing exact-cluster membership, not global geNomad/CheckV/Sourmash reporting metadata. Recreating deleted residual FASTQs at this boundary uses saved read evidence, not earlier alignments. Global clustering, annotation and combined reports retain their genuine shared dependencies.

Existing projects need a one-time refresh to create per-sample selection/host-reference outputs, save missing own-assembly filtered BAMs and adopt the revised accounting/extraction rules. This can still schedule a large initial run; the narrower invalidation applies once those outputs are current. Do not use `--touch` or disable code/input/parameter tracking to hide this migration. Biological selection thresholds, final table paths and stage-summary paths are unchanged; previous all-host pooled references are not deleted automatically.

## Reading reports and figures

The main sample-level table is `FIGURES_AND_TABLES/08_phage_isolates_summary.tot.csv`: `decision`, `recovery_state`, `warnings`, and `failure_reasons` explain recovery/purity; `qc_status` and `qc_warnings` are independent read-QC assessments. The mapping heatmap `08_phage_isolates_mapping_heatmap.tot.png/.svg` shows all seven priority fractions with PASS/REVIEW/FAIL at the right. Cells show percent of all QC-passed reads and `n`, the number of distinct reference sequences receiving assigned reads at that stage, not all available contigs or alternative placements. Host columns count scaffolds/viral regions, and other-sample columns count pooled exact representatives. Unexplained reads have no reference count. Counts are saved as `{stage}_mapped_references` in the main CSV and `mapped_references` in the read-accounting TSV; they are not additive biological genome counts. Stacked-bar figures remain available. The heatmap is displayed in the notebook and included in the HTML report.

Review notes give the measured percentages and configured limits, including the chromosome/host-viral split and retained-contig/cluster counts. The main CSV also reports the dominant **retained** contig's length, geNomad classification/evidence scope and CheckV quality/completeness. These are its own annotations, never inherited from the cluster representative. If there is no mapped retained contig, text fields say `not reported` and numeric annotations are empty. `retained_genomad_virus_calls` and `retained_genomad_plasmid_calls` count original retained contigs; the counts can overlap and are not additional filtering criteria.

When host activity is enabled, the summary reuses its existing TSV for the assigned host and adds `host_prophage_activity_assessment`, `active_host_prophage_count`, `active_host_prophage_regions` and `active_host_prophage_details` (region ID with mean depth and breadth). Only embedded regions called `active` are counted; whole-contig viral candidates are not treated as embedded prophages. Disabled, unavailable and partially assessed activity is labelled explicitly; an unassessed count is empty, not zero. A host with no predicted embedded prophages has a zero count and the status `no_predicted_embedded_prophages`. Activity is supporting evidence, not confirmed induction, and never changes retention or the sample decision. These fields and the decision notes appear in both the notebook and HTML report without new mapping or activity calculations.

`08_phage_isolates_sample_clusters.tot.tsv` reports each sample/cluster combination, member counts split into retained/excluded, representative origin/length/classification/scores/CheckV quality, own-member quality, host-prophage evidence and assigned retained-read percentages. Representative quality is never inherited as the quality of every member. The cluster-summary CSV additionally flags `mixed_genomad_calls`; this means a cluster has both virus and plasmid calls, not that every member is a phage-plasmid. The existing individual `phage_plasmid_candidate` flag remains separate.

The `08_phage_isolates_cluster_coverage.tot` heatmap shows representative versus sample, with original-contig counts in the cells. Coverage is the union of the existing clustering BLAST subject intervals, using exact-representative membership to recover duplicate origins. It is not summed member length, a completeness estimate or read-mapping breadth. Bold labels mean the representative itself is CheckV Complete; red labels mean at least one member matches its assigned host prophage. All clusters are plotted in pages under `08_phage_isolates_clusters.tot/`; `isolate_plot_max_clusters` sets page size. Every Complete representative also gets an individual retained/excluded member-alignment plot and an entry in `complete_cluster_figures.tsv`. These views do not require the virome filtering or legacy AAI jobs.

## Host prophage activity

`host_prophage_activity=True` (default) adds `NOTEBOOKS/08_host_prophage_activity.tot.ipynb` to isolate runs with host references. Set it to `False` to skip this optional report. It reuses full-read `map_to_host` for each sample's assigned host only; it does not enable the all-host identification test. The existing rule produces both masked and unmasked host mappings in the same job; only the unmasked BAM contributes to activity.

The activity depth rule excludes unmapped, secondary and supplementary alignment records, then writes compressed per-base coverage under `06_MAPPING/HOST/ACTIVITY/`. All QC-passed paired and orphan reads are used. The notebook restores uncovered positions to zero, masks scaffold ends and compares each embedded prophage with its own scaffold outside the union of all embedded prophages. Masked positions do not contribute to means or breadth denominators; assessed lengths are included. Display binning does not affect statistics.

Coverage profiles show the prophage plus 20,000 bp to the left and right by default (`host_prophage_activity_plot_flank_bp`), clipped at scaffold boundaries. The displayed dashed baseline and activity statistics still use the whole parent scaffold background, not just those local flanks and not the entire multi-scaffold host genome. Profile coordinates are recorded as 1-based inclusive `plot_start`/`plot_end` in `coverage_profiles.tsv`.

Named activity thresholds use [PropagAtE v1.1 defaults](https://github.com/AnantharamanLab/PropagAtE#flag-descriptions): ratio >=2, Cohen's d >=0.70, prophage mean depth >=1x, breadth >=50% at >=1x, scaffold-end mask 150 bp and minimum prophage length 1,000 bp. The effect uses the unweighted pooled variance in the original implementation. MOSAIC retains its 95% identity/85% aligned-read-length policy, rather than the tool's 97% identity default: this is PropagAtE-style analysis, not an identical stock run.

The TSV `08_host_prophage_activity.tot.tsv` contains all assigned-host sample/region combinations, depth, breadth, background support, ratios, effects and individual threshold checks. Core calls are `active`, `dormant`, `ambiguous` and `not_present`. Zero host background, short regions, fully masked regions and whole-contig viral candidates receive explicit `not_assessed_*` calls. `active` means coverage enrichment compatible with replication, not confirmed induction; `dormant` does not prove inactivity. The main heatmap and per-region coverage profiles are saved and displayed inside the notebook. This evidence never changes contig retention or purity accounting. No new environment or external tool installation is required.

The host report includes CheckM, embedded prophages versus whole-contig viral candidates, CheckV and terminal-repeat evidence. The terminal-repeat candidate flag requires explicit geNomad `DTR`/`ITR` topology or the circularity rule's threshold-qualified `terminal_repeat_type` (using `circularity_min_repeat_bp`). `No terminal repeats` and sub-threshold raw overlaps do not qualify. Whole-host-contig repeat evidence is not assigned to an embedded prophage. BACPHLIP is assessed only for CheckV-complete Caudoviricetes genomes; other sequences remain not assessed. Terminal repeats alone do not establish physical circularity, completeness or biological activity.

## Recovery/purity decisions and separate QC warnings

### Configurable cutoffs

The values below are defaults in `config.yaml`; override them with `--config` when needed. Read alignment identity/length and the independent 95/85 clustering criterion are the shared workflow policies, not additional isolate-specific config switches.

| Parameter | Default | Use |
| --- | --- | --- |
| `min_len` | 1000 bp | Minimum saved SPAdes assembly length, before catalogue creation. |
| `isolate_min_depth` | 5x | Minimum own-assembly mean depth for retention. |
| `isolate_min_length_bp` | 4000 bp | Minimum length unless the short-contig depth requirement passes. |
| `isolate_short_min_depth` | 10x | Depth needed to retain a shorter saved contig. |
| `isolate_host_min_identity` | 90% | Inclusive assigned-host BLAST identity threshold. |
| `isolate_host_min_query_coverage` | 90% | Inclusive assigned-host BLAST query-coverage threshold; both host tests must pass. |
| `isolate_min_retained_mapping_percent` | 70% | `REVIEW` below this own-retained read fraction. |
| `isolate_max_host_mapping_percent` | 10% | `REVIEW` above this combined stage-02/03 fraction. |
| `isolate_max_unexplained_percent` | 10% | `REVIEW` above this final unexplained fraction. |
| `sourmash_contig_catalogue_min_length` | 5000 bp | Inclusive representative-screen cutoff; never removes shorter clusters from the catalogue. |
| `sourmash_clean_min_shared_bp` | 50000 bp | Minimum estimated shared bases for the clean-read gather/support screen. |
| `sourmash_clean_min_reference_fraction` | 0.10 | Minimum reference k-mer containment for a supported clean-read match; not mapping breadth. |

`host_prophage_activity` defaults to `True` in isolate runs with hosts. Its ratio/effect/depth/breadth thresholds and plotting flank are described above. `sourmash_clean_reads`, `sourmash_contig_catalogue`, `host_identification_test`, `map_to_RefSeq` and `metavr_blast` default to `False`. The optional host-identification evidence cutoffs are `isolate_host_test_min_reads=100` and `isolate_host_test_min_breadth_percent=1`.

The summary `decision` describes recovery/purity, independently of read-QC warning levels:

- `PASS`: passes the configured recovery/purity checks. This is not proof of a pure single-virus isolate.
- `REVIEW`: retained-contig count !=1, <70% own-retained mapping, >10% host-associated reads, >10% unexplained reads, supported host-prophage exclusions, unassessed/unavailable host references, or an ambiguous/alternative best host in the optional identification test. These concerns do not erase a recovered genome; a supported host-prophage-associated genome excluded by policy remains `REVIEW`.
- `FAIL`: no QC-passed reads, no assembled contigs, or missing, negative/non-finite or inconsistent read accounting. Stage assignments must conserve entering reads, consecutive stages must agree and all stage denominators must match. A failure takes precedence over review warnings.

`warnings` contains recovery/purity review reasons; `failure_reasons` explains failures. `read_accounting_conserved` retains the original total-count check; `read_accounting_valid` additionally checks stage transitions and count availability. Existing biological and QC thresholds are unchanged.

Read-QC levels remain available in `qc_status` and their individual status columns. The `qc_warnings` field describes low-quality filtering, SuperDeduper PCR duplicates and Kraken Eukaryota estimates, including percentages. These appear in a separate notebook/HTML table. A QC `WARN` or `FAIL` does not automatically make the isolate decision `REVIEW` or `FAIL`; missing QC measurements remain unassessed (`INFO`/missing), not evidence that QC passed. A sample can therefore have `decision=PASS` and a PCR-duplicate warning.

Existing project reports are not edited in place: rerun the workflow to regenerate them. Previously product-clustered outputs must be regenerated using `--forcerun vOUTclustering` when migrating. The general `07_Normalise.py.ipynb` and virome-positive filtering remain unchanged; independent original 95/85 clustering is the intentional repository-wide change.

## Technical appendix: Snakemake dependency graphs

These graphs are generated from the actual `phage_isolates` dependencies, unlike the conceptual step diagrams above. Keep them for checking rule relationships and database provisioning. Click the SVG for readable rule names when zooming.

![Phage isolate rule dependencies, including database downloads](figures/phage_isolates_rulegraph.png)

[Zoomable rule overview](figures/phage_isolates_rulegraph.svg) · [Full two-sample job DAG (PNG)](figures/phage_isolates_dag.png) · [Full DAG (SVG)](figures/phage_isolates_dag.svg)

The example has two samples assigned to two hosts, both Sourmash screens enabled, and the default host-activity report enabled. Optional RefSeq read mapping, METAVR BLAST and all-host identification are off. There are no pre-existing results, databases or downloaded tools in the temporary example, so provisioning nodes are included. Pale orange nodes are database/tool download rules. Other node colours are Snakemake's rule colours, not biological classifications.

The full job DAG retains separate sample, host and accounting-stage jobs. The rule overview collapses repeated jobs by rule name; reusable rules can therefore form loops in that overview even though the expanded job DAG is acyclic. Shared edge routes are bundled for readability; all original nodes and dependency edges remain in the DOT sources. The `remove_user_contaminants_PE` node remains as a bypass/link step in isolate mode, not biological read filtering. Conda package installation and Python downloads inside notebooks do not appear as separate rule nodes.

| Download rule | Used for |
| --- | --- |
| `downloadKrakenDB`, `getKrakenTools` | Read-level contamination reporting. |
| `downloadGenomadDB` | Assembly and host classification evidence. |
| `downloadCheckvDB` | Isolate-contig and host viral-region quality assessment. |
| `downloadCheckMDB` | Host assembly completeness/contamination. |
| `downloadRefSeqViral` | Core cluster-representative BLAST; also reused by optional RefSeq read mapping. The downloader keeps a dated snapshot and undated links. |
| `downloadSourmashRocksDB`, `downloadSourmashTaxonomy` | GTDB RS226 references and matching taxonomy, shared by both optional Sourmash screens. |

Regenerate the technical graphs from the repository root:

```bash
conda activate Mosaic_workflow
python mosaic/docs/generate_isolate_dag.py
```

Graphviz (`dot`) must be available. This generator asks Snakemake for `--rulegraph` and `--dag` in an isolated temporary example. It does not execute jobs, download databases, create Conda environments or alter project results. The DOT sources, SVGs and configuration/version manifest are saved beside the PNGs. To draw the optional branches without replacing the documented baseline, use:

```bash
python mosaic/docs/generate_isolate_dag.py \
  --with-refseq --with-metavr --with-host-test \
  --output-dir /tmp/mosaic-isolate-extended-graphs
```
