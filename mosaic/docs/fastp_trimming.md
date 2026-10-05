# Illumina read trimming

MOSAIC uses fastp for paired-end Illumina reads in both metagenome and phage-isolate workflows. This applies to DNA and RNA-enriched runs; Nanopore and PacBio processing are unchanged.

The command is explicit in `trim_adapters_quality_illumina_PE` in `rules/01_quality_control_short.smk`. It uses R1/R2 overlap analysis and additional paired-end adapter detection (`--detect_adapter_for_pe`). No adapter FASTA is supplied; the existing `db/adapters/adapters.fa` is left untouched and unused.

## Settings

Thresholds are configured in `config.yaml`:

The previous `trimmomatic_*` settings are replaced by the `fastp_*` settings below. Update any old command-line or project-config overrides; supplying a `trimmomatic_*` key no longer changes trimming.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `fastp_cut_front_window_size` | 1 | Window for trimming low-quality leading bases |
| `fastp_cut_front_mean_quality` | 20 | Minimum mean quality in the leading window |
| `fastp_cut_right_window_size` | 4 | Sliding-window size for cutting the low-quality suffix |
| `fastp_cut_right_mean_quality` | 20 | Minimum mean quality in the sliding window |
| `fastp_length_required` | 50 | Minimum retained read length, in bases |
| `fastp_poly_g_min_len` | 10 | Minimum poly-G tail length |
| `fastp_qualified_quality_phred` | 15 | Quality threshold for the separate whole-read filter |
| `fastp_unqualified_percent_limit` | 40 | Maximum percentage of bases below that threshold |
| `fastp_n_base_limit` | 5 | Maximum number of N bases per retained read |

Poly-G tail trimming is explicitly enabled, independent of instrument identifiers in the FASTQ headers. Merging, deduplication, base correction, poly-X trimming and low-complexity filtering remain disabled. The Q20 cutting and 50-bp minimum approximate the previous Trimmomatic settings; algorithms and the additional whole-read filters can produce different retained reads.

## Outputs and reports

The paired, unpaired and concatenated-unpaired FASTQ filenames in `02_CLEAN_DATA` are unchanged. Surviving single mates are saved rather than discarded. The separate R1/R2 unpaired intermediates remain temporary.

When `remove_euk=False` and `contaminants_list` is empty, the final-clean R1, R2 and unpaired FASTQs are relative symlinks to the trimmed reads instead of duplicate copies. In this mode, the concatenated trimmed-unpaired file is retained so its clean-read symlink cannot become dangling. When filtering is enabled, final-clean reads remain real files and the concatenated trimmed-unpaired intermediate stays temporary. With eukaryote removal enabled but no user contaminants, the plain no-eukaryote FASTQs are compressed into real clean-read files.

Each sample gets `01_QC/{sample}_fastp.html`, `.json` and `.log`. The JSON reports are included in the existing post-QC MultiQC report. The `01_QC` notebook displays a fastp summary, links to the individual reports and writes `FIGURES_AND_TABLES/01_fastp_summary.{sampling}.csv`. The fastp after-filtering summary describes retained paired reads; the existing multistep counts also include surviving unpaired reads. The old `trimmomatic` table column is now named `trimmed`.

To request trimming and its reports alone, use the `fastp` target:

```bash
snakemake --use-conda -p fastp --config input_dir=/path/to/00_RAW_DATA -j 8
```

The existing `env1.yaml` now pins fastp instead of Trimmomatic, without removing BBMap or the other dependencies. Normal workflow commands do not need an extra flag.

## Existing projects

An existing cleaned FASTQ is not converted into a fastp result by changing the workflow. Dry-run first and review which downstream steps are scheduled:

```bash
snakemake --use-conda -p runWorkflow --config input_dir=/path/to/00_RAW_DATA -j 8 -n
```

Missing fastp reports and the changed trimming command/environment will normally schedule trimming again. If needed, explicitly request `--forcerun trim_adapters_quality_illumina_PE` with your usual workflow target. This can rerun affected downstream analyses; do not apply it to a running workflow.

fastp's options and paired/unpaired behavior are described in its [official documentation](https://github.com/OpenGene/fastp).

Existing copied clean-read files are not replaced automatically. Symlinks are created the next time `remove_user_contaminants_PE` runs in passthrough mode. Do not replace FASTQs or force their producing rules while a workflow using them is running.
