# Workflow benchmarks

Every computational rule and checkpoint now declares a benchmark. Input-only
targets such as `all` and `runWorkflow` do not run computations themselves.
Each benchmark has its own rule directory under `results_dir/BENCHMARK/`, with
all output wildcards in the filename. User-reference rules also include the
reference basename. For example:

```text
BENCHMARK/rna_assemble_trinity/sample=S1.tsv
BENCHMARK/microdiversity_map/sample=S1.tsv
BENCHMARK/normalise_reads_RefSeq/tot.tsv
```

This replaces shared paths that previously allowed different rules to overwrite
one another's measurements. Existing benchmark files are not moved or deleted.
DRAM-v annotation now passes `{threads}` to its command, respecting Snakemake's
thread allocation.

## Generate the summary

After any workflow stage completes, run this target from the `mosaic/` workflow
directory, using the same results directory and configuration as that stage:

```bash
snakemake --cores 1 runBenchmarkSummary --config results_dir=/path/to/results
```

The report runs in the existing workflow notebook environment. It reads existing
benchmark files without requesting their producing jobs, so this target does not
start missing assemblies or database downloads. Its parameter snapshot tracks
benchmark paths, sizes and modification times; additions, changes and removals
cause it to refresh under Snakemake's normal parameter rerun trigger. If rerun
triggers were customized to ignore parameters, use `--forcerun benchmark_summary`.
Run the report after the analysis finishes so files are not changing during the
snapshot. Report generation is explicitly requested, not automatically appended
to every analysis target.

Outputs:

- `FIGURES_AND_TABLES/09_benchmark_jobs.csv`: each measurement, its rule/wildcards,
  elapsed time, CPU time, sampled memory, I/O and available current input sizes.
- `FIGURES_AND_TABLES/09_benchmark_rules.csv`: measurement counts, summed job time,
  median/maximum runtime, largest sampled RSS, summed CPU time and I/O per rule.
- `FIGURES_AND_TABLES/09_benchmark_summary.html`: portable report with embedded plot.
- `FIGURES_AND_TABLES/09_benchmark_summary.png`: runtime and memory comparison.
- `NOTEBOOKS/09_benchmark_summary.ipynb`: executed notebook.

The summary excludes its own benchmark to avoid repeated self-triggered refreshes.
Old paths that cannot be matched to a current rule are marked `legacy_unattributed`:
their measurements are retained, but shared old filenames cannot identify which
rule last wrote them. Malformed files are listed explicitly as `invalid`.

## Interpretation

These are the measurements present in this results directory, potentially from
multiple executions. They are not a complete execution history. Reruns reuse
the same benchmark path, and absent benchmarks do not identify failed jobs.
Repeated measurements occupy separate rows and are counted in the summary.

Summed job runtime differs from workflow elapsed time because jobs overlap.
The largest sampled per-job RSS differs from total peak workflow memory.
`mean_cpu_cores` is measured CPU seconds divided by elapsed seconds; it is not
allocated threads. Standard Snakemake 7 TSVs do not record allocated threads or
historical input sizes. Therefore `threads_declared_current` and
`input_bytes_existing_current` describe the current configuration and accessible
files only. Deleted temporary inputs, directories and dynamic input functions
can make sizes incomplete. Unknown values stay missing rather than becoming zero.
The report retains native benchmark I/O fields without changing their units.

New benchmark paths can cause Snakemake to schedule completed jobs whose new
benchmark files are missing. Inspect a dry-run before resuming a large analysis.
Historical measurements cannot be recreated by adding a benchmark declaration.

## Batched counting and host mapping

Paired-end read counts are measured once per sample and processing stage:
`countReads_raw`, `countReads_trimmed`, `countReads_noEuk`, `countReads_clean`
and `countReads_norm`. Each job writes the same individual read-count files
used by the existing reports. The single-file counters remain available for
other read layouts. Each stage has its own benchmark; normalized counts also
include the sampling type in the benchmark filename.

`map_to_host` now runs unmasked and prophage-masked host mapping sequentially
in one job per sample/host, retaining all existing report filenames and mapping
settings. Its benchmark measures both mappings together, using eight threads
unless overridden. Earlier measurements under `map_to_host` measured only the
unmasked run, so their timings are not directly comparable. Existing masked
benchmark files are left untouched and appear as legacy measurements.

If a batched job needs to run, all of its outputs are regenerated. This reduces
local job scheduling overhead but does not reduce the number of alignments or
the amount of read data counted. No outputs or completion metadata are touched
automatically to bypass genuine reruns.

## Verification

Run `python -m unittest discover -s mosaic/tests -p test_benchmarks.py -v` in the
workflow environment. Checks load the short-read, RNA, nanopore, pooled-nanopore
and PacBio rules, verify unique paths and matching wildcard sets, execute the
notebook on known/legacy/invalid/empty data, and check the standalone report DAG.
