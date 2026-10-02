# Assembly identifiers

One final FASTA is retained per assembly, in its established location. There
is no `RENAMED` directory and no retained original-ID copy of that FASTA.

The assembler produces a temporary `*.unrenamed.fasta`. The shared
`rename_assembled_contigs` rule in `02_assembly_short.smk` writes the final
FASTA and a neighboring `*.ids.tsv` table. Snakemake removes the temporary
input after successful processing. The helper uses Biopython already available
in the workflow environment; no new environment is needed.

```text
03_CONTIGS/{sample}_spades_filtered_scaffolds.tot.fasta
03_CONTIGS/{sample}_spades_filtered_scaffolds.tot.ids.tsv
03_CONTIGS/RNA/{sample}/rnaviralspades.fasta
03_CONTIGS/RNA/{sample}/rnaviralspades.ids.tsv
03_CONTIGS/RNA/{sample}/megahit.fasta
03_CONTIGS/RNA/{sample}/megahit.ids.tsv
03_CONTIGS/RNA/{sample}/trinity.fasta
03_CONTIGS/RNA/{sample}/trinity.ids.tsv
```

Identifiers follow `{sample}_{assembler}_{number:0{width}d}_len_{length}`, for
example `NT193_MSV_OY_01_metaspades_00001_len_1000`. The counter width is
`max(5, len(str(total_contigs)))`, calculated separately for each assembly.
An assembly of 200,000 contigs therefore starts at `000001` and ends at
`200000`; an assembly of 2,000,000 starts at `0000001`. An initial header-count
pass determines the width without loading all sequences into memory. Counters
follow FASTA order and restart at 1 per assembly. Length is calculated from the
actual sequence, including Ns. Renaming preserves sequence content and case. Existing DNA length/coverage
filtering happens first and is unchanged.

SPAdes assemblies use `metaspades` when the actual command uses `--meta`,
otherwise `spades`. Subsampled assemblies add `_sub` to the sample prefix;
assembly-depth tests add `_{percentage}pct`. Long-read assembly copies use the
actual assembler, and polished stages include their stage in the prefix.

ID tables record `contig_id`, `sample`, `assembler`, `original_id`, `length_bp`
and `original_description`. Assembly graphs keep their original identifiers;
the ID table provides the link to the old FASTA headers. No separate retained
FASTA is required to recover those headers.

The `.ids.tsv` files are lookup-only sidecars. No downstream rule reads or
requires them. Combination, classification, clustering, mapping and metadata
joins use the renamed FASTA identifiers directly. Original assembler headers
are kept only in the sidecars, not copied into the combined provenance, top
metadata or cluster membership tables.

RNA combination and all-assembled concatenation retain the named contig IDs.
Existing RNA exact-sequence/reverse-complement deduplication keeps the first
encountered named ID. Its provenance table records both member and representative
IDs (`contig_id` and `representative`). All-assembled provenance records
`contig_id`, `sample`, `assembler` and `length_bp`, derived from the named FASTAs
and their source paths. MMseqs/CheckV clustering and the final filtered vOTU catalogue filename
remain unchanged.
If a tool extracts or trims a region, the retained identifier refers to its
parent assembly; use the measured sequence length in that tool's output for
the processed sequence.

## Full-assembly geNomad evidence

With `map_to_all_assembled=True`, the existing `genomad_viral_id` rule also runs
on `03_CONTIGS/ALL_ASSEMBLED/all_assembled_contigs.tot.fasta`: **all retained
DNA and selected RNA assembly contigs, before dereplication or top selection**.
Results are under `04_VIRAL_ID/all_assembled_geNomad_tot/`. This is not a
top-100-only classification run. Per-sample DNA and per-assembler RNA
classification still run on their named final FASTAs.

Top-contig metadata subsets the full geNomad results and terminal-repeat reports
directly by `contig_id`, without translating original headers or switching
between per-assembly geNomad summaries. Missing classifications are recorded
as `not reported`. VIBRANT, VirSorter, CheckV, Pharokka and Cenote retain their
existing selected-top-contig scope. geNomad uses `--restart` to regenerate its
results when an input assembly changes.

## Existing projects

Existing assembly FASTAs require a one-time conversion with neighboring ID
tables. Convert only approved projects, while no workflow is using those
files. Conversion uses a same-directory temporary file followed by atomic
replacement, verifies sequence content/order/length, and retains no backup
FASTA. The saved original descriptions allow original headers to be restored
if required.

After conversion, rerun downstream classification, combined catalogues,
mapping, normalization and reports. Previous results contain old IDs and
must not be mixed with renamed assemblies. Assemblies themselves need not
be rebuilt solely for conversion.
Clear the old assembler's Snakemake metadata for each converted final FASTA
with `snakemake --cleanup-metadata <converted FASTA paths>`; otherwise its
old producer/code record can cause an unnecessary assembly rerun. Leave
preprocessing and downstream metadata intact.
