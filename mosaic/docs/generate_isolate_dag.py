#!/usr/bin/env python3
"""Draw the real isolate workflow in a fresh, temporary two-sample project."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


WORKFLOW = Path(__file__).resolve().parents[1]
SETTINGS = {
    "isolates": True,
    "metagenome": False,
    "sourmash_clean_reads": True,
    "sourmash_contig_catalogue": True,
    "sourmash_contig_catalogue_min_length": 5000,
    "host_prophage_activity": True,
    "host_identification_test": False,
    "map_to_RefSeq": False,
    "metavr_blast": False,
    "RNA_enriched": False,
    "map_to_all_assembled": False,
}


def collect_graph(snakemake, directory, config, kind, environment):
    command = [snakemake, "--snakefile", "Snakefile", "phage_isolates",
               "--nolock", "--cores", "1", "--" + kind, "--config"]
    command.extend(f"{key}={value}" for key, value in config.items())
    result = subprocess.run(command, cwd=directory, env=environment, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    # The workflow prints context before Snakemake's DOT document.
    start = result.stdout.index("digraph ")
    graph = result.stdout[start:].strip() + "\n"
    nodes = re.findall(r'^\s*\d+\[label = "([^"\\]+)', graph, re.MULTILINE)
    downloads = sorted({name for name in nodes if name.startswith(("download", "get", "install"))})
    if not downloads:
        raise RuntimeError("The example graph must contain missing database/tool dependencies")
    return graph, nodes, downloads


def render_graph(graph, prefix, dot, title):
    # Keep Snakemake's nodes/edges; change layout and highlight provisioning only.
    lines = []
    for line in graph.splitlines():
        match = re.search(r'^\s*\d+\[label = "([^"\\]+)', line)
        if match and match.group(1).startswith(("download", "get", "install")):
            line = line.replace('style="rounded"', 'style="rounded,filled", fillcolor="#fff0dc"')
        lines.append(line)
    styled = "\n".join(lines) + "\n"
    prefix.with_suffix(".dot").write_text(styled)
    for extension in ["png", "svg"]:
        subprocess.run([dot, "-T" + extension, "-Gdpi=144", "-Grankdir=TB",
                        "-Gconcentrate=true", "-Gsplines=polyline",
                        "-Granksep=0.35", "-Gnodesep=0.18", "-Gbgcolor=white",
                        "-Glabelloc=t", "-Glabel=" + title, "-Gfontname=DejaVu Sans",
                        "-Gfontsize=18", "-Nfontname=DejaVu Sans", "-Nfontsize=10",
                        "-Ecolor=#94a3b8", "-Earrowsize=0.6",
                        "-o", str(prefix.with_suffix("." + extension))],
                       input=styled, text=True, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snakemake", default=shutil.which("snakemake"))
    parser.add_argument("--dot", default=shutil.which("dot"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "figures")
    parser.add_argument("--with-refseq", action="store_true", help="Include optional RefSeq read mapping/detection")
    parser.add_argument("--with-metavr", action="store_true", help="Include optional METAVR BLAST")
    parser.add_argument("--with-host-test", action="store_true", help="Include optional all-host identification mappings")
    args = parser.parse_args()
    if not args.snakemake or not args.dot:
        parser.error("Activate the MOSAIC workflow environment and install Graphviz (dot)")
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    settings = dict(SETTINGS, map_to_RefSeq=args.with_refseq,
                    metavr_blast=args.with_metavr, host_identification_test=args.with_host_test)
    manifest = {"target": "phage_isolates", "settings": settings,
                "samples": ["IsolateA", "IsolateB"], "hosts": ["HostA", "HostB"],
                "snakemake_version": subprocess.check_output([args.snakemake, "--version"], text=True).strip(),
                "graph_scope": "Fresh example project with no installed databases, tools or previous outputs",
                "source_sha256": {str(path.relative_to(WORKFLOW)): hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in [WORKFLOW / "Snakefile", WORKFLOW / "config.yaml",
                                               *sorted((WORKFLOW / "rules").glob("*.smk"))]},
                "graphs": {}}
    with tempfile.TemporaryDirectory(prefix="mosaic_isolate_dag_") as temporary:
        root = Path(temporary)
        workflow = root / "workflow"
        workflow.mkdir()
        for name in ["Snakefile", "config.yaml", "config.py"]:
            shutil.copy2(WORKFLOW / name, workflow / name)
        for name in ["rules", "scripts", "notebooks", "envs"]:
            (workflow / name).symlink_to(WORKFLOW / name, target_is_directory=True)
        project = root / "project"
        reads = project / "00_RAW_DATA"
        reads.mkdir(parents=True)
        hosts = project / "HOST"
        hosts.mkdir()
        for sample, host in zip(manifest["samples"], manifest["hosts"]):
            for mate in ["R1", "R2"]:
                with gzip.open(reads / f"{sample}_{mate}.fastq.gz", "wt") as handle:
                    handle.write(f"@example/{mate[-1]}\n" + "ACGT" * 25 + "\n+\n" + "I" * 100 + "\n")
            (hosts / f"{host}.fasta").write_text(f">{host}_chromosome\n" + "ACGT" * 2500 + "\n")
        (project / "host_mapping_file.tsv").write_text("sample\thost\nIsolateA\tHostA\nIsolateB\tHostB\n")
        config = dict(settings, input_dir=str(reads), results_dir=str(project))
        environment = dict(os.environ, XDG_CACHE_HOME=str(root / "cache"))
        for kind, name, description in [
            ("rulegraph", "phage_isolates_rulegraph", "MOSAIC phage isolates: rule-level dependency overview"),
            ("dag", "phage_isolates_dag", "MOSAIC phage isolates: complete two-sample job DAG"),
        ]:
            graph, nodes, downloads = collect_graph(args.snakemake, workflow, config, kind, environment)
            render_graph(graph, args.output_dir / name, args.dot, description)
            manifest["graphs"][kind] = {"nodes": len(nodes), "rules": sorted(set(nodes)),
                                        "download_and_tool_rules": downloads}
            print(f"{name}: {len(nodes)} nodes; provisioning: {', '.join(downloads)}")
    (args.output_dir / "phage_isolates_graph_info.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
