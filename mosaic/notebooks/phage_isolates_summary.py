from __future__ import annotations

import html
import math
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


STATUS_ORDER = {"PASS": 0, "OK": 0, "INFO": 0, "WARN": 1, "REVIEW": 1, "FAIL": 2}
CHECKV_GOOD = {"Complete", "High-quality", "Medium-quality"}


def _setup_plots() -> None:
    sns.set_style("ticks", {"axes.grid": True})
    sns.set_palette("colorblind")
    plt.rcParams["axes.linewidth"] = 1.5
    plt.rcParams["xtick.major.width"] = 1.5
    plt.rcParams["ytick.major.width"] = 1.5
    plt.rcParams["xtick.major.size"] = 8
    plt.rcParams["ytick.major.size"] = 8
    plt.rcParams["axes.titlepad"] = 16
    plt.rcParams["axes.titlesize"] = 18
    plt.rcParams["axes.labelsize"] = 14
    plt.rcParams["xtick.labelsize"] = 10
    plt.rcParams["ytick.labelsize"] = 10
    plt.rcParams["legend.fontsize"] = 10
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Liberation Sans", "DejaVu Sans"]
    plt.rcParams["text.usetex"] = False
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["savefig.dpi"] = 300


def _as_list(value) -> list:
    if value is None:
        return []
    if isinstance(value, (str, os.PathLike)):
        return [str(value)]
    try:
        return [str(x) for x in value]
    except TypeError:
        return [str(value)]


def _get_named(container, name, default=None):
    try:
        return getattr(container, name)
    except Exception:
        return default


def _ensure_parent(*paths) -> None:
    for path in paths:
        if path:
            Path(str(path)).parent.mkdir(parents=True, exist_ok=True)


def _read_table(path, sep="\t", names=None, header="infer") -> pd.DataFrame:
    path = Path(str(path))
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=names or [])
    try:
        return pd.read_csv(path, sep=sep, names=names, header=header)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=names or [])


def _to_numeric(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    for column in columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def _first_present(df: pd.DataFrame, candidates: list[str]):
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    return None


def _clean_sample_name(value: str, sampling: str = "tot") -> str:
    sample = os.path.basename(str(value))
    replacements = [
        ".fastq.gz",
        ".fq.gz",
        ".fastq",
        ".fq",
        ".fasta",
        ".fa",
        ".fna",
        ".csv",
        ".tsv",
        ".txt",
        f"_spades_filtered_scaffolds.{sampling}",
        f"_spades_scaffolds.{sampling}",
        f"_assembled_contigs.{sampling}",
        f"_viral_contigs.{sampling}",
        f"_unfiltered_contigs.{sampling}",
        f"_{sampling}",
    ]
    for suffix in replacements:
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
    return sample


def sample_from_contig(contig: str) -> str:
    contig = str(contig)
    if "_NODE_" in contig:
        return contig.split("_NODE_")[0]
    if "_provirus_" in contig:
        return contig.split("_provirus_")[0]
    if "-provirus_" in contig:
        return contig.split("-provirus_")[0]
    if "_contig_" in contig:
        return contig.split("_contig_")[0]
    return contig.split("|")[0]


def _status_from_percent(value, warn: float, fail: float, high_bad: bool = True) -> str:
    if pd.isna(value):
        return "INFO"
    if high_bad:
        if value >= fail:
            return "FAIL"
        if value >= warn:
            return "WARN"
        return "PASS"
    if value < fail:
        return "FAIL"
    if value < warn:
        return "WARN"
    return "PASS"


def _worst_status(values) -> str:
    values = [str(value) for value in values if str(value) and str(value) != "nan"]
    if not values:
        return "INFO"
    return max(values, key=lambda value: STATUS_ORDER.get(value, 1))


def _save_table(df: pd.DataFrame, path: str) -> None:
    _ensure_parent(path)
    df.to_csv(path, index=False)


def _concat_nonempty(frames: list[pd.DataFrame], columns: list[str] | None = None) -> pd.DataFrame:
    frames = [frame for frame in frames if frame is not None and not frame.empty]
    if not frames:
        return pd.DataFrame(columns=columns or [])
    return pd.concat(frames, ignore_index=True)


def _placeholder_plot(path_png: str, path_svg: str, title: str, message: str) -> None:
    _ensure_parent(path_png, path_svg)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    ax.set_title(title)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path_png)
    fig.savefig(path_svg)
    plt.close(fig)


def _save_fig(fig, path_png: str, path_svg: str) -> None:
    _ensure_parent(path_png, path_svg)
    fig.tight_layout()
    fig.savefig(path_png)
    fig.savefig(path_svg)
    plt.close(fig)


def _read_id_list(path) -> set[str]:
    ids = set()
    for item in _as_list(path):
        item_path = Path(item)
        if not item_path.exists() or item_path.stat().st_size == 0:
            continue
        with item_path.open() as handle:
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#"):
                    ids.add(line.split()[0].split(",")[0])
    return ids


def load_qc(path: str) -> pd.DataFrame:
    columns = [
        "sample",
        "raw_reads",
        "clean_reads",
        "low_quality_percent",
        "eukaryotic_percent",
        "phix_contaminant_percent",
        "pcr_duplicate_percent",
        "assembly_reads_percent",
    ]
    path = Path(str(path))
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=columns)

    df = pd.read_csv(path)
    if "sample" not in df.columns:
        return pd.DataFrame(columns=columns)
    if "type" in df.columns:
        df = df[df["type"].astype(str).isin(["R1", "merged", "single", "nanopore", "pacbio"])].copy()
        df = df.drop_duplicates("sample", keep="first")

    out = pd.DataFrame({"sample": df["sample"].astype(str)})
    out["raw_reads"] = pd.to_numeric(df.get("raw", np.nan), errors="coerce")
    out["clean_reads"] = pd.to_numeric(df.get("bbduk", df.get("kraken", np.nan)), errors="coerce")
    percent_map = {
        "low_quality_percent": "low_QC_reads_p",
        "eukaryotic_percent": "eukaryotic_reads_p",
        "phix_contaminant_percent": "bbduk_phix174_reads_p",
        "pcr_duplicate_percent": "duplicate_reads_p",
        "assembly_reads_percent": "assembly_reads_p",
    }
    for out_col, in_col in percent_map.items():
        out[out_col] = pd.to_numeric(df.get(in_col, np.nan), errors="coerce") * 100.0

    out["low_quality_status"] = out["low_quality_percent"].apply(
        lambda value: _status_from_percent(value, warn=10, fail=30)
    )
    out["eukaryotic_status"] = out["eukaryotic_percent"].apply(
        lambda value: _status_from_percent(value, warn=10, fail=30)
    )
    out["phix_contaminant_status"] = out["phix_contaminant_percent"].apply(
        lambda value: _status_from_percent(value, warn=10, fail=30)
    )
    out["pcr_duplicate_status"] = out["pcr_duplicate_percent"].apply(
        lambda value: _status_from_percent(value, warn=40, fail=60)
    )
    out["qc_status"] = out[
        [
            "low_quality_status",
            "eukaryotic_status",
            "phix_contaminant_status",
            "pcr_duplicate_status",
        ]
    ].apply(_worst_status, axis=1)
    return out


def load_quast(path: str, sampling: str) -> pd.DataFrame:
    path = Path(str(path))
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=["sample"])
    df = _read_table(path, sep="\t")
    if "Assembly" not in df.columns:
        return pd.DataFrame(columns=["sample"])
    df["sample"] = df["Assembly"].apply(lambda value: _clean_sample_name(value, sampling=sampling))
    numeric_cols = [
        "# contigs (>= 0 bp)",
        "# contigs (>= 1000 bp)",
        "# contigs (>= 5000 bp)",
        "Total length (>= 0 bp)",
        "Total length (>= 1000 bp)",
        "Total length (>= 5000 bp)",
        "Largest contig",
        "N50",
        "GC (%)",
    ]
    return _to_numeric(df, numeric_cols)


def load_checkv(path: str) -> pd.DataFrame:
    names = [
        "contig",
        "contig_length",
        "provirus",
        "proviral_length",
        "gene_count",
        "viral_genes",
        "host_genes",
        "checkv_quality",
        "miuvig_quality",
        "completeness",
        "completeness_method",
        "contamination",
        "kmer_freq",
        "warnings",
    ]
    df = _read_table(path, sep="\t", names=names, header=None)
    if df.empty:
        return df
    if str(df.iloc[0]["contig"]).lower() in {"contig_id", "contig"}:
        df = df.iloc[1:].copy()
    df["sample"] = df["contig"].apply(sample_from_contig)
    return _to_numeric(
        df,
        [
            "contig_length",
            "proviral_length",
            "gene_count",
            "viral_genes",
            "host_genes",
            "completeness",
            "contamination",
            "kmer_freq",
        ],
    )


def load_nucleotide_content(path: str) -> pd.DataFrame:
    df = _read_table(path, sep="\t", header=None)
    if df.empty or df.shape[1] < 6:
        return pd.DataFrame(columns=["contig", "length_nt", "gc_percent"])
    out = df.iloc[:, :6].copy()
    out.columns = ["contig", "length_nt", "A", "C", "G", "T"]
    out = _to_numeric(out, ["length_nt", "A", "C", "G", "T"])
    total = out[["A", "C", "G", "T"]].sum(axis=1).replace(0, np.nan)
    out["gc_percent"] = (out["G"] + out["C"]) * 100.0 / total
    return out[["contig", "length_nt", "gc_percent"]]


def _parse_covstats_name(path: Path, reference_set: str, sampling: str) -> tuple[str, str | None]:
    name = path.name
    host_match = re.match(r"bowtie2(?:_filtered)?_(.+?)_vs_(.+?)(?:_masked_prophages)?_covstats\.txt$", name)
    if host_match and "/HOST/" in str(path):
        return host_match.group(1), host_match.group(2)

    patterns = [
        rf"bowtie2_(.+?)_assembled_contigs\.{re.escape(sampling)}_covstats\.txt$",
        rf"bowtie2_(.+?)_viral_contigs\.{re.escape(sampling)}_covstats\.txt$",
        rf"bowtie2_(.+?)_unfiltered_contigs\.{re.escape(sampling)}_covstats\.txt$",
        rf"bowtie2_(.+?)_{re.escape(sampling)}_covstats\.txt$",
    ]
    for pattern in patterns:
        match = re.match(pattern, name)
        if match:
            return match.group(1), None
    return _clean_sample_name(name, sampling=sampling), None


def parse_covstats_file(path: str, reference_set: str, sampling: str) -> pd.DataFrame:
    path = Path(str(path))
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    df = _read_table(path, sep="\t")
    if df.empty:
        return pd.DataFrame()

    contig_col = _first_present(df, ["Contig", "contig", "Sequence"])
    if not contig_col:
        return pd.DataFrame()

    mean_col = _first_present(df, [col for col in df.columns if str(col).endswith(" Mean")])
    length_col = _first_present(df, ["Length", "length"])
    covered_col = _first_present(df, ["Covered Bases", "covered_bases"])
    read_count_col = _first_present(df, [col for col in df.columns if str(col).endswith(" Read Count")] + ["Read Count"])
    trimmed_mean_col = _first_present(df, [col for col in df.columns if str(col).endswith(" Trimmed Mean")] + ["Trimmed Mean"])
    rpkm_col = _first_present(df, [col for col in df.columns if str(col).endswith(" RPKM")] + ["RPKM"])

    out = pd.DataFrame({"contig": df[contig_col].astype(str)})
    out["sample"], host = _parse_covstats_name(path, reference_set=reference_set, sampling=sampling)
    out["host"] = host
    out["reference_set"] = reference_set
    out["mean_coverage"] = pd.to_numeric(df[mean_col], errors="coerce") if mean_col else np.nan
    out["length"] = pd.to_numeric(df[length_col], errors="coerce") if length_col else np.nan
    out["covered_bases"] = pd.to_numeric(df[covered_col], errors="coerce") if covered_col else np.nan
    out["read_count"] = pd.to_numeric(df[read_count_col], errors="coerce") if read_count_col else np.nan
    out["trimmed_mean"] = pd.to_numeric(df[trimmed_mean_col], errors="coerce") if trimmed_mean_col else np.nan
    out["rpkm"] = pd.to_numeric(df[rpkm_col], errors="coerce") if rpkm_col else np.nan
    out["breadth_percent"] = out["covered_bases"] * 100.0 / out["length"].replace(0, np.nan)
    return out


def load_covstats(paths, reference_set: str, sampling: str) -> pd.DataFrame:
    frames = [parse_covstats_file(path, reference_set=reference_set, sampling=sampling) for path in _as_list(paths)]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame(
            columns=[
                "contig",
                "sample",
                "host",
                "reference_set",
                "mean_coverage",
                "length",
                "covered_bases",
                "read_count",
                "trimmed_mean",
                "rpkm",
                "breadth_percent",
            ]
        )
    return pd.concat(frames, ignore_index=True)


def parse_flagstat_file(path: str, mapping_type: str, sampling: str) -> pd.DataFrame:
    path = Path(str(path))
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    text = path.read_text(errors="replace").splitlines()
    total_reads = np.nan
    mapped_reads = np.nan
    mapped_percent = np.nan
    for line in text:
        if " in total " in line or line.endswith(" in total"):
            match = re.match(r"\s*(\d+)", line)
            if match:
                total_reads = float(match.group(1))
        if " mapped (" in line and "mate" not in line and "singleton" not in line:
            match = re.match(r"\s*(\d+).*\(([\d.]+)%", line)
            if match:
                mapped_reads = float(match.group(1))
                mapped_percent = float(match.group(2))
                break
    name = path.name
    sample = _clean_sample_name(name, sampling=sampling)
    sample = re.sub(r"^bowtie2_flagstats(?:_filtered)?_", "", sample)
    for suffix in ["_assembled_contigs", "_viral_contigs", "_unfiltered_contigs"]:
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
    return pd.DataFrame(
        [
            {
                "sample": sample,
                "mapping_type": mapping_type,
                "flagstat_total_reads": total_reads,
                "flagstat_mapped_reads": mapped_reads,
                "flagstat_mapped_percent": mapped_percent,
            }
        ]
    )


def load_flagstats(paths, mapping_type: str, sampling: str) -> pd.DataFrame:
    frames = [parse_flagstat_file(path, mapping_type=mapping_type, sampling=sampling) for path in _as_list(paths)]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame(columns=["sample", "mapping_type"])
    return pd.concat(frames, ignore_index=True)


def load_relative_hits(path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    names = [
        "query",
        "subject",
        "subject_title",
        "q_start",
        "q_end",
        "q_len",
        "s_len",
        "qcovs",
        "evalue",
        "align_len",
        "pident",
    ]
    df = _read_table(path, sep="\t", names=names, header=None)
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()
    df = _to_numeric(df, ["q_len", "s_len", "qcovs", "evalue", "align_len", "pident"])
    df = df.dropna(subset=["query", "subject"])
    if df.empty:
        return df, pd.DataFrame()
    df["weighted_identity"] = df["pident"] * df["align_len"].fillna(0)
    grouped = (
        df.groupby(["query", "subject", "subject_title", "q_len", "s_len"], dropna=False)
        .agg(
            aligned_bases=("align_len", "sum"),
            weighted_identity=("weighted_identity", "sum"),
            best_evalue=("evalue", "min"),
        )
        .reset_index()
    )
    grouped["pident"] = grouped["weighted_identity"] / grouped["aligned_bases"].replace(0, np.nan)
    grouped["query_coverage_percent"] = (grouped["aligned_bases"] * 100.0 / grouped["q_len"].replace(0, np.nan)).clip(upper=100)
    grouped["subject_coverage_percent"] = (grouped["aligned_bases"] * 100.0 / grouped["s_len"].replace(0, np.nan)).clip(upper=100)
    grouped["length_ratio_query_subject"] = grouped["q_len"] / grouped["s_len"].replace(0, np.nan)
    grouped["good_relative"] = grouped["query_coverage_percent"].ge(90)
    grouped = grouped.sort_values(
        ["query", "good_relative", "query_coverage_percent", "pident", "aligned_bases"],
        ascending=[True, False, False, False, False],
    )
    best = grouped.drop_duplicates("query", keep="first").copy()
    counts = grouped[grouped["good_relative"]].groupby("query").size().rename("n_good_relatives").reset_index()
    best = best.merge(counts, on="query", how="left")
    best["n_good_relatives"] = best["n_good_relatives"].fillna(0).astype(int)
    best = best.rename(
        columns={
            "query": "contig",
            "subject": "best_relative_id",
            "subject_title": "best_relative_title",
            "pident": "best_relative_identity_percent",
            "query_coverage_percent": "best_relative_query_coverage_percent",
            "subject_coverage_percent": "best_relative_subject_coverage_percent",
            "s_len": "best_relative_length",
            "length_ratio_query_subject": "best_relative_length_ratio",
        }
    )
    return grouped, best


def load_host_blast(paths) -> pd.DataFrame:
    frames = []
    names = [
        "contig",
        "host_contig",
        "host_title",
        "q_start",
        "q_end",
        "q_len",
        "s_len",
        "qcovs",
        "evalue",
        "align_len",
        "pident",
    ]
    for item in _as_list(paths):
        path = Path(item)
        if not path.exists() or path.stat().st_size == 0:
            continue
        df = _read_table(path, sep="\t", names=names, header=None)
        if df.empty:
            continue
        host_match = re.search(r"blastn_out_assembly_(.+?)\.tot\.csv$", path.name)
        df["host"] = host_match.group(1) if host_match else path.stem
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["contig", "host", "host_blast_like"])
    out = pd.concat(frames, ignore_index=True)
    out = _to_numeric(out, ["q_len", "s_len", "qcovs", "evalue", "align_len", "pident"])
    out["query_coverage_percent"] = out["qcovs"]
    needs_qcov = out["query_coverage_percent"].isna()
    out.loc[needs_qcov, "query_coverage_percent"] = (
        out.loc[needs_qcov, "align_len"] * 100.0 / out.loc[needs_qcov, "q_len"].replace(0, np.nan)
    )
    out["host_blast_like"] = out["pident"].ge(90) & out["query_coverage_percent"].ge(90)
    out = out.sort_values(["contig", "host_blast_like", "query_coverage_percent", "pident"], ascending=[True, False, False, False])
    return out


def load_clusters(path: str) -> pd.DataFrame:
    path = Path(str(path))
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=["rep", "member"])
    df = _read_table(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["rep", "member"])
    lower = {str(col).lower(): col for col in df.columns}
    rep_col = lower.get("rep") or lower.get("representative") or df.columns[0]
    mem_col = lower.get("mem_d") or lower.get("member") or (df.columns[1] if len(df.columns) > 1 else df.columns[0])
    out = df[[rep_col, mem_col]].copy()
    out.columns = ["rep", "member"]
    return out.dropna()


def load_rpkm(path: str) -> pd.DataFrame:
    path = Path(str(path))
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_csv(path)
    if df.empty:
        return df
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "contig"})
    numeric_cols = [column for column in df.columns if column != "contig"]
    for column in df.columns:
        if column != "contig":
            df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0)
    if numeric_cols:
        df = df.groupby("contig", as_index=False)[numeric_cols].max()
    return df


def _build_contig_table(
    checkv: pd.DataFrame,
    nucleotide: pd.DataFrame,
    covstats: pd.DataFrame,
    host_blast: pd.DataFrame,
    relative_best: pd.DataFrame,
    clusters: pd.DataFrame,
    vibrant_ids: set[str],
    virsorter_ids: set[str],
) -> pd.DataFrame:
    if checkv.empty:
        contigs = pd.DataFrame(columns=["contig", "sample", "contig_length"])
    else:
        contigs = checkv.copy()

    if not nucleotide.empty:
        contigs = contigs.merge(nucleotide, on="contig", how="left")

    if not clusters.empty:
        cluster_counts = (
            clusters.assign(member_sample=clusters["member"].apply(sample_from_contig))
            .groupby("rep")
            .agg(cluster_size=("member", "nunique"), cluster_sample_count=("member_sample", "nunique"))
            .reset_index()
        )
        contigs = contigs.merge(clusters.rename(columns={"member": "contig", "rep": "cluster_rep"}), on="contig", how="left")
        contigs = contigs.merge(cluster_counts.rename(columns={"rep": "cluster_rep"}), on="cluster_rep", how="left")

    if not covstats.empty:
        cov = covstats[covstats["reference_set"].isin(["unfiltered_votu", "assembled_contigs", "viral_contigs", "filtered_votu"])].copy()
        cov["contig_sample"] = cov["contig"].apply(sample_from_contig)
        cov = cov[(cov["sample"] == cov["contig_sample"]) | cov["reference_set"].isin(["filtered_votu"])]
        priority = {"assembled_contigs": 0, "viral_contigs": 1, "unfiltered_votu": 2, "filtered_votu": 3}
        cov["priority"] = cov["reference_set"].map(priority).fillna(9)
        cov = cov.sort_values(["contig", "priority"]).drop_duplicates("contig", keep="first")
        cov = cov[
            [
                "contig",
                "mean_coverage",
                "breadth_percent",
                "read_count",
                "rpkm",
                "covered_bases",
                "reference_set",
            ]
        ].rename(columns={"reference_set": "coverage_source"})
        contigs = contigs.merge(cov, on="contig", how="left")

    if not relative_best.empty:
        keep = [
            "contig",
            "best_relative_id",
            "best_relative_title",
            "best_relative_identity_percent",
            "best_relative_query_coverage_percent",
            "best_relative_subject_coverage_percent",
            "best_relative_length",
            "best_relative_length_ratio",
            "n_good_relatives",
        ]
        contigs = contigs.merge(relative_best[[col for col in keep if col in relative_best.columns]], on="contig", how="left")

    if not host_blast.empty:
        best_host = host_blast.sort_values(
            ["contig", "host_blast_like", "query_coverage_percent", "pident"],
            ascending=[True, False, False, False],
        ).drop_duplicates("contig", keep="first")
        best_host = best_host[["contig", "host", "pident", "query_coverage_percent", "host_blast_like"]].rename(
            columns={
                "host": "best_host_blast",
                "pident": "best_host_identity_percent",
                "query_coverage_percent": "best_host_query_coverage_percent",
            }
        )
        contigs = contigs.merge(best_host, on="contig", how="left")
    else:
        contigs["host_blast_like"] = False

    if "checkv_quality" not in contigs.columns:
        contigs["checkv_quality"] = np.nan
    if "n_good_relatives" not in contigs.columns:
        contigs["n_good_relatives"] = 0

    contigs["sample"] = contigs.get("sample", contigs["contig"].apply(sample_from_contig)).fillna(contigs["contig"].apply(sample_from_contig))
    contigs["length_for_filter"] = contigs.get("contig_length", np.nan)
    if "length_nt" in contigs.columns:
        contigs["length_for_filter"] = contigs["length_for_filter"].fillna(contigs["length_nt"])
    if "length_for_filter" in contigs.columns and "mean_coverage" in contigs.columns:
        contigs["candidate_pass"] = (
            (contigs["length_for_filter"].ge(4000) & contigs["mean_coverage"].ge(10))
            | contigs["mean_coverage"].ge(5)
        )
    else:
        contigs["candidate_pass"] = False
    contigs["vibrant_positive"] = contigs["contig"].isin(vibrant_ids)
    contigs["virsorter_positive"] = contigs["contig"].isin(virsorter_ids)
    contigs["checkv_supported"] = contigs.get("checkv_quality", "").isin(CHECKV_GOOD)
    contigs["relative_supported"] = pd.to_numeric(contigs.get("n_good_relatives", 0), errors="coerce").fillna(0).gt(0)
    contigs["host_blast_like"] = contigs["host_blast_like"].fillna(False).astype(bool)
    evidence_cols = ["vibrant_positive", "virsorter_positive", "checkv_supported", "relative_supported"]
    contigs["viral_evidence_count"] = contigs[evidence_cols].sum(axis=1)

    def label(row) -> str:
        if bool(row.get("host_blast_like", False)):
            return "HOST-LIKE"
        if bool(row.get("candidate_pass", False)) and row.get("viral_evidence_count", 0) >= 1:
            return "VIRAL-CANDIDATE"
        if bool(row.get("candidate_pass", False)):
            return "REVIEW"
        return "LOW-SUPPORT"

    contigs["contig_call"] = contigs.apply(label, axis=1)
    contigs["relative_warning"] = pd.to_numeric(contigs.get("n_good_relatives", 0), errors="coerce").fillna(0).lt(3)
    return contigs


def _build_host_mapping(host_cov: pd.DataFrame, host_masked_cov: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for df, mapping_type in [(host_cov, "host"), (host_masked_cov, "host_masked_prophages")]:
        if df.empty:
            continue
        grouped = (
            df.groupby(["sample", "host"], dropna=False)
            .agg(
                host_mean_coverage=("mean_coverage", "mean"),
                host_covered_bases=("covered_bases", "sum"),
                host_length=("length", "sum"),
                host_read_count=("read_count", "sum"),
            )
            .reset_index()
        )
        grouped["mapping_type"] = mapping_type
        grouped["host_breadth_percent"] = grouped["host_covered_bases"] * 100.0 / grouped["host_length"].replace(0, np.nan)
        frames.append(grouped)

    if not frames:
        return pd.DataFrame(columns=["sample", "host", "mapping_type"])
    out = pd.concat(frames, ignore_index=True)
    clean_reads = qc[["sample", "clean_reads"]].drop_duplicates("sample") if "clean_reads" in qc.columns else pd.DataFrame(columns=["sample", "clean_reads"])
    out = out.merge(clean_reads, on="sample", how="left")
    out["host_reads_percent_of_clean"] = out["host_read_count"] * 100.0 / out["clean_reads"].replace(0, np.nan)
    out["host_mapping_status"] = out["host_reads_percent_of_clean"].apply(lambda value: _status_from_percent(value, warn=30, fail=50))

    host_raw = out[out["mapping_type"] == "host"][
        ["sample", "host", "host_reads_percent_of_clean", "host_mean_coverage", "host_breadth_percent"]
    ].rename(
        columns={
            "host_reads_percent_of_clean": "unmasked_host_reads_percent",
            "host_mean_coverage": "unmasked_host_mean_coverage",
            "host_breadth_percent": "unmasked_host_breadth_percent",
        }
    )
    host_masked = out[out["mapping_type"] == "host_masked_prophages"][
        ["sample", "host", "host_reads_percent_of_clean", "host_mean_coverage", "host_breadth_percent"]
    ].rename(
        columns={
            "host_reads_percent_of_clean": "masked_host_reads_percent",
            "host_mean_coverage": "masked_host_mean_coverage",
            "host_breadth_percent": "masked_host_breadth_percent",
        }
    )
    activity = host_raw.merge(host_masked, on=["sample", "host"], how="outer")
    if not activity.empty:
        activity["prophage_signal_delta_percent"] = activity["unmasked_host_reads_percent"] - activity["masked_host_reads_percent"]
        activity["possible_prophage_signal"] = activity["prophage_signal_delta_percent"].fillna(0).ge(5)
        out = out.merge(activity[["sample", "host", "prophage_signal_delta_percent", "possible_prophage_signal"]], on=["sample", "host"], how="left")
    return out


def _dominant_signal(rpkm: pd.DataFrame, covstats: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    rows = []
    if not rpkm.empty:
        for sample in samples:
            if sample not in rpkm.columns:
                continue
            series = rpkm[["contig", sample]].copy()
            series[sample] = pd.to_numeric(series[sample], errors="coerce").fillna(0)
            series = series.sort_values(sample, ascending=False)
            top = series.iloc[0] if len(series) else None
            second = series.iloc[1] if len(series) > 1 else None
            total = series[sample].sum()
            rows.append(
                {
                    "sample": sample,
                    "top_contig": top["contig"] if top is not None else np.nan,
                    "top_rpkm": top[sample] if top is not None else np.nan,
                    "second_contig": second["contig"] if second is not None else np.nan,
                    "second_rpkm": second[sample] if second is not None else 0,
                    "top_second_ratio": (top[sample] / second[sample]) if second is not None and second[sample] > 0 else np.inf,
                    "top_fraction": (top[sample] / total) if top is not None and total > 0 else np.nan,
                }
            )
    elif not covstats.empty:
        cov = covstats[covstats["reference_set"].isin(["assembled_contigs", "unfiltered_votu", "filtered_votu"])].copy()
        cov = cov.dropna(subset=["contig"])
        cov = cov.groupby(["sample", "contig"], as_index=False).agg(rpkm=("rpkm", "max"))
        cov = cov.sort_values(["sample", "rpkm"], ascending=[True, False])
        for sample, sub in cov.groupby("sample"):
            sub = sub.reset_index(drop=True)
            top = sub.iloc[0] if len(sub) else None
            second = sub.iloc[1] if len(sub) > 1 else None
            total = sub["rpkm"].fillna(0).sum()
            rows.append(
                {
                    "sample": sample,
                    "top_contig": top["contig"] if top is not None else np.nan,
                    "top_rpkm": top["rpkm"] if top is not None else np.nan,
                    "second_contig": second["contig"] if second is not None else np.nan,
                    "second_rpkm": second["rpkm"] if second is not None else 0,
                    "top_second_ratio": (top["rpkm"] / second["rpkm"]) if second is not None and second["rpkm"] > 0 else np.inf,
                    "top_fraction": (top["rpkm"] / total) if top is not None and total > 0 else np.nan,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["sample", "top_contig", "top_rpkm", "second_contig", "second_rpkm", "top_second_ratio", "top_fraction"])
    out["dominant_signal_status"] = out.apply(
        lambda row: "WARN" if pd.notna(row.get("second_rpkm")) and row.get("second_rpkm", 0) > 0 and row.get("top_second_ratio", np.inf) < 3 else "PASS",
        axis=1,
    )
    return out


def _build_sample_summary(
    samples: list[str],
    qc: pd.DataFrame,
    quast: pd.DataFrame,
    contigs: pd.DataFrame,
    flagstats: pd.DataFrame,
    host_mapping: pd.DataFrame,
    dominant: pd.DataFrame,
    config_flags: dict,
) -> pd.DataFrame:
    summary = pd.DataFrame({"sample": samples})
    for df in [qc, quast]:
        if not df.empty and "sample" in df.columns:
            summary = summary.merge(df.drop_duplicates("sample"), on="sample", how="left")

    if not flagstats.empty:
        pivot = flagstats.pivot_table(
            index="sample",
            columns="mapping_type",
            values=["flagstat_mapped_reads", "flagstat_mapped_percent"],
            aggfunc="first",
        )
        pivot.columns = ["_".join(col).strip("_") for col in pivot.columns]
        pivot = pivot.reset_index()
        summary = summary.merge(pivot, on="sample", how="left")

    for mapping_type in ["assembled_contigs", "viral_contigs", "unfiltered_contigs", "filtered_votu"]:
        col = f"flagstat_mapped_reads_{mapping_type}"
        if col in summary.columns:
            percent_col = f"{mapping_type}_reads_percent_of_clean"
            summary[percent_col] = (summary[col] / 2.0) * 100.0 / summary["clean_reads"].replace(0, np.nan)

    if "assembled_contigs_reads_percent_of_clean" in summary.columns:
        summary["assembly_mapping_status"] = summary["assembled_contigs_reads_percent_of_clean"].apply(
            lambda value: _status_from_percent(value, warn=70, fail=50, high_bad=False)
        )
    else:
        summary["assembly_mapping_status"] = "INFO"

    if not contigs.empty:
        grouped = (
            contigs.groupby("sample")
            .agg(
                n_contigs=("contig", "nunique"),
                n_candidate_pass=("candidate_pass", "sum"),
                n_host_like_contigs=("host_blast_like", "sum"),
                n_viral_candidate_contigs=("contig_call", lambda values: (values == "VIRAL-CANDIDATE").sum()),
                n_review_contigs=("contig_call", lambda values: (values == "REVIEW").sum()),
                n_good_checkv=("checkv_supported", "sum"),
                max_completeness=("completeness", "max"),
                median_contig_coverage=("mean_coverage", "median"),
            )
            .reset_index()
        )
        grouped["n_remaining_candidate_contigs"] = (
            contigs[contigs["candidate_pass"] & ~contigs["host_blast_like"]].groupby("sample").size()
        ).reindex(grouped["sample"]).fillna(0).astype(int).values
        summary = summary.merge(grouped, on="sample", how="left")

    if not host_mapping.empty:
        host_best = (
            host_mapping[host_mapping["mapping_type"] == "host"]
            .sort_values(["sample", "host_reads_percent_of_clean", "host_breadth_percent"], ascending=[True, False, False])
            .drop_duplicates("sample", keep="first")
        )
        host_best = host_best[
            [
                "sample",
                "host",
                "host_reads_percent_of_clean",
                "host_mean_coverage",
                "host_breadth_percent",
                "host_mapping_status",
                "possible_prophage_signal",
            ]
        ].rename(columns={"host": "top_host"})
        summary = summary.merge(host_best, on="sample", how="left")

    if not dominant.empty:
        summary = summary.merge(dominant, on="sample", how="left")

    count_cols = [
        "n_contigs",
        "n_candidate_pass",
        "n_host_like_contigs",
        "n_viral_candidate_contigs",
        "n_review_contigs",
        "n_good_checkv",
        "n_remaining_candidate_contigs",
    ]
    for col in count_cols:
        if col not in summary.columns:
            summary[col] = 0
        summary[col] = summary[col].fillna(0).astype(int)

    summary["viral_candidate_status"] = summary.get("n_viral_candidate_contigs", 0).apply(lambda value: "FAIL" if value == 0 else "PASS")
    summary["candidate_count_status"] = summary.get("n_remaining_candidate_contigs", 0).apply(
        lambda value: "PASS" if value == 1 else ("FAIL" if value == 0 else "WARN")
    )
    summary["phage_isolates_flag_status"] = "PASS" if config_flags.get("isolates") and not config_flags.get("metagenome") else "WARN"
    summary["remove_euk_status"] = "PASS" if not config_flags.get("remove_euk") else "WARN"

    status_cols = [
        "qc_status",
        "assembly_mapping_status",
        "host_mapping_status",
        "dominant_signal_status",
        "viral_candidate_status",
        "candidate_count_status",
        "phage_isolates_flag_status",
        "remove_euk_status",
    ]
    for col in status_cols:
        if col not in summary.columns:
            summary[col] = "INFO"
    summary["decision"] = summary[status_cols].apply(_worst_status, axis=1)
    return summary


def _plot_decisions(summary: pd.DataFrame, output_png: str, output_svg: str) -> None:
    if summary.empty or "decision" not in summary.columns:
        _placeholder_plot(output_png, output_svg, "Sample Decisions", "No sample decision data available")
        return
    order = ["PASS", "WARN", "REVIEW", "FAIL", "INFO"]
    counts = summary["decision"].value_counts().reindex(order).dropna().reset_index()
    counts.columns = ["decision", "samples"]
    fig, ax = plt.subplots(figsize=(6.5, 4))
    sns.barplot(data=counts, x="decision", y="samples", ax=ax, order=[x for x in order if x in counts["decision"].tolist()])
    ax.set_title("Phage Isolate Decisions")
    ax.set_xlabel("")
    ax.set_ylabel("Samples")
    _save_fig(fig, output_png, output_svg)


def _plot_remaining(summary: pd.DataFrame, output_png: str, output_svg: str) -> None:
    if summary.empty or "n_remaining_candidate_contigs" not in summary.columns:
        _placeholder_plot(output_png, output_svg, "Remaining Candidate Contigs", "No candidate contig data available")
        return
    df = summary.sort_values("n_remaining_candidate_contigs", ascending=False)
    height = max(4, min(16, 0.25 * len(df) + 2))
    fig, ax = plt.subplots(figsize=(9, height))
    sns.barplot(data=df, y="sample", x="n_remaining_candidate_contigs", ax=ax, color=sns.color_palette("colorblind")[0])
    ax.axvline(1, color="black", linestyle="--", linewidth=1)
    ax.set_title("Candidate Contigs After Notebook Filters")
    ax.set_xlabel("Candidate contigs not host-like")
    ax.set_ylabel("")
    _save_fig(fig, output_png, output_svg)


def _plot_host_viral(contigs: pd.DataFrame, output_png: str, output_svg: str) -> None:
    if contigs.empty or "contig_call" not in contigs.columns:
        _placeholder_plot(output_png, output_svg, "Contig Calls", "No contig call data available")
        return
    counts = contigs.groupby(["sample", "contig_call"]).size().unstack(fill_value=0)
    for col in ["VIRAL-CANDIDATE", "HOST-LIKE", "REVIEW", "LOW-SUPPORT"]:
        if col not in counts.columns:
            counts[col] = 0
    counts = counts[["VIRAL-CANDIDATE", "HOST-LIKE", "REVIEW", "LOW-SUPPORT"]]
    counts = counts.loc[counts.sum(axis=1).sort_values(ascending=False).index]
    height = max(4, min(16, 0.25 * len(counts) + 2))
    fig, ax = plt.subplots(figsize=(9, height))
    counts.plot(kind="barh", stacked=True, ax=ax, color=sns.color_palette("colorblind", n_colors=4))
    ax.set_title("Viral-Like, Host-Like, And Low-Support Contigs")
    ax.set_xlabel("Contigs")
    ax.set_ylabel("")
    ax.legend(title="")
    _save_fig(fig, output_png, output_svg)


def _plot_completeness_coverage(contigs: pd.DataFrame, output_png: str, output_svg: str) -> None:
    needed = {"mean_coverage", "completeness", "length_for_filter"}
    if contigs.empty or not needed.issubset(contigs.columns):
        _placeholder_plot(output_png, output_svg, "Completeness vs Coverage", "No CheckV or coverage data available")
        return
    df = contigs.dropna(subset=["mean_coverage", "completeness"]).copy()
    if df.empty:
        _placeholder_plot(output_png, output_svg, "Completeness vs Coverage", "No contigs had both CheckV completeness and coverage")
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    size = np.sqrt(df["length_for_filter"].fillna(0).clip(lower=1)) / 6
    sns.scatterplot(
        data=df,
        x="mean_coverage",
        y="completeness",
        hue="checkv_quality" if "checkv_quality" in df.columns else None,
        size=size,
        sizes=(20, 220),
        ax=ax,
        legend="brief",
    )
    ax.axvline(5, color="grey", linestyle="--", linewidth=1)
    ax.axvline(10, color="grey", linestyle=":", linewidth=1)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_title("CheckV Completeness vs Coverage")
    ax.set_xlabel("Mean coverage")
    ax.set_ylabel("Completeness (%)")
    _save_fig(fig, output_png, output_svg)


def _plot_dominant(dominant: pd.DataFrame, output_png: str, output_svg: str) -> None:
    if dominant.empty:
        _placeholder_plot(output_png, output_svg, "Dominant Phage Signal", "No RPKM or coverage signal available")
        return
    df = dominant.copy()
    df["top_second_ratio_plot"] = df["top_second_ratio"].replace(np.inf, np.nan).fillna(df["top_second_ratio"].replace(np.inf, np.nan).max())
    df["top_second_ratio_plot"] = df["top_second_ratio_plot"].fillna(10).clip(upper=100)
    df = df.sort_values("top_second_ratio_plot")
    height = max(4, min(16, 0.25 * len(df) + 2))
    fig, ax = plt.subplots(figsize=(9, height))
    sns.barplot(data=df, y="sample", x="top_second_ratio_plot", ax=ax, color=sns.color_palette("colorblind")[2])
    ax.axvline(3, color="black", linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_title("Dominant Phage Signal")
    ax.set_xlabel("Top / second RPKM ratio")
    ax.set_ylabel("")
    _save_fig(fig, output_png, output_svg)


def _plot_rpkm_heatmap(rpkm: pd.DataFrame, output_png: str, output_svg: str) -> None:
    if rpkm.empty or rpkm.shape[1] < 3:
        _placeholder_plot(output_png, output_svg, "vOTU RPKM Heatmap", "No RPKM matrix available")
        return
    matrix = rpkm.set_index("contig")
    matrix = matrix.loc[matrix.max(axis=1).sort_values(ascending=False).head(50).index]
    matrix = np.log10(matrix + 1)
    height = max(5, min(18, 0.2 * len(matrix) + 2))
    width = max(7, min(18, 0.25 * matrix.shape[1] + 4))
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(matrix, cmap="viridis", ax=ax)
    ax.set_title("Top vOTU RPKM Signals")
    ax.set_xlabel("Sample")
    ax.set_ylabel("vOTU")
    _save_fig(fig, output_png, output_svg)


def _plot_fragmentation(quast: pd.DataFrame, output_png: str, output_svg: str) -> None:
    cols = ["# contigs (>= 1000 bp)", "# contigs (>= 5000 bp)", "N50", "Largest contig"]
    if quast.empty or not any(col in quast.columns for col in cols):
        _placeholder_plot(output_png, output_svg, "Assembly Fragmentation", "No QUAST assembly table available")
        return
    df = quast.copy()
    contig_col = _first_present(df, ["# contigs (>= 5000 bp)", "# contigs (>= 1000 bp)", "# contigs (>= 0 bp)"])
    df = df.sort_values(contig_col, ascending=False)
    height = max(4, min(16, 0.25 * len(df) + 2))
    fig, ax = plt.subplots(figsize=(9, height))
    sns.barplot(data=df, y="sample", x=contig_col, ax=ax, color=sns.color_palette("colorblind")[4])
    ax.set_title("Assembly Fragmentation")
    ax.set_xlabel(contig_col)
    ax.set_ylabel("")
    _save_fig(fig, output_png, output_svg)


def _render_html(
    output_html: str,
    summary: pd.DataFrame,
    contigs: pd.DataFrame,
    relatives: pd.DataFrame,
    host_mapping: pd.DataFrame,
    config_flags: dict,
    plot_paths: dict[str, str],
) -> None:
    _ensure_parent(output_html)
    mean_remaining = np.nan
    if "n_remaining_candidate_contigs" in summary.columns and not summary.empty:
        mean_remaining = summary["n_remaining_candidate_contigs"].mean()
    flag_rows = "".join(
        f"<tr><th>{html.escape(str(key))}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in config_flags.items()
    )
    plot_links = "".join(
        f'<li><a href="{html.escape(os.path.basename(path))}">{html.escape(label)}</a></li>'
        for label, path in plot_paths.items()
        if path
    )
    summary_view = summary.copy()
    contig_view = contigs.copy()
    if not contig_view.empty:
        sort_cols = [col for col in ["sample", "contig_call", "length_for_filter", "mean_coverage"] if col in contig_view.columns]
        contig_view = contig_view.sort_values(sort_cols, ascending=[True, True, False, False][: len(sort_cols)])
    relatives_view = relatives.copy()
    if not relatives_view.empty and "good_relative" in relatives_view.columns:
        relatives_view = relatives_view[relatives_view["good_relative"]].head(200)
    host_view = host_mapping.copy()

    html_text = f"""
<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Phage Isolates Summary</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 28px; color: #222; }}
    h1, h2 {{ color: #17202a; }}
    table {{ border-collapse: collapse; width: 100%; margin-bottom: 28px; font-size: 13px; }}
    th, td {{ border: 1px solid #ddd; padding: 6px 8px; text-align: left; }}
    th {{ background: #f4f6f8; position: sticky; top: 0; }}
    .note {{ background: #f7fbff; border-left: 4px solid #2c7fb8; padding: 12px 16px; margin-bottom: 18px; }}
    .warn {{ background: #fff7e6; border-left: 4px solid #e69f00; padding: 12px 16px; margin-bottom: 18px; }}
    .tablewrap {{ overflow-x: auto; }}
  </style>
</head>
<body>
  <h1>Phage Isolates Summary</h1>
  <div class="note">
    This notebook summarizes isolate-level QC, assembly support, viral evidence, host-like contigs, host read mapping,
    closest relatives, and dominant phage signal. It only flags eukaryotic/PhiX/host signals; it does not remove reads.
  </div>
  <div class="warn">
    Mean remaining candidate contigs per sample:
    <strong>{'NA' if pd.isna(mean_remaining) else round(float(mean_remaining), 2)}</strong>.
    Candidate rule: length >= 4 kb and coverage >= 10x, or coverage >= 5x.
  </div>
  <h2>Run Flags</h2>
  <div class="tablewrap"><table>{flag_rows}</table></div>
  <h2>Plots</h2>
  <ul>{plot_links}</ul>
  <h2>Sample Decisions</h2>
  <div class="tablewrap">{summary_view.to_html(index=False, escape=True)}</div>
  <h2>Contig Triage</h2>
  <div class="tablewrap">{contig_view.head(500).to_html(index=False, escape=True)}</div>
  <h2>Closest Relatives With >=90% Query Coverage</h2>
  <div class="tablewrap">{relatives_view.to_html(index=False, escape=True)}</div>
  <h2>Host Mapping</h2>
  <div class="tablewrap">{host_view.to_html(index=False, escape=True)}</div>
</body>
</html>
"""
    Path(output_html).write_text(html_text)


def build_phage_isolates_summary(snakemake):
    _setup_plots()
    sampling = str(_get_named(snakemake.params, "sampling", "tot"))
    samples = list(_get_named(snakemake.params, "samples", []))
    config_flags = {
        "isolates": bool(_get_named(snakemake.params, "isolates", False)),
        "metagenome": bool(_get_named(snakemake.params, "metagenome", False)),
        "microbial": bool(_get_named(snakemake.params, "microbial", False)),
        "remove_euk": bool(_get_named(snakemake.params, "remove_euk", False)),
        "sourmash": bool(_get_named(snakemake.params, "sourmash", False)),
        "sampling": sampling,
    }

    qc = load_qc(_get_named(snakemake.input, "df_counts_paired", ""))
    if not samples and not qc.empty:
        samples = qc["sample"].astype(str).tolist()

    quast = load_quast(_get_named(snakemake.input, "quast", ""), sampling=sampling)
    checkv = load_checkv(_get_named(snakemake.input, "checkv", ""))
    nucleotide = load_nucleotide_content(_get_named(snakemake.input, "nucleotide_content", ""))
    clusters = load_clusters(_get_named(snakemake.input, "clusters", ""))
    vibrant_ids = _read_id_list(_get_named(snakemake.input, "vibrant_positive", []))
    virsorter_ids = _read_id_list(_get_named(snakemake.input, "virsorter_positive", []))
    relatives, relative_best = load_relative_hits(_get_named(snakemake.input, "viral_refseq_blast", ""))
    host_blast = load_host_blast(_get_named(snakemake.input, "host_blast", []))
    rpkm = load_rpkm(_get_named(snakemake.input, "rpkm", ""))

    covstats = _concat_nonempty(
        [
            load_covstats(_get_named(snakemake.input, "filtered_covstats", []), "filtered_votu", sampling),
            load_covstats(_get_named(snakemake.input, "assembled_covstats", []), "assembled_contigs", sampling),
            load_covstats(_get_named(snakemake.input, "viral_covstats", []), "viral_contigs", sampling),
            load_covstats(_get_named(snakemake.input, "unfiltered_covstats", []), "unfiltered_votu", sampling),
        ],
    )
    host_cov = load_covstats(_get_named(snakemake.input, "host_covstats", []), "host", sampling)
    host_masked_cov = load_covstats(_get_named(snakemake.input, "host_masked_covstats", []), "host_masked_prophages", sampling)
    host_mapping = _build_host_mapping(host_cov, host_masked_cov, qc)

    flagstats = _concat_nonempty(
        [
            load_flagstats(_get_named(snakemake.input, "filtered_flagstats", []), "filtered_votu", sampling),
            load_flagstats(_get_named(snakemake.input, "assembled_flagstats", []), "assembled_contigs", sampling),
            load_flagstats(_get_named(snakemake.input, "viral_flagstats", []), "viral_contigs", sampling),
            load_flagstats(_get_named(snakemake.input, "unfiltered_flagstats", []), "unfiltered_contigs", sampling),
        ],
    )

    contigs = _build_contig_table(
        checkv=checkv,
        nucleotide=nucleotide,
        covstats=covstats,
        host_blast=host_blast,
        relative_best=relative_best,
        clusters=clusters,
        vibrant_ids=vibrant_ids,
        virsorter_ids=virsorter_ids,
    )
    if not samples and not contigs.empty:
        samples = sorted(contigs["sample"].dropna().astype(str).unique())

    dominant = _dominant_signal(rpkm, covstats, samples)
    summary = _build_sample_summary(
        samples=samples,
        qc=qc,
        quast=quast,
        contigs=contigs,
        flagstats=flagstats,
        host_mapping=host_mapping,
        dominant=dominant,
        config_flags=config_flags,
    )

    outputs = snakemake.output
    _save_table(summary, outputs.summary_csv)
    _save_table(contigs, outputs.contig_csv)
    _save_table(relatives, outputs.closest_relatives_csv)
    _save_table(host_mapping, outputs.host_mapping_csv)

    _plot_decisions(summary, outputs.decisions_png, outputs.decisions_svg)
    _plot_remaining(summary, outputs.remaining_png, outputs.remaining_svg)
    _plot_host_viral(contigs, outputs.host_viral_png, outputs.host_viral_svg)
    _plot_completeness_coverage(contigs, outputs.completeness_coverage_png, outputs.completeness_coverage_svg)
    _plot_dominant(dominant, outputs.dominant_signal_png, outputs.dominant_signal_svg)
    _plot_rpkm_heatmap(rpkm, outputs.rpkm_heatmap_png, outputs.rpkm_heatmap_svg)
    _plot_fragmentation(quast, outputs.assembly_fragmentation_png, outputs.assembly_fragmentation_svg)

    plot_paths = {
        "Sample decisions": outputs.decisions_png,
        "Remaining candidate contigs": outputs.remaining_png,
        "Host-like and viral-like contigs": outputs.host_viral_png,
        "Completeness vs coverage": outputs.completeness_coverage_png,
        "Dominant phage signal": outputs.dominant_signal_png,
        "vOTU RPKM heatmap": outputs.rpkm_heatmap_png,
        "Assembly fragmentation": outputs.assembly_fragmentation_png,
    }
    _render_html(
        output_html=outputs.summary_html,
        summary=summary,
        contigs=contigs,
        relatives=relatives,
        host_mapping=host_mapping,
        config_flags=config_flags,
        plot_paths=plot_paths,
    )

    return {
        "summary": summary,
        "contigs": contigs,
        "closest_relatives": relatives,
        "host_mapping": host_mapping,
        "config_flags": config_flags,
    }
