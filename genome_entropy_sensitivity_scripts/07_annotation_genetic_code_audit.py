#!/usr/bin/env python3
"""Re-evaluate exact native CDS matches, coverage and genetic codes from phage JSON/GenBank pairs."""

import gzip
from collections import defaultdict
from _shared import *
from sklearn.metrics import roc_auc_score, average_precision_score


def main():
    p = parser(__doc__)
    p.add_argument("--repo", required=True, type=Path)
    p.add_argument("--json-dir", required=True, type=Path)
    p.add_argument("--gbk-dir", required=True, type=Path)
    p.add_argument("--orfs", required=True, help="Native orfs.parquet")
    p.add_argument(
        "--oof", help="Native corrected-label oof.parquet (row,label,p_intrinsic)"
    )
    p.add_argument("--coverage", type=float, default=0.8)
    p.add_argument("--audit-n", type=int, default=100)
    p.add_argument(
        "--refit",
        action="store_true",
        help="Run native 03_models.py on all four filtered scenarios",
    )
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args()
    start(a)
    sys.path.insert(0, str(a.repo / "src"))
    from genome_entropy.io.genbank import (
        extract_cds_features,
        evaluate_orf_genbank_cds_match,
        normalise_orf_interval,
    )
    from genome_entropy.orf.types import OrfRecord

    rows = []
    codes = defaultdict(set)
    files = sorted(list(a.json_dir.glob("*.json")) + list(a.json_dir.glob("*.json.gz")))
    if not files:
        raise FileNotFoundError("No JSON files")
    for jp in files:
        part = jp.name.removesuffix(".gz").removesuffix(".json")
        gp = a.gbk_dir / (part + ".gbk.gz")
        if not gp.exists():
            gp = a.gbk_dir / (part + ".gbk")
        by = defaultdict(list)
        for cds in extract_cds_features(gp):
            by[cds.parent_id].append(cds)
            codes[cds.parent_id].add(cds.translation_table)
        # Native CDS table records actual transl_table. Parse it independently for exclusions.
        from Bio import SeqIO

        with gzip.open(gp, "rt") if str(gp).endswith(".gz") else open(gp) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    if feature.type == "CDS":
                        codes[rec.id].add(
                            int(feature.qualifiers.get("transl_table", ["11"])[0])
                        )
        with gzip.open(jp, "rt") if str(jp).endswith(".gz") else open(jp) as f:
            data = json.load(f)
        for rec in data:
            iid = str(rec["input_id"])
            for oid, feature in (rec.get("features") or {}).items():
                loc, md = feature["location"], feature["metadata"]
                aa = feature["protein"]["aa_sequence"]
                orf = OrfRecord(
                    parent_id=iid,
                    orf_id=oid,
                    start=int(loc["start"]),
                    end=int(loc["end"]),
                    strand=loc["strand"],
                    frame=int(loc["frame"]),
                    nt_sequence="",
                    aa_sequence=aa,
                    table_id=int(md.get("table_id", 11)),
                    has_start_codon=bool(md["has_start_codon"]),
                    has_stop_codon=bool(md["has_stop_codon"]),
                )
                covs = []
                for cds in by[iid]:
                    if cds.strand != orf.strand:
                        continue
                    result = evaluate_orf_genbank_cds_match(orf, cds)
                    if result.matched:
                        iv = normalise_orf_interval(
                            orf.start, orf.end, orf.strand, cds.record_length
                        )
                        covs.append(
                            (
                                result.overlap_nt / iv.length,
                                str(getattr(cds, "locus_tag", "")),
                            )
                        )
                rows.append(
                    dict(
                        part=part,
                        input_id=iid,
                        orf_id=oid,
                        matched=bool(covs),
                        max_matched_cds_orf_coverage=max(
                            (x[0] for x in covs), default=0
                        ),
                        matching_CDS=len(covs),
                    )
                )
    audit = pd.DataFrame(rows)
    unique(audit, ["input_id", "orf_id"], "native ORF audit")
    audit["alternate_code"] = audit.input_id.map(
        lambda x: any(t != 11 for t in codes[x])
    )
    orfs = read(a.orfs)
    unique(orfs, ["row"], "orfs")
    require(
        orfs,
        [
            "row",
            "input_id",
            "orf_id",
            "in_genbank",
            "three_di_entropy",
            "protein_length",
        ],
    )
    d = orfs.merge(audit, on=["part", "input_id", "orf_id"], validate="one_to_one")
    if len(d) != len(orfs):
        raise ValueError("JSON and native ORF table universes differ")
    save(d, a.out / "match_audit.tsv.gz")
    pd.DataFrame(
        [
            {"input_id": x, "codes": ",".join(map(str, sorted(v)))}
            for x, v in codes.items()
            if any(t != 11 for t in v)
        ]
    ).to_csv(a.out / "alternate_code_genomes.tsv", sep="\t", index=False)
    # Exclude ambiguous positive matches, never relabel them as non-coding.
    d["low_coverage_positive"] = d.matched & d.max_matched_cds_orf_coverage.lt(
        a.coverage
    )
    if a.oof:
        oof = read(a.oof)
        unique(oof, ["row"], "OOF")
        d = d.merge(
            oof[["row"] + [c for c in oof if c.startswith("p_")]],
            on="row",
            how="left",
            validate="one_to_one",
        )
    masks = {
        "all": np.ones(len(d), bool),
        "coverage_filtered": ~d.low_coverage_positive,
        "code11_only": ~d.alternate_code,
        "both": ~d.low_coverage_positive & ~d.alternate_code,
    }
    rows = []
    counts = []
    for scenario, mask in masks.items():
        z = d[mask].copy()
        labels = z[["row", "orf_id", "matched"]].rename(
            columns={"matched": "in_genbank"}
        )
        labels.to_parquet(a.out / f"{scenario}.labels.parquet", index=False)
        native = orfs[orfs.row.isin(z.row)]
        native.to_parquet(a.out / f"{scenario}.orfs.parquet", index=False)
        for metric in ["three_di_entropy"] + [c for c in z if c.startswith("p_")]:
            zz = z.dropna(subset=[metric])
            y = zz.matched.astype(int)
            rows.append(
                dict(
                    scenario=scenario,
                    score=metric,
                    n=len(zz),
                    positives=y.sum(),
                    auc=roc_auc_score(y, zz[metric]) if y.nunique() == 2 else np.nan,
                    ap=(
                        average_precision_score(y, zz[metric])
                        if y.nunique() == 2
                        else np.nan
                    ),
                )
            )
        if a.refit:
            subprocess.run(
                [
                    sys.executable,
                    str(a.repo / "Validation/Inphared/scripts/03_models.py"),
                    "--orfs",
                    str(a.out / f"{scenario}.orfs.parquet"),
                    "--labels",
                    str(a.out / f"{scenario}.labels.parquet"),
                    "--sets",
                    "intrinsic",
                    "--out",
                    str(a.out / f"{scenario}_refit"),
                    "--threads",
                    str(a.threads),
                ],
                check=True,
            )
        source = z
        if a.refit:
            refit = read(a.out / f"{scenario}_refit/oof.parquet")
            source = z[[c for c in ["row", "matched", "context"] if c in z]].merge(
                refit[["row", "p_intrinsic"]], on="row", validate="one_to_one"
            )
        for col in [x for x in source if x.startswith("p_")]:
            eligible = ~source.matched
            if "context" in source:
                eligible &= ~source.context.isin(["same_stop", "same_frame"])
            counts.append(
                dict(
                    scenario=scenario,
                    score=col,
                    n=len(source),
                    positive_labels=int(source.matched.sum()),
                    high_scoring_unmatched=int((eligible & source[col].ge(0.8)).sum()),
                    context_available="context" in source,
                    refitted=a.refit,
                )
            )
    save(pd.DataFrame(counts), a.out / "candidate_counts.tsv")
    save(pd.DataFrame(rows), a.out / "metrics.tsv")
    sample = d.groupby(
        ["matched", "alternate_code", "low_coverage_positive"], group_keys=False
    ).apply(lambda g: g.sample(min(len(g), a.audit_n), random_state=a.seed))
    save(sample, a.out / "manual_review_panel.tsv")
    print(
        "Without --refit, OOF metrics and candidate counts use fixed existing scores; --refit re-estimates scores and counts. Exact annotation matches use the checked-out repository matcher."
    )


if __name__ == "__main__":
    main()
