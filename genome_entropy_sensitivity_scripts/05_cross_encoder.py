#!/usr/bin/env python3
"""Stratified ModernProst/ProstT5 3Di benchmark; optional native repository GPU inference."""

from _shared import *
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold


def main():
    p = parser(__doc__)
    p.add_argument(
        "--metadata",
        required=True,
        help="TSV/parquet qid,dataset,group_id,label,length,three_di_entropy; optional arm,score,weight",
    )
    p.add_argument(
        "--sequences",
        required=True,
        help="Native LOGAN seqs/*.tsv or joined phage table: query_id/qid,aa,three_di",
    )
    p.add_argument(
        "--dataset",
        default="panel",
        help="Dataset label when metadata has no dataset column",
    )
    p.add_argument("--orfs", help="Native phage orfs.parquet to enrich selection rows")
    p.add_argument("--labels", help="Native phage corrected labels.parquet")
    p.add_argument(
        "--encoded",
        help="Precomputed TSV qid,prostt5_three_di for selected panel; otherwise run inference",
    )
    p.add_argument("--repo", type=Path)
    p.add_argument(
        "--revision",
        help="Pinned Hugging Face model revision; required for new inference",
    )
    p.add_argument("--model", default="Rostlab/ProstT5_fp16")
    p.add_argument("--device", default="cuda")
    p.add_argument("--encoding-size", type=int, default=4000)
    p.add_argument("--n-per-stratum", type=int, default=300)
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--threshold", type=float, default=2.5)
    a = p.parse_args()

    start(a)
    d = read(a.metadata)

    if a.orfs:
        o = read(a.orfs)
        unique(o, ["row"], "orfs")
        d = d.merge(
            o[[x for x in o if x not in d or x == "row"]],
            on="row",
            validate="one_to_one",
        )

    if a.labels:
        lab = read(a.labels)
        unique(lab, ["row"], "labels")
        d = d.drop(columns=["in_genbank"], errors="ignore").merge(
            lab[["row", "in_genbank"]],
            on="row",
            validate="one_to_one",
        )

    d = norm(d)

    if "dataset" not in d:
        d["dataset"] = a.dataset

    if "group_id" not in d:
        gc = next(
            (x for x in ["accession", "genome", "input_id"] if x in d),
            None,
        )
        if gc is None:
            raise ValueError("Need group_id or accession/genome/input_id")
        d["group_id"] = d[gc]

    if "label" not in d:
        if "in_genbank" in d:
            d["label"] = boolean(d.in_genbank)
        elif "arm" in d:
            d["label"] = d.arm.isin(
                ["annotated", "annotated_cds", "prodigal"]
            )
        else:
            raise ValueError("Need operational label or native arm/in_genbank")
    elif not d.label.astype(str).str.lower().isin(
        ["true", "false", "1", "0"]
    ).all():
        known = {
            "prodigal_complete",
            "prodigal_partial",
            "antisense_overlap",
            "frameshift_overlap",
            "intergenic",
            "marginal_overlap",
            "same_frame_unmatched",
            "frame_undefined_overlap",
        }
        if not set(d.label) <= known:
            raise ValueError(
                "Unknown operational label values; provide explicit Boolean labels"
            )
        d["label"] = d.label.isin(
            ["prodigal_complete", "prodigal_partial"]
        )

    require(
        d,
        ["dataset", "group_id", "label", "length", "three_di_entropy"],
    )
    unique(d, ["qid"], "metadata")
    d["label"] = boolean(d.label).astype(int)

    d["length_bin"] = pd.cut(
        d.length,
        [0, 50, 90, 100, 150, 300, 500, 1000, np.inf],
        right=False,
    ).astype(str)

    d["near_threshold"] = (
        d.three_di_entropy - a.threshold
    ).abs() <= 0.15

    strata = [
        "dataset",
        "label",
        "length_bin",
        "near_threshold",
    ] + (["arm"] if "arm" in d else [])

    # Iterate over groups explicitly so grouping columns are retained
    # under both pandas 2 and pandas 3.
    groups = [
        g.sample(
            n=min(len(g), a.n_per_stratum),
            random_state=a.seed,
        )
        for _, g in d.groupby(
            strata,
            dropna=False,
            observed=True,
        )
    ]

    if not groups:
        raise ValueError("No ORFs available for the cross-encoder panel")

    panel = pd.concat(groups, ignore_index=True)

    seq = read(a.sequences)

    if "qid" not in seq:
        if "query_id" in seq:
            seq["qid"] = seq.query_id
        elif "row" in seq and "row" in panel:
            seq = seq.merge(
                panel[["row", "qid"]],
                on="row",
                validate="one_to_one",
            )
        else:
            raise ValueError(
                "Sequence table needs qid/query_id, or native phage row"
            )

    require(seq, ["qid", "aa", "three_di"])
    unique(seq, ["qid"], "sequences")

    panel = panel.merge(
        seq[["qid", "aa", "three_di"]],
        on="qid",
        how="left",
        validate="one_to_one",
    )

    if panel.aa.isna().any():
        raise ValueError(
            "Selected panel sequences missing; extract them from native caches first"
        )

    if not (panel.aa.str.len() == panel.length).all():
        raise ValueError("AA sequence lengths differ from metadata")

    modern = panel.three_di.map(entropy)

    if np.max(abs(modern - panel.three_di_entropy)) > 0.001:
        raise ValueError(
            "ModernProst token entropy differs from stored values/configuration"
        )

    save(panel, a.out / "panel.tsv.gz")

    if a.encoded:
        encoded = read(a.encoded)
        require(encoded, ["qid", "prostt5_three_di"])
        unique(encoded, ["qid"], "ProstT5 results")

        panel = panel.merge(
            encoded[["qid", "prostt5_three_di"]],
            on="qid",
            how="left",
            validate="one_to_one",
        )

        if panel.prostt5_three_di.isna().any():
            raise ValueError("ProstT5 results missing selected panel IDs")
    else:
        if not a.repo or not a.revision:
            raise ValueError(
                "Inference requires --repo and --revision (a pinned model snapshot)"
            )

        from huggingface_hub import snapshot_download

        snapshot = snapshot_download(a.model, revision=a.revision)

        sys.path.insert(0, str(a.repo / "src"))
        from genome_entropy.encode3di.prostt5 import ProstT5ThreeDiEncoder

        encoder = ProstT5ThreeDiEncoder(
            model_name=snapshot,
            device=a.device,
        )

        panel["prostt5_three_di"] = encoder.encode(
            panel.aa.tolist(),
            encoding_size=a.encoding_size,
        )

        save(
            panel[["qid", "prostt5_three_di"]],
            a.out / "prostt5_tokens.tsv.gz",
        )

    if not (panel.prostt5_three_di.str.len() == panel.length).all():
        raise ValueError(
            "ProstT5 output must have one token per AA; no silent truncation"
        )

    panel["prostt5_entropy"] = panel.prostt5_three_di.map(entropy)
    panel["entropy_difference"] = (
        panel.prostt5_entropy - panel.three_di_entropy
    )
    panel["modern_fixed"] = panel.three_di_entropy >= a.threshold
    panel["prost_fixed"] = panel.prostt5_entropy >= a.threshold

    rows = []
    cal = []

    for (ds, lb), g in panel.groupby(["dataset", "length_bin"]):
        if g.label.nunique() < 2:
            continue

        rows.append(
            dict(
                dataset=ds,
                length_bin=lb,
                n=len(g),
                spearman=g.three_di_entropy.corr(
                    g.prostt5_entropy,
                    method="spearman",
                ),
                modern_auc=roc_auc_score(
                    g.label,
                    g.three_di_entropy,
                ),
                prost_auc=roc_auc_score(
                    g.label,
                    g.prostt5_entropy,
                ),
                modern_ap=average_precision_score(
                    g.label,
                    g.three_di_entropy,
                ),
                prost_ap=average_precision_score(
                    g.label,
                    g.prostt5_entropy,
                ),
                fixed_agreement=(
                    g.modern_fixed == g.prost_fixed
                ).mean(),
                fixed_jaccard=(
                    (g.modern_fixed & g.prost_fixed).sum()
                    / max(
                        1,
                        (g.modern_fixed | g.prost_fixed).sum(),
                    )
                ),
            )
        )

    for ds, g in panel.groupby("dataset"):
        if g.group_id.nunique() < a.folds:
            continue

        for fold, (tr, te) in enumerate(
            GroupKFold(a.folds).split(g, groups=g.group_id)
        ):
            train, test = g.iloc[tr], g.iloc[te]

            if train.label.nunique() < 2:
                continue

            # Equal 95% positive-label sensitivity, learned on
            # training groups only.
            for col in ["three_di_entropy", "prostt5_entropy"]:
                threshold = train.loc[
                    train.label == 1, col
                ].quantile(0.05)

                for _, r in test.iterrows():
                    cal.append(
                        dict(
                            qid=r.qid,
                            dataset=ds,
                            fold=fold,
                            encoder=col,
                            threshold=threshold,
                            selected=r[col] >= threshold,
                            label=r.label,
                        )
                    )

    save(
        pd.DataFrame(rows),
        a.out / "within_length_metrics.tsv",
    )
    save(
        pd.DataFrame(cal),
        a.out / "group_heldout_thresholds.tsv.gz",
    )
    save(
        panel.sort_values(
            "entropy_difference",
            key=abs,
            ascending=False,
        ),
        a.out / "paired_results.tsv.gz",
    )

    print(
        "AUC/AP refer to this deliberately stratified panel and operational "
        "labels; threshold overlap is an entropy-screen diagnostic, not "
        "reapplication of the original 5-feature candidate classifier."
    )


if __name__ == "__main__":
    main()
