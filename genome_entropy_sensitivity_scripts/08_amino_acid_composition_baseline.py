#!/usr/bin/env python3
"""Grouped CV comparing AA entropy, 20-AA frequencies and predicted 3Di; native LOGAN/phage adapters."""

from _shared import *
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.ensemble import HistGradientBoostingClassifier

AA = "ACDEFGHIKLMNPQRSTVWY"


def main():
    p = parser(__doc__)
    p.add_argument("--metadata", required=True)
    p.add_argument(
        "--sequences",
        required=True,
        help="Native LOGAN sequence extracts or joined phage sequences",
    )
    p.add_argument(
        "--dataset", choices=["logan", "phage", "canonical"], default="canonical"
    )
    p.add_argument(
        "--labels", help="Native phage corrected labels.parquet, join on row"
    )
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--max-iter", type=int, default=300)
    p.add_argument("--group-column", help="Override accession/input_id/group_id")
    a = p.parse_args()
    start(a)
    d = read(a.metadata)
    seq = read(a.sequences)
    if a.dataset == "logan":
        require(
            d,
            [
                "accession",
                "contig_n",
                "start",
                "aa_len",
                "stop",
                "strand",
                "label",
                "trunc",
                "weight",
            ],
        )
        d = d[d.trunc == "complete"].copy()
        d["orf_id"] = (
            d.accession
            + "_"
            + d.contig_n.astype(str)
            + ":"
            + d.start.astype(str)
            + "-"
            + (d.start + 3 * (d.aa_len + d.stop)).astype(str)
            + ":"
            + d.strand.astype(str)
        )
        d["qid"] = d.orf_id
        d["y"] = (d.label == "prodigal_complete").astype(int)
        d["group_id"] = d.accession
        if "orf_id" not in seq:
            require(seq, ["query_id"])
            seq["orf_id"] = seq.query_id.str.split("|", n=1).str[-1]
        seq["qid"] = seq.orf_id
    elif a.dataset == "phage":
        require(d, ["row", "input_id", "orf_id", "in_genbank"])
        if a.labels:
            lab = read(a.labels)
            unique(lab, ["row"], "corrected labels")
            d = d.drop(columns="in_genbank").merge(
                lab[["row", "in_genbank"]], on="row", validate="one_to_one"
            )
        d["y"] = boolean(d.in_genbank).astype(int)
        d["group_id"] = d.input_id
        d["qid"] = d.row.astype(str)
        require(seq, ["row"])
        seq["qid"] = seq.row.astype(str)
    else:
        d = norm(d)
        require(d, ["label", "group_id"])
        d["y"] = boolean(d.label).astype(int)
        if "qid" not in seq:
            seq["qid"] = seq.query_id
    unique(d, ["qid"], "metadata")
    unique(seq, ["qid"], "sequences")
    require(seq, ["aa"])
    d = norm(d)
    d = d.merge(seq[["qid", "aa"]], on="qid", how="left", validate="one_to_one")
    if d.aa.isna().any():
        raise ValueError(
            "Sequences must cover EVERY selected training/test ORF; Foldseek selection alone is not the native model sample"
        )
    if not (d.aa.str.len() == d.length).all():
        raise ValueError("AA lengths disagree")
    if "aa_entropy" not in d:
        d["aa_entropy"] = d.aa.map(entropy)
    require(d, ["three_di_entropy"])
    group = a.group_column or "group_id"
    require(d, [group])
    # Frequencies use full length as denominator; unknown residues reported separately, not silently discarded.
    freq = np.asarray([[s.upper().count(x) / len(s) for x in AA] for s in d.aa])
    unknown = 1 - freq.sum(axis=1)
    length = np.log(d.length.to_numpy(float)).reshape(-1, 1)
    ent = d.aa_entropy.to_numpy().reshape(-1, 1)
    struct = d.three_di_entropy.to_numpy().reshape(-1, 1)
    sets = {
        "length": length,
        "length_aa_entropy": np.c_[length, ent],
        "length_aa_frequencies": np.c_[length, freq],
        "length_aa_entropy_3di": np.c_[length, ent, struct],
        "length_aa_frequencies_3di": np.c_[length, freq, struct],
    }
    y = d.y.to_numpy(int)
    weights = d.weight.to_numpy(float) if "weight" in d else np.ones(len(d))
    if d[group].nunique() < a.folds:
        raise ValueError("Too few independent groups")
    splits = list(GroupKFold(a.folds).split(d, y, d[group]))
    rows = []
    pred = d[["qid", group, "y", "length"]].copy()
    pred["weight"] = weights
    for name, X in sets.items():
        oof = np.full(len(d), np.nan)
        for fold, (tr, te) in enumerate(splits):
            if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
                raise ValueError(
                    "A grouped fold has only one class; revise grouping/fold count"
                )
            # Training remains unweighted to match the LOGAN analysis. No random validation split inside groups.
            model = HistGradientBoostingClassifier(
                max_iter=a.max_iter,
                max_leaf_nodes=63,
                learning_rate=0.1,
                early_stopping=False,
                random_state=a.seed,
            )
            model.fit(X[tr], y[tr])
            oof[te] = model.predict_proba(X[te])[:, 1]
            rows.append(
                dict(
                    features=name,
                    fold=str(fold),
                    n=len(te),
                    auc=roc_auc_score(y[te], oof[te], sample_weight=weights[te]),
                    ap=average_precision_score(
                        y[te], oof[te], sample_weight=weights[te]
                    ),
                )
            )
        rows.append(
            dict(
                features=name,
                fold="pooled",
                n=len(d),
                auc=roc_auc_score(y, oof, sample_weight=weights),
                ap=average_precision_score(y, oof, sample_weight=weights),
            )
        )
        pred[name] = oof
        bins = pd.cut(
            d.length, [0, 50, 90, 100, 150, 200, 300, 500, np.inf], right=False
        )
        for lb in bins.unique():
            ix = (bins == lb).to_numpy()
            if len(np.unique(y[ix])) == 2:
                rows.append(
                    dict(
                        features=name,
                        fold="length:" + str(lb),
                        n=ix.sum(),
                        auc=roc_auc_score(y[ix], oof[ix], sample_weight=weights[ix]),
                        ap=average_precision_score(
                            y[ix], oof[ix], sample_weight=weights[ix]
                        ),
                    )
                )
    save(pd.DataFrame(rows), a.out / "metrics.tsv")
    pred.to_parquet(a.out / "oof.parquet", index=False)
    (a.out / "sequence_qc.json").write_text(
        json.dumps(
            {
                "n": len(d),
                "unknown_fraction_mean": float(unknown.mean()),
                "unknown_fraction_max": float(unknown.max()),
                "positive_fraction_unweighted": float(y.mean()),
                "groups": int(d[group].nunique()),
                "early_stopping": False,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
