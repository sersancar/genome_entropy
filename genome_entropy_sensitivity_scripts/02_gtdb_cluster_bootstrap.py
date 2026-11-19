#!/usr/bin/env python3
"""Genome bootstrap retaining comparator draws; paired and own-genome uncertainty, influence exclusions."""

from _shared import *


def main():
    p = parser(__doc__)
    gtdb_args(p)
    p.add_argument("--boots", type=int, default=1000)
    p.add_argument(
        "--exclude-genomes", help="TSV genome; prespecified poorly annotated assemblies"
    )
    p.add_argument(
        "--top-k",
        type=int,
        default=13,
        help="Also run a transparently data-driven top-C7-excess influence analysis",
    )
    a = p.parse_args()
    start(a)
    c, s, ann, pairs, _ = gtdb_load(a)
    require(c, ["genome"])
    ci = c.set_index("qid")
    si = s.set_index("qid")
    effect = (
        pd.DataFrame(
            {
                "genome": ci.loc[pairs.candidate_qid, "genome"].to_numpy(),
                "delta_C7": ci.loc[pairs.candidate_qid, "C7"].astype(int).to_numpy()
                - si.loc[pairs.comparator_qid, "C7"].astype(int).to_numpy(),
            }
        )
        .groupby("genome")
        .delta_C7.sum()
        .sort_values(ascending=False)
    )
    save(effect.rename("C7_excess").reset_index(), a.out / "genome_influence.tsv")
    exclusions = {
        "full": set(),
        "without_top_C7_excess": set(effect.head(a.top_k).index),
    }
    if a.exclude_genomes:
        exclusions["without_prespecified"] = set(read(a.exclude_genomes).genome)
    allrows = []
    bootrows = []
    rng = np.random.default_rng(a.seed)
    for scenario, excluded in exclusions.items():
        cc = c[~c.genome.isin(excluded)]
        ss = s[~s.genome.isin(excluded)]
        aa = ann[~ann.genome.isin(excluded)]
        pp = pairs[
            pairs.candidate_qid.isin(cc.qid) & pairs.comparator_qid.isin(ss.qid)
        ].copy()
        if pp.empty or aa.empty:
            continue
        point = pd.DataFrame(gtdb_summary(cc, ss, aa, pp, scenario))
        allrows.append(point)
        genomes = np.union1d(
            np.union1d(cc.genome.unique(), ss.genome.unique()), aa.genome.unique()
        )
        for assignment in ["paired", "own_genome"]:
            for b in range(a.boots):
                counts = rng.multinomial(
                    len(genomes), np.ones(len(genomes)) / len(genomes)
                )
                weights = dict(zip(genomes, counts))
                for row in gtdb_summary(cc, ss, aa, pp, scenario, weights, assignment):
                    row.update(replicate=b, assignment=assignment)
                    bootrows.append(row)
    point = pd.concat(allrows, ignore_index=True)
    boot = pd.DataFrame(bootrows)
    save(point, a.out / "point_estimates.tsv")
    save(boot, a.out / "replicates.tsv.gz")
    rows = []
    for (label, ass, metric), g in boot.groupby(["label", "assignment", "metric"]):
        for field in [
            "candidate_rate",
            "comparator_rate",
            "annotated_rate",
            "difference",
            "mixture",
            "implied_count",
            "excess",
        ]:
            if field not in g:
                continue
            v = g[field].dropna()
            if len(v):
                rows.append(
                    dict(
                        scenario=label,
                        assignment=ass,
                        metric=metric,
                        statistic=field,
                        lo=v.quantile(0.025),
                        hi=v.quantile(0.975),
                        valid_replicates=len(v),
                        total_replicates=a.boots,
                    )
                )
    save(pd.DataFrame(rows), a.out / "bootstrap_intervals.tsv")
    # Conditional independent-binomial intervals for direct comparison; comparator Kish ESS.
    baseline = point[point.label == "full"]
    w = pairs.comparator_qid.value_counts()
    neff = w.sum() ** 2 / (w * w).sum()
    delta = []
    for _, r in baseline[baseline.metric.isin(DBS)].iterrows():
        cr, sr, ar, f = r.candidate_rate, r.comparator_rate, r.annotated_rate, r.mixture
        var = (
            (
                cr * (1 - cr) / len(pairs)
                + (1 - f) ** 2 * sr * (1 - sr) / neff
                + f * f * ar * (1 - ar) / len(ann)
            )
            / (ar - sr) ** 2
            if ar > sr
            else np.nan
        )
        delta.append(
            dict(
                db=r.metric,
                mixture=f,
                lo=f - 1.96 * np.sqrt(var),
                hi=f + 1.96 * np.sqrt(var),
                comparator_kish_ess=neff,
            )
        )
    save(pd.DataFrame(delta), a.out / "independent_delta_intervals.tsv")
    total = effect.sum()
    (a.out / "influence.json").write_text(
        json.dumps(
            {
                "total_C7_excess": float(total),
                "top_k_excess": float(effect.head(a.top_k).sum()),
                "top_k_share": (
                    float(effect.head(a.top_k).sum() / total) if total else None
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
