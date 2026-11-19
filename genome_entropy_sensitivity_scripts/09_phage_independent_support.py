#!/usr/bin/env python3
"""Exclude same-genome/ANI-cluster targets before reducing raw viral hits; alternative-caller overlap and review panel."""

from _shared import *


def main():
    p = parser(__doc__)
    p.add_argument(
        "--selection",
        required=True,
        help="Native phage selection.parquet enriched with input_id OR canonical qid,arm,input_id,length",
    )
    p.add_argument(
        "--orfs",
        help="Native orfs.parquet; joins input_id by row when selection lacks it",
    )
    p.add_argument(
        "--hits",
        nargs="+",
        required=True,
        help="RAW multi-hit qdb.DB.MODE.m8, including completed empty files",
    )
    p.add_argument(
        "--target-provenance",
        required=True,
        help="TSV db,target,target_genome; multiple source genomes per target allowed",
    )
    p.add_argument(
        "--genome-groups",
        required=True,
        help="TSV input_id,group: 95% ANI clusters covering query AND viral source genomes",
    )
    p.add_argument("--search-evalue", required=True, type=float)
    p.add_argument(
        "--phanotate-bed",
        help="TSV input_id,start,end,strand, 0-based half-open; OR use --phanotate-output",
    )
    p.add_argument(
        "--phanotate-output",
        help="Standard PHANOTATE tabular output (# headers ignored), START STOP FRAME CONTIG SCORE; 1-based inclusive",
    )
    p.add_argument(
        "--candidate-pairs",
        help="Native/canonical replacement pairs table, retained for review alongside evidence",
    )
    p.add_argument("--boots", type=int, default=500)
    p.add_argument("--review-n", type=int, default=100)
    a = p.parse_args()
    start(a)
    d = norm(read(a.selection))
    if a.orfs:
        o = read(a.orfs)
        unique(o, ["row"], "orfs")
        d = d.merge(
            o[
                [
                    c
                    for c in ["row", "input_id", "g_start", "g_end", "strand"]
                    if c not in d or c == "row"
                ]
            ],
            on="row",
            validate="many_to_one",
        )
    if "input_id" not in d:
        raise ValueError("--orfs needed to map native selection row -> genome")
    require(d, ["qid", "arm", "input_id"])
    unique(d, ["qid"], "selection")
    if any("/best/" in x or "/fullcov/" in x for x in a.hits):
        raise ValueError(
            "Target exclusion requires raw multi-hit search output, not reduced best hits"
        )
    h, seen = hitfiles(a.hits)
    check_searches(seen, {"afdb_swissprot", "pdb100", "cath50", "bfvd", "phold"})
    if a.search_evalue < 1e-3:
        raise ValueError("Original search must retain hits through E=1e-3")
    prov = read(a.target_provenance)
    require(prov, ["db", "target", "target_genome"])
    groups = read(a.genome_groups)
    require(groups, ["input_id", "group"])
    unique(groups, ["input_id"], "genome groups")
    gm = groups.set_index("input_id")["group"]
    if set(d.input_id) - set(gm.index):
        raise ValueError("Query genomes missing from ANI groups")
    # Conservatively exclude a deduplicated target if ANY known source belongs to the query cluster.
    pg = prov.groupby(["db", "target"]).target_genome.agg(set).to_dict()
    qg = d.set_index("qid").input_id.to_dict()
    viral = h.db.isin(["bfvd", "phold"])
    same = []
    related = []
    unknown = []
    for r in h.itertuples():
        if r.query not in qg:
            raise ValueError(f"Hit query absent from selection: {r.query}")
        sources = pg.get((r.db, r.target), set())
        query = qg[r.query]
        unk = (r.db in ["bfvd", "phold"]) and (
            not sources or any(x not in gm.index for x in sources)
        )
        same.append(query in sources)
        related.append(any(x in gm.index and gm[x] == gm[query] for x in sources))
        unknown.append(unk)
    h["same_genome"] = same
    h["related"] = related
    h["unknown_provenance"] = unknown
    filtered = {
        "all": h,
        "non_viral": h[~viral],
        "exclude_same_genome_known": h[
            ~viral | (~h.same_genome & ~h.unknown_provenance)
        ],
        "exclude_ANI_cluster_known": h[~viral | (~h.related & ~h.unknown_provenance)],
    }
    rows = []
    detail = []
    for name, hh in filtered.items():
        full = coverage(hh, {db for db, m in seen})
        x = d.copy()
        x["full_struct"] = x.qid.isin(full)
        x["evidence_filter"] = name
        detail.append(x)
        for arm, g in x.groupby("arm"):
            r, lo, hi = rate_ci(g, "full_struct", "input_id", a.boots, a.seed)
            rows.append(dict(filter=name, arm=arm, n=len(g), rate=r, lo=lo, hi=hi))
    save(pd.DataFrame(rows), a.out / "independence_adjusted_rates.tsv")
    pd.concat(detail).to_parquet(a.out / "per_query.parquet", index=False)
    save(
        h.groupby(["db", "mode"])[["same_genome", "related", "unknown_provenance"]]
        .sum()
        .reset_index(),
        a.out / "provenance_exclusions.tsv",
    )
    if a.phanotate_bed and a.phanotate_output:
        raise ValueError("Choose one alternative-caller input")
    if a.phanotate_bed or a.phanotate_output:
        if a.phanotate_bed:
            b = read(a.phanotate_bed)
        else:
            b = pd.read_csv(
                a.phanotate_output,
                sep="\t",
                comment="#",
                header=None,
                names=["start", "end", "strand", "input_id", "score"],
            )
            lo = b[["start", "end"]].min(axis=1)
            hi = b[["start", "end"]].max(axis=1)
            b["start"] = lo - 1
            b["end"] = hi
        require(b, ["input_id", "start", "end", "strand"])
        require(d, ["g_start", "g_end", "strand"])
        agreement = []
        for r in d.itertuples():
            z = b[(b.input_id == r.input_id) & (b.strand == r.strand)]
            stop = (z.end == r.g_end) if r.strand == "+" else (z.start == r.g_start)
            ov = np.maximum(
                0, np.minimum(z.end, r.g_end) - np.maximum(z.start, r.g_start)
            )
            # Stop agreement PLUS >=80% mutual coding-interval coverage (start changes reported separately).
            cov = (ov / (r.g_end - r.g_start) >= 0.8) & (ov / (z.end - z.start) >= 0.8)
            agreement.append(
                dict(
                    qid=r.qid,
                    arm=r.arm,
                    same_stop=bool(stop.any()),
                    mutual80_same_stop=bool((stop & cov).any()),
                )
            )
        ag = pd.DataFrame(agreement)
        save(ag, a.out / "alternative_caller_agreement.tsv")
        save(
            ag.groupby("arm")[["same_stop", "mutual80_same_stop"]].mean().reset_index(),
            a.out / "alternative_caller_rates.tsv",
        )
    cand = pd.concat(detail)
    cand = cand[cand.arm.eq("candidate")]
    review = cand.groupby(["evidence_filter", "full_struct"], group_keys=False).apply(
        lambda g: g.sample(min(len(g), a.review_n), random_state=a.seed)
    )
    save(review, a.out / "manual_review_panel.tsv")
    if a.candidate_pairs:
        save(read(a.candidate_pairs), a.out / "replacement_pairs_for_review.tsv")
    print(
        "Known-provenance exclusion is conservative; report unknown target provenance. Rates are not biological FDRs. No new genome ANI or PHANOTATE runs are inferred from missing inputs."
    )


if __name__ == "__main__":
    main()
