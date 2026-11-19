#!/usr/bin/env python3
"""Select strict entropy/length/truncation-matched LOGAN shadows, then evaluate searched pairs."""

from _shared import *


def main():
    p = parser(__doc__)
    p.add_argument("--selection", required=True)
    p.add_argument(
        "--pool",
        required=True,
        help="Glob of LOGAN candidate/control scan tables or metadata-enriched shadow pool",
    )
    p.add_argument("--match-by", choices=["accession", "taxonomy"], default="accession")
    p.add_argument(
        "--hits",
        nargs="+",
        help="Completed raw or fullcov qdb.*.m8 searches for both arms",
    )
    p.add_argument(
        "--searched-queries",
        help="TSV/parquet listing ALL query_id values submitted, including no-hit queries",
    )
    p.add_argument("--boots", type=int, default=500)
    a = p.parse_args()
    start(a)
    c = norm(read(a.selection))
    c = c[c.arm == "candidate"].copy()
    pool = read(a.pool)
    if "query_id" not in pool and "qid" not in pool:
        if "orf_id" not in pool:
            require(
                pool, ["accession", "contig_n", "start", "aa_len", "stop", "strand"]
            )
            pool["end"] = pool.start + 3 * (pool.aa_len + pool.stop)
            pool["orf_id"] = (
                pool.accession
                + "_"
                + pool.contig_n.astype(str)
                + ":"
                + pool.start.astype(str)
                + "-"
                + pool.end.astype(str)
                + ":"
                + pool.strand
            )
        pool["query_id"] = "entropy_shadow|" + pool.orf_id.astype(str)
    s = norm(pool)
    if "prodigal_class" in s:
        s = s[s.prodigal_class.isin(["antisense_overlap", "frameshift_overlap"])].copy()
    elif "context" in s:
        s = s[s.context.isin(["antisense", "frameshift"])].copy()
    elif "label" in s:
        s = s[s.label.isin(["antisense_overlap", "frameshift_overlap"])].copy()
    else:
        raise ValueError(
            "Pool needs prodigal_class or context to select actual shadows"
        )
    if "trunc" not in s:
        raise ValueError(
            "Pool needs trunc: use repository candidate-scan tables or explicitly annotate truncation"
        )
    require(c, ["accession"])
    require(s, ["accession"])
    pairs = match(
        c,
        s,
        a.match_by,
        a.seed,
        prefer="accession" if a.match_by == "taxonomy" else None,
    )
    save(pairs, a.out / "pairs.tsv")
    save(c[~c.qid.isin(pairs.candidate_qid)], a.out / "unmatched_candidates.tsv")
    selected = s[s.qid.isin(pairs.comparator_qid)].copy()
    selected["query_id"] = selected.qid
    selected["arm"] = "entropy_shadow"
    if "chunk" in selected:
        selected["chunk"] = pd.to_numeric(selected.chunk)
    save(selected, a.out / "shadow_selection.tsv")
    selected.to_parquet(a.out / "shadow_selection.parquet", index=False)
    (a.out / "balance.json").write_text(
        json.dumps(
            {
                "candidates": len(c),
                "matched": len(pairs),
                "coverage": len(pairs) / len(c),
                "median_abs_entropy_difference": float(pairs.d_entropy.abs().median()),
                "median_abs_relative_length_difference": float(
                    pairs.relative_d_length.abs().median()
                ),
            },
            indent=2,
        )
    )
    if not a.hits:
        print(
            "Selection complete. Extract sequences with native LOGAN 22_extract_sequences.py, search the selected shadows, then rerun with --hits and --searched-queries."
        )
        return
    if not a.searched_queries:
        raise ValueError(
            "--searched-queries required: no-hit and unsearched are different"
        )
    searched = norm(read(a.searched_queries))
    needed = set(pairs.candidate_qid) | set(pairs.comparator_qid)
    if needed - set(searched.qid):
        raise ValueError(
            f"{len(needed-set(searched.qid))} matched queries not documented as searched"
        )
    h, seen = hitfiles(a.hits)
    check_searches(seen, {"afdb50", "afdb_swissprot", "pdb100", "cath50", "bfvd"})
    rows = []
    ci = c.set_index("qid")
    si = s.set_index("qid")
    for name, dbs in groups_of_db(seen).items():
        if not dbs:
            continue
        full = coverage(h, dbs)
        x = pairs.copy()
        x["accession"] = ci.loc[x.candidate_qid, "accession"].to_numpy()
        x["c"] = x.candidate_qid.isin(full)
        x["s"] = x.comparator_qid.isin(full)
        x["difference"] = x.c.astype(int) - x.s.astype(int)
        x["length_bin"] = pd.cut(
            ci.loc[x.candidate_qid, "length"].to_numpy(),
            [100, 150, 200, 300, 500, 1000, np.inf],
            right=False,
        ).astype(str)
        for label, g in [("all", x)] + list(x.groupby("length_bin")):
            diff, lo, hi = rate_ci(g, "difference", "accession", a.boots, a.seed)
            rows.append(
                dict(
                    db_set=name,
                    length_bin=str(label),
                    pairs=len(g),
                    candidate_rate=g.c.mean(),
                    shadow_rate=g.s.mean(),
                    difference=diff,
                    lo=lo,
                    hi=hi,
                )
            )
    save(pd.DataFrame(rows), a.out / "matched_evidence.tsv")


if __name__ == "__main__":
    main()
