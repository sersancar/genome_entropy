#!/usr/bin/env python3
"""Threshold grid from RAW multi-hit output; optional matched-control mixture estimates."""

from itertools import product
from _shared import *


def main():
    p = parser(__doc__)
    p.add_argument(
        "--queries",
        required=True,
        help="TSV/parquet qid,arm,length,score or three_di_entropy; includes every searched query",
    )
    p.add_argument(
        "--hits",
        nargs="+",
        required=True,
        help="RAW qdb.DB.MODE.m8 globs, not best/ or fullcov/",
    )
    p.add_argument(
        "--search-evalue",
        type=float,
        required=True,
        help="Maximum E-value retained by original search",
    )
    p.add_argument("--score-column", default="three_di_entropy")
    p.add_argument("--thresholds", default="2.3,2.5,2.7")
    p.add_argument("--coverages", default=".7,.8,.9")
    p.add_argument("--evalues", default="1e-5,1e-3")
    p.add_argument(
        "--pairs",
        help="TSV candidate_qid,comparator_qid; preserve one row per candidate draw",
    )
    p.add_argument(
        "--descriptions",
        help="Native GTDB description directory; enables C3/C7 at each cutoff (requires trunc and length)",
    )
    p.add_argument("--positive-arm", default="annotated_cds")
    p.add_argument("--candidate-arm", default="candidate")
    p.add_argument("--query-tag", default="qdb")
    p.add_argument(
        "--ascertainment-floor",
        type=float,
        required=True,
        help="Lowest selection score searched; prevents extrapolation below query sampling threshold",
    )
    a = p.parse_args()
    start(a)
    d = norm(read(a.queries))
    unique(d, ["qid"], "queries")
    require(d, ["arm", a.score_column])
    th = list(map(float, a.thresholds.split(",")))
    cov = list(map(float, a.coverages.split(",")))
    ev = list(map(float, a.evalues.split(",")))
    if min(th) < a.ascertainment_floor:
        raise ValueError(
            "Threshold below searched universe: extract/search new queries before lowering selection cutoff"
        )
    if max(ev) > a.search_evalue:
        raise ValueError(
            "Cannot relax E beyond the original search cutoff; rerun search"
        )
    if any("/best/" in x or "/fullcov/" in x for x in a.hits):
        raise ValueError(
            "Supply raw multi-hit search output; reduced hits are not valid for this grid"
        )
    h, seen = hitfiles(a.hits, a.query_tag)
    dbs = {db for db, m in seen}
    check_searches(seen, dbs)
    pairs = read(a.pairs) if a.pairs else None
    descriptions = {}
    if a.descriptions:
        require(d, ["length", "trunc"])
        for db in ["afdb_swissprot", "pdb100"]:
            desc = read(str(Path(a.descriptions) / (db + ".desc.tsv")))
            require(desc, ["target", "description"])
            descriptions[db] = set(
                desc.loc[desc.description.astype(str).str.len().gt(0), "target"]
            )
    rows = []
    for t, v, e in product(th, cov, ev):
        selected = d[(d.arm == a.candidate_arm) & d[a.score_column].ge(t)]
        for dbset, db in {"all": dbs, **{x: {x} for x in sorted(dbs)}}.items():
            f = coverage(h, db, v, e)
            q = coverage(h, db, v, e, mode="seq")
            r = dict(
                threshold=t,
                coverage=v,
                evalue=e,
                db=dbset,
                searched_candidates=len(selected),
                candidate_hits=selected.qid.isin(f).sum(),
                candidate_rate=selected.qid.isin(f).mean(),
                structure_only_rate=(
                    selected.qid.isin(f) & ~selected.qid.isin(q)
                ).mean(),
            )
            if pairs is not None:
                pp = pairs[pairs.candidate_qid.isin(selected.qid)]
                if set(pp.comparator_qid) - set(d.qid):
                    raise ValueError(
                        "Comparator queries absent from searched-query manifest"
                    )
                pos = d[d.arm == a.positive_arm]
                cr = pp.candidate_qid.isin(f).mean()
                sr = pp.comparator_qid.isin(f).mean()
                ar = pos.qid.isin(f).mean()
                r.update(
                    paired_candidates=len(pp),
                    paired_candidate_rate=cr,
                    comparator_rate=sr,
                    positive_rate=ar,
                    difference=cr - sr,
                    comparator_excess=len(pp) * (cr - sr),
                    mixture=mixture(cr, sr, ar),
                )
            if a.descriptions and dbset == "all":
                st = h[h["mode"].eq("struct") & h.evalue.le(e)]
                full = st[st.qcov.ge(v) & st.tcov.ge(v)]
                ndb = full.groupby("query").db.nunique()
                curated = set(
                    full.loc[
                        full.db.isin(["afdb_swissprot", "pdb100", "cath50"]), "query"
                    ]
                )
                named = set()
                # Match native convention: product attached to top-by-bits structural target per DB.
                top = st.sort_values("bits", ascending=False).drop_duplicates(
                    ["query", "db"]
                )
                for db, targets in descriptions.items():
                    named.update(
                        top.loc[
                            top.db.eq(db)
                            & top.target.isin(targets)
                            & top.qcov.ge(v)
                            & top.tcov.ge(v),
                            "query",
                        ]
                    )
                strong = set(top.loc[top.evalue.lt(1e-10), "query"])
                complete = (
                    d.trunc.astype(str).str.lower().isin(["false", "complete", "0"])
                )
                c3 = set(d.loc[complete & d.qid.isin(set(full["query"])), "qid"])
                c7 = set(
                    d.loc[
                        complete
                        & d.length.ge(100)
                        & d.qid.isin(curated)
                        & d.qid.isin(set(ndb[ndb >= 2].index))
                        & d.qid.isin(named)
                        & d.qid.isin(strong),
                        "qid",
                    ]
                )
                for key, ids in [("C3", c3), ("C7", c7)]:
                    r[key + "_candidates"] = selected.qid.isin(ids).sum()
                    if pairs is not None:
                        pp = pairs[pairs.candidate_qid.isin(selected.qid)]
                        r[key + "_paired_excess"] = int(
                            pp.candidate_qid.isin(ids).sum()
                        ) - int(pp.comparator_qid.isin(ids).sum())
            rows.append(r)
    save(pd.DataFrame(rows), a.out / "threshold_grid.tsv")
    print(
        "Counts describe supplied searched queries, not all genome ORFs. Changing the score threshold requires rematching controls to claim a fully recalibrated estimate. Search output remains subject to upstream max-seqs/prefilter sensitivity."
    )


if __name__ == "__main__":
    main()
