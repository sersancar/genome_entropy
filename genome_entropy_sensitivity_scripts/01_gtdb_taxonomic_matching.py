#!/usr/bin/env python3
"""Strict same-genome/species/genus controls; evaluate old/new controls on identical candidates."""

import pandas as pd
from _shared import *


def main():
    p = parser(__doc__)
    gtdb_args(p)
    p.add_argument(
        "--pool", required=True, help="Glob of native 39_shadow_pool clean-pool tables"
    )
    p.add_argument(
        "--taxonomy",
        required=True,
        help="TSV genome,species,genus; exact assembly IDs, no silent prefix stripping",
    )
    p.add_argument("--levels", default="genome,species,genus")
    p.add_argument("--repeats", type=int, default=3)
    a = p.parse_args()
    start(a)
    c, s, ann, old, best = gtdb_load(a)
    pool = read(a.pool)
    pool["group"] = "shadow_hi"
    pool = norm(pool)
    require(pool, ["cds_frame_class"])
    pool = pool[
        pool.cds_frame_class.isin(["opposite strand", "same strand, frameshift"])
    ].copy()
    pool["qid"] = qid(pool, "shadow_hi")
    unique(pool, ["qid"], "clean pool")
    tax = read(a.taxonomy)
    require(tax, ["genome", "species", "genus"])
    unique(tax, ["genome"], "taxonomy")
    for col in ["species", "genus"]:
        mp = tax.set_index("genome")[col]
        c[col] = c.genome.map(mp)
        pool[col] = pool.genome.map(mp)
    for d in [c, pool]:
        d["length"] = pd.to_numeric(d.length)
        d["three_di_entropy"] = pd.to_numeric(d.three_di_entropy)
    results = []
    todo = []
    balance = []
    unmatched = []
    for level in a.levels.split(","):
        if level not in ["genome", "species", "genus"]:
            raise ValueError(level)
        for rep in range(a.repeats):
            label = f"{level}_seed{a.seed+rep}"
            pairs = match(c, pool, level, a.seed + rep)
            save(pairs, a.out / f"{label}.pairs.tsv")
            ids = set(pairs.candidate_qid)
            cc = c[c.qid.isin(ids)]
            missing = c[~c.qid.isin(ids)].copy()
            missing["restriction"] = label
            unmatched.append(missing)
            balance.append(
                dict(
                    label=label,
                    all_candidates=len(c),
                    matched=len(pairs),
                    coverage=len(pairs) / len(c),
                    median_abs_d_entropy=pairs.d_entropy.abs().median(),
                    median_abs_relative_d_length=pairs.relative_d_length.abs().median(),
                )
            )
            if pairs.empty:
                continue
            oldsub = old[old.candidate_qid.isin(ids)]
            results += gtdb_summary(cc, s, ann, oldsub, label + ":original")
            # New control evidence may not exist yet. Never treat an unsearched ORF as a no-hit.
            searched_ids = set(s.qid)
            new_ids = set(pairs.comparator_qid) - searched_ids
            if new_ids:
                t = pool[pool.qid.isin(new_ids)].copy()
                t["group"] = "shadow_hi"
                t["restriction"] = label
                todo.append(t)
                continue
            results += gtdb_summary(cc, s, ann, pairs, label + ":restricted")
    save(pd.DataFrame(balance), a.out / "matching_coverage.tsv")
    save(pd.DataFrame(results), a.out / "estimates.tsv")
    if unmatched:
        save(pd.concat(unmatched), a.out / "unmatched_candidates.tsv.gz")
    if todo:
        t = pd.concat(todo).drop_duplicates("qid")
        save(t, a.out / "to_search.tsv")
        print(
            f"{len(t)} new controls need extraction, searches and context annotation; estimates.tsv contains ORIGINAL-arm references only where restricted hits are pending. See README."
        )
    print(
        "Annotated positive arm is held fixed: taxonomic restriction controls the negative baseline, not all possible reference-group bias."
    )


if __name__ == "__main__":
    main()
