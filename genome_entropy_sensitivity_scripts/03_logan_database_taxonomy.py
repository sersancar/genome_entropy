#!/usr/bin/env python3
"""Reaggregate native LOGAN search results by database set, taxonomy and length; accession bootstrap."""

from _shared import *


def main():
    p = parser(__doc__)
    p.add_argument("--selection", required=True)
    p.add_argument("--taxonomy", required=True)
    p.add_argument(
        "--hits",
        nargs="+",
        required=True,
        help="Quoted globs of native search/fullcov/qdb.*.m8",
    )
    p.add_argument("--boots", type=int, default=500)
    p.add_argument(
        "--population-counts",
        help="TSV length_bin,trunc,count: full candidate population for extrapolation",
    )
    a = p.parse_args()
    start(a)
    d = norm(read(a.selection))
    unique(d, ["qid"], "selection")
    require(d, ["arm", "accession", "contig_n", "length", "trunc"])
    tax = read(a.taxonomy)
    require(tax, ["contig", "contig_group"])
    unique(tax, ["contig"], "taxonomy")
    d["contig"] = d.accession + "_" + d.contig_n.astype(str)
    d = d.merge(
        tax[["contig", "contig_group"]], on="contig", how="left", validate="many_to_one"
    )
    # Not fetched is unknown ascertainment, not a confidently unclassified contig.
    d["taxonomy"] = d.contig_group.fillna("not_fetched")
    d["taxonomy"] = d.taxonomy.replace(
        {
            "bacterial": "Bacteria",
            "archaeal": "Archaea",
            "eukaryotic": "Eukaryota",
            "plants": "Eukaryota",
            "animals": "Eukaryota",
            "fungi": "Eukaryota",
            "other eukaryotes": "Eukaryota",
        }
    )
    d["length_bin"] = pd.cut(
        d.length, [100, 150, 200, 300, 500, 1000, np.inf], right=False
    ).astype(str)
    h, seen = hitfiles(a.hits)
    sets = groups_of_db(seen)
    check_searches(seen, {"afdb50", "afdb_swissprot", "pdb100", "cath50", "bfvd"})
    rows = []
    details = []
    for name, dbs in sets.items():
        if not dbs:
            continue
        x = d.copy()
        x["struct"] = x.qid.isin(coverage(h, dbs))
        x["seq"] = x.qid.isin(coverage(h, dbs, mode="seq"))
        x["both"] = x.struct & x.seq
        x["structure_only"] = x.struct & ~x.seq
        x["sequence_only"] = x.seq & ~x.struct
        x["db_set"] = name
        details.append(x)
        for scope, base in [("all", x)] + [
            (str(t), g) for t, g in x.groupby("taxonomy")
        ]:
            for arm, g in base.groupby("arm"):
                for field in [
                    "struct",
                    "seq",
                    "both",
                    "structure_only",
                    "sequence_only",
                ]:
                    rate, lo, hi = rate_ci(g, field, "accession", a.boots, a.seed)
                    rows.append(
                        dict(
                            db_set=name,
                            taxonomy=scope,
                            arm=arm,
                            evidence=field,
                            n=len(g),
                            accessions=g.accession.nunique(),
                            rate=rate,
                            lo=lo,
                            hi=hi,
                        )
                    )
    save(pd.DataFrame(rows), a.out / "rates.tsv")
    detail = pd.concat(details)
    detail.to_parquet(a.out / "per_query.parquet", index=False)
    # Common candidate length/truncation distribution within each taxonomic stratum.
    std = []
    for (ds, tx), g in detail.groupby(["db_set", "taxonomy"]):
        target = g[g.arm == "candidate"].groupby(["length_bin", "trunc"]).size()
        if target.empty:
            continue
        for arm, z in g.groupby("arm"):
            stats = z.groupby(["length_bin", "trunc"]).struct.agg(["mean", "size"])
            support = target.index.intersection(stats.index)
            mass = target.loc[support].sum() / target.sum()
            row = dict(
                db_set=ds,
                taxonomy=tx,
                arm=arm,
                candidate_mass_covered=mass,
                standardised_rate=np.nan,
            )
            if mass >= 1 - 1e-12:
                z = z.copy()
                cells = list(zip(z.length_bin, z.trunc))
                z["std_weight"] = [
                    target.get(k, 0) / stats.loc[k, "size"] for k in cells
                ]
                rate, lo, hi = rate_ci(
                    z, "struct", "accession", a.boots, a.seed, "std_weight"
                )
                row.update(standardised_rate=rate, lo=lo, hi=hi)
            std.append(row)
    save(pd.DataFrame(std), a.out / "length_truncation_standardised.tsv")
    if a.population_counts:
        pop = read(a.population_counts)
        require(pop, ["length_bin", "trunc", "count"])
        unique(pop, ["length_bin", "trunc"], "population counts")
        pop["trunc"] = pop.trunc.astype(str)
        output = []
        for ds, g in detail[detail.arm == "candidate"].groupby("db_set"):
            n = g.groupby(["length_bin", "trunc"]).size()
            weights = pop.set_index(["length_bin", "trunc"])["count"]
            if set(weights.index) - set(n.index):
                raise ValueError(
                    "Population cells lack candidate validation queries; cannot extrapolate"
                )
            g = g.copy()
            g["w"] = [weights.get(k, 0) / n[k] for k in zip(g.length_bin, g.trunc)]
            rate, lo, hi = rate_ci(g, "struct", "accession", a.boots, a.seed, "w")
            output.append(
                dict(
                    db_set=ds,
                    rate=rate,
                    lo=lo,
                    hi=hi,
                    population=weights.sum(),
                    expected_matches=rate * weights.sum(),
                    lo_matches=lo * weights.sum(),
                    hi_matches=hi * weights.sum(),
                )
            )
        save(pd.DataFrame(output), a.out / "population_extrapolation.tsv")


if __name__ == "__main__":
    main()
