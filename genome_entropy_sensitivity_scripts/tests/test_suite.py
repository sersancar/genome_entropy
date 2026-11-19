"""Small deterministic end-to-end fixtures; no GPU, network, or reference DB needed.
Run: python tests/test_suite.py --repo /path/to/genome_entropy
"""

import argparse, sys, subprocess, tempfile, json
from pathlib import Path
import pandas as pd, numpy as np

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P))
from _shared import entropy

ap = argparse.ArgumentParser()
ap.add_argument("--repo", required=True)
args = ap.parse_args()
repo = str(Path(args.repo).resolve())
with tempfile.TemporaryDirectory() as td:
    t = Path(td)

    def write(d, name):
        p = t / name
        p.parent.mkdir(parents=True, exist_ok=True)
        if name.endswith(".parquet"):
            pd.DataFrame(d).to_parquet(p, index=False)
        else:
            pd.DataFrame(d).to_csv(p, sep="\t", index=False)
        return str(p)

    def run(name, *opts, fail=False):
        cmd = [sys.executable, str(next(P.glob(name + "*.py"))), *map(str, opts)]
        r = subprocess.run(
            cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        if fail:
            assert r.returncode != 0, (name, "expected rejection")
            return
        assert r.returncode == 0, (name, r.stdout, r.stderr)
        print("PASS", name)

    # Native GTDB inputs, 4 genomes, paired shadows; all native DB/mode outputs populated.
    wanted = []
    ctx = []
    draws = []
    tax = []
    pool = []
    for g in range(4):
        genome = "G" + str(g)
        tax.append(
            dict(genome=genome, species="sp" + str(g), genus="gen" + str(g // 2))
        )
        for arm in ["candidate", "shadow_hi", "annotated_cds"]:
            for k in range(2):
                r = dict(
                    domain="bac",
                    chunk="001",
                    genome=genome,
                    input_id=genome + "_ctg",
                    orf_id=str(k),
                    group=arm,
                    aa_length=150 + k,
                    three_di_entropy=3.0 + k * 0.01,
                    protein_entropy=4.0,
                    truncated=False,
                    truncated_calc=False,
                    cds_frame_class="opposite strand",
                )
                wanted.append(r)
                if arm != "annotated_cds":
                    ctx.append(r)
                if arm == "shadow_hi":
                    pool.append(r)
                    draws.append(
                        dict(
                            **{x: r[x] for x in ["genome", "input_id", "orf_id"]},
                            cand_genome=genome,
                            cand_input_id=r["input_id"],
                            cand_orf_id=r["orf_id"],
                            match_level="same genome",
                        )
                    )
    w = write(wanted, "wanted.tsv")
    write(ctx, "context/test.context.tsv.gz")
    write(pool, "pool.tsv")
    comp = write(draws, "comparator.tsv")
    taxp = write(tax, "taxonomy.tsv")
    # write() doesn't gzip plain TSV named gz; native reader requires actual gzip.
    pd.DataFrame(ctx).to_csv(
        t / "context/test.context.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    for db in ["afdb_swissprot", "pdb100", "cath50", "bfvd"]:
        write(
            [dict(target="T", description="ATPase", cath_superfamily="x")],
            f"desc/{db}.desc.tsv",
        )
        for mode in ["struct", "seq"]:
            rows = []
            for r in wanted:
                hit = (
                    r["group"] == "annotated_cds"
                    or (mode == "struct" and r["group"] == "candidate")
                    or (r["group"] == "shadow_hi" and r["orf_id"] == "0")
                )
                if hit:
                    rows.append(
                        [
                            "|".join(
                                [r["genome"], r["input_id"], r["orf_id"], r["group"]]
                            ),
                            "T",
                            0.3,
                            150,
                            150,
                            150,
                            0.9,
                            0.9,
                            1e-20,
                            100,
                            2,
                            "Bacteria",
                        ]
                    )
            (t / "best").mkdir(exist_ok=True)
            pd.DataFrame(rows).to_csv(
                t / f"best/qdb.{db}.{mode}.m8", sep="\t", header=False, index=False
            )
    base = [
        "--repo",
        repo,
        "--wanted",
        w,
        "--context",
        t / "context",
        "--comparator",
        comp,
        "--search-dirs",
        t / "best",
        "--descriptions",
        t / "desc",
    ]
    run(
        "01",
        *base,
        "--pool",
        t / "pool.tsv",
        "--taxonomy",
        taxp,
        "--repeats",
        1,
        "--out",
        t / "01",
    )
    run("02", *base, "--boots", 10, "--top-k", 1, "--out", t / "02")
    # LOGAN/native evidence fixture.
    q = []
    seq = []
    md = []
    for i in range(24):
        arm = "candidate" if i % 2 else "prodigal"
        seqaa = ("ACDEFGHIKLMNPQRSTVWY" * 6) if i % 2 == 0 else ("AAAACCCCDDDD" * 10)
        tokens = ("ACDEFGHIKL" * 12) if i % 2 == 0 else ("ACD" * 40)
        q.append(
            dict(
                query_id="q" + str(i),
                arm=arm,
                accession="A" + str(i // 4),
                contig_n=i,
                length=120,
                trunc="complete",
                three_di_entropy=entropy(tokens),
                label=i % 2 == 0,
                group_id="A" + str(i // 4),
                dataset="test",
                weight=1.0,
                score=0.9,
            )
        )
        seq.append(dict(query_id="q" + str(i), aa=seqaa, three_di=tokens))
    selection = write(q, "selection.parquet")
    seqpath = write(seq, "seq.tsv")
    kt = write(
        [
            dict(
                contig=r["accession"] + "_" + str(r["contig_n"]),
                contig_group="Bacteria" if i % 2 else "plants",
            )
            for i, r in enumerate(q)
        ],
        "kt.tsv",
    )
    (t / "raw").mkdir()
    for db in ["afdb_swissprot", "pdb100", "cath50", "bfvd", "afdb50"]:
        for mode in ["struct", "seq"]:
            h = [
                [
                    r["query_id"],
                    "T",
                    0.4,
                    120,
                    120,
                    120,
                    0.9,
                    0.9,
                    1e-12,
                    100,
                    2,
                    "Bacteria",
                ]
                for i, r in enumerate(q)
                if i % 3 == 0
            ]
            pd.DataFrame(h).to_csv(
                t / f"raw/qdb.{db}.{mode}.m8", sep="\t", header=False, index=False
            )
    hits = str(t / "raw/qdb.*.m8")
    run(
        "03",
        "--selection",
        selection,
        "--taxonomy",
        kt,
        "--hits",
        hits,
        "--boots",
        5,
        "--out",
        t / "03",
    )
    sp = []
    for r in q:
        if r["arm"] == "candidate":
            sp.append(
                dict(
                    query_id="s" + r["query_id"],
                    accession=r["accession"],
                    length=r["length"],
                    trunc=r["trunc"],
                    three_di_entropy=r["three_di_entropy"],
                    context="antisense",
                )
            )
    poolpath = write(sp, "shadowpool.parquet")
    run("04", "--selection", selection, "--pool", poolpath, "--out", t / "04")
    pairs4 = pd.read_csv(t / "04/pairs.tsv", sep="\t")
    assert len(pairs4) == 12
    assert pairs4.comparator_qid.is_unique
    assert pairs4.d_entropy.abs().max() <= 0.025
    assert pairs4.relative_d_length.abs().max() <= 0.1
    # Missing completed searches must fail, rather than silently look like no-hits.
    run(
        "03",
        "--selection",
        selection,
        "--taxonomy",
        kt,
        "--hits",
        str(t / "raw/qdb.pdb100.*.m8"),
        "--out",
        t / "03bad",
        fail=True,
    )
    encoded = write(
        [dict(qid=r["query_id"], prostt5_three_di=r["three_di"]) for r in seq],
        "encoded.tsv",
    )
    run(
        "05",
        "--metadata",
        selection,
        "--sequences",
        seqpath,
        "--encoded",
        encoded,
        "--folds",
        3,
        "--out",
        t / "05",
    )
    run(
        "06",
        "--queries",
        selection,
        "--hits",
        hits,
        "--search-evalue",
        0.001,
        "--ascertainment-floor",
        0,
        "--thresholds",
        "1.0,2.5",
        "--out",
        t / "06",
    )
    run(
        "06",
        "--queries",
        selection,
        "--hits",
        hits,
        "--search-evalue",
        0.001,
        "--ascertainment-floor",
        2.5,
        "--thresholds",
        "2.3",
        "--out",
        t / "06bad",
        fail=True,
    )
    run(
        "08",
        "--metadata",
        selection,
        "--sequences",
        seqpath,
        "--folds",
        3,
        "--max-iter",
        5,
        "--out",
        t / "08",
    )
    # Provenance: two hits for candidate q1, related higher-scoring target excluded, remote hit retained.
    ph = []
    for i, r in enumerate(q):
        ph.append(dict(**r, input_id="V" + str(i // 4)))
    phpath = write(ph, "phage.parquet")
    prov = []
    for db in ["afdb_swissprot", "pdb100", "cath50", "bfvd", "afdb50"]:
        prov.append(dict(db=db, target="T", target_genome="remote"))
    prov.append(dict(db="phold", target="T", target_genome="remote"))
    for mode in ["struct", "seq"]:
        (t / f"raw/qdb.phold.{mode}.m8").write_text(
            (t / f"raw/qdb.bfvd.{mode}.m8").read_text()
        )
    provenance = write(prov, "provenance.tsv")
    groups = write(
        [dict(input_id="V" + str(i), group="cluster" + str(i)) for i in range(6)]
        + [dict(input_id="remote", group="remote")],
        "groups.tsv",
    )
    run(
        "09",
        "--selection",
        phpath,
        "--hits",
        hits,
        "--target-provenance",
        provenance,
        "--genome-groups",
        groups,
        "--search-evalue",
        0.001,
        "--boots",
        5,
        "--out",
        t / "09",
    )
    # Exact native CDS matcher with matching and nonmatching ORFs, plus an alternative-code record.
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO

    rec = SeqRecord(
        Seq("GCT" * 240 + "TAA" + "GCT" * 10),
        id="virus",
        name="virus",
        description="fixture",
    )
    rec.annotations["molecule_type"] = "DNA"
    rec.features = [
        SeqFeature(
            FeatureLocation(0, 363, strand=1),
            type="CDS",
            qualifiers={
                "translation": ["A" * 120],
                "transl_table": ["11"],
                "locus_tag": ["cds1"],
            },
        )
    ]
    (t / "gbk").mkdir()
    SeqIO.write([rec], t / "gbk/part.gbk", "genbank")
    (t / "json").mkdir()
    features = {}
    for k in range(3):
        features[str(k)] = {
            "location": dict(start=1, end=723 if k == 2 else 363, strand="+", frame=0),
            "metadata": dict(table_id=11, has_start_codon=False, has_stop_codon=True),
            "protein": {
                "aa_sequence": ("C" if k == 1 else "A") * (240 if k == 2 else 120)
            },
        }
    (t / "json/part.json").write_text(
        json.dumps([dict(input_id="virus", features=features)])
    )
    orfs = write(
        [
            dict(
                row=k,
                part="part",
                input_id="virus",
                orf_id=str(k),
                in_genbank=k != 1,
                three_di_entropy=float(k),
                protein_length=240 if k == 2 else 120,
            )
            for k in range(3)
        ],
        "orfs.parquet",
    )
    run(
        "07",
        "--repo",
        repo,
        "--json-dir",
        t / "json",
        "--gbk-dir",
        t / "gbk",
        "--orfs",
        orfs,
        "--out",
        t / "07",
    )
    audit = pd.read_csv(t / "07/match_audit.tsv.gz", sep="\t")
    assert audit.matched.sum() == 2
    strict = pd.read_parquet(t / "07/coverage_filtered.labels.parquet")
    assert len(strict) == 2 and strict.in_genbank.sum() == 1
    # Exact numeric regression: c=1, s=.5, a=1 -> mixture=1; C7 excess=4.
    point = pd.read_csv(t / "02/point_estimates.tsv", sep="\t")
    r = point[(point.label == "full") & (point.metric == "pdb100")].iloc[0]
    assert np.isclose(r.mixture, 1.0) and np.isclose(r.difference, 0.5)
    assert point[(point.label == "full") & (point.metric == "C7")].iloc[0].excess == 4
    # Coverage is applied before hit-set reduction: a lower-bit full-length match must survive.
    from _shared import coverage

    raw = pd.DataFrame(
        [
            dict(
                query="x",
                target="partial",
                db="bfvd",
                mode="struct",
                qcov=0.4,
                tcov=0.4,
                evalue=1e-12,
                bits=200,
            ),
            dict(
                query="x",
                target="remote_full",
                db="bfvd",
                mode="struct",
                qcov=0.9,
                tcov=0.9,
                evalue=1e-10,
                bits=100,
            ),
        ]
    )
    assert coverage(raw, {"bfvd"}) == {"x"}
    print(
        "All nine entry points passed synthetic end-to-end checks. GPU inference and full native model refits are not tested."
    )
