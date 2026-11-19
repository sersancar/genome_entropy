"""Shared readers/statistics for the nine Genome Entropy sensitivity entry points."""

from pathlib import Path
import argparse, glob, importlib.util, json, sys, hashlib, subprocess
import numpy as np
import pandas as pd

COMMIT = "09968f5eccf6076892ea5f48feff94580999cd1b"
KEY = ["genome", "input_id", "orf_id"]
DBS = ["afdb_swissprot", "pdb100", "cath50", "bfvd"]
M8 = [
    "query",
    "target",
    "fident",
    "alnlen",
    "qlen",
    "tlen",
    "qcov",
    "tcov",
    "evalue",
    "bits",
    "taxid",
    "taxname",
    "theader",
]


def require(d, cols, name="table"):
    missing = set(cols) - set(d.columns)
    if missing:
        raise ValueError(f"{name}: missing columns {sorted(missing)}")


def read(path):
    paths = sorted(glob.glob(str(path)))
    if not paths:
        raise FileNotFoundError(path)
    frames = []
    for p in paths:
        x = (
            pd.read_parquet(p)
            if p.endswith(".parquet")
            else pd.read_csv(p, sep="\t", keep_default_na=False)
        )
        for c in KEY + ["chunk", "query_id", "accession"]:
            if c in x:
                x[c] = x[c].astype(str)
        frames.append(x)
    return pd.concat(frames, ignore_index=True)


def boolean(s):
    z = (
        s.astype(str)
        .str.lower()
        .map(
            {
                "true": True,
                "false": False,
                "1": True,
                "0": False,
                "1.0": True,
                "0.0": False,
            }
        )
    )
    if z.isna().any():
        raise ValueError(f"Invalid Boolean values: {s[z.isna()].unique()[:10]}")
    return z.astype(bool)


def unique(d, key, name):
    require(d, key, name)
    if d.duplicated(key).any():
        raise ValueError(f"{name}: duplicate keys {key}")


def qid(d, group=None):
    require(d, KEY)
    return (
        d[KEY].astype(str).agg("|".join, axis=1)
        + "|"
        + (d["group"].astype(str) if group is None else group)
    )


def module(repo, rel, name):
    p = Path(repo) / rel
    if not p.exists():
        raise FileNotFoundError(p)
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def parser(doc):
    p = argparse.ArgumentParser(description=doc)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--seed", type=int, default=102)
    return p


def start(a):
    a.out.mkdir(parents=True, exist_ok=True)
    meta = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    meta["reference_commit"] = COMMIT
    meta["argv"] = sys.argv
    if getattr(a, "repo", None) is not None:
        meta["checkout_commit"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=a.repo, text=True
        ).strip()
    (a.out / "run.json").write_text(json.dumps(meta, indent=2) + "\n")


def save(d, path):
    d.to_csv(path, sep="\t", index=False)


def entropy(seq):
    _, n = np.unique(list(seq), return_counts=True)
    p = n / n.sum() if len(n) else n
    return float(-(p * np.log2(p)).sum())


def norm(d):
    d = d.copy()
    aliases = {
        "aa_len": "length",
        "aa_length": "length",
        "protein_length": "length",
        "query_id": "qid",
        "group": "arm",
    }
    for old, new in aliases.items():
        if old in d and new not in d:
            d[new] = d[old]
    if "qid" not in d:
        if set(KEY) <= set(d) and "group" in d:
            d["qid"] = qid(d)
        else:
            raise ValueError("Need query_id/qid or GTDB genome,input_id,orf_id,group")
    for x in ["three_di_entropy", "aa_entropy"]:
        if x not in d and x + "_q4" in d:
            d[x] = pd.to_numeric(d[x + "_q4"]) / 1e4
    if "truncated" in d:
        d["trunc"] = boolean(d.truncated).astype(str)
    elif "truncated_calc" in d:
        d["trunc"] = boolean(d.truncated_calc).astype(str)
    elif "contig_edge" in d:
        d["trunc"] = boolean(d.contig_edge).astype(str)
    if "trunc" in d:
        d["trunc"] = d.trunc.astype(str)
    return d


def hitfiles(paths, tag="qdb"):
    """Headerless repository 12/13-field M8 files; filter requested query tag.
    Caller must supply raw multi-hit output for variable cutoffs or target exclusion.
    Files must be named TAG.DB.MODE.m8, including completed EMPTY outputs.
    """
    out = []
    seen = set()
    for pattern in paths:
        files = sorted(glob.glob(pattern))
        if not files:
            raise FileNotFoundError(pattern)
        for path in files:
            stem = Path(path).name.removesuffix(".gz").removesuffix(".m8")
            bits = stem.split(".")
            if len(bits) != 3:
                raise ValueError(f"Expected TAG.DB.MODE.m8: {path}")
            t, db, mode = bits
            if t != tag:
                continue
            if mode not in ["struct", "seq"]:
                raise ValueError(f"Unknown mode {mode}")
            seen.add((db, mode))
            try:
                d = pd.read_csv(
                    path, sep="\t", header=None, quoting=3, keep_default_na=False
                )
            except pd.errors.EmptyDataError:
                continue
            if d.shape[1] not in [12, 13]:
                raise ValueError(f"{path}: expected 12 or 13 columns, got {d.shape[1]}")
            d.columns = M8[: d.shape[1]]
            d["db"] = db
            d["mode"] = mode
            for c in ["qcov", "tcov", "evalue", "bits", "fident"]:
                d[c] = pd.to_numeric(d[c], errors="raise")
            if not d.qcov.between(0, 1).all() or not d.tcov.between(0, 1).all():
                raise ValueError("Coverage must be fractions")
            out.append(d)
    if not seen:
        raise ValueError(f"No search files with tag {tag}")
    return (
        pd.concat(out, ignore_index=True)
        if out
        else pd.DataFrame(columns=M8 + ["db", "mode"])
    ), seen


def coverage(h, dbs, cov=0.8, e=1e-3, mode="struct"):
    return set(
        h.loc[
            h.db.isin(dbs)
            & h["mode"].eq(mode)
            & h.qcov.ge(cov)
            & h.tcov.ge(cov)
            & h.evalue.le(e),
            "query",
        ]
    )


def groups_of_db(seen):
    dbs = {d for d, m in seen}
    cur = {"pdb100", "cath50"} & dbs
    return {
        "all": dbs,
        "no_afdb50": dbs - {"afdb50"},
        "pdb_cath": cur,
        "swissprot": {"afdb_swissprot"} & dbs,
    }


def check_searches(seen, dbs, modes=("struct", "seq")):
    missing = {(d, m) for d in dbs for m in modes} - seen
    if missing:
        raise ValueError(
            f"Missing completed search files (empty is valid): {sorted(missing)}"
        )


def rate_ci(d, value, cluster, boots=500, seed=102, weight=None):
    require(d, [value, cluster])
    w = np.ones(len(d)) if weight is None else d[weight].to_numpy(float)
    if not len(d) or w.sum() == 0:
        return np.nan, np.nan, np.nan
    p = float(np.average(d[value], weights=w))
    g = (
        d.assign(_num=d[value].to_numpy(float) * w, _den=w)
        .groupby(cluster)[["_num", "_den"]]
        .sum()
    )
    if len(g) < 2 or boots < 1:
        return p, np.nan, np.nan
    rng = np.random.default_rng(seed)
    v = g.to_numpy()
    samples = []
    for _ in range(boots):
        x = v[rng.integers(len(v), size=len(v))].sum(axis=0)
        samples.append(x[0] / x[1])
    return p, *np.quantile(samples, [0.025, 0.975])


def match(
    c, p, group, seed=102, entropy_caliper=0.025, length_caliper=0.1, prefer=None
):
    """Strict relative length +/-10%, no widening, no replacement; entropy-indexed groups."""
    require(c, ["qid", "length", "three_di_entropy", "trunc", group])
    require(p, ["qid", "length", "three_di_entropy", "trunc", group])
    unique(c, ["qid"], "candidates")
    unique(p, ["qid"], "pool")
    if set(c.qid) & set(p.qid):
        raise ValueError("Candidate and control pool overlap")
    p = p.reset_index(drop=True)
    rng = np.random.default_rng(seed)
    p = p.iloc[rng.permutation(len(p))].reset_index(drop=True)
    indices = {}
    for key, g in p.groupby([group, "trunc"], dropna=False):
        ix = g.sort_values("three_di_entropy", kind="stable").index.to_numpy()
        indices[key] = (ix, p.loc[ix, "three_di_entropy"].to_numpy(float))
    used = set()
    rows = []
    order = c.sample(frac=1, random_state=seed).sort_values(
        "three_di_entropy", ascending=False, kind="stable"
    )
    for _, x in order.iterrows():
        if pd.isna(x[group]) or str(x[group]) == "":
            continue
        ix, ent = indices.get(
            (x[group], x.trunc), (np.array([], dtype=int), np.array([]))
        )
        ix = ix[
            np.searchsorted(
                ent, x.three_di_entropy - entropy_caliper
            ) : np.searchsorted(ent, x.three_di_entropy + entropy_caliper, side="right")
        ]
        ix = [
            j
            for j in ix
            if j not in used
            and abs(p.loc[j, "length"] / x.length - 1) <= length_caliper
        ]
        if prefer and ix:
            same = [j for j in ix if p.loc[j, prefer] == x[prefer]]
            if same:
                ix = same
        if not ix:
            continue
        j = min(
            ix,
            key=lambda j: abs(p.loc[j, "three_di_entropy"] - x.three_di_entropy)
            / entropy_caliper
            + abs(np.log(p.loc[j, "length"] / x.length)) / length_caliper,
        )
        used.add(j)
        y = p.loc[j]
        rows.append(
            {
                "candidate_qid": x.qid,
                "comparator_qid": y.qid,
                "candidate_genome": x.get("genome", x.get("accession", "")),
                "comparator_genome": y.get("genome", y.get("accession", "")),
                "d_entropy": y.three_di_entropy - x.three_di_entropy,
                "relative_d_length": y.length / x.length - 1,
                "match_group": x[group],
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "candidate_qid",
            "comparator_qid",
            "candidate_genome",
            "comparator_genome",
            "d_entropy",
            "relative_d_length",
            "match_group",
        ],
    )


def gtdb_args(p):
    p.add_argument("--repo", required=True, type=Path)
    p.add_argument("--wanted", required=True)
    p.add_argument("--context", required=True)
    p.add_argument("--comparator", required=True)
    p.add_argument(
        "--search-dirs",
        required=True,
        help="Comma-separated repository best/ directories",
    )
    p.add_argument("--descriptions", required=True)


def gtdb_load(a):
    m = module(a.repo, "Validation/gtdb/scripts/22_rank_candidates.py", "rank_native")
    w = read(a.wanted)
    w["qid"] = qid(w)
    unique(w, ["qid"], "wanted")
    ctx = m.load_context(a.context, {"candidate", "shadow_hi"})
    unique(ctx, ["qid"], "context")
    require(ctx, ["truncated_calc"])
    ctx["truncated_calc"] = boolean(ctx.truncated_calc)
    best = m.load_hits(a.search_dirs.split(","), m.DBS)
    for db in m.DBS:
        for mode in ["struct", "seq"]:
            if (db, mode) not in best:
                raise ValueError(
                    f"GTDB native best input missing/empty {db},{mode}; supply completed populated baseline or adapt native widening for empty database"
                )
    wide = m.widen(ctx, best, m.load_descriptions(a.descriptions, m.DBS), m.DBS)
    tc = m.tiers(wide)
    wide["C3"] = tc["C3 C2 + any full-length structural hit"]
    wide["C7"] = tc["C7 strict combined rule"]
    c = wide[wide.group == "candidate"].copy()
    p = wide[wide.group == "shadow_hi"].copy()
    draws = read(a.comparator)
    require(draws, KEY + ["cand_genome", "cand_input_id", "cand_orf_id"])
    draws["candidate_qid"] = (
        draws[["cand_genome", "cand_input_id", "cand_orf_id"]]
        .astype(str)
        .agg("|".join, axis=1)
        + "|candidate"
    )
    draws["comparator_qid"] = qid(draws, "shadow_hi")
    unique(draws, ["candidate_qid"], "comparator draws")
    if set(draws.candidate_qid) != set(c.qid):
        raise ValueError(
            "Comparator candidate universe differs from context candidate universe"
        )
    if set(draws.comparator_qid) - set(p.qid):
        raise ValueError("Comparator context is incomplete")
    ann = w[w.group == "annotated_cds"].copy()
    for db in DBS:
        h = best[(db, "struct")]
        full = set(h.loc[h.qcov.ge(0.8) & h.tcov.ge(0.8) & h.evalue.le(1e-3), "query"])
        for d in [c, p, ann]:
            d[db] = d.qid.isin(full)
    return norm(c), norm(p), norm(ann), draws, best


def mixture(c, s, a):
    return (c - s) / (a - s) if a > s + 1e-12 else np.nan


def gtdb_summary(c, p, ann, pairs, label, weights=None, assignment="paired"):
    ci = c.set_index("qid").loc[pairs.candidate_qid].reset_index()
    si = p.set_index("qid").loc[pairs.comparator_qid].reset_index()
    cw = np.ones(len(ci))
    sw = cw
    aw = np.ones(len(ann))
    if weights is not None:
        cw = ci.genome.map(weights).fillna(0).to_numpy(float)
        sw = (
            cw
            if assignment == "paired"
            else si.genome.map(weights).fillna(0).to_numpy(float)
        )
        aw = ann.genome.map(weights).fillna(0).to_numpy(float)
    if min(cw.sum(), sw.sum(), aw.sum()) == 0:
        return []
    rows = []
    for db in DBS:
        cr = np.average(ci[db], weights=cw)
        sr = np.average(si[db], weights=sw)
        ar = np.average(ann[db], weights=aw)
        f = mixture(cr, sr, ar)
        rows.append(
            dict(
                label=label,
                metric=db,
                n_candidate=cw.sum(),
                candidate_rate=cr,
                comparator_rate=sr,
                annotated_rate=ar,
                difference=cr - sr,
                mixture=f,
                implied_count=f * cw.sum(),
            )
        )
    for k in ["C3", "C7"]:
        cr = np.average(ci[k], weights=cw)
        sr = np.average(si[k], weights=sw)
        rows.append(
            dict(
                label=label,
                metric=k,
                n_candidate=cw.sum(),
                candidate_rate=cr,
                comparator_rate=sr,
                difference=cr - sr,
                excess=(cr - sr) * cw.sum(),
            )
        )
    return rows
