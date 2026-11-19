# Genome Entropy: nine sensitivity analyses

These are nine separate analysis entry points plus the required shared reader/statistics
module `_shared.py`. Keep them in the same directory. Analysis outputs are written
under `--out`; input annotation, search, cache and manuscript files are read-only.
Use a fresh output directory for each analysis configuration: rerunning with the
same output path can replace generated results. Optional model download/inference
can also populate the model cache outside `--out`.

## Repository inspected

https://github.com/sersancar/genome_entropy

Branch available from the fork at inspection: `main`.
Reference commit: `09968f5eccf6076892ea5f48feff94580999cd1b`.
No manuscript branch was present in the fetched branch inventory. This package uses
the deposited `Validation/` implementation, not unobserved manuscript-only changes.
Each run records the actual local checkout commit and CLI arguments in run.json
when `--repo` is supplied. Core numerical versions are not pinned: record your
installed environment (`python -m pip freeze > environment.txt`) with the results.
The original large per-ORF tables, caches and structure databases are NOT in this
package and were not available for validation here.

Inspected code:
- GTDB: 22_rank_candidates.py, 39_shadow_pool.py, 40_match_clean_comparator.py,
  41_comparator_summary.py, 43_assembly_strata_v7.py.
- LOGAN: 09_build_analysis_tables.py, 12_sample_orfs.py, 13_models.py,
  18_score_candidates.py, 21_candidates_select.py, 22_extract_sequences.py,
  23_foldseek_search.slurm, 24_candidates_evidence.py, 29_contig_taxonomy.py.
- Phage: 01_build_orf_table.py, 02_relabel.py, 03_models.py,
  06_foldseek_search.slurm, 07_evidence.py.
- Package: native CDS matcher and ProstT5ThreeDiEncoder.encode().

## Setup

Use Python >=3.10 (tested with 3.12); put the package beside or inside your checkout.

```bash
python -m pip install -r requirements.txt
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1
python sensitivity_scripts/01_gtdb_taxonomic_matching.py --help
```

For new GPU encoding (05), install the repository and its documented encoding
requirements, including torch, transformers, sentencepiece and huggingface_hub.
For native phage model refits (07), also install xgboost and the native analysis
requirements. Install optional dependencies in a dedicated environment. No GPU or
reference database is needed to run the synthetic checks:

```bash
python tests/test_suite.py --repo /path/to/genome_entropy
```

## Important statistical contracts

- Annotation/caller labels are operational labels, not biological ground truth.
- Strict matching uses entropy +/-0.025 bits, relative length +/-10% and identical
  truncation. This deliberately uses relative length, whereas native script 40
  implements its length caliper on the log scale. The original comparator is kept
  unchanged; the small difference in the new matching rule is explicit.
- Native C3/C7 and descriptions are reused in 01/02. In 06, cutoff-specific C3/C7
  use raw hits with native product-name conventions. These are evidence criteria,
  not validated gene counts. The native interpretable-product convention is a
  nonempty Swiss-Prot/PDB description; no new functional-name blacklist is invented.
- Mixture fractions are NEVER clipped to [0,1]. Values outside that range can reflect sampling variability or incompatibility
  with reference assumptions; nonpositive or near-zero reference denominators
  produce NaN.
- Bootstrap intervals address the specified genome/accession grouping conditional
  on data, selection and matching. They do not include predictor training
  uncertainty, cover every cross-group dependency or prove mixture assumptions.
- Existing GTDB summaries use best-by-bits hits. 01/02 preserve that convention to
  reproduce the baseline. LOGAN supports best *qualifying* full-length hits in
  fullcov/. For variable coverage or target exclusion, use RAW multi-hit output:
  a lower-scoring valid hit must survive when a top hit fails the new rule.
- For the new raw/fullcov readers, empty completed search files are valid no-hit
  results and required missing database/mode outputs cause failure. The native
  GTDB loader reused in 01/02 requires populated best-hit tables for every expected
  database/mode; an entirely empty GTDB arm needs a loader adaptation. A file's
  existence alone does not verify search completion. Use a complete query manifest
  and job/completion records to distinguish no-hit from unsearched queries.
- Do not call raw output exhaustive: prefiltering and --max-seqs can still censor
  hits. Record search parameters. A relaxed E cutoff requires rerunning the search.
- Run bacteria and archaea separately; do not pool database ceilings across domains.

## 01 — GTDB taxonomically restricted comparators

### Question

Does the GTDB candidate excess persist when comparator shadows are restricted to
the same genome, species or genus?

### Design

- Construct new strict matches separately within each taxonomic restriction,
  without replacement or caliper widening.
- Match truncation status, predicted 3Di entropy within 0.025 bits and relative
  length within 10%. The new relative-length rule differs slightly from the
  original log-length caliper; retain that distinction when interpreting changes.
- For each restriction and seed, compare the new and original comparator schemes
  on exactly the candidates for which restricted matching succeeds.
- Hold the annotated-CDS positive reference arm fixed. This tests the comparator
  baseline; it does not fully control taxonomic differences in the positive arm.
- Repeat matching across seeds and inspect unmatched candidates. Same-species
  matching may add little beyond same-genome matching in species representatives.
- Search and annotate newly selected controls before reporting restricted-arm
  evidence. Selecting a control is not evidence that it has no hits.

### Report

- Numbers and fractions of candidates matched under each restriction and seed.
- Entropy and relative-length differences within pairs; characteristics of the
  unmatched population, summarised from the exported metadata.
- Original and restricted comparator full-length hit rates, candidate rates and
  their absolute differences on the same candidate subset.
- Database-specific mixture estimates and C3/C7 excesses where evidence is
  available; clearly identify scenarios pending new searches.

### Interpretation

Persistent excess supports candidate enrichment under a more taxonomically
restricted negative baseline. Similar mixture estimates additionally support the
numerical estimate on the matched subset. Neither result establishes an unbiased
coding fraction for all candidates: matching coverage, the fixed positive arm and
structural-database representation remain relevant. Do not extrapolate an estimate
from the matchable subset to unmatched candidates without further assumptions.

```bash
python 01_gtdb_taxonomic_matching.py \
  --repo /path/to/genome_entropy \
  --wanted /data/matched_clean/bac/wanted_v7.tsv \
  --context /data/matched_clean/bac/context_v7 \
  --comparator /data/matched_clean/bac/comparator.tsv.gz \
  --search-dirs /data/old_search/best,/data/new_search/best \
  --descriptions /data/foldseek_db/descriptions \
  --pool '/data/matched_clean/bac/pool/bac_*.pool.tsv.gz' \
  --taxonomy /data/bac_genome_taxonomy.tsv \
  --levels genome,species,genus --repeats 3 \
  --out results/01_bac
```

Native candidate and pool IDs are `genome|input_id|orf_id|group`. Verify exact
identifier equality across all inputs, including leading zeros in numeric-looking
columns; TSV type inference must not be assumed to preserve them automatically. Taxonomy TSV requires exact `genome,species,genus`. Convert GTDB
metadata accession prefixes explicitly and validate the mapping before running;
the script intentionally does not guess whether RS_/GB_ prefixes should be removed.

Outputs: per-seed pairs, matching_coverage.tsv, unmatched_candidates.tsv.gz,
estimates.tsv, and to_search.tsv when newly selected controls lack context/evidence.
Within-species matching is often redundant because GTDB has one representative
per species; it is included as a diagnostic, not expected to increase coverage.
Original and restricted comparator estimates use IDENTICAL candidate subsets and
hold the annotated-CDS reference arm fixed. Thus this is primarily a sensitivity
of the negative baseline, not a completely taxonomically balanced mixture model.

If new controls need searching, the script writes only the original-arm reference
estimates for those scenarios. It never invents no-hit evidence. `to_search.tsv`
retains the native wanted-list columns. Use native GTDB `14b_extract_orf_seqs.py`
(`--archive --wanted --out-dir`), `15_build_query_db.py` (`--aa --3di --db-dir`),
and your native search workflow. Annotate their contexts with `20_orf_context.py`
and add new `*.context.tsv.gz` and best-hit outputs to a separate combined input
location. Rerun 01 with the enriched context/search paths and the ORIGINAL
comparator.tsv.gz. Original comparator rows remain the reference for comparison.

## 02 — Genome bootstrap and influential assemblies

### Question

How sensitive are the GTDB uncertainty intervals and estimated excesses to
within-genome dependence and to a small number of influential assemblies?

### Design

- Bootstrap genomes represented in the supplied query arms and retain comparator
  draw multiplicity. This replaces ORF-independent resampling with genome clusters.
- Evaluate two assignments: comparator draws follow the candidate genome in the
  paired analysis, or their own source genome in the alternative analysis.
- Recompute hit-rate differences, mixture fractions, implied counts and C3/C7
  excesses in each replicate.
- Compare full data with exclusion of a prespecified set of poorly annotated
  genomes and, separately, a data-driven exclusion of the largest C7 contributors.
- Keep uncertainty from matching seeds separate: script 01 varies matching, while
  script 02 conditions on the supplied comparator draws.

### Report

- Genome-bootstrap percentile intervals alongside the independent-binomial
  delta-method reference intervals.
- Valid replicate counts, the resampling assignment and the remaining paired
  candidate count in each exclusion scenario.
- Full-data and exclusion-scenario point estimates.
- Genome contributions to C7 excess and the fraction contributed by the top-k set.
  Summarise any prespecified exclusion set separately from the data-driven set.

### Interpretation

Wider intervals indicate reduced precision under genome clustering; changes in the
point estimate after exclusion indicate sensitivity to influential assemblies.
A concentrated excess can support recovery of unannotated regions rather than a
uniform rate of isolated missed genes. Genome resampling does not fully account
for protein-family or taxonomic dependence across genomes, and neither bootstrap
assignment captures every dependency introduced by cross-genome matching.

```bash
python 02_gtdb_cluster_bootstrap.py \
  --repo /path/to/genome_entropy \
  --wanted /data/matched_clean/bac/wanted_v7.tsv \
  --context /data/matched_clean/bac/context_v7 \
  --comparator /data/matched_clean/bac/comparator.tsv.gz \
  --search-dirs /data/old_search/best,/data/new_search/best \
  --descriptions /data/foldseek_db/descriptions \
  --exclude-genomes /data/prespecified_poorly_annotated.tsv \
  --boots 1000 --top-k 13 --out results/02_bac
```

`--exclude-genomes` is optional, one column `genome`. Full, prespecified exclusion
and a DATA-DRIVEN top-C7-excess exclusion are reported separately. Define the
prespecified set using annotation completeness, not structural results.

Paired bootstrap assigns a comparator draw its candidate genome's multiplier;
own_genome assigns it its source genome's multiplier, testing dependence across
cross-genome draws. Annotated positives follow their own genomes. Multiplicity of
reused shadows is preserved. Influential-genome exclusions remove pairs touching
excluded genomes on either side, so report remaining paired candidate coverage.
The resampling universe is the genomes represented in the supplied query arms.
Genome influence is reported directly for C7. Claims about the contribution to
another evidence criterion require a corresponding calculation.

Outputs: point_estimates.tsv, replicates.tsv.gz, bootstrap_intervals.tsv,
independent_delta_intervals.tsv (native-style independent-binomial reference),
genome_influence.tsv and influence.json. In 01, repeat seeds assess matching
variability; 02 intervals condition on the supplied matching, not on re-matching.

## 03 — LOGAN database/taxonomy reaggregation

### Question

Does metagenomic homology enrichment persist without AFDB50 and among ORFs on
contigs classified as bacterial or archaeal?

### Design

- Reaggregate full-length support across all databases, all except AFDB50,
  PDB100/CATH50, and AFDB/Swiss-Prot.
- Stratify by contig taxonomy, retaining eukaryotic, viral, unclassified and
  not-fetched categories separately. Classify using the query contig, not its hit.
- Retain the original query arms. Within each taxonomic stratum, standardise arm
  rates to the candidate length/truncation distribution where common support exists.
- Resample accessions for uncertainty. Use explicit population cell counts only
  when extrapolating beyond the validation sample.

### Report

- Sample sizes, accession counts, rates and bootstrap intervals for candidates
  and controls in each database set and taxonomy stratum.
- Structural support, sequence support, support from both searches and support
  unique to either search mode.
- Standardised structural-support rates and the candidate distribution covered by
  control cells; withhold estimates where required cells are missing.
- When population counts are supplied, population-weighted support rates and
  extrapolated numbers of expected matches, distinct from sample proportions.

### Interpretation

Persistence without AFDB50 and within classified prokaryotic contigs strengthens
the interpretation that enrichment is not solely driven by that database or by
eukaryotic sequences outside Prodigal's scope. It does not establish biological
coding status or quantify all missed prokaryotic genes. Classification uncertainty,
missing contigs and selection into the searched sample limit generalisation.

```bash
python 03_logan_database_taxonomy.py \
  --selection /data/logan/qc/candidates/selection.parquet \
  --taxonomy /data/logan/taxonomy/contig_taxonomy.tsv \
  --hits '/data/logan/qc/candidates/search/fullcov/qdb.*.m8' \
  --boots 1000 --out results/03_logan
```

Taxonomy file is the native output from `29_contig_taxonomy.py`. Plants, animals,
fungi and other eukaryotes are combined into Eukaryota. Viruses remain separate.
Unfetched contigs remain `not_fetched`, not `unclassified`. The original contig
screen did not necessarily cover every frameshift control: inspect that stratum.

Expected databases: afdb50, afdb_swissprot, pdb100, cath50, bfvd; struct AND seq
files for each. Output sets: all, no_afdb50, pdb_cath, swissprot. `no_afdb50`
still contains BFVD; it does not mean experimentally validated support.

Outputs: rates.tsv, per_query.parquet, length_truncation_standardised.tsv.
Stratified rates and accession-bootstrap intervals include structure, sequence,
both, structure-only, sequence-only. Standardised rates are withheld if controls
lack a candidate length/truncation cell; candidate_mass_covered reports overlap.
Weights are held fixed within each bootstrap; intervals are conditional on the
chosen standardisation distribution, not on estimating that distribution anew.

Optional `--population-counts` TSV has `length_bin,trunc,count`, using EXACT bin
labels from per_query.parquet (`[100.0, 150.0)` etc.; copy them programmatically).
Counts must cover the complete/5'-open candidate population intended for
extrapolation. This produces population_extrapolation.tsv; do not use the search
sample counts as the population. Taxonomy-stratified extrapolation would require
population taxonomy counts, which are not inferred from this selected query set.

## 04 — Entropy-matched LOGAN shadows

### Question

Do high-scoring intergenic ORFs retain excess homology support relative to shadows
with comparable length, truncation and predicted structural-token entropy?

### Design

- Match candidates to antisense or frameshift shadows with strict entropy and
  relative-length calipers and identical truncation status.
- Restrict matching to the same accession by default. Taxonomy matching is an
  alternative only when taxonomy has been supplied for both inputs.
- Keep the original controls as separate comparisons and record matching coverage
  and residual imbalance.
- Extract and search selected shadows that lack evidence. Document all submitted
  query IDs, including no-hit queries, before evaluating support.
- Compare the candidate and matched-shadow arms on the matched subset and use
  candidate-accession bootstrap intervals for paired differences.

### Report

- Pair counts, matching coverage, entropy/length differences and unmatched metadata.
- Candidate and matched-shadow full-length structural hit rates, with and without
  AFDB50 and in the other supported database sets.
- Absolute paired rate differences and accession-bootstrap intervals, overall and
  within length bins.

### Interpretation

A persistent excess supports enrichment beyond matching these particular
structural-complexity and length summaries. An attenuated contrast is compatible
with selection on those features, but does not causally partition the original
contrast: the comparator population and matchable candidate subset also change.
Low coverage limits the scope of any conclusion. The bootstrap assignment does
not eliminate dependence across accessions when taxonomy matching crosses them.

```bash
python 04_logan_entropy_matched_shadows.py \
  --selection /data/logan/qc/candidates/selection.parquet \
  --pool '/data/logan/analysis/v1/candidates/*.parquet' \
  --match-by accession --out results/04_logan
```

The pool is native `18_score_candidates.py` output: the script reconstructs
orf_id/query_id from accession,contig_n,start,aa_len,stop,strand. It accepts native
`label` or `prodigal_class` shadows. The existing sampled control pool may be thin
at high entropy. Do not relax calipers silently: report coverage, or build a larger
shadow pool from the per-ORF analysis cache. The code does NOT scan 90 billion
rows to enlarge that pool automatically.

Outputs: pairs.tsv, shadow_selection.tsv/parquet, unmatched_candidates.tsv,
balance.json. Keep original controls; these are ADDITIONAL selected shadows.
Use native LOGAN `22_extract_sequences.py --selection ... --shard ...` on the
output shadow_selection.parquet, with LOGAN_ROOT pointing at the original cache.
It retains chunk and orf_id so native extraction works. Build/search a separate
query DB with these query IDs; do not replace the original arm silently.

After searching, rerun with the SAME selection/pool/seed plus:

```bash
  --hits '/data/combined_search/fullcov/qdb.*.m8' \
  --searched-queries /data/combined_query_manifest.parquet
```

The combined manifest must list every original and shadow query_id, including
no-hit queries. Combined search globs must include ALL databases being compared.
Outputs matched_evidence.tsv with paired accession-bootstrap differences, overall
and by length. `--match-by taxonomy` requires an explicitly supplied `taxonomy`
column on BOTH metadata inputs; accession preference is retained within taxonomy.
No taxonomy is guessed from gene-hit targets.

## 05 — Alternative encoder

### Question

Are 3Di-entropy discrimination and entropy-screen selections robust to replacing
ModernProst with ProstT5?

### Design

- Build a stratified panel containing annotation/caller-positive ORFs, shadows,
  intergenic candidates and low-scoring annotated ORFs across length bins.
- Include threshold-adjacent and long sequences. Record the panel composition;
  a candidate-search selection alone is a pilot rather than a broad benchmark.
- Encode the same AA sequences using a pinned ProstT5 revision and fixed settings,
  or supply precomputed ProstT5 tokens for exactly those selected IDs.
- Compare 3Di entropy only; do not assume equivalent 12-state outputs.
- Separate reuse of the unchanged 2.5-bit cutoff from recalibration at 95%
  positive-label sensitivity using training groups only, evaluated in held-out groups.

### Report

- Paired entropies and entropy differences; within-length-bin Spearman correlation,
  ROC AUC and AP on the selected panel.
- Agreement and Jaccard overlap for the fixed entropy cutoff.
- Held-out selection decisions and training-derived thresholds for each encoder;
  summarise their agreement from the exported decisions if required.
- Discordant sequences and their length, arm and annotation context, for review.

### Interpretation

Fixed-cutoff agreement tests threshold portability; held-out recalibration tests
whether discrimination transfers when score scales differ. These are checks of
an entropy screen, not a rerun of the original multi-feature candidate classifier.
Panel AP and overlap depend on deliberate sampling and should not be presented
as population precision or catalogue-wide overlap. ProstT5 provides an alternative
predictor, but its AlphaFold-derived training means agreement is not fully
independent biological validation ([Heinzinger et al., 2024](https://doi.org/10.1093/nargab/lqae150)).

```bash
python 05_cross_encoder.py \
  --metadata /data/logan/qc/candidates/selection.parquet \
  --sequences '/data/logan/qc/candidates/seqs/*.tsv' \
  --dataset logan --repo /path/to/genome_entropy \
  --model Rostlab/ProstT5_fp16 --revision FULL_HUGGINGFACE_COMMIT \
  --device cuda --encoding-size 4000 --n-per-stratum 300 \
  --out results/05_logan
```

Supply a real pinned model revision, not the literal example value. Native phage
selection works with `--orfs tables/orfs.parquet --labels tables/labels.parquet`
and `--sequences tables/seqs.parquet`. The sequence table is joined by native row.
A GTDB wanted table can provide metadata, but its AA and ModernProst token FASTAs
must first be joined by exact query ID into qid,aa,three_di TSV/parquet.

Canonical metadata: qid,dataset,group_id,label,length,three_di_entropy, optionally
arm and score. Native aliases aa_len/aa_length/protein_length, query_id/group,
*_q4 are accepted. Groups infer from accession/genome/input_id. If labels are not
provided, native annotated/annotated_cds/prodigal arms are positive and other arms
are negative CALLER labels. A panel of only gene-caller disagreements is not a
representative coding-status benchmark. Build a broad metadata/sequence panel
including known positives, general shadows, threshold-adjacent candidates and
low-scoring annotations; a LOGAN candidate selection alone is only a pilot.
Sampling strata: dataset,label,length bin,near-2.5-bit flag,arm when present.
Keep deliberate oversampling explicit; panel AP is NOT population AP.

Inference calls the repository ProstT5 encoder on a pinned downloaded snapshot.
AA/ModernProst token lengths and stored ModernProst entropies are checked before
encoding. The script refuses token-length mismatch rather than silently trimming.
For CPU testing or precomputed output, replace the inference options with
`--encoded prostt5.tsv` (qid,prostt5_three_di).

Outputs: panel.tsv.gz, prostt5_tokens.tsv.gz (new inference), paired_results.tsv.gz,
within_length_metrics.tsv, group_heldout_thresholds.tsv.gz. The unchanged 2.5-bit
cutoff is compared directly. Recalibration uses a positive-label 95% sensitivity
cutoff learned separately for each encoder in training groups, evaluated in held-out
groups. The original multi-feature classifier is NOT applied to incompatible
ProstT5 feature vectors: overlap here means an entropy-screen overlap. Rebuilding
that classifier needs new training encodings and no equivalent 12-state output
is assumed. Structural model predictions are corroboration, not experimental truth.

## 06 — Threshold / coverage / E sensitivity

### Question

Do candidate enrichment and conditional mixture estimates remain similar across a
small prespecified range of selection and homology-support thresholds?

### Design

- Vary the candidate entropy or classifier-score threshold, mutual alignment
  coverage and search E-value cutoff.
- Use raw multi-hit search output and apply coverage/significance criteria before
  calculating support. Reduced best-hit tables are insufficient for this grid.
- Restrict selection thresholds to the population that was actually searched.
  Lower thresholds require new queries; relaxed E-values may require new searches.
- Preserve supplied pairs for a fixed-matching diagnostic. Rematch controls in a
  separate analysis if the selected candidate population changes materially.
- Enable C3/C7 calculations only with the required GTDB metadata and descriptions.

### Report

- Threshold settings, searched-query counts, structural support rates and
  structure-only support.
- When pairs and a positive reference arm are supplied, paired candidate and
  comparator rates, their difference, comparator excess and mixture fraction.
- With GTDB descriptions, C3/C7 candidate counts and paired excesses.
- Identify all quantities as conditional on the supplied searched universe;
  catalogue-wide candidate counts require the complete scan, not the query sample.

### Interpretation

Stability supports robustness within the settings and searched population tested.
Sensitivity identifies dependence on operational choices; it is not a reason to
select whichever setting gives the largest excess. This grid does not include
new uncertainty intervals, re-match controls automatically or establish unbiased
population counts from a selected query sample.

```bash
python 06_threshold_coverage_sensitivity.py \
  --queries /data/search_query_metadata.tsv \
  --hits '/data/raw_search/qdb.*.m8' \
  --search-evalue 0.001 --ascertainment-floor 2.3 \
  --score-column three_di_entropy --thresholds 2.3,2.5,2.7 \
  --coverages .7,.8,.9 --evalues 1e-5,1e-3 \
  --pairs results/01_bac/genome_seed102.pairs.tsv \
  --descriptions /data/foldseek_db/descriptions \
  --out results/06_gtdb
```

The query metadata must include ALL searched arms and their qid,arm,score/entropy,
length and truncation. Use native GTDB wanted/context rows with the correct qid.
`--query-tag` handles GTDB files not named qdb.*.m8. Supplying --descriptions adds
C3/C7 candidate counts and paired excesses at each setting. Counts are searched-
universe counts. Optional --pairs preserves matching draws; --positive-arm defaults
to annotated_cds. The command example assumes the expanded query population has already been
searched down to 2.3 bits. It is not valid for the original >=2.5-bit query set.
Also, the pairs in the example must describe the same candidate universe as the
query metadata; a restricted-subset pair file cannot represent all candidates.

For LOGAN use --positive-arm prodigal, --score-column score,
--thresholds .633,.784,.874,.963 and omit GTDB descriptions unless explicitly
preparing a separate GTDB-style criterion analysis.

CRITICAL: native GTDB candidates were searched only above 2.5 bits. To evaluate
2.3 you must extract/search the additional 2.3–2.5 population and an appropriate
control pool first. With existing outputs use --ascertainment-floor 2.5 and a grid
such as 2.5,2.6,2.7. Likewise, the LOGAN 20k high-score sample does NOT support
population candidate counts or lower thresholds. Native candidate scan tables
can give full score counts; they do not give unperformed homology searches.
Input --ascertainment-floor and --search-evalue are explicit assertions by the
operator; the file contents alone cannot establish their truth.

Raw multi-hit output is mandatory. Combining best/ and fullcov/ is insufficient.
A fixed matching is a conditional diagnostic; re-run 01/04 to rematch when changing
candidate selection materially. No entire-universe counts or new coverage are
imputed. Output: threshold_grid.tsv.

## 07 — Exact annotation/coverage/genetic-code audit

### Question

How much do annotation-match coverage and alternative genetic codes affect the
phage labels, discrimination metrics and high-scoring unmatched-ORF counts?

### Design

- Recompute exact ORF/CDS matches with the native strand, phase, overlap and
  translation-identity matcher, then measure coverage of the ORF by qualifying CDS.
- Compare all ORFs, exclusion of low-coverage matched positives, exclusion of
  alternative-code genomes, and both exclusions together.
- Exclude ambiguous low-coverage positives rather than relabelling them negative.
- Identify alternative-code genomes from the supplied annotation; the draft's 104
  genomes are an expected reference count to check, not a hard-coded exclusion list.
- Use existing OOF scores for fixed-score comparisons, or request native model
  refits for newly trained scores and candidate counts.
- Export a stratified metadata panel for manual review. Same-frame unmatched ORFs
  and replacement pairs need focused follow-up; they are not automatically verified.

### Report

- Match and label counts, qualifying-CDS coverage and identified genetic codes.
- Discrimination metrics and high-scoring unmatched counts for each scenario,
  clearly separating fixed-score results from refitted-model results.
- Whether contextual exclusions were available for the candidate count.
- Review-panel metadata and findings from subsequent manual inspection, if done.

### Interpretation

Changes identify dependence on label construction or on genomes translated under
a different code. Filtering changes the evaluation population, while refitting also
changes the learned model; these effects should be distinguished. This audit does
not implement alternative-code ORF calling or validate every proposed replacement.
The historical GTDB 30-aa/90-aa discrepancy is already explained in the draft and
is not an unresolved result to be reconciled by this phage audit.

```bash
python 07_annotation_genetic_code_audit.py \
  --repo /path/to/genome_entropy \
  --json-dir /data/phage/genbank_split_json \
  --gbk-dir /data/phage/genbank_split \
  --orfs /data/phage/tables/orfs.parquet \
  --oof /data/phage/rank_intrinsic/oof.parquet \
  --coverage .8 --audit-n 100 --out results/07_phage
```

Native JSON/GBK part basenames must match. Calls the actual native
`evaluate_orf_genbank_cds_match`, including strand, phase, overlap and translation
identity, then measures the largest qualifying CDS coverage of each ORF. This
avoids mistaking a merely overlapping contextual CDS for the annotation-matched
CDS. Any non-11 CDS translation table flags the whole genome.

Four scenarios: all, coverage_filtered, code11_only, both. Low-coverage matched
positives are EXCLUDED, never relabelled negative. Outputs corrected filtered
labels/orfs parquet tables, match_audit.tsv.gz, alternate_code_genomes.tsv,
metrics.tsv, candidate_counts.tsv and manual_review_panel.tsv. Without --refit,
classifier metrics/counts use existing OOF scores. With --refit, native 03_models.py
is called for all four scenarios with the intrinsic feature set and 5-fold genome
holdouts, producing new scores and candidate counts (requires enough genomes/classes
per fold). Native refit directories retain their own full metrics. Compare them,
not the fixed-score metrics.tsv, for model-retraining effects.

GTDB's 30-aa/90-aa minimum was already explained and fixed in the text: this
package does not reinterpret historical run minima. Candidate counts are not
unique proteins after dereplication. The review panel is metadata for manual
inspection, not automated confirmation of replacements.

## 08 — Stronger AA composition baseline

### Question

Does predicted 3Di entropy improve discrimination beyond the full amino-acid
frequency vector, rather than only beyond amino-acid entropy?

### Design

- Use the same ORFs, operational labels, grouped folds and learner for all feature
  comparisons: length; length plus AA entropy; length plus 20-AA frequencies;
  and each sequence-feature set with 3Di entropy added.
- Require AA sequences for every selected training/evaluation ORF and report
  unknown-residue content.
- Retain native population evaluation weights when available; train unweighted,
  following the LOGAN sampling analysis.
- Disable random early-stopping validation so it does not introduce an additional
  within-group validation split. Grouped outer folds remain identical across models.

### Report

- Pooled and fold-specific ROC AUC/AP for all five feature sets.
- Within-length-bin ROC AUC/AP and paired differences between relevant feature sets,
  calculated from the common-fold metrics or exported OOF predictions.
- Sample sizes, operational label prevalence, evaluation weighting, independent
  group counts and sequence-quality diagnostics.
- Distinguish descriptive fold variation from formal uncertainty intervals, which
  this script does not calculate.

### Interpretation

An improvement over the frequency-vector baseline supports added predictive
information beyond that representation of AA composition under this learner.
Strong composition performance would narrow the claim of structural-feature
superiority, rather than negate the usefulness of 3Di entropy. Because the draft
now limits its claim to AA entropy, this remains an informative extension rather
than a prerequisite for that narrower wording. The controlled learner differs
from the original phage XGBRF model, so absolute metrics are not a direct replication.

Canonical metadata columns are `qid,label,group_id,length,three_di_entropy`,
with optional `weight`. The sequence table needs `qid,aa`.

```bash
python 08_amino_acid_composition_baseline.py \
  --metadata /data/model_metadata.parquet --sequences /data/model_sequences.parquet \
  --dataset canonical --folds 5 --max-iter 300 --out results/08_baseline
```

Native LOGAN: `--dataset logan --metadata 'analysis/v1/sample/*.parquet'` accepts
native label,trunc,weight and q4 features, reconstructs ORF IDs and keeps complete
ORFs. Sequence metadata must cover EVERY ORF in that sample: Foldseek query
extracts alone are too small. Extract those sequences from original chunk caches
before running; this script refuses incomplete joins.
Native phage: `--dataset phage --metadata tables/orfs.parquet --sequences
tables/seqs.parquet --labels tables/labels.parquet` uses corrected labels and row
joins. The latter is already the repository's standard sequence table.

Identical grouped folds for all five baselines; no random split leaking accessions
across folds. HistGradientBoosting follows native LOGAN settings except early
stopping is disabled (the native random validation subset is not group-aware).
Full 20-AA frequency vector uses total sequence length as denominator; unknown
residues are reported in sequence_qc.json. Training unweighted, evaluation weighted
when native population weights exist. This is not the native phage XGBRF model;
it provides a controlled sequence-feature comparison under one fixed learner.
Outputs pooled/fold/length metrics and OOF scores. Do not claim identical old and
new AUCs without reproducing the original baseline on the same sample first.

## 09 — Selective phage corroboration

### Question

How much phage homology support persists after excluding same-genome or
close-relative viral targets, and how do the candidate ORFs agree with an
alternative phage gene caller?

### Design

- Use raw multi-hit results and a target-to-source-genome provenance table.
- Compare all databases, non-viral databases, exclusion of same-genome viral targets
  and exclusion of viral targets in the query's ANI cluster.
- Cover queries and viral target source genomes with the same ANI-grouping rule.
  For deduplicated targets, exclude a target if any known source is related.
- Exclude viral hits of unknown provenance from the known-provenance analyses and
  quantify those exclusions separately from related-target exclusions.
- Optionally compare separately generated PHANOTATE calls using same-stop and
  same-stop plus mutual-coverage criteria on consistent coordinates.
- Export a review panel and, when supplied, the replacement-pair table for focused
  examination of domains, genomic context, remote homologues and possible overlaps.

### Report

- Candidate/control support rates and genome-bootstrap intervals for each evidence
  filter, plus per-query support.
- Counts of hits with same-genome, related-genome or unknown provenance, noting
  that these flags can overlap and are hit counts rather than unique-query counts.
- If alternative calls are supplied, per-query and per-arm same-stop and
  mutual-coverage agreement.
- Findings from subsequent manual review, reported separately from the generated
  metadata panel. Protein-family holdouts are not implemented by this script.

### Interpretation

Support surviving close-relative exclusion is less dependent on direct annotation
recurrence, but remains subject to database overlap and predictor dependence.
A decrease can reflect removal of genuine conserved homologues as well as
annotation circularity. Unknown-provenance exclusions reduce ascertainment and
must not be described as disproving those ORFs. Agreement with another caller is
computational corroboration, not verified translation or evidence that an
overlapping annotated gene must be removed.

```bash
python 09_phage_independent_support.py \
  --selection /data/phage/foldseek_intrinsic/selection.parquet \
  --orfs /data/phage/tables/orfs.parquet \
  --hits '/data/phage/foldseek_intrinsic/search/raw/qdb.*.m8' \
  --target-provenance /data/viral_target_provenance.tsv \
  --genome-groups /data/query_and_target_ANI_groups.tsv \
  --search-evalue .001 --boots 1000 --out results/09_phage
```

Provenance schema: db,target,target_genome; duplicate db/target rows are permitted
for a deduplicated target with several source genomes. Genome-groups schema:
input_id,group. Must cover queries AND viral target source genomes, using the same
95%-ANI rule for all. Native query-only ANI groups are insufficient for new viral
sources: assign targets using their genome sequences before this analysis.
This package does not invent a target -> genome mapping from ambiguous headers.

RAW multi-hit input is mandatory: exclude targets before calculating support.
If ANY known source is related, a deduplicated target is conservatively excluded.
Viral hits with unknown provenance are removed from known-provenance exclusion
arms and counted explicitly. Non-viral hits are retained. Outputs all-database,
non-viral, same-genome-excluded and ANI-cluster-excluded rates, accession/genome
bootstrap intervals and per-query support. These are not biological false-discovery
rates or completely independent evidence.

Optional alternative caller input:
- --phanotate-bed TSV input_id,start,end,strand (0-based, half-open); OR
- --phanotate-output standard PHANOTATE tabular START STOP FRAME CONTIG SCORE
  (1-based inclusive, # comment lines ignored).
Run PHANOTATE separately on exactly the analysed genome sequences using the
version/options you report. The script compares same-stop and >=80% mutual
coverage + same-stop support, not merely any genomic overlap. Preserve record IDs.
--candidate-pairs can attach a replacement table for review. The manual panel is
stratified by evidence filter/support and provides a starting point for domain,
context, remote-homology and genuine-overlap assessment. It does not perform wet-lab
validation. Protein-family holdouts remain an optional follow-up, not an implemented
part of these nine primary analyses.

## Validation and scope

All nine entry points were exercised against deterministic synthetic native-format
fixtures, including exact native C3/C7 and CDS matching; disallowed relaxed cutoffs
and incomplete search sets are rejected. Explicit checks include comparator
uniqueness/calipers, known mixture/C7 values, exclusion (not negative relabelling)
of low-coverage positive matches, and preservation of valid lower-scoring hits.
GPU model inference, full-scale model refits and biological outcomes were not run.

Before interpreting manuscript results: reproduce the original native GTDB summary
using 41_comparator_summary.py, compare with the full baseline output from 02,
then investigate any difference before interpreting sensitivity estimates. On the
real dataset also inspect common support, missing taxonomy, source version,
search completion and valid bootstrap replicate counts. None of the scripts makes
an unperformed search look like a negative result.
