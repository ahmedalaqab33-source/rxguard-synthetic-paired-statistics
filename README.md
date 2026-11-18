# RxGuard Synthetic Paired Statistics

Python and independent base-R arithmetic for paired binary comparisons on deliberately
invented rows. Author: **Ahmed Alaqab**. **Experimental research software; no clinical results.**

## Purpose and methods

The toolkit separates unknown reference labels from binary labels, preserves within-row
pairing, and resamples complete clusters with multiplicity. The estimand is the
**output-weighted paired difference**, not the unweighted mean of cluster means.
Percentile intervals use linear/type-7 quantiles. A two-sided exact discordance test is
restricted to independent single-item clusters. It rejects clustered items for that test.

`paired_stats_v1.py` provides validation, paired summaries, a deterministic shared
resampling plan, cluster bootstrap and the exact test. `paired_stats_v1.R` independently
constructs the same eight-row toy table and checks a supplied synthetic schedule.
The examples are separate from clinical benchmark labels, reviewer cases and study results.

## Requirements and setup

Python 3.10+ syntax; local verification used Python 3.13.16, NumPy 2.5.4 and base R 4.6.1.
R is optional for Python tests. Other runtime/platform combinations are unverified.

```text
python -m pip install -r requirements.txt
python -B run_tests.py
```

Use a separate environment if installing. The code never installs dependencies itself.
See [verification](docs/VERIFICATION.md) for actual observed results and their limits.

## Invented Python/R example

Choose a new output directory outside the repository and supply the actual Rscript path:

```text
python -B examples/python_r_parity.py --output ../paired-statistics-example --rscript /path/to/Rscript
```

The example constructs eight invented rows, generates a 2,000-replicate shared schedule,
compares all replicate differences and interval endpoints with base R, and writes a
`PARITY_SUMMARY.json` explicitly marked as invented arithmetic. No external data is read.
On Windows use Rscript.exe and quote paths containing spaces. Output logs remain external.

## Input/output contract

In-memory rows declare `scope="SYNTHETIC"`, `invented=True`, unique item IDs, cluster IDs,
binary candidate/comparator flags and reference 0/1/null. Missing method outputs are
rejected rather than silently treated as negative. Callers must explicitly pass
`synthetic=True`. This is a caller assertion, not a privacy/security boundary.

Summaries expose paired rates, gains/losses, reference-stratum denominators and unknown
reference counts. Bootstrap output includes every replicate and interval bounds.
The toy positive-stratum difference is 0.25; it is a constructed arithmetic example,
not a clinical effect, power assumption or benchmark performance estimate.

## Limitations

Passing unit/parity checks does not establish interval coverage, sample-size adequacy,
clinical validity, independent adjudication, causal validity or external validation.
Percentile intervals may be degenerate; only two eligible clusters are mechanically
required, which is not an adequacy recommendation. No clinical data adapter, patient
model, reference-standard promotion or automatic reviewer labeling is included.
No hosted-CI or independent-human replication is claimed.

## Data, attribution and maintenance

No real data or clinical labels are distributed. [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md)
and [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json) define the release and source lineage.
MIT license, copyright Ahmed Alaqab. Cite the repository and exact commit with
[CITATION.cff](CITATION.cff); no DOI or publication is claimed. Use GitHub Issues for
maintenance questions and reproducible software bug reports.
