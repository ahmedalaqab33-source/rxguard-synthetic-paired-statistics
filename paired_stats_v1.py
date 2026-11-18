"""Paired statistical mechanics on explicitly invented binary rows only.

Not connected to benchmark hazard or reviewer labels. No dataset loader,
clinical model, reference-standard promotion or automatic adjudication.
"""
from collections import defaultdict
from math import comb
import random
import numpy as np

VERSION = "invented-paired-mechanics-0.1.0"


def validate_rows(rows, *, synthetic=False):
    if synthetic is not True or not isinstance(rows, list) or not rows:
        raise ValueError("Nonempty invented row list and explicit synthetic=True required")
    seen = set()
    for row in rows:
        if not isinstance(row, dict) or row.get("scope") != "SYNTHETIC" or row.get("invented") is not True:
            raise ValueError("Invented-row declaration required")
        for field in ("item_id", "cluster_id"):
            if not isinstance(row.get(field), str) or not row[field].strip():
                raise ValueError("Nonempty item/cluster IDs required")
        if row["item_id"] in seen:
            raise ValueError("Item IDs must be unique; methods must be paired within row")
        seen.add(row["item_id"])
        ref = row.get("reference")
        if ref is not None and (type(ref) is not int or ref not in (0, 1)):
            raise ValueError("Invented reference must be 0/1/None")
        for field in ("candidate", "comparator"):
            if type(row.get(field)) is not int or row[field] not in (0, 1):
                raise ValueError("Missing/failed/nonbinary method outputs require a prespecified policy; not negatives")
    return rows


def _stratum(rows, reference):
    if type(reference) is not int or reference not in (0, 1):
        raise ValueError("Reference stratum must be 0 or 1")
    return [r for r in rows if r["reference"] == reference]


def paired_summary(rows, *, reference=1, synthetic=False):
    validate_rows(rows, synthetic=synthetic)
    selected = _stratum(rows, reference)
    if not selected:
        raise ValueError("No evaluable rows in requested reference stratum")
    n = len(selected)
    gains = sum(r["candidate"] == 1 and r["comparator"] == 0 for r in selected)
    losses = sum(r["candidate"] == 0 and r["comparator"] == 1 for r in selected)
    return {"scope": "SYNTHETIC_MECHANICS_ONLY", "version": VERSION,
            "reference_stratum": reference, "n_total": len(rows),
            "n_unknown_reference": sum(r["reference"] is None for r in rows),
            "n_stratum": n, "n_clusters": len({r["cluster_id"] for r in selected}),
            "candidate_rate": sum(r["candidate"] for r in selected) / n,
            "comparator_rate": sum(r["comparator"] for r in selected) / n,
            "paired_difference": (gains - losses) / n,
            "gains": gains, "losses": losses,
            "clinical_interpretation": "PROHIBITED_FOR_THIS_INVENTED_FIXTURE"}


def resampling_plan(rows, *, reference=1, draws=2000, seed=20261023, synthetic=False):
    validate_rows(rows, synthetic=synthetic)
    clusters = sorted({r["cluster_id"] for r in _stratum(rows, reference)})
    if not clusters or type(draws) is not int or draws < 2 or type(seed) is not int:
        raise ValueError("Eligible clusters, >=2 draws and integer seed required")
    rng = random.Random(seed)
    return [[rng.choice(clusters) for _ in clusters] for _ in range(draws)]


def cluster_bootstrap(rows, plan, *, reference=1, synthetic=False):
    validate_rows(rows, synthetic=synthetic)
    groups = defaultdict(list)
    for row in _stratum(rows, reference):
        groups[row["cluster_id"]].append(row)
    if len(groups) < 2:
        raise ValueError("At least two eligible clusters required; no clinical adequacy implied")
    if not isinstance(plan, list) or len(plan) < 2:
        raise ValueError("At least two shared resampling replicates required")
    values = []
    for replicate in plan:
        if not isinstance(replicate, list) or len(replicate) != len(groups) or any(c not in groups for c in replicate):
            raise ValueError("Each replicate must draw K eligible cluster IDs with replacement")
        # Repeated sampled clusters contribute repeatedly; never collapse to a set.
        sampled = [row for cluster in replicate for row in groups[cluster]]
        values.append(sum(r["candidate"] - r["comparator"] for r in sampled) / len(sampled))
    bounds = np.quantile(values, [.025, .975], method="linear")
    return {"scope": "SYNTHETIC_MECHANICS_ONLY", "replicates": len(plan),
            "interval": "percentile-95-type7", "lower": float(bounds[0]), "upper": float(bounds[1]),
            "bootstrap_unit": "eligible_cluster", "estimand": "output_weighted_paired_difference",
            "values": values}


def exact_discordance_p(rows, *, reference=1, synthetic=False):
    validate_rows(rows, synthetic=synthetic)
    selected = _stratum(rows, reference)
    summary = paired_summary(rows, reference=reference, synthetic=True)
    if len(selected) != len({r["cluster_id"] for r in selected}):
        raise ValueError("Exact discordance test requires independent single-item clusters")
    n = summary["gains"] + summary["losses"]
    if not n:
        return 1.0
    k = min(summary["gains"], summary["losses"])
    return min(1.0, 2 * sum(comb(n, j) for j in range(k + 1)) / (2 ** n))


def invented_rows():
    # Labels invented solely to test arithmetic; never reviewer judgments.
    values = [("i1", "c1", 1, 1, 0), ("i2", "c1", 1, 1, 1),
              ("i3", "c2", 1, 0, 1), ("i4", "c3", 1, 1, 0),
              ("i5", "c3", 0, 1, 0), ("i6", "c4", 0, 0, 1),
              ("i7", "c4", 0, 0, 0), ("i8", "c5", None, 1, 1)]
    return [dict(scope="SYNTHETIC", invented=True, item_id=i, cluster_id=c,
                 reference=ref, candidate=a, comparator=b) for i, c, ref, a, b in values]
