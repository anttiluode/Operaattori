from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

PRIMARY_Q = 1.05
ATTACKER_Q = 1.10
DRIVES = (0.5, 1.0, 2.0)
FINITE_SIGN_EPS_MV = 1.0e-5


def load_shards(root: Path) -> list[dict[str, Any]]:
    shards = []
    for path in sorted(root.rglob("*.json")):
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if obj.get("protocol") == "FINITE_DECOMPOSITION_PROTOCOL.md":
            shards.append(obj)
    return shards


def shard(shards: list[dict[str, Any]], drive: float, q: float) -> dict[str, Any]:
    rows = [
        x
        for x in shards
        if abs(float(x["drive_scale"]) - drive) < 1e-12
        and abs(float(x["metric_scale"]) - q) < 1e-12
    ]
    if len(rows) != 1:
        raise RuntimeError(f"need one shard drive={drive} q={q}; got {len(rows)}")
    return rows[0]


def percentile(values: list[float], p: float) -> float | None:
    a = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    return float(np.quantile(a, p)) if len(a) else None


def classify(records: list[dict[str, Any]], operating_points: list[dict[str, Any]]) -> dict[str, Any]:
    full_eval = [r for r in records if int(r["finite_full_sign"]) != 0]
    trans_eval = [r for r in records if int(r["finite_transport_sign"]) != 0]

    full_acc = (
        float(np.mean([int(r["tangent_full_sign"]) == int(r["finite_full_sign"]) for r in full_eval]))
        if full_eval
        else None
    )
    trans_acc = (
        float(np.mean([int(r["tangent_transport_sign"]) == int(r["finite_transport_sign"]) for r in trans_eval]))
        if trans_eval
        else None
    )

    predicted_eval = [
        r
        for r in records
        if bool(r["predicted_reversal"])
        and int(r["finite_full_sign"]) != 0
        and int(r["finite_transport_sign"]) != 0
    ]
    actual = [r for r in records if bool(r["actual_finite_reversal"])]
    tp = [r for r in actual if bool(r["predicted_reversal"])]
    precision = float(len(tp) / len(predicted_eval)) if predicted_eval else None
    recall = float(len(tp) / len(actual)) if actual else None

    replay_max = float(
        max(float(x["base_replay_relative_rms_error"]) for x in operating_points)
    )
    all_converged = bool(all(bool(r["full_all_converged"]) for r in records))

    full_abs = [float(r["full_abs_prediction_error_mV"]) for r in records]
    trans_abs = [float(r["transport_abs_prediction_error_mV"]) for r in records]
    full_rel = [r.get("full_relative_prediction_error") for r in records]
    trans_rel = [r.get("transport_relative_prediction_error") for r in records]

    return {
        "records": len(records),
        "operating_points": len(operating_points),
        "all_finite_full_solves_converged": all_converged,
        "max_base_replay_relative_rms_error": replay_max,
        "internal_replay_identity_pass": replay_max <= 1.0e-10,
        "full_sign_evaluable": len(full_eval),
        "transport_sign_evaluable": len(trans_eval),
        "full_tangent_finite_sign_accuracy": full_acc,
        "transport_tangent_finite_sign_accuracy": trans_acc,
        "predicted_reversals": int(sum(bool(r["predicted_reversal"]) for r in records)),
        "predicted_reversals_finite_evaluable": len(predicted_eval),
        "actual_finite_reversals": len(actual),
        "reversal_true_positive": len(tp),
        "predicted_reversal_precision": precision,
        "predicted_reversal_recall": recall,
        "full_abs_prediction_error_mV_median": percentile(full_abs, 0.5),
        "full_abs_prediction_error_mV_p90": percentile(full_abs, 0.9),
        "transport_abs_prediction_error_mV_median": percentile(trans_abs, 0.5),
        "transport_abs_prediction_error_mV_p90": percentile(trans_abs, 0.9),
        "full_relative_prediction_error_median": percentile(full_rel, 0.5),
        "full_relative_prediction_error_p90": percentile(full_rel, 0.9),
        "transport_relative_prediction_error_median": percentile(trans_rel, 0.5),
        "transport_relative_prediction_error_p90": percentile(trans_rel, 0.9),
    }


def subgroup(records: list[dict[str, Any]], key: str, value: Any) -> dict[str, Any]:
    rows = [r for r in records if r[key] == value]
    # no replay criterion inside subgroups; pass empty dummy operating points is not useful
    full_eval = [r for r in rows if int(r["finite_full_sign"]) != 0]
    trans_eval = [r for r in rows if int(r["finite_transport_sign"]) != 0]
    pred_eval = [
        r
        for r in rows
        if bool(r["predicted_reversal"])
        and int(r["finite_full_sign"]) != 0
        and int(r["finite_transport_sign"]) != 0
    ]
    actual = [r for r in rows if bool(r["actual_finite_reversal"])]
    tp = [r for r in actual if bool(r["predicted_reversal"])]
    return {
        "records": len(rows),
        "full_sign_accuracy": (
            float(np.mean([int(r["tangent_full_sign"]) == int(r["finite_full_sign"]) for r in full_eval]))
            if full_eval else None
        ),
        "transport_sign_accuracy": (
            float(np.mean([int(r["tangent_transport_sign"]) == int(r["finite_transport_sign"]) for r in trans_eval]))
            if trans_eval else None
        ),
        "predicted_reversals": int(sum(bool(r["predicted_reversal"]) for r in rows)),
        "actual_finite_reversals": len(actual),
        "reversal_precision": float(len(tp) / len(pred_eval)) if pred_eval else None,
        "reversal_recall": float(len(tp) / len(actual)) if actual else None,
    }


def aggregate(shards: list[dict[str, Any]], q: float) -> dict[str, Any]:
    selected = [shard(shards, drive, q) for drive in DRIVES]
    records = [r for s in selected for r in s["records"]]
    operating = [x for s in selected for x in s["operating_points"]]
    out = classify(records, operating)
    out["metric_scale"] = q
    out["by_drive"] = {
        str(d): subgroup(records, "drive_scale", d) for d in DRIVES
    }
    out["by_kind"] = {
        kind: subgroup(records, "kind", kind) for kind in ("length", "diameter")
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    shards = load_shards(args.root)
    for q in (PRIMARY_Q, ATTACKER_Q):
        for drive in DRIVES:
            shard(shards, drive, q)

    primary = aggregate(shards, PRIMARY_Q)
    attacker = aggregate(shards, ATTACKER_Q)

    def enough(v: float | None, threshold: float) -> bool:
        return v is not None and np.isfinite(v) and v >= threshold

    criteria = {
        "all_finite_full_solves_converged": bool(primary["all_finite_full_solves_converged"]),
        "internal_replay_identity_le_1e-10": bool(primary["internal_replay_identity_pass"]),
        "full_sign_accuracy_ge_0p90": enough(primary["full_tangent_finite_sign_accuracy"], 0.90),
        "transport_sign_accuracy_ge_0p90": enough(primary["transport_tangent_finite_sign_accuracy"], 0.90),
        "reversal_precision_ge_0p75": enough(primary["predicted_reversal_precision"], 0.75),
        "reversal_recall_ge_0p75": enough(primary["predicted_reversal_recall"], 0.75),
    }
    passed = all(criteria.values())

    result = {
        "object": "finite prediction audit of transport versus NMDA-feedback derivative decomposition",
        "protocol": "FINITE_DECOMPOSITION_PROTOCOL.md",
        "receipt_counts": {"shards": len(shards), "expected": 6},
        "primary_5_percent": primary,
        "curvature_attacker_10_percent": attacker,
        "primary_criteria": criteria,
        "classification": (
            "FINITE_DECOMPOSITION_PREDICTIVE"
            if passed
            else "FINITE_DECOMPOSITION_EXPLANATORY_ONLY_OR_LOCAL"
        ),
        "claim_boundary": (
            "This is an internal finite-intervention test of the released reduced model. "
            "It is not external validation in Jaxley or direct NEURON and does not establish novelty."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"classification": result["classification"], "criteria": criteria, "primary": primary, "attacker": attacker}, indent=2))


if __name__ == "__main__":
    main()
