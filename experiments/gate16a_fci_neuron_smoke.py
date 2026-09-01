from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np


FCI_COMMIT = "75ad8b4d81a7f51bf888b30650c543592340db06"
MODEL_REL = Path(
    "simulating_neurons/neuron_models/human/eyal/"
    "Human_L23_PC_0603_11_937_Eyal_passive_dends_simple_soma"
)
MORPHOLOGY = "2013_03_06_cell11_1125_H41_06.asc"


def git_head(root: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        text=True,
    ).strip()


def import_builder(fci_root: Path):
    sys.path.insert(0, str(fci_root))
    path = fci_root / MODEL_REL / "get_standard_model.py"
    spec = importlib.util.spec_from_file_location("fci_cell1125_builder", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def choose_dendritic_synapse(syn_df):
    rows = []
    for i, row in syn_df.iterrows():
        seg = row["segments"]
        name = seg.sec.name()
        if ("dend" in name or "apic" in name) and "axon" not in name:
            rows.append((int(i), row))
    if not rows:
        raise RuntimeError("released model exposed no dendritic synapses")
    return rows[len(rows) // 2]


def run_one_event(cell, syn_df, *, event_ms: float = 20.0, tstop_ms: float = 120.0):
    from neuron import h

    row_index, row = choose_dendritic_synapse(syn_df)
    seg = row["segments"]
    syn = row["exc_synapses"]
    netcon = row["exc_netcons"]

    soma_v = h.Vector().record(cell.soma[0](0.5)._ref_v)
    local_v = h.Vector().record(seg._ref_v)
    i_ampa = h.Vector().record(syn._ref_i_AMPA)
    i_nmda = h.Vector().record(syn._ref_i_NMDA)
    g_ampa = h.Vector().record(syn._ref_g_AMPA)
    g_nmda = h.Vector().record(syn._ref_g_NMDA)
    time = h.Vector().record(h._ref_t)

    h.dt = 0.025
    h.finitialize(-70.0)
    h.fcurrent()
    netcon.event(float(event_ms))
    h.continuerun(float(tstop_ms))

    traces = {
        "time_ms": np.asarray(time, dtype=float),
        "soma_v_mV": np.asarray(soma_v, dtype=float),
        "local_v_mV": np.asarray(local_v, dtype=float),
        "i_AMPA_nA": np.asarray(i_ampa, dtype=float),
        "i_NMDA_nA": np.asarray(i_nmda, dtype=float),
        "g_AMPA_uS": np.asarray(g_ampa, dtype=float),
        "g_NMDA_uS": np.asarray(g_nmda, dtype=float),
    }
    return row_index, seg, syn, netcon, traces


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("results/gate16a/fci_neuron_smoke.json"),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    head = git_head(fci_root)
    if head != FCI_COMMIT:
        raise RuntimeError(f"FCI checkout is {head}, expected pinned {FCI_COMMIT}")

    model_folder = fci_root / MODEL_REL
    if not (model_folder / "morphologies" / MORPHOLOGY).exists():
        raise FileNotFoundError("cell-1125 morphology missing from pinned FCI checkout")

    compiled_dirs = [
        p for p in model_folder.iterdir()
        if p.is_dir() and p.name.startswith(("x86_64", "arm64", "aarch64"))
    ]
    if not compiled_dirs:
        raise RuntimeError("model mechanisms were not compiled with nrnivmodl")

    builder = import_builder(fci_root)
    cell, syn_df = builder.create_cell(path=str(model_folder) + "/")

    row_index, seg, syn, netcon, traces = run_one_event(cell, syn_df)

    for name, values in traces.items():
        if values.size == 0 or not np.all(np.isfinite(values)):
            raise FloatingPointError(f"non-finite or empty trace: {name}")

    t = traces["time_ms"]
    pre = t < 19.0
    post = (t >= 20.0) & (t <= 100.0)
    if not np.any(pre) or not np.any(post):
        raise RuntimeError("trace windows missing")

    soma_baseline = float(np.median(traces["soma_v_mV"][pre]))
    local_baseline = float(np.median(traces["local_v_mV"][pre]))
    soma_peak_depol = float(np.max(traces["soma_v_mV"][post]) - soma_baseline)
    local_peak_depol = float(np.max(traces["local_v_mV"][post]) - local_baseline)
    peak_abs_ampa = float(np.max(np.abs(traces["i_AMPA_nA"][post])))
    peak_abs_nmda = float(np.max(np.abs(traces["i_NMDA_nA"][post])))
    peak_g_ampa = float(np.max(traces["g_AMPA_uS"][post]))
    peak_g_nmda = float(np.max(traces["g_NMDA_uS"][post]))

    summary = {
        "gate": "16a",
        "purpose": "pinned released-model time-domain smoke test",
        "fci_commit": head,
        "model_folder": str(MODEL_REL),
        "morphology": MORPHOLOGY,
        "released_builder": "get_standard_model.py",
        "released_mechanism": "AMPANMDA_EMS.mod",
        "dendritic_synapse_row": int(row_index),
        "section": str(seg.sec.name()),
        "segment_x": float(seg.x),
        "synapse_parameters": {
            "gamma_per_mV": float(syn.gamma),
            "NMDA_ratio": float(syn.NMDA_ratio),
            "netcon_weight_uS": float(netcon.weight[0]),
            "tau_r_AMPA_ms": float(syn.tau_r_AMPA),
            "tau_d_AMPA_ms": float(syn.tau_d_AMPA),
            "tau_r_NMDA_ms": float(syn.tau_r_NMDA),
            "tau_d_NMDA_ms": float(syn.tau_d_NMDA),
        },
        "trace": {
            "samples": int(len(t)),
            "dt_ms": float(t[1] - t[0]) if len(t) > 1 else None,
            "soma_baseline_mV": soma_baseline,
            "local_baseline_mV": local_baseline,
            "soma_peak_depolarization_mV": soma_peak_depol,
            "local_peak_depolarization_mV": local_peak_depol,
            "peak_abs_AMPA_current_nA": peak_abs_ampa,
            "peak_abs_NMDA_current_nA": peak_abs_nmda,
            "peak_AMPA_conductance_uS": peak_g_ampa,
            "peak_NMDA_conductance_uS": peak_g_nmda,
        },
        "classification": "RELEASED_FCI_CELL1125_TIME_DOMAIN_RUNS",
        "stopping_line": (
            "This is only an execution/representation receipt. It earns a bounded "
            "Gate-16 dynamic clustered-vs-dispersed assay; it is not an FCI result."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("Operaattori Gate 16a — released FCI cell 1125 time-domain smoke")
    print()
    print(f"FCI commit:                 {head}")
    print(f"section / x:                {seg.sec.name()} / {seg.x:.4f}")
    print(f"released gamma:             {syn.gamma:.3f} /mV")
    print(f"released weight:            {netcon.weight[0]:.6f} uS")
    print(f"tau AMPA rise/decay:        {syn.tau_r_AMPA:.3f} / {syn.tau_d_AMPA:.3f} ms")
    print(f"tau NMDA rise/decay:        {syn.tau_r_NMDA:.3f} / {syn.tau_d_NMDA:.3f} ms")
    print(f"soma peak depolarization:   {soma_peak_depol:.6g} mV")
    print(f"local peak depolarization:  {local_peak_depol:.6g} mV")
    print(f"peak |I_AMPA|:              {peak_abs_ampa:.6g} nA")
    print(f"peak |I_NMDA|:              {peak_abs_nmda:.6g} nA")
    print(f"peak g_AMPA:                {peak_g_ampa:.6g} uS")
    print(f"peak g_NMDA:                {peak_g_nmda:.6g} uS")
    print()
    print("classification: RELEASED_FCI_CELL1125_TIME_DOMAIN_RUNS")

    assert len(t) > 1000
    assert abs(float(syn.gamma) - 0.078) < 1e-12
    assert abs(float(netcon.weight[0]) - 0.00088) < 1e-12
    assert peak_abs_ampa > 1e-8
    assert peak_abs_nmda > 1e-8
    assert peak_g_ampa > 1e-8
    assert peak_g_nmda > 1e-8
    assert local_peak_depol > 1e-5
    assert soma_peak_depol > 1e-8


if __name__ == "__main__":
    main()
