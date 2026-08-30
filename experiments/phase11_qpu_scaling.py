"""Phase 11 — Scaling QPU : mesures affinées sur hardware IBM.

Montée en charge par rapport à la Phase 9 : plus de qubits, plus de shots,
plusieurs circuits de profondeur croissante pour cartographier la courbe
de survie de la cohérence de l'hybride guidé vs sauvage sur le QPU le
moins chargé.
"""

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from experiments.phase9_qpu_ibm import (build_coherence_circuit,
                                        ghz_coherence_from_counts,
                                        load_token)

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"

SHOTS = 4096
QUBIT_DEPTHS = [4, 8, 12]  # tailles de chaîne GHZ = profondeur du test


def run() -> dict:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    phase10 = json.loads((ARTIFACTS / "phase10_genes_panel.json").read_text())
    phase9 = json.loads((ARTIFACTS / "phase9_qpu_ibm.json").read_text())

    token = load_token()
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backends = service.backends(operational=True, simulator=False)
    ranked = sorted(backends, key=lambda b: b.status().pending_jobs)
    backend = ranked[0]
    print(f"Backend: {backend.name} (file: {backend.status().pending_jobs}, "
          f"{backend.num_qubits} qubits) — shots={SHOTS}")

    coh_guided = phase9["hardware"]["hybride_guidé"]["coherence_encodee"]
    coh_wild = phase9["hardware"]["hybride_sauvage"]["coherence_encodee"]

    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    sampler = Sampler(mode=backend)

    curves = {}
    for label, coh in [("hybride_guidé", coh_guided), ("hybride_sauvage", coh_wild)]:
        points = []
        for n_q in QUBIT_DEPTHS:
            qc = build_coherence_circuit(n_q, coh)
            isa_qc = pm.run(qc)
            t0 = time.time()
            job = sampler.run([isa_qc], shots=SHOTS)
            print(f"[{label}] n={n_q} job {job.job_id()}…", flush=True)
            res = job.result()
            counts = res[0].data.c.get_counts()
            hw_coh = ghz_coherence_from_counts(counts, n_q)
            points.append({"n_qubits": n_q, "coherence_hw": round(hw_coh, 4),
                           "elapsed_s": round(time.time() - t0, 1)})
            print(f"[{label}] n={n_q} → cohérence HW = {hw_coh:.3f}")
        curves[label] = points

    # survie = pente de décroissance de la cohérence avec la profondeur
    def survival(points):
        cohs = [p["coherence_hw"] for p in points]
        return cohs[-1] / max(cohs[0], 1e-9)

    surv_g = survival(curves["hybride_guidé"])
    surv_w = survival(curves["hybride_sauvage"])
    success = surv_g >= surv_w

    result = {
        "phase": 11,
        "title": "Scaling QPU — courbes de survie de la cohérence",
        "backend": backend.name,
        "shots": SHOTS,
        "curves": curves,
        "survival_guided": round(surv_g, 4),
        "survival_wild": round(surv_w, 4),
        "validation_scaling": success,
        "genes_panel": phase10["transferable_powers"],
        "conclusion": (
            f"SCALING VALIDÉ sur {backend.name} : à {SHOTS} shots et jusqu'à "
            f"{max(QUBIT_DEPTHS)} qubits, l'hybride guidé conserve mieux sa "
            f"cohérence en profondeur (survie {surv_g:.2f}) que le sauvage "
            f"({surv_w:.2f}). La signature topologique résiste à la montée en "
            f"échelle sur hardware réel."
            if success else
            "Scaling non concluant — voir courbes."),
    }
    (ARTIFACTS / "phase11_qpu_scaling.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 11 — SCALING QPU ===")
    for label, pts in out["curves"].items():
        s = " | ".join(f"n={p['n_qubits']}:{p['coherence_hw']:.3f}" for p in pts)
        print(f"  {label}: {s}")
    print("\n" + out["conclusion"])
