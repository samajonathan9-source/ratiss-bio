"""Phase 14 — Stress thermique sur structures RÉELLES (robustesse).

Les vraies structures C-alpha (AlphaFold DB) sont soumises à un stress
thermique : perturbation brownienne croissante (37°C corps humain → 41°C
fièvre), corrélée au pLDDT réel de chaque résidu — les zones désordonnées
(pLDDT < 70) absorbent plus d'énergie, comme en dynamique moléculaire.

Métrique : ΔP_sig sous stress = dérive topologique réelle. Une structure
qui s'effondre sous 41°C est un candidat dangereux (dénaturation) ; une
structure stable est apte au transfert.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from bio114.topology import compute_signature
from experiments.phase13_alphafold_probe import pdb_to_ca_cloud

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"
STRUCTS = ARTIFACTS / "real_structures"

# énergie thermique relative : kT à 310K (37°C) vs 314K (41°C)
TEMP_STRESS = {"37C": 1.0, "41C": 1.8}
RNG = np.random.default_rng(114)


def thermal_perturbation(cloud: np.ndarray, conf: list[float],
                         stress: float) -> np.ndarray:
    """Bruit brownien modulé par le désordre local réel (pLDDT)."""
    conf_arr = np.array(conf)
    disorder = 1.0 - conf_arr / 100.0            # 0 = ordonné, 1 = désordonné
    sigma = 0.5 * stress * (0.3 + disorder)      # amplitude en Å
    noise = RNG.normal(0, sigma[:, None], cloud.shape)
    return cloud + noise


def run() -> dict:
    phase13 = json.loads((ARTIFACTS / "phase13_alphafold.json").read_text())
    results = {}
    for name, info in phase13["results"].items():
        if not info.get("found"):
            continue
        pdb = (STRUCTS / f"{name}.pdb").read_text()
        cloud, conf = pdb_to_ca_cloud(pdb)
        base_sig = phase13["results"][name]["psig_real"]["p_sig"]

        row = {"p_sig_37C": base_sig}
        for temp, stress in TEMP_STRESS.items():
            if temp == "37C":
                continue
            hot = thermal_perturbation(cloud, conf, stress)
            sig_hot = compute_signature(hot, f"{name}_{temp}")
            drift = abs(sig_hot.p_sig - base_sig)
            row[f"p_sig_{temp}"] = round(sig_hot.p_sig, 4)
            row["drift_41C"] = round(drift, 4)
            row["denatured"] = sig_hot.p_sig < 0.003
        row["stable_under_fever"] = not row.get("denatured", False)
        results[name] = row
        mark = "✅" if row["stable_under_fever"] else "🛑"
        print(f"[{name}] {mark} P_sig 37°C={row['p_sig_37C']:.4f} → "
              f"41°C={row.get('p_sig_41C', 0):.4f} (Δ={row.get('drift_41C', 0):.4f})")

    stable = [n for n, r in results.items() if r["stable_under_fever"]]
    result = {
        "phase": 14,
        "title": "Stress thermique sur structures réelles (37°C → 41°C)",
        "method": "bruit brownien modulé par pLDDT réel (désordre local)",
        "results": results,
        "stable_count": len(stable),
        "conclusion": (
            f"{len(stable)}/{len(results)} structures réelles restent "
            f"topologiquement stables sous fièvre (41°C). Les dérives ΔP_sig "
            f"mesurées sur les vraies coordonnées servent de tension "
            f"topologique physique pour le garde-fou."
            if results else "Aucune structure à tester."),
        "phase14_success": bool(stable),
    }
    (ARTIFACTS / "phase14_thermal_stress.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 14 ===")
    print(out["conclusion"])
