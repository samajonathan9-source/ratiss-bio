"""Phase 5 — L'œil quantique : cohérence de l'hybride.

Mesure la décohérence quantique (tunneling des protons) de l'hybride
guidé vs l'hybride sauvage. Si le QPU IBM est accessible (clé dans
secrets/), signale le pont vers le calcul quantique réel.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bio114.sequences import SEQUENCES
from bio114.decoherence import measure_decoherence, qpu_available, qpu_coupling_estimate
from bio114.eth import ETH_HUMAN

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"


def run() -> dict:
    guided_seq = (ARTIFACTS / "guided_hybrid_seq.txt").read_text().strip()
    wild_seq = (ARTIFACTS / "wild_hybrid_seq.txt").read_text().strip()

    guided = measure_decoherence(guided_seq, ETH_HUMAN.temperature, "hybride_guidé")
    wild = measure_decoherence(wild_seq, ETH_HUMAN.temperature, "hybride_sauvage")
    human_ref = measure_decoherence(SEQUENCES["human_tp53"], ETH_HUMAN.temperature, "humain_référence")

    for r in (guided, wild, human_ref):
        print(f"[{r.name}] GC={r.n_gc} AT={r.n_at} | "
              f"tunneling={r.tunneling_rate:.2e} s⁻¹ | "
              f"cohérence={r.coherence_score:.3f} | {r.verdict}")

    qpu = qpu_available()
    coupling = qpu_coupling_estimate(guided_seq)

    result = {
        "phase": 5,
        "title": "L'œil quantique — décohérence",
        "reports": [guided.to_dict(), wild.to_dict(), human_ref.to_dict()],
        "qpu_ibm_accessible": qpu,
        "guided_coupling_estimate": round(coupling, 4),
        "conclusion": (
            f"L'hybride guidé maintient une cohérence de {guided.coherence_score:.3f} "
            f"(référence humaine {human_ref.coherence_score:.3f}) — "
            + ("fusion viable au niveau quantique. " if guided.coherence_score > 0.5
               else "décohérence trop rapide — rejet. ")
            + (f"Pont QPU IBM disponible (couplage estimé {coupling:.2f})."
               if qpu else "QPU IBM non configuré — estimation analytique utilisée.")
        ),
    }

    (ARTIFACTS / "phase5_decoherence.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 5 ===")
    print(out["conclusion"])
