"""Phase 10 — Panneau multi-gènes : longévité, écholocation, hibernation.

L'Exp. 114 a prouvé le protocole sur le gène immunitaire. Cette phase
étend le pipeline complet (ponts → synthèse guidée → garde-fou) aux trois
autres pouvoirs de la chauve-souris :

- SIRT6 (longévité / réparation ADN)
- Prestin (écholocation — moteur cochléaire ultra-rapide)
- CIRBP (hibernation — résistance au froid / choc métabolique)

Pour chaque gène : pont de réciprocité → hybride guidé → verdict garde-fou.
Un gène est 'transférable' si son hybride guidé est maintenu par le
garde-fou. Sinon : repliement KTN:Li (comme l'hybride sauvage).
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bio114.sequences import SEQUENCES, HUMAN_TP53
from bio114.bridges import find_bridges
from bio114.synthesis import synthesize
from bio114.guard import TopologicalGuard

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"
FIGURES = Path(__file__).resolve().parent.parent / "figures"

# pouvoirs à transférer : (gène chauve-souris, rôle)
POWERS = {
    "bat_longevity_sirt6": "longévité / réparation ADN",
    "bat_echolocation_prestin": "écholocation (moteur cochléaire)",
    "bat_hibernation_cirbp": "hibernation / choc métabolique",
}


def run() -> dict:
    guard = TopologicalGuard(HUMAN_TP53, SEQUENCES["bat_immune_ifit"])

    results = {}
    for gene, role in POWERS.items():
        bat_seq = SEQUENCES[gene]
        # ponts de réciprocité avec le contexte humain TP53
        bridges = find_bridges(HUMAN_TP53, bat_seq, top_k=1)
        if not bridges:
            results[gene] = {"role": role, "transferable": False,
                             "reason": "aucun pont de réciprocité trouvé"}
            print(f"[{gene}] 🛑 aucun pont")
            continue

        best = bridges[0]
        hybrid = synthesize(HUMAN_TP53, best)
        verdict = guard.inspect(hybrid.sequence, f"hybride_{gene}")
        transferable = verdict.action == "MAINTIEN"

        results[gene] = {
            "role": role,
            "bridge_compat": best.to_dict()["topo_compat"],
            "hybrid": hybrid.to_dict(),
            "guard": verdict.to_dict(),
            "transferable": transferable,
        }
        mark = "✅" if transferable else "🛑"
        print(f"[{gene}] {mark} compat={best.topo_compat:.3f} "
              f"P_sig={verdict.psig:.4f} → {verdict.action} ({role})")

    transferable_powers = [g for g, r in results.items() if r.get("transferable")]
    result = {
        "phase": 10,
        "title": "Panneau multi-gènes — pouvoirs de la chauve-souris",
        "results": results,
        "transferable_powers": transferable_powers,
        "conclusion": (
            f"{len(transferable_powers)}/{len(POWERS)} pouvoirs transférables "
            f"par synthèse guidée : "
            + (", ".join(POWERS[g] for g in transferable_powers)
               if transferable_powers else "aucun")
            + ". Le garde-fou topologique filtre chaque gène individuellement."),
        "phase10_success": len(transferable_powers) > 0,
    }

    (ARTIFACTS / "phase10_genes_panel.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))
    _make_figure(results)
    return result


def _make_figure(results: dict) -> None:
    FIGURES.mkdir(exist_ok=True)
    genes = list(results.keys())
    psigs = [r["guard"]["psig"] if "guard" in r else 0 for r in results.values()]
    colors = ["#2ca02c" if r.get("transferable") else "#d62728"
              for r in results.values()]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar([g.replace("bat_", "") for g in genes], psigs, color=colors)
    ax.axhline(0.003, color="orange", linestyle="--",
               label="seuil d'effondrement (0.003)")
    ax.set_ylabel("P_sig de l'hybride guidé")
    ax.set_title("Phase 10 — Transférabilité des pouvoirs chauve-souris")
    ax.tick_params(axis="x", rotation=15)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "phase10_genes_panel.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 10 ===")
    print(out["conclusion"])
