"""Phase 2 — Le bain thermodynamique : l'ADN ne vit pas dans le vide.

Teste la stabilité thermodynamique de chaque séquence dans son bain
naturel (humain / chauve-souris) et dans le bain de l'autre espèce.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bio114.sequences import SEQUENCES
from bio114.eth import (ETH_HUMAN, ETH_BAT, sequence_stability,
                        thermodynamic_collapse)

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"
FIGURES = Path(__file__).resolve().parent.parent / "figures"


def run() -> dict:
    results = {}
    for name, seq in SEQUENCES.items():
        own_bath = ETH_HUMAN if name.startswith("human") else ETH_BAT
        other_bath = ETH_BAT if name.startswith("human") else ETH_HUMAN

        in_own = sequence_stability(seq, own_bath)
        in_other = sequence_stability(seq, other_bath)
        results[name] = {"own_bath": in_own, "foreign_bath": in_other}
        print(f"[{name}] stabilité dans {in_own['bath']}={in_own['stability']:.3f} | "
              f"dans {in_other['bath']}={in_other['stability']:.3f} "
              f"{'⚠️ EFFONDREMENT' if thermodynamic_collapse(in_other['stability']) else ''}")

    # clé de l'expérience : l'ADN chauve-souris reste stable dans eth_humain ?
    bat_in_human = results["bat_immune_ifit"]["foreign_bath"]["stability"]
    human_in_own = results["human_tp53"]["own_bath"]["stability"]

    result = {
        "phase": 2,
        "title": "Le bain thermodynamique virtuel (eth)",
        "stabilities": results,
        "human_stability_own": human_in_own,
        "bat_stability_in_human_bath": bat_in_human,
        "bat_survives_human_bath": not thermodynamic_collapse(bat_in_human),
        "conclusion": (
            f"L'ADN immunitaire de chauve-souris conserve une stabilité de "
            f"{bat_in_human:.3f} dans le bain humain (seuil 0.35) — "
            + ("il survit à l'environnement humain. La synthèse est envisageable."
               if not thermodynamic_collapse(bat_in_human) else
               "il s'effondre dans l'environnement humain. Rejet.")
        ),
    }

    ARTIFACTS.mkdir(exist_ok=True)
    (ARTIFACTS / "phase2_eth.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))

    _make_figure(results)
    return result


def _make_figure(results: dict) -> None:
    FIGURES.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5))
    names = list(results)
    own = [results[n]["own_bath"]["stability"] for n in names]
    foreign = [results[n]["foreign_bath"]["stability"] for n in names]
    x = range(len(names))
    ax.bar([i - 0.2 for i in x], own, width=0.4, label="bain naturel", color="#2ca02c")
    ax.bar([i + 0.2 for i in x], foreign, width=0.4, label="bain étranger", color="#ff7f0e")
    ax.axhline(0.35, color="red", linestyle="--", label="seuil d'effondrement")
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, rotation=20, ha="right")
    ax.set_ylabel("stabilité (0–1)")
    ax.set_title("Phase 2 — Stabilité thermodynamique dans les bains eth")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "phase2_eth_stability.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 2 ===")
    print(out["conclusion"])
