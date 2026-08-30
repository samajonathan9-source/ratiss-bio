"""Phase 1 — L'Invariant topologique : observer avant de toucher.

Modélise l'ADN humain et chauve-souris en nuages de points 3D,
calcule les signatures P_sig de référence (les deux invariants).
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from bio114.sequences import SEQUENCES
from bio114.topology import sequence_to_point_cloud, compute_signature, signature_distance

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"
FIGURES = Path(__file__).resolve().parent.parent / "figures"


def run() -> dict:
    signatures = {}
    clouds = {}

    for name, seq in SEQUENCES.items():
        cloud = sequence_to_point_cloud(seq)
        clouds[name] = cloud
        sig = compute_signature(cloud, name)
        signatures[name] = sig
        print(f"[{name}] n={sig.n_points} pts | H1={sig.h1_count} boucles "
              f"({sig.h1_persistence:.1f}) | H2={sig.h2_count} cavités "
              f"({sig.h2_persistence:.1f}) | P_sig={sig.p_sig:.4f}")

    # Invariants de référence : moyenne par espèce
    human = [signatures[k] for k in signatures if k.startswith("human")]
    bat = [signatures[k] for k in signatures if k.startswith("bat")]
    human_psig = sum(s.p_sig for s in human) / len(human)
    bat_psig = sum(s.p_sig for s in bat) / len(bat)

    # Distance inter-espèces vs intra-espèces
    intra_human = signature_distance(human[0], human[1])
    intra_bat = signature_distance(bat[0], bat[1])
    inter = sum(signature_distance(h, b) for h in human for b in bat) / (len(human) * len(bat))

    result = {
        "phase": 1,
        "title": "L'Invariant topologique",
        "signatures": {k: s.to_dict() for k, s in signatures.items()},
        "human_invariant_psig": round(human_psig, 4),
        "bat_invariant_psig": round(bat_psig, 4),
        "intra_species_distance": {
            "human": round(intra_human, 4),
            "bat": round(intra_bat, 4),
        },
        "inter_species_distance": round(inter, 4),
        "conclusion": (
            f"Distance inter-espèces ({inter:.3f}) vs intra-humain ({intra_human:.3f}) : "
            + ("les deux espèces ont des signatures topologiques distinctes — "
               "la topologie de l'information diffère au-delà des lettres."
               if inter > max(intra_human, intra_bat) else
               "les signatures se recouvrent partiellement — ponts possibles.")
        ),
    }

    ARTIFACTS.mkdir(exist_ok=True)
    (ARTIFACTS / "phase1_invariant.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))

    _make_figures(clouds, signatures)
    return result


def _make_figures(clouds: dict, signatures: dict) -> None:
    FIGURES.mkdir(exist_ok=True)

    fig = plt.figure(figsize=(14, 10))
    for idx, name in enumerate(clouds, 1):
        ax = fig.add_subplot(2, 2, idx, projection="3d")
        pts = clouds[name]
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=pts[:, 2],
                   cmap="viridis", s=3, alpha=0.7)
        ax.set_title(f"{name} — P_sig={signatures[name].p_sig:.3f}", fontsize=9)
        ax.set_xlabel("x (Å)", fontsize=7)
        ax.set_ylabel("y (Å)", fontsize=7)
        ax.set_zlabel("z (Å)", fontsize=7)
    fig.suptitle("Phase 1 — Nuages de points ADN (double hélice 3D)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURES / "phase1_point_clouds.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    names = list(signatures)
    vals = [signatures[n].p_sig for n in names]
    colors = ["#1f77b4" if n.startswith("human") else "#d62728" for n in names]
    ax.bar(names, vals, color=colors)
    ax.set_ylabel("P_sig")
    ax.set_title("Phase 1 — Signatures topologiques de référence")
    ax.tick_params(axis="x", rotation=20)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / "phase1_psig_bar.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 1 ===")
    print(out["conclusion"])
    print(f"Invariant humain : P_sig = {out['human_invariant_psig']}")
    print(f"Invariant chauve-souris : P_sig = {out['bat_invariant_psig']}")
