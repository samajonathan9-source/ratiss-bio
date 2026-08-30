"""Phase 4 — La synthèse par petits brins : pas de monstre.

Insère le meilleur pont chauve-souris dans le génome humain virtuel,
soumet l'hybride au bain humain et vérifie qu'il ne s'effondre pas.
Compare aussi un hybride sauvage (fusion brutale) pour la démonstration.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bio114.sequences import SEQUENCES
from bio114.bridges import find_bridges, Bridge
from bio114.synthesis import synthesize
from bio114.eth import ETH_HUMAN

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"


def run() -> dict:
    human = SEQUENCES["human_tp53"]
    bat = SEQUENCES["bat_immune_ifit"]

    # 1. Hybride guidé : insertion au meilleur pont de réciprocité
    bridges = find_bridges(human, bat, window=120, bath=ETH_HUMAN, top_k=3)
    guided = synthesize(human, bridges[0], ETH_HUMAN)

    # 2. Hybride sauvage : fusion brutale bout-à-bout (pour comparaison)
    wild_bridge = Bridge(bat_fragment=bat, human_context="", position=len(human) // 2,
                         topo_compat=0.0, bat_psig=0.0, human_local_psig=0.0,
                         functional_gain=0.0)
    wild_bridge.insertion_point = len(human) // 2
    wild_seq = human[:len(human) // 2] + bat + human[len(human) // 2:]
    from bio114.topology import sequence_to_point_cloud, compute_signature
    from bio114.eth import apply_bath_to_cloud, sequence_stability, thermodynamic_collapse
    cloud = apply_bath_to_cloud(sequence_to_point_cloud(wild_seq), ETH_HUMAN)
    wild_sig = compute_signature(cloud, "hybrid_wild")
    wild_stab = sequence_stability(wild_seq, ETH_HUMAN)

    print(f"[hybride guidé ] P_sig={guided.psig:.4f} stabilité={guided.stability_in_human_bath:.3f} "
          f"thermo_ok={guided.thermodynamic_ok}")
    print(f"[hybride sauvage] P_sig={wild_sig.p_sig:.4f} stabilité={wild_stab['stability']:.3f} "
          f"thermo_ok={not thermodynamic_collapse(wild_stab['stability'])}")

    result = {
        "phase": 4,
        "title": "La synthèse par petits brins",
        "guided_hybrid": guided.to_dict(),
        "wild_hybrid": {
            "psig": round(wild_sig.p_sig, 4),
            "stability_in_human_bath": round(wild_stab["stability"], 4),
            "thermodynamic_ok": not thermodynamic_collapse(wild_stab["stability"]),
        },
        "conclusion": (
            f"L'hybride guidé (P_sig={guided.psig:.4f}, stabilité "
            f"{guided.stability_in_human_bath:.3f}) "
            + ("survit au bain humain — synthèse viable. "
               if guided.thermodynamic_ok else "s'effondre — rejet. ")
            + f"L'hybride sauvage sert de témoin (P_sig={wild_sig.p_sig:.4f})."
        ),
    }

    ARTIFACTS.mkdir(exist_ok=True)
    (ARTIFACTS / "phase4_synthesis.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))

    # sauvegarder la séquence hybride guidée pour la Phase 5
    (ARTIFACTS / "guided_hybrid_seq.txt").write_text(guided.sequence)
    (ARTIFACTS / "wild_hybrid_seq.txt").write_text(wild_seq)
    return result


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 4 ===")
    print(out["conclusion"])
