"""Phase 6 — Le garde-fou : guérir sans perdre son humanité.

Le garde-fou topologique inspecte les quatre structures : humain pur,
chauve-souris pur, hybride guidé, hybride sauvage. Il décide pour chacune :
MAINTIEN ou REPLIEMENT KTN:Li (arrêt d'urgence).
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bio114.sequences import SEQUENCES
from bio114.guard import TopologicalGuard
from bio114.eth import ETH_HUMAN

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"


def run() -> dict:
    guard = TopologicalGuard(
        human_ref_seq=SEQUENCES["human_tp53"],
        bat_ref_seq=SEQUENCES["bat_immune_ifit"],
        bath=ETH_HUMAN,
    )

    guided_seq = (ARTIFACTS / "guided_hybrid_seq.txt").read_text().strip()
    wild_seq = (ARTIFACTS / "wild_hybrid_seq.txt").read_text().strip()

    structures = {
        "humain_pur": SEQUENCES["human_tp53"],
        "chauve_souris_pure": SEQUENCES["bat_immune_ifit"],
        "hybride_guidé": guided_seq,
        "hybride_sauvage": wild_seq,
    }

    verdicts = {}
    for name, seq in structures.items():
        v = guard.inspect(seq, name)
        verdicts[name] = v
        emoji = "✅" if v.action == "MAINTIEN" else "🛑"
        print(f"{emoji} [{name}] P_sig={v.psig:.4f} dist_humain={v.distance_to_human:.4f} "
              f"→ {v.action}")

    result = {
        "phase": 6,
        "title": "Le garde-fou topologique",
        "human_invariant_psig": round(guard.human_sig.p_sig, 4),
        "bat_invariant_psig": round(guard.bat_sig.p_sig, 4),
        "verdicts": {k: v.to_dict() for k, v in verdicts.items()},
        "conclusion": (
            f"Garde-fou actif. Hybride guidé : {verdicts['hybride_guidé'].action} — "
            f"hybride sauvage : {verdicts['hybride_sauvage'].action}. "
            "Le système protège l'identité humaine tout en permettant la régénération."
        ),
    }

    (ARTIFACTS / "phase6_guard.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 6 ===")
    for name, v in out["verdicts"].items():
        print(f"  {name}: {v['action']}")
        print(f"    → {v['reason']}")
