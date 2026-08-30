"""Phase 3 — Les ponts de réciprocité : les invariants communs.

Cherche, pour chaque fenêtre du génome humain, le fragment chauve-souris
topologiquement compatible (même signature locale) mais fonctionnellement
supérieur (plus riche structurellement).
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bio114.sequences import SEQUENCES
from bio114.bridges import find_bridges
from bio114.eth import ETH_HUMAN

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"


def run() -> dict:
    # contexte humain : TP53 (réparation) — cible de la régénération
    human = SEQUENCES["human_tp53"]
    # source chauve-souris : gène immunitaire IFIT (résistance virale)
    bat = SEQUENCES["bat_immune_ifit"]

    bridges = find_bridges(human, bat, window=120, bath=ETH_HUMAN, top_k=5)

    for b in bridges:
        print(f"[pont @{b.position}] compat={b.topo_compat:.3f} "
              f"gain={b.functional_gain:.4f} "
              f"(humain P_sig={b.human_local_psig:.4f} → chauve-souris {b.bat_psig:.4f})")

    best = bridges[0] if bridges else None
    result = {
        "phase": 3,
        "title": "Les ponts de réciprocité",
        "human_context": "human_tp53 (réparation)",
        "bat_source": "bat_immune_ifit (immunité/longévité)",
        "n_bridges_found": len(bridges),
        "bridges": [b.to_dict() for b in bridges],
        "best_bridge": best.to_dict() if best else None,
        "conclusion": (
            f"{len(bridges)} ponts de réciprocité identifiés. Le meilleur "
            f"(position {best.position}) offre une compatibilité topologique de "
            f"{best.topo_compat:.2f} avec un gain fonctionnel de {best.functional_gain:.4f} — "
            "le fragment chauve-souris peut s'insérer sans briser la forme humaine."
            if best else "Aucun pont trouvé.")
    }

    ARTIFACTS.mkdir(exist_ok=True)
    (ARTIFACTS / "phase3_bridges.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 3 ===")
    print(out["conclusion"])
