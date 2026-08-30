"""Phase 16 — Mutagenèse dirigée in silico (ordres RATISS Mère).

Applique les corrections sur les vraies séquences Exp. 115 :
- SIRT6 : re-humanisation du loop L7-L8 (résidus 124/127/131 → K/G/V)
- Prestin : restauration Lys742 (ancrage PIP2 du domaine STAS)
- CIRBP : construction de la chimère CIRBP_Hybrid_Safe
  (RRM chauve-souris 1-95 + queue IDP humaine 96+)

Chaque séquence corrigée est re-passée au garde-fou topologique.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bio114.mutagenesis import (apply_mutations, build_cirbp_hybrid_safe,
                                SIRT6_REHUMANIZATION, PRESTIN_STAS_REPAIR,
                                CIRBP_RRM_END)
from bio114.ncbi import protein_to_dna_proxy
from bio114.guard import TopologicalGuard
from bio114.sequences import HUMAN_TP53, BAT_IMMUNE_IFIT

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"
REAL = ARTIFACTS / "real_sequences"
CORRECTED = ARTIFACTS / "corrected_sequences"


def read_fasta(path: Path) -> str:
    return "".join(l for l in path.read_text().splitlines()
                   if not l.startswith(">"))


def run() -> dict:
    CORRECTED.mkdir(parents=True, exist_ok=True)
    guard = TopologicalGuard(HUMAN_TP53, BAT_IMMUNE_IFIT)
    reports = {}

    # --- SIRT6 : re-humanisation du loop L7-L8 -------------------------
    bat_sirt6 = read_fasta(REAL / "bat_sirt6.fasta")
    r = apply_mutations(bat_sirt6, SIRT6_REHUMANIZATION, "sirt6",
                        notes="Re-humanisation loop L7-L8 — débloque k_on "
                              "(fixation enzymatique) en gardant le cœur "
                              "catalytique chauve-souris.")
    (CORRECTED / "sirt6_rehumanized.fasta").write_text(
        f">SIRT6_ReHumanized_L7L8 (Exp. 116)\n{r.corrected_sequence}\n")
    v = guard.inspect(protein_to_dna_proxy(r.corrected_sequence),
                      "sirt6_rehumanized")
    reports["sirt6"] = {**r.to_dict(), "guard": v.to_dict()}
    print(f"[SIRT6] mutations {r.to_dict()['mutations']} → "
          f"garde-fou: {v.action} (P_sig={v.psig:.4f})")

    # --- Prestin : restauration Lys742 ---------------------------------
    bat_prestin = read_fasta(REAL / "bat_prestin.fasta")
    r = apply_mutations(bat_prestin, PRESTIN_STAS_REPAIR, "prestin",
                        notes="Restauration Lys742 — ré-ancrage aux lipides "
                              "PIP2, stoppe la fuite énergétique de 18%.")
    (CORRECTED / "prestin_stas_humanized.fasta").write_text(
        f">Prestin_STAS_Humanized_K742 (Exp. 116)\n{r.corrected_sequence}\n")
    v = guard.inspect(protein_to_dna_proxy(r.corrected_sequence),
                      "prestin_stas_humanized")
    reports["prestin"] = {**r.to_dict(), "guard": v.to_dict()}
    print(f"[PRESTIN] mutation {r.to_dict()['mutations']} → "
          f"garde-fou: {v.action} (P_sig={v.psig:.4f})")

    # --- CIRBP : chimère Hybrid_Safe -----------------------------------
    bat_cirbp = read_fasta(REAL / "bat_cirbp.fasta")
    human_cirbp = read_fasta(REAL / "human_cirbp.fasta")
    chimera = build_cirbp_hybrid_safe(bat_cirbp, human_cirbp)
    (CORRECTED / "cirbp_hybrid_safe.fasta").write_text(
        f">CIRBP_Hybrid_Safe RRM_bat(1-{CIRBP_RRM_END})+IDP_human (Exp. 116)\n"
        f"{chimera}\n")
    v = guard.inspect(protein_to_dna_proxy(chimera), "cirbp_hybrid_safe")
    reports["cirbp_hybrid_safe"] = {
        "target": "cirbp", "composition":
            f"RRM chauve-souris (1-{CIRBP_RRM_END}) + queue IDP humaine",
        "length": len(chimera), "sequence": chimera,
        "guard": v.to_dict(),
        "notes": "Zéro nouvelle interface de surface (Condition 10). "
                 "Validation par Xi_IDP en Phase 17.",
    }
    print(f"[CIRBP] chimère Hybrid_Safe ({len(chimera)} aa) → "
          f"garde-fou: {v.action} (P_sig={v.psig:.4f})")

    n_ok = sum(1 for r in reports.values()
               if r["guard"]["action"] == "MAINTIEN")
    result = {
        "phase": 16,
        "title": "Mutagenèse dirigée — corrections RATISS Mère",
        "reports": reports,
        "n_validated": n_ok,
        "conclusion": (
            f"{n_ok}/3 constructions corrigées validées par le garde-fou. "
            f"Séquences déposées dans artifacts/corrected_sequences/ — "
            f"prêtes pour la validation Xi_IDP (Phase 17) et la dynamique "
            f"moléculaire (Phase 18)."),
        "phase16_success": n_ok > 0,
    }
    (ARTIFACTS / "phase16_mutagenesis.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 16 ===")
    print(out["conclusion"])
