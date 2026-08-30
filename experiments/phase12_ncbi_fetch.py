"""Phase 12 — Extraction des séquences biologiques RÉELLES (NCBI).

Télécharge les vraies séquences protéiques RefSeq de SIRT6, Prestin
(SLC26A5) et CIRBP chez l'humain et la chauve-souris (Myotis lucifugus).
Fin du synthétique : c'est le vrai code source du vivant.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bio114.ncbi import TARGETS, fetch_target, protein_to_dna_proxy

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"
REAL = ARTIFACTS / "real_sequences"


def run() -> dict:
    REAL.mkdir(parents=True, exist_ok=True)
    results = {}
    for name, (gene, organism) in TARGETS.items():
        print(f"[NCBI] {name} ({gene} @ {organism})…", flush=True)
        t = fetch_target(name, gene, organism)
        if t["found"]:
            # FASTA de la protéine réelle
            (REAL / f"{name}.fasta").write_text(
                f">{t['header']}\n{t['sequence']}\n")
            # proxy ADN pour le moteur topologique
            dna = protein_to_dna_proxy(t["sequence"])
            (REAL / f"{name}_dna.fasta").write_text(
                f">{name} (proxy ADN de {t['accession']})\n{dna}\n")
            t["dna_proxy_length"] = len(dna)
            print(f"  ✅ {t['accession']} — {t['length']} aa "
                  f"({len(dna)} pb proxy)")
        else:
            print("  ⚠️ non trouvé sur RefSeq")
        results[name] = t

    found = [n for n, r in results.items() if r["found"]]
    pairs = all(results[n]["found"] for n in results)

    result = {
        "phase": 12,
        "title": "Extraction NCBI — vraies séquences SIRT6/Prestin/CIRBP",
        "source": "NCBI RefSeq protein (E-utilities)",
        "results": results,
        "all_found": pairs,
        "conclusion": (
            f"{len(found)}/6 séquences réelles récupérées depuis NCBI "
            f"RefSeq. Le code source du vivant est dans artifacts/real_sequences/. "
            f"Les proxies ADN sont prêts pour le moteur topologique."
            if found else "Aucune séquence trouvée — vérifier l'accès NCBI."),
        "phase12_success": bool(found),
    }
    (ARTIFACTS / "phase12_ncbi.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 12 ===")
    print(out["conclusion"])
