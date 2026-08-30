"""Phase 8 — Cascades organiques : du gène à l'organe.

Compare la propagation multi-échelle (gène → protéine → cellule → tissu →
organe) de l'hybride guidé (stable) vs l'hybride sauvage (effondré).
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bio114.cascade import (translate_protein, protein_folding_coherence,
                            run_cascade, organ_transformation_risk)

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"
FIGURES = Path(__file__).resolve().parent.parent / "figures"


def run() -> dict:
    guided_seq = (ARTIFACTS / "guided_hybrid_seq.txt").read_text().strip()
    wild_seq = (ARTIFACTS / "wild_hybrid_seq.txt").read_text().strip()

    phase4 = json.loads((ARTIFACTS / "phase4_synthesis.json").read_text())
    guided_psig = phase4["guided_hybrid"]["psig"]
    wild_psig = phase4["wild_hybrid"]["psig"]

    results = {}
    for name, seq, psig in [("hybride_guidé", guided_seq, guided_psig),
                            ("hybride_sauvage", wild_seq, wild_psig)]:
        protein = translate_protein(seq)
        prot_coh = protein_folding_coherence(protein)
        cascade = run_cascade(psig, prot_coh, name)
        risk = organ_transformation_risk(cascade)

        results[name] = {
            "protein_length": len(protein),
            "protein_folding_coherence": round(prot_coh, 4),
            "cascade": [lv.to_dict() for lv in cascade],
            "organ_risk": risk,
        }
        print(f"\n[{name}] protéine={len(protein)} aa, repliement={prot_coh:.3f}")
        for lv in cascade:
            mark = "✅" if lv.stable else "🛑"
            print(f"  {mark} {lv.name:18s} ({lv.scale}): cohérence={lv.coherence:.3f}")
        print(f"  → {risk['verdict']}")

    guided_stable = results["hybride_guidé"]["organ_risk"]["risk_ratio"] == 0
    wild_unstable = results["hybride_sauvage"]["organ_risk"]["risk_ratio"] > 0

    result = {
        "phase": 8,
        "title": "Cascades organiques — du gène à l'organe",
        "cascades": results,
        "conclusion": (
            "L'hybride guidé propage une cascade STABLE jusqu'à l'organe : "
            "la régénération fonctionne sans transformation sauvage. "
            if guided_stable else "")
        + (
            "L'hybride sauvage déclenche une cascade CHAOTIQUE — le garde-fou "
            "de l'Exp. 114 (repliement KTN:Li) l'avait déjà rejeté au niveau du "
            "gène. Double validation."
            if wild_unstable else ""),
        "phase8_success": guided_stable and wild_unstable,
    }

    (ARTIFACTS / "phase8_cascade.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))

    _make_figure(results)
    return result


def _make_figure(results: dict) -> None:
    FIGURES.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for name, data in results.items():
        levels = data["cascade"]
        names = [lv["level"] for lv in levels]
        coh = [lv["coherence"] for lv in levels]
        color = "#2ca02c" if "guidé" in name else "#d62728"
        ax.plot(names, coh, marker="o", label=name, color=color, linewidth=2)
    ax.axhline(0.4, color="orange", linestyle="--", label="seuil de stabilité (0.4)")
    ax.set_ylabel("cohérence structurelle")
    ax.set_title("Phase 8 — Cascade organique : gène → protéine → cellule → tissu → organe")
    ax.legend()
    ax.set_ylim(0, 1.05)
    fig.tight_layout()
    fig.savefig(FIGURES / "phase8_cascade.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT PHASE 8 ===")
    print(out["conclusion"])
