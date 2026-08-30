"""Phase 7 — Preuve et verdict : démonstration par le fonctionnement.

Agrège toutes les phases, génère le tableau de verdict final, les figures
comparatives, et produit le rapport scientifique.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARTIFACTS = Path(__file__).resolve().parent.parent / "artifacts"
FIGURES = Path(__file__).resolve().parent.parent / "figures"
DOCS = Path(__file__).resolve().parent.parent / "docs"
NAMES = {1: "invariant", 2: "eth", 3: "bridges", 4: "synthesis", 5: "decoherence", 6: "guard"}


def load(phase: int) -> dict:
    return json.loads((ARTIFACTS / f"phase{phase}_{NAMES[phase]}.json").read_text())


def run() -> dict:
    p6 = load(6)
    v = p6["verdicts"]

    tableau = [
        {"etat": "ADN humain pur", "p_sig": v["humain_pur"]["psig"],
         "verdict": v["humain_pur"]["action"], "role": "Invariant de référence"},
        {"etat": "ADN chauve-souris", "p_sig": v["chauve_souris_pure"]["psig"],
         "verdict": v["chauve_souris_pure"]["action"], "role": "Invariant animal"},
        {"etat": "Hybride sauvage", "p_sig": v["hybride_sauvage"]["psig"],
         "verdict": v["hybride_sauvage"]["action"], "role": "Rejet — effondrement"},
        {"etat": "Hybride réussi (guidé)", "p_sig": v["hybride_guidé"]["psig"],
         "verdict": v["hybride_guidé"]["action"], "role": "Régénération, humanité intacte"},
    ]

    success = (v["hybride_guidé"]["action"] == "MAINTIEN"
               and v["hybride_sauvage"]["action"] == "REPLIEMENT_KTN_LI")

    result = {
        "phase": 7,
        "title": "Preuve et verdict final — Expérience n°114",
        "loi": "R = P_sig ; ΔW = η·φ·P_sig·C",
        "tableau_verdict": tableau,
        "phases_resume": {f"phase{i}_{NAMES[i]}": load(i)["conclusion"] for i in range(1, 7)},
        "experience_reussie": success,
        "conclusion_finale": (
            "EXPÉRIENCE N°114 RÉUSSIE. La synthèse guidée par les ponts de "
            "réciprocité produit un hybride topologiquement stable qui survit au "
            "bain humain, maintient sa cohérence quantique, et est validé par le "
            "garde-fou. L'hybride sauvage (fusion brutale) est rejeté par "
            "repliement KTN:Li. Guérir sans perdre son humanité."
            if success else
            "EXPÉRIENCE N°114 INCONCLUSIVE. Voir les verdicts de phase."
        ),
    }

    (ARTIFACTS / "phase7_verdict.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False))

    _make_figure(tableau)
    _write_report(result, tableau)
    return result


def _make_figure(tableau: list) -> None:
    FIGURES.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    etats = [t["etat"] for t in tableau]
    psigs = [t["p_sig"] for t in tableau]
    colors = ["#2ca02c" if t["verdict"] == "MAINTIEN" else "#d62728" for t in tableau]
    bars = ax.bar(etats, psigs, color=colors)
    ax.axhline(0.003, color="orange", linestyle="--", label="seuil d'effondrement (0.003)")
    ax.set_ylabel("P_sig")
    ax.set_title("Expérience n°114 — Verdict : signature topologique des 4 états")
    ax.tick_params(axis="x", rotation=15)
    ax.legend()
    for bar, t in zip(bars, tableau):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.0002,
                f"{t['p_sig']:.4f}\n{t['verdict']}", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "phase7_verdict.png", dpi=130)
    plt.close(fig)


def _write_report(result: dict, tableau: list) -> None:
    lines = [
        "# RAPPORT — EXPÉRIENCE N°114",
        "## Bio-Topologie Quantique Computationnelle — Protocole Morbius Contrôlé",
        "",
        "**Loi fondatrice :** R = P_sig · ΔW = η·φ·P_sig·C",
        "**Auteurs :** Jonathan Evina & JOHNKING0 — RATIS Labs",
        "",
        "---",
        "",
        "## Verdict final",
        "",
        "| État | P_sig | Verdict | Rôle |",
        "|---|---|---|---|",
    ]
    for t in tableau:
        mark = "✅" if t["verdict"] == "MAINTIEN" else "🛑"
        lines.append(f"| {t['etat']} | {t['p_sig']:.4f} | {mark} {t['verdict']} | {t['role']} |")
    lines += [
        "",
        f"**Résultat :** {result['conclusion_finale']}",
        "",
        "---",
        "",
        "## Résumé des phases",
        "",
    ]
    for phase, concl in result["phases_resume"].items():
        lines.append(f"- **{phase}** : {concl}")
    lines += [
        "",
        "## Méthode transdisciplinaire",
        "",
        "Chaque phase transfère une loi validée de l'arsenal RATIS vers la biologie :",
        "- Tissage 3D KTN:Li → repliement de l'ADN",
        "- Corps thermodynamique HYBRID MIND → bain eth",
        "- Crypto-net-veo (furtivité) → compatibilité thérapeutique",
        "- SUPER ACTIVATION Skynet → sélection par P_sig",
        "- Moteur de décohérence → stabilité quantique des liaisons",
        "- Repliement KTN:Li → anti-mutation (Conditions 9 & 10 AGI)",
        "",
        "**Guérir sans perdre son humanité.**",
    ]
    DOCS.mkdir(exist_ok=True)
    (DOCS / "RAPPORT_EXPERIENCE_114.md").write_text("\n".join(lines))


if __name__ == "__main__":
    out = run()
    print("\n=== VERDICT FINAL — EXPÉRIENCE N°114 ===")
    for t in out["tableau_verdict"]:
        m = "✅" if t["verdict"] == "MAINTIEN" else "🛑"
        print(f"  {m} {t['etat']}: P_sig={t['p_sig']:.4f} → {t['verdict']}")
    print("\n" + out["conclusion_finale"])
