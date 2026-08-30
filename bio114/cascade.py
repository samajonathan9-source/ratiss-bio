"""Phase 8 — Modélisation des cascades organiques (multi-échelle).

L'ADN hybride est topologiquement stable (Exp. 114). Mais une modification
génétique se propage en cascade : gène → protéine → cellule → tissu →
organe. C'est là que naît le syndrome Morbius : une mutation réussie au
nanomètre peut déclencher une transformation macroscopique sauvage.

Modèle : chaque échelle est une couche de propagation. La cohérence
topologique (P_sig) de l'étage inférieur module la stabilité de l'étage
supérieur. Si la cohérence s'effondre à un étage, la cascade devient
chaotique → transformation organique sauvage.

Transfert transdisciplinaire : propagation multi-échelle inspirée du
moteur de décohérence (ratiss-topological-decoherence-engine) — une
perturbation locale se diffuse et s'amplifie ou s'éteint selon la
cohérence du milieu. Même logique que la propagation d'une attaque dans
un réseau (Crypto-net-veo-) : l'organe est le réseau, la mutation est
l'intrus.
"""

from dataclasses import dataclass, field

import numpy as np

# Codons → acides aminés (table standard simplifiée, codons de départ/stop
# + principaux). Pour la simulation du recrutement protéique.
CODON_TABLE = {
    "ATG": "Met", "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "GTT": "Val", "GTC": "Val",
    "GTA": "Val", "GTG": "Val", "TCT": "Ser", "TCC": "Ser", "TCA": "Ser",
    "TCG": "Ser", "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "CAT": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln", "AAT": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys", "GAT": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu", "TGT": "Cys", "TGC": "Cys",
    "TGG": "Trp", "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    "TAA": "STOP", "TAG": "STOP", "TGA": "STOP",
}

# Hydrophobie des acides aminés (échelle de Kyte-Doolittle simplifiée).
# Détermine le repliement protéique : les hydrophobes se cachent au cœur.
HYDROPHOBICITY = {
    "Ile": 4.5, "Val": 4.2, "Leu": 3.8, "Phe": 2.8, "Cys": 2.5,
    "Met": 1.9, "Ala": 1.8, "Gly": -0.4, "Thr": -0.7, "Ser": -0.8,
    "Trp": -0.9, "Tyr": -1.3, "Pro": -1.6, "His": -3.2, "Gln": -3.5,
    "Asn": -3.5, "Glu": -3.5, "Asp": -3.5, "Lys": -3.9, "Arg": -4.5,
    "STOP": 0.0,
}


@dataclass
class CascadeLevel:
    """Un étage de la cascade organique."""

    name: str
    scale: str               # échelle (nm, µm, mm, cm)
    coherence: float         # cohérence structurelle 0–1
    amplification: float     # facteur d'amplification de la perturbation
    stable: bool
    detail: str

    def to_dict(self) -> dict:
        return {
            "level": self.name, "scale": self.scale,
            "coherence": round(self.coherence, 4),
            "amplification": round(self.amplification, 4),
            "stable": self.stable, "detail": self.detail,
        }


def translate_protein(seq: str) -> list[str]:
    """Traduit une séquence ADN en chaîne d'acides aminés."""
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = CODON_TABLE.get(codon, "Gly")
        if aa == "STOP":
            break
        protein.append(aa)
    return protein


def protein_folding_coherence(protein: list[str]) -> float:
    """Cohérence du repliement protéique.

    Une protéine fonctionnelle a un cœur hydrophobe stable entouré d'une
    enveloppe hydrophile. La cohérence mesure l'alternance équilibrée
    hydrophobe/hydrophile (proxy du repliement correct).
    """
    if len(protein) < 3:
        return 0.0
    hydro = [HYDROPHOBICITY.get(aa, 0.0) for aa in protein]
    # variance de l'hydrophobie : trop uniforme = pas de structure,
    # équilibrée = repliement stable (cœur + surface)
    arr = np.array(hydro)
    # fenêtre glissante : alternance locale = signal de repliement
    diffs = np.abs(np.diff(arr))
    fold_signal = float(np.mean(diffs) / 4.0)  # normalisé
    return float(np.clip(fold_signal, 0.0, 1.0))


def run_cascade(gene_psig: float, protein_coherence: float,
                name: str = "hybride") -> list[CascadeLevel]:
    """Propage la perturbation du gène vers l'organe, niveau par niveau.

    La cohérence de chaque étage dépend de l'étage précédent : une perte
    de cohérence s'amplifie (cascade chaotique) ou s'atténue (absorption).
    """
    levels = []
    # État initial : cohérence du gène (P_sig normalisé × repliement protéique)
    coherence = float(np.clip(gene_psig * 100 * 0.6 + protein_coherence * 0.4, 0.0, 1.0))

    cascade_def = [
        ("Gène hybride", "nm", coherence, 1.0),
        ("Protéine", "nm", protein_coherence, 1.1),
        ("Cellule", "µm", None, 1.3),      # cohérence héritée
        ("Tissu", "mm", None, 1.6),
        ("Organe", "cm", None, 2.0),
    ]

    current = coherence
    for level_name, scale, base_coh, amp in cascade_def:
        if base_coh is not None:
            current = float(np.clip(base_coh * 0.7 + current * 0.3, 0.0, 1.0))
        else:
            # propagation : la cohérence de l'étage inférieur module l'étage.
            # Si cohérence > seuil (0.4) → la structure absorbe la perturbation
            # (atténuation). Sinon → amplification chaotique.
            if current > 0.4:
                current = float(np.clip(current * (1 + 0.05 / amp), 0.0, 1.0))
            else:
                current = float(np.clip(current * (1 - 0.20 * amp), 0.0, 1.0))

        stable = current > 0.4
        detail = ("structure absorbée, fonction maintenue" if stable else
                  "cascade chaotique — transformation sauvage")
        levels.append(CascadeLevel(
            name=level_name, scale=scale, coherence=current,
            amplification=amp, stable=stable, detail=detail))

    return levels


def organ_transformation_risk(levels: list[CascadeLevel]) -> dict:
    """Risque de transformation organique : proportion d'étages instables
    et profondeur atteinte par la cascade chaotique."""
    unstable = [lv for lv in levels if not lv.stable]
    depth = next((lv.name for lv in levels if not lv.stable), "aucun")
    risk = len(unstable) / len(levels)
    verdict = ("organe stable — pas de transformation sauvage"
               if risk == 0 else
               f"⚠️ cascade instable à partir de l'étage '{depth}' — "
               f"risque de transformation {risk * 100:.0f}%")
    return {
        "unstable_levels": [lv.name for lv in unstable],
        "first_failure": depth,
        "risk_ratio": round(risk, 3),
        "verdict": verdict,
    }
