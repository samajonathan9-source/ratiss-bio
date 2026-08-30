"""Phase 2 — Le bain thermodynamique virtuel (eth).

L'ADN ne vit pas dans le vide. Il baigne dans un milieu : température,
pH, ions, tension mécanique. Ce module recrée ces conditions et applique
la contrainte physique aux brins d'ADN virtuels.

Transfert transdisciplinaire : hérite du corps thermodynamique de
ratiss-Skynet/thermo_emotions.py (HYBRID MIND) — l'état thermique module
la structure comme l'émotion module la génération dans le MCT.
"""

from dataclasses import dataclass

import numpy as np

# Constantes physiques
R_GAS = 1.987e-3  # kcal/(mol·K), constante des gaz

# ΔG d'empilement nearest-neighbor (kcal/mol à 37°C, SantaLucia 1998).
# Clé = dinucléotide du brin sens. Plus négatif = empilement plus stable.
NN_DG37 = {
    "AA": -1.00, "AT": -0.88, "AG": -1.28, "AC": -1.44,
    "TA": -0.58, "TT": -1.00, "TG": -1.45, "TC": -1.28,
    "GA": -1.30, "GT": -1.44, "GG": -1.84, "GC": -2.17,
    "CA": -1.45, "CT": -1.30, "CG": -2.24, "CC": -1.84,
}
DEFAULT_DG = -1.2


@dataclass
class ThermalBath:
    """Un environnement thermodynamique (bain) pour l'ADN."""

    name: str
    temperature: float        # K
    ph: float
    ionic_strength: float     # mol/L (force ionique)
    supercoil_stress: float   # rad/pb — tension mécanique du super-enroulement

    def dg_correction(self) -> float:
        """Correction du ΔG selon la température (van 't Hoff simplifié).

        À température élevée, l'empilement se déstabilise (dénaturation).
        Référence : 310 K (37°C, corps humain).
        """
        return (self.temperature - 310.0) * 0.02  # kcal/mol par K d'écart


# Les deux bains de référence
ETH_HUMAN = ThermalBath(
    name="eth_humain",
    temperature=310.0,        # 37°C
    ph=7.4,
    ionic_strength=0.15,      # ~physiologique
    supercoil_stress=-0.05,   # super-enroulement négatif physiologique
)

ETH_BAT = ThermalBath(
    name="eth_chauve_souris",
    temperature=313.0,        # 40°C en vol (métabolisme extrême)
    ph=7.3,
    ionic_strength=0.16,
    supercoil_stress=-0.08,   # tension plus forte (cellules très actives)
)


def sequence_stability(seq: str, bath: ThermalBath) -> dict:
    """Stabilité thermodynamique d'une séquence dans un bain donné.

    Retourne ΔG total (kcal/mol), le ΔG moyen par paire, et un score
    de stabilité normalisé 0–1. Plus le ΔG est négatif, plus la double
    hélice est stable dans ce bain.
    """
    corr = bath.dg_correction()
    n = len(seq)
    total = 0.0
    for i in range(n - 1):
        dinuc = seq[i:i + 2]
        total += NN_DG37.get(dinuc, DEFAULT_DG) + corr

    # pH extrême déstabilise (déprotonation des bases)
    ph_penalty = abs(bath.ph - 7.35) * 1.5
    # force ionique stabilise (écrantage des phosphates négatifs)
    ionic_bonus = -bath.ionic_strength * 2.0

    dg_total = total + ph_penalty + ionic_bonus
    dg_per_bp = dg_total / max(1, n - 1)

    # Score de stabilité : 1 = très stable, 0 = dénaturé.
    # ΔG/pb typique entre -2.5 (très stable) et -0.3 (dénaturé).
    stability = float(np.clip((-dg_per_bp - 0.3) / 2.2, 0.0, 1.0))

    return {
        "bath": bath.name,
        "temperature": bath.temperature,
        "dg_total": round(dg_total, 2),
        "dg_per_bp": round(dg_per_bp, 3),
        "stability": round(stability, 4),
    }


def thermodynamic_collapse(stability: float, threshold: float = 0.35) -> bool:
    """Effondrement thermodynamique : la structure se dénature si la
    stabilité chute sous le seuil (P_sig s'effondre → échec)."""
    return stability < threshold


def apply_bath_to_cloud(points: np.ndarray, bath: ThermalBath) -> np.ndarray:
    """Applique la contrainte du bain au nuage de points.

    La température agite thermiquement (bruit brownien), la tension de
    super-enroulement tord la structure. Retourne le nuage déformé.
    """
    rng = np.random.default_rng(seed=int(bath.temperature * 100))
    # agitation thermique proportionnelle à l'écart de température
    thermal_noise = (bath.temperature - 280.0) * 0.03  # Å
    noisy = points + rng.normal(0, thermal_noise, points.shape)

    # torsion par super-enroulement autour de l'axe z du centre de masse
    center = noisy.mean(axis=0)
    rel = noisy - center
    z = rel[:, 2]
    angle = bath.supercoil_stress * (z - z.min()) / max(1.0, (z.max() - z.min()))
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    x = rel[:, 0] * cos_a - rel[:, 1] * sin_a
    y = rel[:, 0] * sin_a + rel[:, 1] * cos_a
    twisted = np.stack([x, y, rel[:, 2]], axis=1) + center
    return twisted
