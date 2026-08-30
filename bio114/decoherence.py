"""Phase 5 — L'œil quantique : la décohérence de la structure hybride.

À l'échelle atomique, la physique classique ne suffit plus. Les liaisons
hydrogène qui maintiennent les deux brins d'ADN obéissent à la mécanique
quantique : tunneling des protons (mutations spontanées), forces de
Van der Waals (compaction).

Ce module mesure à quelle vitesse la structure hybride perd sa cohérence
quantique. Décohérence rapide → structure instable → rejet. Cohérence
tenue → fusion viable.

Transfert transdisciplinaire : moteur hérité de
ratiss-topological-decoherence-engine et des protocoles RATISS-ODV-AEON.
Pont vers le QPU IBM (QPU-Ratiss-COSMOS) pour les calculs au-delà du
classique — clé dans secrets/.
"""

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Constantes quantiques
HBAR = 1.054e-34          # J·s
K_B = 1.381e-23           # J/K
PROTON_MASS = 1.673e-27   # kg

# Barrière de tunneling pour le transfert de proton dans une liaison H
# (A-T : 2 liaisons, G-C : 3 liaisons) — hauteur effective en eV.
TUNNEL_BARRIER_EV = {"AT": 0.6, "GC": 0.9}
BARRIER_WIDTH_A = 0.7     # largeur de barrière en Å


@dataclass
class DecoherenceReport:
    """Rapport de décohérence d'une structure."""

    name: str
    n_gc: int
    n_at: int
    tunneling_rate: float        # s⁻¹ — taux de tunneling des protons
    decoherence_time_ns: float   # ns — temps de cohérence
    coherence_score: float       # 0–1, 1 = cohérence maintenue
    qpu_used: bool
    verdict: str

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "n_gc": self.n_gc,
            "n_at": self.n_at,
            "tunneling_rate": f"{self.tunneling_rate:.3e}",
            "decoherence_time_ns": round(self.decoherence_time_ns, 3),
            "coherence_score": round(self.coherence_score, 4),
            "qpu_used": self.qpu_used,
            "verdict": self.verdict,
        }


def _tunneling_rate(barrier_ev: float, temperature: float) -> float:
    """Taux de tunneling WKB d'un proton à travers une barrière.

    κ = exp(-2·d·√(2m·V)/ħ) — approximation Wentzel-Kramers-Brillouin.
    """
    v_joule = barrier_ev * 1.602e-19
    d = BARRIER_WIDTH_A * 1e-10
    kappa = 2 * d * np.sqrt(2 * PROTON_MASS * v_joule) / HBAR
    attempt_freq = K_B * temperature / HBAR  # fréquence de tentative
    return float(attempt_freq * np.exp(-kappa))


def measure_decoherence(seq: str, temperature: float = 310.0,
                        name: str = "structure") -> DecoherenceReport:
    """Mesure la décohérence quantique d'une séquence.

    Compte les liaisons GC (3 ponts H, plus stables quantiquement) et AT
    (2 ponts H, tunneling plus probable). Plus il y a de GC, plus la
    cohérence est maintenue longtemps.
    """
    n_gc = sum(1 for b in seq if b in "GC")
    n_at = sum(1 for b in seq if b in "AT")

    rate_gc = _tunneling_rate(TUNNEL_BARRIER_EV["GC"], temperature)
    rate_at = _tunneling_rate(TUNNEL_BARRIER_EV["AT"], temperature)

    # taux de tunneling moyen pondéré par la composition
    n = max(1, n_gc + n_at)
    avg_rate = (n_gc * rate_gc + n_at * rate_at) / n

    # temps de décohérence : inverse du taux total de tunneling
    total_rate = avg_rate * n
    decoherence_time_s = 1.0 / max(total_rate, 1e-30)
    decoherence_time_ns = decoherence_time_s * 1e9

    # score de cohérence : les liaisons GC maintiennent la cohérence.
    # Normalisé : >50% GC → cohérence élevée.
    gc_fraction = n_gc / n
    coherence_score = float(np.clip(0.3 + 1.2 * gc_fraction, 0.0, 1.0))

    verdict = ("cohérence maintenue — fusion viable"
               if coherence_score > 0.5 else
               "décohérence rapide — structure instable, rejet")

    return DecoherenceReport(
        name=name, n_gc=n_gc, n_at=n_at,
        tunneling_rate=avg_rate,
        decoherence_time_ns=decoherence_time_ns,
        coherence_score=coherence_score,
        qpu_used=False,
        verdict=verdict,
    )


def qpu_available() -> bool:
    """Le QPU IBM est-il accessible ? (clé dans secrets/)."""
    token_file = Path(__file__).resolve().parent.parent / "secrets" / "IBM_QUANTUM_TOKEN"
    if not token_file.exists():
        # chercher n'importe quel fichier dans secrets/
        secrets_dir = Path(__file__).resolve().parent.parent / "secrets"
        if secrets_dir.exists():
            files = [f for f in secrets_dir.iterdir() if f.is_file()]
            return len(files) > 0
        return False
    content = token_file.read_text().strip()
    return len(content) > 0 and "placeholder" not in content.lower()


def qpu_coupling_estimate(seq: str) -> float:
    """Estimation du couplage quantique via le QPU IBM.

    Si la clé est disponible, cette fonction serait le point d'entrée
    vers IBM Quantum (qiskit) pour simuler les interactions électroniques
    réelles des liaisons hybrides — au-delà du calcul classique.

    Retourne ici une estimation analytique du couplage de Van der Waals.
    """
    n = len(seq)
    gc = sum(1 for b in seq if b in "GC")
    # couplage proportionnel à la densité de GC (dipôles plus forts)
    return float(np.clip(gc / max(1, n) * 1.4, 0.0, 1.0))
