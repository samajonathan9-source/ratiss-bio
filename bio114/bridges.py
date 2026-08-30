"""Phase 3 — Les ponts de réciprocité.

À l'échelle invisible (topologique), chercher les invariants communs entre
les deux espèces. Une boucle H1 chez la chauve-souris peut correspondre à
une boucle H1 humaine sous la même tension thermique → motif compatible.

Transfert transdisciplinaire : méthode issue de Crypto-net-veo- (détection
de régimes et anomalies furtives). Une séquence étrangère compatible est
comme une attaque furtive indétectable dans le trafic — ici on inverse
l'arme : la furtivité topologique devient compatibilité thérapeutique.
"""

from dataclasses import dataclass

import numpy as np

from .topology import sequence_to_point_cloud, compute_signature, TopologicalSignature
from .eth import ThermalBath, ETH_HUMAN, apply_bath_to_cloud


@dataclass
class Bridge:
    """Un pont de réciprocité : un fragment chauve-souris topologiquement
    compatible avec un contexte humain."""

    bat_fragment: str
    human_context: str
    position: int
    topo_compat: float       # proximité topologique (0–1, 1=identique)
    bat_psig: float
    human_local_psig: float
    functional_gain: float   # gain fonctionnel estimé (immunité/réparation)

    def to_dict(self) -> dict:
        return {
            "bat_fragment": self.bat_fragment,
            "human_context": self.human_context,
            "position": self.position,
            "topo_compat": round(self.topo_compat, 4),
            "bat_psig": round(self.bat_psig, 4),
            "human_local_psig": round(self.human_local_psig, 4),
            "functional_gain": round(self.functional_gain, 4),
        }


def _local_psig(seq: str, bath: ThermalBath) -> float:
    """P_sig d'un fragment dans un bain donné (sans sous-échantillonnage)."""
    from .topology import compute_persistence, _total_persistence
    cloud = sequence_to_point_cloud(seq)
    cloud = apply_bath_to_cloud(cloud, bath)
    dgms = compute_persistence(cloud)
    h1_p, _ = _total_persistence(dgms["h1"])
    h2_p, _ = _total_persistence(dgms["h2"])
    raw = h1_p + 1.5 * h2_p
    return float(1.0 - np.exp(-raw / (15.0 * max(1, len(cloud)))))


def find_bridges(human_seq: str, bat_seq: str, window: int = 120,
                 bath: ThermalBath = ETH_HUMAN, top_k: int = 5) -> list[Bridge]:
    """Cherche les ponts de réciprocité entre génome humain et chauve-souris.

    Pour chaque fenêtre du génome humain, cherche le fragment chauve-souris
    dont la signature topologique locale (dans le même bain) est la plus
    proche. Le gain fonctionnel est estimé par la richesse en boucles H1
    du fragment chauve-souris (proxy de complexité fonctionnelle).
    """
    bridges = []
    n_h, n_b = len(human_seq), len(bat_seq)

    for pos in range(0, n_h - window, window):
        human_frag = human_seq[pos:pos + window]
        h_psig = _local_psig(human_frag, bath)

        best = None
        for bpos in range(0, n_b - window, window):
            bat_frag = bat_seq[bpos:bpos + window]
            b_psig = _local_psig(bat_frag, bath)
            compat = 1.0 - min(1.0, abs(b_psig - h_psig) / max(h_psig, 1e-6))
            # gain fonctionnel : le fragment chauve-souris apporte plus de
            # structure (P_sig supérieur) tout en restant compatible
            gain = max(0.0, b_psig - h_psig)
            score = compat + gain
            if best is None or score > best[0]:
                best = (score, compat, gain, bpos, bat_frag, b_psig)

        if best is not None:
            _, compat, gain, bpos, bat_frag, b_psig = best
            bridges.append(Bridge(
                bat_fragment=bat_frag,
                human_context=human_frag,
                position=pos,
                topo_compat=compat,
                bat_psig=b_psig,
                human_local_psig=h_psig,
                functional_gain=gain,
            ))

    # trier par compatibilité × gain, garder les meilleurs
    bridges.sort(key=lambda b: b.topo_compat * (1 + b.functional_gain), reverse=True)
    return bridges[:top_k]
