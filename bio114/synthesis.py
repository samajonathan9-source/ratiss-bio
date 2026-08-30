"""Phase 4 — La synthèse par petits brins.

Pas de copier-coller brutal. On insère le fragment chauve-souris dans le
contexte humain au point de plus forte compatibilité topologique, puis on
soumet l'hybride au bain humain.

Transfert transdisciplinaire : sélection du meilleur candidat par
maximisation de P_sig — le même principe que la SUPER ACTIVATION de
ratiss-Skynet (8 candidats générés, la topologie choisit, pas le hasard).
"""

from dataclasses import dataclass

from .topology import (sequence_to_point_cloud, compute_signature,
                       TopologicalSignature)
from .eth import ThermalBath, ETH_HUMAN, apply_bath_to_cloud, sequence_stability, thermodynamic_collapse
from .bridges import Bridge


@dataclass
class HybridCandidate:
    """Une instance hybride candidate."""

    sequence: str
    insertion_point: int
    psig: float
    stability_in_human_bath: float
    thermodynamic_ok: bool
    signature: TopologicalSignature

    def to_dict(self) -> dict:
        return {
            "length": len(self.sequence),
            "insertion_point": self.insertion_point,
            "psig": round(self.psig, 4),
            "stability_in_human_bath": round(self.stability_in_human_bath, 4),
            "thermodynamic_ok": self.thermodynamic_ok,
        }


def synthesize(human_seq: str, bridge: Bridge,
               bath: ThermalBath = ETH_HUMAN) -> HybridCandidate:
    """Synthétise un hybride en insérant le fragment chauve-souris au pont
    de réciprocité (remplacement local, pas ajout brutal)."""
    pos = bridge.insertion_point if hasattr(bridge, "insertion_point") else bridge.position
    frag = bridge.bat_fragment
    w = len(frag)
    hybrid = human_seq[:pos] + frag + human_seq[pos + w:]

    # signature topologique de l'hybride dans le bain humain
    cloud = sequence_to_point_cloud(hybrid)
    cloud = apply_bath_to_cloud(cloud, bath)
    sig = compute_signature(cloud, "hybrid_candidate")

    # stabilité thermodynamique
    stab = sequence_stability(hybrid, bath)
    ok = not thermodynamic_collapse(stab["stability"])

    return HybridCandidate(
        sequence=hybrid,
        insertion_point=pos,
        psig=sig.p_sig,
        stability_in_human_bath=stab["stability"],
        thermodynamic_ok=ok,
        signature=sig,
    )
