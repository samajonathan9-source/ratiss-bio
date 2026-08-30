"""Phase 6 — Le garde-fou topologique : guérir sans perdre son humanité.

Surveille en continu la signature P_sig de l'hybride vs l'invariant
humain. Si la structure dérive vers la signature animale (ou s'effondre),
ALERTE ROUGE → repliement KTN:Li : retour forcé à la structure humaine
stable.

Transfert transdisciplinaire : Conditions 9 & 10 de la fiche AGI —
honnêteté mesurable sur l'incertitude (« je ne peux pas intégrer ce gène
sans détruire la forme humaine, j'arrête ») + mécanisme d'arrêt fiable.
Le même repliement qui interrompt l'hallucination dans Skynet interrompt
ici la mutation sauvage.
"""

from dataclasses import dataclass

from .topology import (sequence_to_point_cloud, compute_signature,
                       TopologicalSignature, signature_distance)
from .eth import ThermalBath, ETH_HUMAN, apply_bath_to_cloud

# Frontière absolue : au-delà de cette distance de l'invariant humain,
# la structure n'est plus humaine.
DRIFT_THRESHOLD = 0.15
# Seuil d'effondrement : P_sig sous cette valeur = structure dénaturée.
COLLAPSE_PSIG = 0.003


@dataclass
class GuardVerdict:
    """Le verdict du garde-fou topologique."""

    structure_name: str
    psig: float
    distance_to_human: float
    collapsed: bool
    drifted_to_animal: bool
    action: str          # "MAINTIEN" | "REPLIEMENT_KTN_LI"
    reason: str

    def to_dict(self) -> dict:
        return {
            "structure": self.structure_name,
            "psig": round(self.psig, 4),
            "distance_to_human": round(self.distance_to_human, 4),
            "collapsed": self.collapsed,
            "drifted_to_animal": self.drifted_to_animal,
            "action": self.action,
            "reason": self.reason,
        }


class TopologicalGuard:
    """Le garde-fou : surveille et, si nécessaire, replie (KTN:Li)."""

    def __init__(self, human_ref_seq: str, bat_ref_seq: str,
                 bath: ThermalBath = ETH_HUMAN):
        self.bath = bath
        # invariant humain de référence
        h_cloud = apply_bath_to_cloud(sequence_to_point_cloud(human_ref_seq), bath)
        self.human_sig = compute_signature(h_cloud, "human_invariant")
        # signature animale (pour mesurer la direction de la dérive)
        b_cloud = apply_bath_to_cloud(sequence_to_point_cloud(bat_ref_seq), bath)
        self.bat_sig = compute_signature(b_cloud, "bat_invariant")

    def inspect(self, seq: str, name: str) -> GuardVerdict:
        """Inspecte une structure et décide : maintien ou repliement."""
        cloud = apply_bath_to_cloud(sequence_to_point_cloud(seq), self.bath)
        sig = compute_signature(cloud, name)

        dist_human = signature_distance(sig, self.human_sig)
        dist_bat = signature_distance(sig, self.bat_sig)

        collapsed = sig.p_sig < COLLAPSE_PSIG
        # dérive vers l'animal : plus proche de la chauve-souris que de l'humain
        # ET au-delà du seuil de dérive autorisé
        drifted = (dist_bat < dist_human) and (dist_human > DRIFT_THRESHOLD)

        if collapsed:
            action = "REPLIEMENT_KTN_LI"
            reason = (f"Effondrement thermodynamique (P_sig={sig.p_sig:.4f} < "
                      f"{COLLAPSE_PSIG}). Structure dénaturée — repliement forcé "
                      "vers la forme humaine stable.")
        elif drifted:
            action = "REPLIEMENT_KTN_LI"
            reason = (f"Dérive vers la signature animale (distance humain="
                      f"{dist_human:.3f} > seuil {DRIFT_THRESHOLD}, plus proche de "
                      "la chauve-souris). Mutation sauvage détectée — repliement "
                      "forcé. Honnêteté sur l'incertitude : je ne peux pas intégrer "
                      "ce gène sans détruire la forme humaine, j'arrête.")
        else:
            action = "MAINTIEN"
            reason = (f"Structure stable et proche de l'invariant humain "
                      f"(P_sig={sig.p_sig:.4f}, distance={dist_human:.3f}). "
                      "Régénération sans perte d'humanité.")

        return GuardVerdict(
            structure_name=name,
            psig=sig.p_sig,
            distance_to_human=dist_human,
            collapsed=collapsed,
            drifted_to_animal=drifted,
            action=action,
            reason=reason,
        )
