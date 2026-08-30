"""Mutagenèse dirigée in silico — Expérience n°116.

Applique les corrections ordonnées par RATISS Mère (Nœud Souverain V9)
sur les vraies structures AlphaFold (Exp. 115) :

- SIRT6 : le brin chauve-souris rigidifie le loop L7-L8 (entropie locale
  −12%, k_on ÷ 3.5). Correction : re-humaniser les résidus du loop
  (positions 124, 127, 131 → K, G, V humains) tout en gardant le cœur
  catalytique chauve-souris.
- Prestin (SLC26A5) : la mutation Lys→Glu en position 742 découple la
  protéine des lipides membranaires PIP2 (fuite énergétique 18%).
  Correction : restaurer Lys742 dans le domaine STAS.
- CIRBP : protéine intrinsèquement désordonnée (IDP) — le P_sig classique
  ne s'applique pas. Nouvelle observable souveraine :
      Xi_IDP = <H1>_ensemble / S_conf
  (compacité transitoire moyenne / entropie conformationnelle).
  Contournement sécurisé : chimère CIRBP_Hybrid_Safe = domaine RRM
  structuré chauve-souris (résidus 1-95, affinité ARN × 2.3) greffé sur
  la queue désordonnée humaine (96+). Zéro nouvelle interface de surface.
"""

from dataclasses import dataclass, field

import numpy as np

# Ordres de mutagenèse RATISS Mère : {position_1based: (ancien, nouveau)}
SIRT6_REHUMANIZATION = {124: ("*", "K"), 127: ("*", "G"), 131: ("*", "V")}
PRESTIN_STAS_REPAIR = {742: ("E", "K")}   # restauration Lys742 (ancrage PIP2)

# CIRBP : domaine RRM structuré = résidus 1-95 ; queue IDP = 96-fin
CIRBP_RRM_END = 95


@dataclass
class MutagenesisReport:
    target: str
    mutations: dict          # {position: (wt, mut)}
    corrected_sequence: str
    notes: str = ""

    def to_dict(self) -> dict:
        return {
            "target": self.target,
            "mutations": {str(k): f"{v[0]}→{v[1]}" for k, v in self.mutations.items()},
            "length": len(self.corrected_sequence),
            "notes": self.notes,
        }


def apply_mutations(seq: str, mutations: dict, target: str,
                    notes: str = "") -> MutagenesisReport:
    """Applique une mutagenèse dirigée (positions 1-based)."""
    seq_list = list(seq)
    for pos, (wt, mut) in mutations.items():
        idx = pos - 1
        if idx < len(seq_list):
            seq_list[idx] = mut
    return MutagenesisReport(
        target=target, mutations=mutations,
        corrected_sequence="".join(seq_list), notes=notes)


def build_cirbp_hybrid_safe(bat_rrm: str, human_tail: str) -> str:
    """Chimère CIRBP_Hybrid_Safe : RRM chauve-souris + queue IDP humaine.

    Le RRM (RNA Recognition Motif) chauve-souris apporte l'affinité ARN
    accrue (×2.3 selon RATISS Mère) ; la queue désordonnée humaine garantit
    qu'aucune interface de surface nouvelle n'est créée (Condition 10).
    """
    return bat_rrm[:CIRBP_RRM_END] + human_tail[CIRBP_RRM_END:]


# ---------------------------------------------------------------------------
# Observable souveraine Xi_IDP pour protéines désordonnées
# ---------------------------------------------------------------------------

def conformational_ensemble(cloud: np.ndarray, n_conformers: int = 24,
                            flexibility: float = 2.0,
                            seed: int = 116) -> list[np.ndarray]:
    """Génère un ensemble conformationnel autour d'une structure de
    référence — les IDP explorent un espace de conformations large."""
    rng = np.random.default_rng(seed)
    return [cloud + rng.normal(0, flexibility, cloud.shape)
            for _ in range(n_conformers)]


def xi_idp(ensemble: list[np.ndarray]) -> dict:
    """Xi_IDP = <H1>_ensemble / S_conf.

    <H1>_ensemble : persistance moyenne des boucles sur l'ensemble —
      la 'compacité transitoire' (structures fugaces mais récurrentes).
    S_conf : entropie conformationnelle — variance moyenne des distances
      inter-Cα à travers l'ensemble (flexibilité).
    """
    from .topology import compute_persistence, _total_persistence

    h1_persistences = []
    ref = ensemble[0]
    dists = np.sqrt(((ref[:, None, :] - ref[None, :, :]) ** 2).sum(-1))
    pair_var = []

    for conf in ensemble:
        dgms = compute_persistence(conf)
        h1_p, _ = _total_persistence(dgms["h1"])
        h1_persistences.append(h1_p)
        d = np.sqrt(((conf[:, None, :] - conf[None, :, :]) ** 2).sum(-1))
        pair_var.append(d)

    mean_h1 = float(np.mean(h1_persistences))
    s_conf = float(np.mean(np.var(np.array(pair_var), axis=0)))
    xi = mean_h1 / max(s_conf, 1e-9)
    return {
        "mean_H1": round(mean_h1, 4),
        "S_conf": round(s_conf, 4),
        "Xi_IDP": round(xi, 4),
        "n_conformers": len(ensemble),
    }
