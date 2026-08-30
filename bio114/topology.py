"""Phase 1 — L'Invariant topologique.

L'ADN modélisé comme nuage de points 3D. La séquence guide le repliement
local via les propriétés physico-chimiques des bases (angles de roll
dinucléotidiques), comme dans la réalité structurale. Analysé par
homologie persistante (Vietoris-Rips) → signature P_sig.

Transfert transdisciplinaire : même traitement que le tissage 3D du
cristal ferroélectrique KTN:Li (repo Travaux). L'ADN possède des
propriétés ferro/piézoélectriques — les bases ont des moments dipolaires.
"""

from dataclasses import dataclass, field

import numpy as np

# Paramètres B-DNA (Angströms)
HELIX_RADIUS = 10.0       # rayon de l'hélice
RISE_PER_BP = 3.4         # élévation par paire de bases
BP_PER_TURN = 10.5        # paires de bases par tour
TWIST = 2 * np.pi / BP_PER_TURN

# Encombrement stérique par base (modifie localement le rayon)
BASE_BULK = {"A": 1.0, "G": 1.0, "C": 0.85, "T": 0.85}  # purines > pyrimidines

# Angles de roll dinucléotidiques (degrés) — flexibilité de l'empilement.
# Inspirés des tables de Olson et al. : la séquence détermine comment
# l'hélice se courbe localement. Clé = dinucléotide (5'->3').
DINUC_ROLL_DEG = {
    "AA": -1.0, "AT": -2.5, "AG": -0.6, "AC": 1.8,
    "TA": 4.5, "TT": -1.0, "TG": 4.8, "TC": 2.9,
    "GA": 5.6, "GT": 1.8, "GG": -1.8, "GC": 0.0,
    "CA": 9.0, "CT": 5.6, "CG": 6.1, "CC": -1.8,
}
DEFAULT_ROLL_DEG = 0.0

# Rigidité locale (3 ponts H pour GC vs 2 pour AT)
GC_STACKING = {"G": 1.3, "C": 1.3, "A": 1.0, "T": 1.0}


@dataclass
class TopologicalSignature:
    """Signature topologique P_sig d'un repliement d'ADN."""

    name: str
    n_points: int
    h0_persistence: float = 0.0   # composantes connexes
    h1_persistence: float = 0.0   # boucles
    h2_persistence: float = 0.0   # cavités
    h1_count: int = 0
    h2_count: int = 0
    p_sig: float = 0.0            # la signature globale (loi R = P_sig)
    diagram: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "n_points": self.n_points,
            "h0_persistence": round(self.h0_persistence, 4),
            "h1_persistence": round(self.h1_persistence, 4),
            "h2_persistence": round(self.h2_persistence, 4),
            "h1_count": self.h1_count,
            "h2_count": self.h2_count,
            "p_sig": round(self.p_sig, 4),
        }


def _rodrigues(vec: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotation de vec autour de axis (unitaire) par angle (Rodrigues)."""
    return (vec * np.cos(angle)
            + np.cross(axis, vec) * np.sin(angle)
            + axis * np.dot(axis, vec) * (1 - np.cos(angle)))


def sequence_to_point_cloud(seq: str, supercoil: float = 0.0) -> np.ndarray:
    """Convertit une séquence en nuage de points 3D (deux brins).

    Le repliement dépend de la séquence : chaque dinucléotide impose un
    angle de roll qui courbe localement l'hélice. Des motifs comme CA/TG
    (fort roll) créent des courbures marquées → boucles et cavités.
    supercoil : torsion supplémentaire (rad/pb) simulant le super-enroulement.
    """
    from .sequences import complement

    comp = complement(seq)
    n = len(seq)

    # Cumul des courbures : chaque dinucléotide penche l'axe de l'hélice,
    # la chaîne se replie progressivement (comme la chromatine).
    axis = np.zeros((n, 3))
    direction = np.array([0.0, 0.0, 1.0])
    pos = np.zeros(3)
    for i in range(n):
        if i > 0:
            dinuc = seq[i - 1:i + 1]
            roll = np.radians(DINUC_ROLL_DEG.get(dinuc, DEFAULT_ROLL_DEG))
            perp = np.cross(direction, [1.0, 0.0, 0.0])
            if np.linalg.norm(perp) < 1e-9:
                perp = np.array([0.0, 1.0, 0.0])
            perp /= np.linalg.norm(perp)
            direction = _rodrigues(direction, perp, roll)
            direction /= np.linalg.norm(direction)
            pos = pos + direction * RISE_PER_BP
        axis[i] = pos

    axis -= axis.mean(axis=0)

    pts = np.zeros((n, 3))
    for i, (b, cb) in enumerate(zip(seq, comp)):
        theta = i * TWIST + i * supercoil
        r1 = HELIX_RADIUS * BASE_BULK[b] * GC_STACKING[b]
        r2 = HELIX_RADIUS * BASE_BULK[cb] * GC_STACKING[cb]
        tangent = axis[min(i + 1, n - 1)] - axis[max(i - 1, 0)]
        if np.linalg.norm(tangent) < 1e-9:
            tangent = np.array([0.0, 0.0, 1.0])
        tangent /= np.linalg.norm(tangent)
        u = np.cross(tangent, [1.0, 0.0, 0.0])
        if np.linalg.norm(u) < 1e-9:
            u = np.cross(tangent, [0.0, 1.0, 0.0])
        u /= np.linalg.norm(u)
        v = np.cross(tangent, u)
        off1 = r1 * (np.cos(theta) * u + np.sin(theta) * v)
        off2 = r2 * (np.cos(theta + np.pi * 0.9) * u + np.sin(theta + np.pi * 0.9) * v)
        pts[i] = axis[i] + off1
    return pts


def compute_persistence(points: np.ndarray) -> dict:
    """Homologie persistante Vietoris-Rips (H0, H1, H2) via ripser."""
    from ripser import ripser

    pts = points
    if len(pts) > 220:
        idx = np.linspace(0, len(pts) - 1, 220).astype(int)
        pts = pts[idx]
    result = ripser(pts, maxdim=2, thresh=30.0)
    return {f"h{k}": result["dgms"][k] for k in range(3)}


def _total_persistence(dgm: np.ndarray) -> tuple[float, int]:
    """Persistance totale et compte des features significatives."""
    if len(dgm) == 0:
        return 0.0, 0
    finite = dgm[np.isfinite(dgm[:, 1])]
    if len(finite) == 0:
        return 0.0, 0
    pers = finite[:, 1] - finite[:, 0]
    significant = pers[pers > 1.0]  # au-delà du bruit (> 1 Å)
    return float(np.sum(significant)), int(len(significant))


def compute_signature(points: np.ndarray, name: str) -> TopologicalSignature:
    """Calcule la signature topologique complète d'un nuage de points.

    P_sig combine les persistances H1 (boucles) et H2 (cavités de
    repliement) — la 'forme' informationnelle de l'ADN. Normalisée 0–1.
    """
    dgms = compute_persistence(points)

    h0_p, _ = _total_persistence(dgms["h0"])
    h1_p, h1_c = _total_persistence(dgms["h1"])
    h2_p, h2_c = _total_persistence(dgms["h2"])

    # Loi R = P_sig : H2 (cavités de repliement) pèse plus.
    raw = h1_p + 1.5 * h2_p
    p_sig = float(1.0 - np.exp(-raw / (15.0 * max(1, len(points)))))

    return TopologicalSignature(
        name=name,
        n_points=len(points),
        h0_persistence=h0_p,
        h1_persistence=h1_p,
        h2_persistence=h2_p,
        h1_count=h1_c,
        h2_count=h2_c,
        p_sig=p_sig,
        diagram={k: v.tolist() for k, v in dgms.items()},
    )


def signature_distance(a: TopologicalSignature, b: TopologicalSignature) -> float:
    """Distance entre deux signatures (pour détecter la dérive)."""
    return float(
        abs(a.p_sig - b.p_sig)
        + 0.05 * abs(a.h1_count - b.h1_count) / max(1, a.h1_count + b.h1_count)
        + 0.05 * abs(a.h2_count - b.h2_count) / max(1, a.h2_count + b.h2_count)
    )
