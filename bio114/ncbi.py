"""Accès aux séquences biologiques RÉELLES via NCBI E-utilities.

Expérience n°115 : on quitte l'in silico pur. Les séquences sont
téléchargées depuis les bases publiques NCBI (RefSeq protein), sans
inscription (3 requêtes/s max, politesse respectée).

Source : https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
"""

import re
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# politesse NCBI : max 3 requêtes/seconde sans clé API
_DELAY = 0.4
_last_call = 0.0

# Cibles de l'Expérience 115 : (nom, gène, organisme)
TARGETS = {
    "human_sirt6":   ("SIRT6",   "Homo sapiens"),
    "human_prestin": ("SLC26A5", "Homo sapiens"),   # Prestin = SLC26A5
    "human_cirbp":   ("CIRBP",   "Homo sapiens"),
    "bat_sirt6":     ("SIRT6",   "Myotis lucifugus"),
    "bat_prestin":   ("SLC26A5", "Myotis lucifugus"),
    "bat_cirbp":     ("CIRBP",   "Myotis lucifugus"),
}


def _get(url: str) -> str:
    global _last_call
    wait = _DELAY - (time.time() - _last_call)
    if wait > 0:
        time.sleep(wait)
    _last_call = time.time()
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode()


def esearch_ids(gene: str, organism: str, db: str = "protein") -> list[str]:
    """Recherche les identifiants RefSeq d'un gène chez un organisme."""
    term = f'{gene}[Gene Name] AND "{organism}"[Organism] AND refseq[filter]'
    url = (f"{EUTILS}/esearch.fcgi?db={db}&retmax=5&term="
           + urllib.parse.quote(term))
    xml = _get(url)
    return re.findall(r"<Id>(\d+)</Id>", xml)


def efetch_fasta(ids: list[str], db: str = "protein") -> str:
    """Télécharge les séquences FASTA (protéines) correspondantes."""
    if not ids:
        return ""
    url = f"{EUTILS}/efetch.fcgi?db={db}&id={','.join(ids)}&rettype=fasta&retmode=text"
    return _get(url)


def parse_fasta(fasta: str) -> list[dict]:
    """Parse un FASTA multi-entrées → [{header, sequence}]."""
    entries = []
    header, seq_lines = None, []
    for line in fasta.splitlines():
        if line.startswith(">"):
            if header:
                entries.append({"header": header,
                                "sequence": "".join(seq_lines)})
            header, seq_lines = line[1:], []
        else:
            seq_lines.append(line.strip())
    if header:
        entries.append({"header": header, "sequence": "".join(seq_lines)})
    return entries


def fetch_target(name: str, gene: str, organism: str) -> dict:
    """Récupère la première séquence RefSeq d'un gène/organisme."""
    ids = esearch_ids(gene, organism)
    if not ids:
        return {"name": name, "gene": gene, "organism": organism,
                "found": False, "sequence": None}
    fasta = efetch_fasta(ids[:1])
    entries = parse_fasta(fasta)
    if not entries:
        return {"name": name, "gene": gene, "organism": organism,
                "found": False, "sequence": None}
    e = entries[0]
    accession = e["header"].split()[0]
    return {"name": name, "gene": gene, "organism": organism,
            "found": True, "accession": accession, "header": e["header"],
            "length": len(e["sequence"]), "sequence": e["sequence"]}


def protein_to_dna_proxy(protein_seq: str) -> str:
    """Convertit une séquence protéique en proxy ADN pour le moteur
    topologique (chaque acide aminé → son codon le plus fréquent humain).

    La topologie travaille sur des chaînes ACGT ; on encode la protéine
    réelle en codons pour conserver l'information de séquence dans le
    formalisme de l'Exp. 114.
    """
    codon_for = {
        "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT",
        "Q": "CAA", "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT",
        "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCT",
        "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT",
        "*": "TAA",
    }
    return "".join(codon_for.get(aa, "GGT") for aa in protein_seq)
