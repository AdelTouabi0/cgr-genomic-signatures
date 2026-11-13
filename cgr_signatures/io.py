"""Lecture de séquences ADN depuis des fichiers GenBank ou FASTA."""

from pathlib import Path

from Bio import SeqIO

# Extensions reconnues pour chaque format Biopython.
FORMATS = {
    ".gb": "genbank", ".gbk": "genbank", ".genbank": "genbank",
    ".fa": "fasta", ".fasta": "fasta", ".fna": "fasta",
}


def read_sequence(filepath, max_length=None):
    """Lit un fichier GenBank ou FASTA et renvoie la séquence concaténée.

    Tous les enregistrements du fichier (chromosomes, plasmides...) sont
    concaténés et mis en majuscules. ``max_length`` permet de tronquer la
    séquence pour accélérer les tests sur de gros génomes.
    """
    path = Path(filepath)
    fmt = FORMATS.get(path.suffix.lower())
    if fmt is None:
        raise ValueError(f"Format non reconnu pour {path.name}")
    parts = [str(record.seq).upper() for record in SeqIO.parse(path, fmt)]
    sequence = "".join(parts)
    if max_length is not None:
        sequence = sequence[:max_length]
    return sequence


def list_genome_files(directory):
    """Renvoie la liste triée des fichiers de séquences d'un dossier."""
    directory = Path(directory)
    return sorted(p for p in directory.iterdir()
                  if p.suffix.lower() in FORMATS)
