"""Télécharge les génomes de référence depuis NCBI (GenBank).

Usage :
    python scripts/download_genomes.py --email prenom.nom@exemple.fr

NCBI demande une adresse e-mail pour identifier les requêtes Entrez.
Les fichiers sont écrits dans data/ au format GenBank.
"""

import argparse
import time
from pathlib import Path

from Bio import Entrez

# Jeu de données par défaut : deux levures proches, une levure plus
# éloignée et deux bactéries. On s'attend à ce que le dendrogramme
# regroupe les levures entre elles et sépare les bactéries.
#
# Chaque entrée est soit un numéro d'accession RefSeq, soit une requête
# Entrez : dans ce cas, le premier résultat est utilisé et affiché, pour
# que l'on sache toujours exactement quelle séquence a été analysée.
GENOMES = {
    # Chromosome IV (le plus long) de chaque levure, pour un
    # téléchargement raisonnable tout en restant représentatif.
    "S_cerevisiae_chrIV": "NC_001136",
    "S_paradoxus_chrIV": ("Saccharomyces paradoxus[Organism] AND "
                          "chromosome IV[Title] AND complete sequence[Title]"),
    "K_lactis_chrD": ("Kluyveromyces lactis[Organism] AND "
                      "chromosome D[Title] AND complete sequence[Title]"),
    "E_coli_K12": "NC_000913",
    "B_subtilis_168": "NC_000964",
}


def resolve(entry):
    """Renvoie un identifiant NCBI à partir d'une accession ou d'une requête."""
    if "[" not in entry:
        return entry
    with Entrez.esearch(db="nuccore", term=entry, retmax=1) as handle:
        ids = Entrez.read(handle)["IdList"]
    if not ids:
        raise RuntimeError(f"Aucun résultat NCBI pour : {entry}")
    return ids[0]


def download(accession, output_path):
    """Télécharge un enregistrement GenBank complet (avec la séquence)."""
    with Entrez.efetch(db="nuccore", id=accession,
                       rettype="gbwithparts", retmode="text") as handle:
        output_path.write_text(handle.read(), encoding="utf-8")


def main():
    """Télécharge chaque génome de la liste s'il n'existe pas déjà."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--email", required=True,
                        help="Adresse e-mail transmise à NCBI Entrez")
    parser.add_argument("-o", "--output", type=Path, default=Path("data"))
    args = parser.parse_args()

    Entrez.email = args.email
    args.output.mkdir(parents=True, exist_ok=True)
    for name, entry in GENOMES.items():
        path = args.output / f"{name}.gb"
        if path.exists():
            print(f"{name} : déjà présent, ignoré")
            continue
        accession = resolve(entry)
        download(accession, path)
        # Première ligne LOCUS : vérifie quelle séquence a été obtenue.
        locus = path.read_text(encoding="utf-8").splitlines()[0]
        print(f"{name} : {locus}")
        # NCBI limite à 3 requêtes par seconde sans clé API.
        time.sleep(0.5)
    print("Terminé.")


if __name__ == "__main__":
    main()
