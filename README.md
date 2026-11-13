# Signatures génomiques par Chaos Game Representation

![tests](https://github.com/AdelTouabi0/cgr-genomic-signatures/actions/workflows/tests.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Licence MIT](https://img.shields.io/badge/licence-MIT-green)

> *Alignment-free genome comparison: Chaos Game Representation, k-mer frequency
> images (FCGR) and UPGMA clustering.*

Comparer des génomes sans les aligner : chaque génome est transformé en une
image de fréquences de k-mers (sa **signature génomique**), puis les images
sont comparées par une distance et regroupées en arbre UPGMA.

## Principe

1. **CGR (Chaos Game Representation)** : chaque base est un coin du carré
   unité (A en haut à gauche, T en haut à droite, C en bas à gauche, G en bas
   à droite). On part du centre et, pour chaque base lue, on avance à
   mi-chemin vers son coin. Le nuage de points obtenu est une fractale
   caractéristique du génome.
2. **FCGR (Frequency CGR)** : après *k* déplacements, la position du point
   encode exactement les *k* dernières bases. En découpant le carré en une
   grille 2^k x 2^k, chaque case correspond donc à un k-mer, et on y compte
   ses occurrences. On obtient une image de taille fixe, comparable entre
   génomes de longueurs différentes.
3. **Comparaison** : les images normalisées sont aplaties en vecteurs, on
   calcule une distance (cosinus par défaut), puis on construit un arbre par
   UPGMA (`average linkage`).

Le lien entre CGR et k-mers est vérifié par un test unitaire
(`test_point_after_k_steps_falls_in_kmer_cell`).

## Installation

```bash
git clone https://github.com/AdelTouabi0/cgr-genomic-signatures.git
cd cgr-genomic-signatures
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Utilisation

```bash
# 1. Télécharger le jeu de données (3 levures, 2 bactéries) depuis NCBI
python scripts/download_genomes.py --email prenom.nom@exemple.fr

# 2. Lancer la pipeline
python -m cgr_signatures.cli data/ -k 6 --metric cosine -o results/
```

Options principales : `-k` (taille des k-mers), `--metric`
(`cosine`, `euclidean`, `cityblock`), `--max-length` (tronque les génomes pour
un test rapide), `--no-cgr` (saute les nuages de points).

Sorties dans `results/` : nuage CGR et heatmap FCGR par génome, matrice de
distances (CSV et image) et dendrogramme UPGMA.

## Résultats

*Section à compléter après exécution sur les génomes NCBI.*

Hypothèse testée : les levures (*S. cerevisiae*, *S. paradoxus*,
*K. lactis*) doivent se regrouper entre elles, et les deux bactéries
(*E. coli*, *B. subtilis*) former un groupe séparé, avec *S. cerevisiae* et
*S. paradoxus* (même genre) les plus proches.

## Structure

```
cgr_signatures/
├── io.py          # lecture GenBank / FASTA
├── cgr.py         # coordonnées CGR et nuage de points
├── fcgr.py        # comptage vectorisé des k-mers et matrice FCGR
├── phylogeny.py   # distances, UPGMA, figures
└── cli.py         # pipeline en ligne de commande
scripts/download_genomes.py
tests/test_fcgr.py
```

## Tests

```bash
pytest -q
```

## Contexte

Projet universitaire de bio-informatique (2026), réécrit et restructuré pour
ce portfolio : comptage des k-mers vectorisé avec NumPy, gestion des bases
ambiguës, tests unitaires et intégration continue. Développé avec
l'assistance d'outils d'IA (Claude).

**Référence** : Jeffrey H. J. (1990). Chaos game representation of gene
structure. *Nucleic Acids Research*, 18(8), 2163-2170.

## Auteur

Adel Touabi, Master 1 Bio-informatique, Université Paris Cité.
