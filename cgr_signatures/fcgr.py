"""Frequency CGR (FCGR) : matrice des fréquences de k-mers.

Le carré CGR est découpé en une grille 2^k x 2^k. Chaque case correspond à
un unique k-mer, et sa valeur est la fréquence de ce k-mer dans le génome.
On obtient ainsi une « image » de taille fixe, quelle que soit la longueur
du génome : c'est la signature génomique.
"""

from functools import lru_cache
from itertools import product

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from .cgr import CORNERS

matplotlib.use("Agg")

BASES = "ACGT"
# Code entier de chaque base, -1 pour toute autre lettre.
_ENCODING = np.full(256, -1, dtype=np.int64)
for _i, _base in enumerate(BASES):
    _ENCODING[ord(_base)] = _i


def kmer_cell(kmer):
    """Renvoie la case (ligne, colonne) d'un k-mer dans la grille FCGR.

    On rejoue la CGR sur le k-mer depuis le centre du carré : après k
    étapes, le point tombe dans la case qui correspond à ce k-mer.
    La ligne 0 est en haut (comme dans une image).
    """
    k = len(kmer)
    size = 2 ** k
    x, y = 0.5, 0.5
    for base in kmer:
        cx, cy = CORNERS[base]
        x = (x + cx) / 2
        y = (y + cy) / 2
    col = min(int(x * size), size - 1)
    row = min(int((1 - y) * size), size - 1)
    return row, col


@lru_cache(maxsize=None)
def _cell_lookup(k):
    """Table code entier de k-mer -> (ligne, colonne), calculée une fois."""
    rows = np.empty(4 ** k, dtype=np.int64)
    cols = np.empty(4 ** k, dtype=np.int64)
    for code, letters in enumerate(product(BASES, repeat=k)):
        rows[code], cols[code] = kmer_cell("".join(letters))
    return rows, cols


def kmer_label_grid(k):
    """Renvoie une grille 2^k x 2^k contenant le k-mer de chaque case."""
    size = 2 ** k
    grid = np.empty((size, size), dtype=object)
    for letters in product(BASES, repeat=k):
        kmer = "".join(letters)
        grid[kmer_cell(kmer)] = kmer
    return grid


def count_kmers(sequence, k):
    """Compte les k-mers de la séquence (fenêtre glissante, pas de 1).

    Renvoie un tableau de longueur 4^k indexé par le code entier du k-mer
    (A=0, C=1, G=2, T=3 en base 4). Les k-mers contenant une base
    ambiguë (N...) sont ignorés.
    """
    if k < 1:
        raise ValueError("k doit être supérieur ou égal à 1")
    codes = _ENCODING[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
    n = len(codes) - k + 1
    if n <= 0:
        return np.zeros(4 ** k, dtype=np.int64)
    # Code de chaque k-mer : somme des bases pondérées par 4^position.
    kmer_codes = np.zeros(n, dtype=np.int64)
    valid = np.ones(n, dtype=bool)
    for i in range(k):
        window = codes[i:i + n]
        kmer_codes = kmer_codes * 4 + window
        valid &= window >= 0
    return np.bincount(kmer_codes[valid], minlength=4 ** k)


def fcgr_matrix(sequence, k, normalize=True):
    """Construit la matrice FCGR (2^k x 2^k) de la séquence.

    Avec ``normalize=True``, la matrice contient des fréquences (somme à 1),
    ce qui permet de comparer des génomes de tailles différentes.
    """
    counts = count_kmers(sequence, k)
    rows, cols = _cell_lookup(k)
    size = 2 ** k
    matrix = np.zeros((size, size), dtype=float)
    matrix[rows, cols] = counts
    total = matrix.sum()
    if normalize and total > 0:
        matrix /= total
    return matrix


def plot_fcgr(matrix, k, output_path, title=""):
    """Enregistre la heatmap FCGR (échelle logarithmique) en PNG.

    Le log(1 + x) rend visibles les k-mers rares, écrasés sinon par les
    plus fréquents. Pour k <= 3, chaque case est annotée de son k-mer.
    """
    # Passage en fréquences « par million » pour un log lisible,
    # que la matrice contienne des comptes ou des fréquences.
    total = matrix.sum()
    values = np.log1p(matrix / total * 1e6) if total > 0 else matrix
    fig, ax = plt.subplots(figsize=(6, 5.2))
    image = ax.imshow(values, cmap="magma", interpolation="nearest")
    fig.colorbar(image, ax=ax, label="log(1 + fréquence x 10^6)")
    if k <= 3:
        labels = kmer_label_grid(k)
        for (row, col), kmer in np.ndenumerate(labels):
            ax.text(col, row, kmer, ha="center", va="center",
                    fontsize=9 if k < 3 else 6, color="white")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"FCGR k={k} {title}".strip())
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
