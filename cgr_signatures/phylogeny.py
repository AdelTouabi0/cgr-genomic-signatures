"""Comparaison de génomes à partir de leurs signatures FCGR."""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform

matplotlib.use("Agg")

METRICS = ("euclidean", "cosine", "cityblock")


def signature_vectors(matrices):
    """Aplatit une liste de matrices FCGR en un tableau (n_genomes, 4^k)."""
    return np.array([matrix.ravel() for matrix in matrices])


def distance_matrix(vectors, metric="cosine"):
    """Renvoie la matrice de distances carrée et sa forme condensée.

    Un vecteur nul (génome vide ou trop court) rendrait la distance cosinus
    indéfinie : on lève une erreur explicite plutôt que de produire des NaN.
    """
    if metric not in METRICS:
        raise ValueError(f"Métrique inconnue : {metric} (choix : {METRICS})")
    empty = np.where(~vectors.any(axis=1))[0]
    if len(empty):
        raise ValueError(f"Signature vide pour les génomes d'indice {list(empty)}")
    condensed = pdist(vectors, metric=metric)
    return squareform(condensed), condensed


def upgma_tree(condensed):
    """Construit l'arbre UPGMA (méthode « average » de scipy)."""
    return linkage(condensed, method="average")


def plot_dendrogram(tree, labels, output_path, metric="cosine", k=None):
    """Enregistre le dendrogramme UPGMA en PNG."""
    fig, ax = plt.subplots(figsize=(8, 0.6 * len(labels) + 2))
    dendrogram(tree, labels=labels, orientation="left",
               leaf_font_size=10, ax=ax, color_threshold=0)
    k_text = f", k={k}" if k else ""
    ax.set_title(f"Dendrogramme UPGMA (distance {metric}{k_text})")
    ax.set_xlabel("Distance")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_distance_heatmap(dist, labels, output_path, metric="cosine"):
    """Enregistre la matrice de distances sous forme de heatmap annotée."""
    fig, ax = plt.subplots(figsize=(1.2 * len(labels) + 2,
                                    1.0 * len(labels) + 1.5))
    image = ax.imshow(dist, cmap="viridis")
    fig.colorbar(image, ax=ax, label=f"Distance {metric}")
    ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), labels)
    for (i, j), value in np.ndenumerate(dist):
        ax.text(j, i, f"{value:.3f}", ha="center", va="center",
                color="white" if value < dist.max() / 2 else "black",
                fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
