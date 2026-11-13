"""Chaos Game Representation (CGR) d'une séquence ADN.

Chaque nucléotide est associé à un coin du carré unité. On part du centre
(0.5, 0.5) et, pour chaque base lue, on se déplace à mi-chemin entre la
position courante et le coin de cette base. Après k étapes, la position du
point encode exactement les k dernières bases lues : c'est ce qui relie la
CGR aux fréquences de k-mers (voir ``fcgr.py``).
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

# Coins du carré unité : A en haut à gauche, T en haut à droite,
# C en bas à gauche, G en bas à droite (convention de Jeffrey, 1990).
CORNERS = {"A": (0.0, 1.0), "T": (1.0, 1.0), "C": (0.0, 0.0), "G": (1.0, 0.0)}


def cgr_points(sequence):
    """Renvoie un tableau (n, 2) des coordonnées CGR de la séquence.

    Les bases autres que A, C, G et T (N, R, Y...) sont ignorées.
    """
    x, y = 0.5, 0.5
    points = []
    for base in sequence:
        corner = CORNERS.get(base)
        if corner is None:
            continue
        x = (x + corner[0]) / 2
        y = (y + corner[1]) / 2
        points.append((x, y))
    return np.array(points)


def plot_cgr(points, output_path, title="CGR", max_points=500_000, seed=0):
    """Trace le nuage de points CGR et l'enregistre en PNG.

    Au-delà de ``max_points`` points, un sous-échantillon aléatoire est
    tracé pour garder un fichier image léger.
    """
    if len(points) > max_points:
        rng = np.random.default_rng(seed)
        points = points[rng.choice(len(points), max_points, replace=False)]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(points[:, 0], points[:, 1], s=0.05, c="black", marker=".")
    # Lettres des coins pour faciliter la lecture.
    for base, (cx, cy) in CORNERS.items():
        ax.text(cx, cy, base, fontsize=14, weight="bold",
                ha="center", va="center", color="tab:red")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
