"""Pipeline en ligne de commande : CGR, FCGR et dendrogramme UPGMA.

Exemple :
    python -m cgr_signatures.cli data/ -k 6 --metric cosine -o results/
"""

import argparse
from pathlib import Path

import numpy as np

from .cgr import cgr_points, plot_cgr
from .fcgr import fcgr_matrix, plot_fcgr
from .io import list_genome_files, read_sequence
from .phylogeny import (
    METRICS,
    distance_matrix,
    plot_dendrogram,
    plot_distance_heatmap,
    signature_vectors,
    upgma_tree,
)


def parse_args(argv=None):
    """Lit les arguments de la ligne de commande."""
    parser = argparse.ArgumentParser(
        description="Signatures génomiques par CGR et clustering UPGMA.")
    parser.add_argument("genome_dir", type=Path,
                        help="Dossier contenant des fichiers GenBank/FASTA")
    parser.add_argument("-k", type=int, default=6,
                        help="Taille des k-mers (défaut : 6)")
    parser.add_argument("--metric", choices=METRICS, default="cosine",
                        help="Distance entre signatures (défaut : cosine)")
    parser.add_argument("-o", "--output", type=Path, default=Path("results"),
                        help="Dossier de sortie (défaut : results/)")
    parser.add_argument("--max-length", type=int, default=None,
                        help="Tronque chaque génome (tests rapides)")
    parser.add_argument("--no-cgr", action="store_true",
                        help="Ne pas tracer les nuages de points CGR")
    return parser.parse_args(argv)


def main(argv=None):
    """Exécute la pipeline complète et écrit les résultats sur disque."""
    args = parse_args(argv)
    files = list_genome_files(args.genome_dir)
    if len(files) < 2:
        raise SystemExit("Il faut au moins 2 génomes pour les comparer.")
    args.output.mkdir(parents=True, exist_ok=True)

    labels, matrices = [], []
    print(f"Pipeline CGR : {len(files)} génomes, k={args.k}, "
          f"distance={args.metric}")
    for path in files:
        name = path.stem
        sequence = read_sequence(path, max_length=args.max_length)
        print(f"  {name:<30} {len(sequence):>12,} pb")
        if not args.no_cgr:
            plot_cgr(cgr_points(sequence), args.output / f"{name}_cgr.png",
                     title=name)
        matrix = fcgr_matrix(sequence, args.k)
        plot_fcgr(matrix, args.k, args.output / f"{name}_fcgr_k{args.k}.png",
                  title=name)
        labels.append(name)
        matrices.append(matrix)

    vectors = signature_vectors(matrices)
    dist, condensed = distance_matrix(vectors, metric=args.metric)
    header = ",".join(["genome"] + labels)
    rows = [",".join([label] + [f"{v:.6f}" for v in dist[i]])
            for i, label in enumerate(labels)]
    (args.output / "distance_matrix.csv").write_text(
        "\n".join([header] + rows) + "\n", encoding="utf-8")
    plot_distance_heatmap(dist, labels, args.output / "distance_matrix.png",
                          metric=args.metric)
    plot_dendrogram(upgma_tree(condensed), labels,
                    args.output / "upgma_dendrogram.png",
                    metric=args.metric, k=args.k)
    np.save(args.output / f"fcgr_k{args.k}.npy", vectors)
    print(f"Résultats écrits dans {args.output}/")


if __name__ == "__main__":
    main()
