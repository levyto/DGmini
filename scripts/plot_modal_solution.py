#!/usr/bin/env python3
#
# This script reads DGmini modal output stored in output/solution.dat
# and plots the reconstructed solution.
# 
# Usage example:
#   python3 scripts/plot_modal_solution.py --show-element-boundaries

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from modal_solution_io import (
    read_modal_solution,
    reconstruct_on_element,
    reconstruct_global_solution
)

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot DGmini modal solution from element-wise Legendre coefficients"
    )
    parser.add_argument(
        "--input_file",
        type=str,
        default='output/solution.dat',
        help="Path to DGmini modal output file",
    )
    parser.add_argument(
        "--points-per-element",
        type=int,
        default=50,
        help="Number of reconstruction points per element",
    )
    parser.add_argument(
        "--show-element-boundaries",
        action="store_true",
        help="Draw vertical lines at element interfaces",
    )

    args = parser.parse_args()

    input_path = Path(args.input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")

    x_left, x_right, coeffs = read_modal_solution(str(input_path))

    x, u = reconstruct_global_solution(
        x_left,
        x_right,
        coeffs,
        args.points_per_element,
    )

    plt.figure()
    plt.plot(x, u, color="red", linewidth=1.5)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$u(x)$")
    plt.grid(True)

    if args.show_element_boundaries:
        interfaces = np.concatenate(([x_left[0]], x_right))
        for xb in interfaces:
            plt.axvline(x=xb, color="blue", linewidth=0.8)

    plt.show()


if __name__ == "__main__":
    main()