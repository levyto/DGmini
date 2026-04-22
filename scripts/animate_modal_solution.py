#!/usr/bin/env python3
#
# Animate DGmini modal output stored in files such as:
#   output/solution_0000.dat
#   output/solution_0001.dat
#   ...
#   output/solution_XXXX.dat
#   output/solution_final.dat
#
# Usage examples:
#   python3 scripts/animate_modal_solution.py
#   python3 scripts/animate_modal_solution.py --show-element-boundaries
#   python3 scripts/animate_modal_solution.py --fps 10 --save output/solution.gif

import argparse
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy.polynomial.legendre import legval

from utils_modal_solution_io import (
    read_modal_solution,
    reconstruct_on_element,
    reconstruct_global_solution,
    read_time
)

def frame_sort_key(path: Path) -> tuple[int, str]:
    """
    Sort files like:
      solution_0000.dat
      solution_0001.dat
      ...
      solution_XXXX.dat
      solution_final.dat

    Numeric frames come first in numeric order, 'final' goes last.
    """
    name = path.stem

    if name.endswith("_final"):
        return (10**12, name)

    match = re.search(r"_(\d+)$", name)
    if match:
        return (int(match.group(1)), name)

    return (10**11, name)


def collect_solution_files(pattern: str) -> list[Path]:
    files = sorted(Path().glob(pattern), key=frame_sort_key)

    if not files:
        raise FileNotFoundError(f"No files found matching pattern: {pattern}")

    return files


def load_all_frames(
    files: list[Path],
    n_plot_points: int,
) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray]:
    """
    Load and reconstruct all frames.

    Returns
    -------
    x_frames : list of x arrays
    u_frames : list of u arrays
    interfaces : array of element interfaces from first frame
    """
    x_frames = []
    u_frames = []
    interfaces = None

    for file in files:
        x_left, x_right, coeffs = read_modal_solution(str(file))
        x, u = reconstruct_global_solution(x_left, x_right, coeffs, n_plot_points)

        x_frames.append(x)
        u_frames.append(u)

        if interfaces is None:
            interfaces = np.concatenate(([x_left[0]], x_right))

    assert interfaces is not None
    return x_frames, u_frames, interfaces


def compute_global_ylim(u_frames: list[np.ndarray], pad_ratio: float = 0.05) -> tuple[float, float]:
    u_min = min(np.min(u) for u in u_frames)
    u_max = max(np.max(u) for u in u_frames)

    if np.isclose(u_min, u_max):
        pad = 1.0
    else:
        pad = pad_ratio * (u_max - u_min)

    return u_min - pad, u_max + pad


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Animate DGmini modal solution from element-wise Legendre coefficients"
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="output/solution_*.dat",
        help="Glob pattern for modal output files",
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
    parser.add_argument(
        "--fps",
        type=int,
        default=8,
        help="Frames per second for animation",
    )
    parser.add_argument(
        "--save",
        type=str,
        default="",
        help="Optional output file for saving animation (.gif or .mp4)",
    )

    args = parser.parse_args()

    files = collect_solution_files(args.pattern)
    x_frames, u_frames, interfaces = load_all_frames(files, args.points_per_element)
    times = [read_time(f) for f in files]

    x_min = np.min(x_frames[0])
    x_max = np.max(x_frames[0])
    y_min, y_max = compute_global_ylim(u_frames)

    fig, ax = plt.subplots()
    (line,) = ax.plot([], [], color="red", linewidth=1.5)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$u(x)$")
    ax.grid(True)

    if args.show_element_boundaries:
        for xb in interfaces:
            ax.axvline(x=xb, color="blue", linewidth=0.8)

    title = ax.set_title(f"t = {times[0]:.4f}")

    def init():
        line.set_data([], [])
        title.set_text(f"t = {times[0]:.4f}")
        return line, title

    def update(frame_idx: int):
        x = x_frames[frame_idx]
        u = u_frames[frame_idx]
        line.set_data(x, u)
        title.set_text(f"t = {times[frame_idx]:.4f}")
        return line, title

    interval_ms = 1000.0 / args.fps

    anim = FuncAnimation(
        fig,
        update,
        frames=len(files),
        init_func=init,
        interval=interval_ms,
        blit=False,
        repeat=True,
    )

    if args.save:
        output_path = Path(args.save)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if output_path.suffix.lower() == ".gif":
            anim.save(output_path, writer="pillow", fps=args.fps)
        elif output_path.suffix.lower() == ".mp4":
            anim.save(output_path, writer="ffmpeg", fps=args.fps)
        else:
            raise ValueError("Supported output formats are .gif and .mp4")

    plt.show()


if __name__ == "__main__":
    main()