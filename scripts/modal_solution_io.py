import numpy as np
from numpy.polynomial.legendre import legval

def read_modal_solution(filename: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read DGmini modal output.

    Expected columns per data row:
        element_index  x_left  x_right  c_0  c_1 ... c_p

    Returns
    -------
    x_left : (n_elements) array
    x_right: (n_elements) array
    coeffs : (n_elements, n_coeffs) array
    """
    rows = []

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()
            values = [float(x) for x in parts]
            rows.append(values)

    if not rows:
        raise ValueError(f"No data rows found in file: {filename}")

    data = np.asarray(rows, dtype=float)

    if data.shape[1] < 4:
        raise ValueError("Expected at least 4 columns: element x_left x_right c_0")

    x_left  = data[:, 1]
    x_right = data[:, 2]
    coeffs  = data[:, 3:]

    return x_left, x_right, coeffs


def reconstruct_on_element(   
    x_left: float,
    x_right: float,
    coeffs: np.ndarray,
    n_plot_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Reconstruct modal polynomial on one physical element [x_left, x_right],
    coeffs are Legendre coefficients on reference element xi in [-1, 1].

    Returns
    -------
    x : (n_plot_points) array
    u : (n_plot_points) array
    """
    xi = np.linspace(-1.0, 1.0, n_plot_points)

    # Map xi -> x
    x = 0.5 * (x_right - x_left) * xi + 0.5 * (x_right + x_left)

    # Evaluate u_h(xi) = sum_j c_j P_j(xi)
    u = legval(xi, coeffs)

    return x, u


def reconstruct_global_solution(
    x_left: np.ndarray,
    x_right: np.ndarray,
    coeffs: np.ndarray,
    n_plot_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Reconstruct the full piecewise-polynomial solution over all elements

    Returns
    -------
    x : (n_elements * n_plot_points) array
    u : (n_elements * n_plot_points) array
    """
    x_all = []
    u_all = []

    n_elements = len(x_left)

    for e in range(n_elements):
        x_e, u_e = reconstruct_on_element(
            x_left[e], x_right[e], coeffs[e], n_plot_points
        )

        x_all.append(x_e)
        u_all.append(u_e)

    return np.concatenate(x_all), np.concatenate(u_all)