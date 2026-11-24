"""Regenerate the single-panel boundary comparison in Figure 1."""

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
from scipy.signal import fftconvolve

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def polynomial_boundary(z, numerators, denominator):
    """Sample the conformal polynomial, scaled to give its image area pi."""
    coefficients = np.asarray(numerators) * (1.0 / denominator)
    powers = 2 * np.arange(len(coefficients)) + 3
    scale = 1 / np.sqrt(1 + np.sum(powers * coefficients**2))
    return scale * (z + sum(a * z**j for j, a in zip(powers, coefficients)))


def corner_boundary(z, beta, coefficients):
    """Integrate the Taylor series of the derivative in the variable z squared."""
    terms = 20000
    orders = 2 * np.arange(terms) + 1
    binomial = np.zeros(terms)
    exponential = np.zeros(terms)
    binomial[0] = exponential[0] = 1
    for k in range(1, terms):
        binomial[k] = binomial[k - 1] * (beta - k + 1) / k
        exponential[k] = sum(
            j * coefficients[j - 1] * exponential[k - j]
            for j in range(1, min(k, len(coefficients)) + 1)
        ) / k
    derivative = fftconvolve(binomial, exponential)[:terms]
    area_factor = np.sum(derivative**2 / orders)
    derivative /= np.sqrt(area_factor)
    return z * np.polynomial.polynomial.polyval(z * z, derivative / orders)


def make_figure(output):
    parameters = json.loads(Path(__file__).with_name("geometry_parameters.json").read_text())
    z = np.exp(2j * np.pi * np.arange(8192) / 8192)
    polynomial = polynomial_boundary(
        z, parameters["polynomial_numerators"], parameters["polynomial_denominator"]
    )
    fitted = corner_boundary(z, parameters["corner_exponent"], parameters["corner_coefficients"])
    theta = np.linspace(0, 2 * np.pi, 4097)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(5.0, 4.7), constrained_layout=True)
    ax.fill(fitted.real, fitted.imag, color="#d9e7f1", alpha=.9)
    ax.plot(np.cos(theta), np.sin(theta), color="#7f8790", ls="--", lw=1.3,
            label="Disk, same area")
    ax.plot(polynomial.real, polynomial.imag, color="#ce7b31", lw=1.4,
            label="Explicit smooth competitor")
    ax.plot(fitted.real, fitted.imag, color="#174b70", lw=1.8,
            label="Numerical cornered candidate")
    tip = fitted[len(fitted) // 4]
    ax.scatter([tip.real, -tip.real], [tip.imag, -tip.imag], s=22, color="#174b70", zorder=5)
    ax.annotate(
        r"$\omega\approx257.45^\circ$", xy=(tip.real, tip.imag), xytext=(.29, 1.12),
        arrowprops={"arrowstyle": "-", "color": "#174b70"}, color="#174b70", fontsize=11,
    )
    ax.set(xlim=(-1.33, 1.33), ylim=(-1.17, 1.26), ylabel="$y$", aspect="equal")
    ax.set_title(r"$n=N\approx5.26894043$, area $\pi$", fontsize=12)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.07), frameon=False, fontsize=10)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path, default=Path("figure1-replay.pdf"),
                        help="destination for the regenerated PDF")
    make_figure(parser.parse_args().output)
