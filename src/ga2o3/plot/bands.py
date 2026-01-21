from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from ga2o3.io.band import read_band
from ga2o3.util import mpl_rc


_DEFAULT_RC = {
    "font.family": "serif",
    "font.size": 11,
    "axes.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
}


def plot_bands_only(
    alpha_band_file: str | Path,
    beta_band_file: str | Path,
    *,
    e_min: float = -1.0,
    e_max: float = 5.0,
    energy_shift: float = 0.0,
    band_legend: bool = True,
    figsize: Tuple[float, float] = (6, 5),
    ax: Optional[plt.Axes] = None,
    rc_params: Optional[dict] = None,
):
    """Plot only spin-polarized bands (alpha solid, beta dashed).

    - All axis labels and tick numbers are removed, but tick marks are kept.
    - High-symmetry vertical lines are drawn if the file contains XAXIS ticks.

    Returns
    -------
    (fig, ax)
    """
    rc = dict(_DEFAULT_RC)
    if rc_params:
        rc.update(rc_params)

    with mpl_rc(rc):
        band_a = read_band(alpha_band_file)
        band_b = read_band(beta_band_file)

        k_a, E_a = band_a.k, band_a.energies_ev
        k_b, E_b = band_b.k, band_b.energies_ev

        if E_a.size:
            E_a = E_a - energy_shift
        if E_b.size:
            E_b = E_b - energy_shift

        if ax is None:
            fig, ax_band = plt.subplots(figsize=figsize)
        else:
            ax_band = ax
            fig = ax.figure

        ax_band.minorticks_on()
        ax_band.tick_params(
            axis="both",
            which="both",
            labelbottom=False,
            labelleft=False,
            labelright=False,
            labeltop=False,
        )

        if E_a.size:
            for ib in range(E_a.shape[1]):
                ax_band.plot(k_a, E_a[:, ib], "-", color="orange", lw=0.9)
        if E_b.size:
            for ib in range(E_b.shape[1]):
                ax_band.plot(k_b, E_b[:, ib], "--", color="k", lw=0.8)

        hs_x = band_a.xaxis_ticks or band_b.xaxis_ticks
        if hs_x:
            ax_band.set_xlim(0, hs_x[-1])
            for x in hs_x:
                ax_band.axvline(x, lw=0.6, alpha=0.35, color="k")
        else:
            if E_a.size:
                ax_band.set_xlim(float(np.min(k_a)), float(np.max(k_a)))
            else:
                ax_band.set_xlim(float(np.min(k_b)), float(np.max(k_b)))

        ax_band.set_ylim(e_min, e_max)

        if band_legend:
            handles = [
                Line2D([0], [0], color="orange", lw=0.9, ls="-", label="Spin up (α)"), 
                Line2D([0], [0], color="k", lw=0.8, ls="--", label="Spin down (β)"), 
            ]
            ax_band.legend(
                handles=handles,
                loc="upper right",
                frameon=True,
                framealpha=0.95,
                fontsize=12,
            )

        fig.tight_layout()
        return fig, ax_band
