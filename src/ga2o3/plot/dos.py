from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from ga2o3.io.band import read_band
from ga2o3.io.doss import read_doss
from ga2o3.util import mpl_rc


_DEFAULT_RC = {
    "font.family": "serif",
    "font.size": 11,
    "axes.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
}


def plot_spin_resolved_dos(
    alpha_path: str | Path,
    beta_path: str | Path,
    *,
    e_range: Tuple[float, float] = (-0.5, 0.0),
    dos_range: Tuple[float, float] = (-1000.0, 1000.0),
    drop_last_projection: bool = True,
    labels_proj: Optional[List[str]] = None,
    colors: Optional[List[str]] = None,
    figsize: Tuple[float, float] = (4, 6),
    show_fermi: bool = True,
    ax: Optional[plt.Axes] = None,
    rc_params: Optional[dict] = None,
):
    """Spin-resolved DOS (ALPHA filled, BETA dashed) for multiple projections.

    Parameters
    ----------
    alpha_path, beta_path:
        Paths to *_ALPHA.DOSS and *_BETA.DOSS, or the same file twice if NSPIN=1.
    drop_last_projection:
        Many CRYSTAL files include a final "total" column; set True to drop it.

    Returns
    -------
    (fig, ax)
    """
    rc = dict(_DEFAULT_RC)
    if rc_params:
        rc.update(rc_params)

    with mpl_rc(rc):
        doss_a = read_doss(alpha_path)
        doss_b = read_doss(beta_path)

        if not np.allclose(doss_a.energies_hartree, doss_b.energies_hartree):
            raise ValueError("Energy grids in ALPHA and BETA files do not match.")

        cols_a = doss_a.dos
        cols_b = doss_b.dos

        if alpha_path == beta_path and (doss_a.nspin == 1 or doss_a.nspin is None):
            # single-spin (unpolarized) file: split equally for visual symmetry
            cols_a = cols_a / 2.0
            cols_b = cols_a.copy()

        if drop_last_projection and cols_a.shape[1] > 1:
            cols_a = cols_a[:, :-1]
            cols_b = cols_b[:, :-1]

        n_proj = cols_a.shape[1]

        energies_ev = doss_a.energies_ev

        if labels_proj is None:
            labels_proj = [f"Proj {i+1}" for i in range(n_proj)]
        if colors is None:
            base = ["blue", "#d62728", "#16d5d5", "#7a8831"]
            colors = (base + [f"C{i}" for i in range(len(base), n_proj)])[:n_proj]

        if ax is None:
            fig, ax1 = plt.subplots(figsize=figsize)
        else:
            ax1 = ax
            fig = ax.figure

        for i in range(n_proj):
            ya = cols_a[:, i]
            yb = cols_b[:, i]
            ax1.fill_betweenx(energies_ev, 0.0, ya, alpha=1.0, color=colors[i], linewidth=0.0)
            ax1.plot(yb, energies_ev, color=colors[i], linewidth=1.0, linestyle="--")

        ax1.set_xlim(dos_range[0], dos_range[1])
        ax1.set_ylim(e_range[0], e_range[1])

        ax1.set_xlabel("")
        ax1.set_ylabel("")
        ax1.tick_params(top=True, right=True, labelbottom=False, labelleft=False)
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])

        if show_fermi:
            ax1.axhline(0.0, linestyle="--", linewidth=1.0, color="black")

        ax1.grid(False)
        fig.tight_layout()
        return fig, ax1


def plot_band_and_dos(
    alpha_band_file: str | Path,
    beta_band_file: str | Path,
    alpha_dos_file: str | Path,
    beta_dos_file: str | Path,
    *,
    e_min: float = -1.0,
    e_max: float = 5.0,
    energy_shift: float = 0.0,
    fermi_energy_ev: float = 0.0,
    dos_xmin: float = -1000.0,
    dos_xmax: float = 1000.0,
    max_projections: int = 4,
    proj_labels: Optional[Sequence[str]] = None,
    proj_colors: Optional[Sequence[str]] = None,
    fill_alpha: float = 0.8,
    band_legend: bool = False,
    figsize: Tuple[float, float] = (8, 6),
    width_ratios: Tuple[float, float] = (4, 4),
    rc_params: Optional[dict] = None,
):
    """Combined band structure + DOS panel (1x2), sharing the energy axis.

    This is a generalized version of your notebook plotting routine:
    - left: bands (alpha solid orange, beta dashed black)
    - right: DOS (alpha filled + line, beta dashed)
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

        doss_a = read_doss(alpha_dos_file)
        doss_b = read_doss(beta_dos_file)

        cols_a = doss_a.dos
        cols_b = doss_b.dos
        energies_ev = doss_a.energies_ev

        if Path(alpha_dos_file) == Path(beta_dos_file) and (doss_a.nspin == 1 or doss_a.nspin is None):
            cols_a = cols_a / 2.0
            cols_b = cols_a.copy()

        # trim projections
        cols_a = cols_a[:, :max_projections]
        cols_b = cols_b[:, :max_projections]
        n_proj = cols_a.shape[1]

        if proj_labels is None:
            proj_labels = [f"Proj {i+1}" for i in range(n_proj)]
        if proj_colors is None:
            base = ["blue", "#d62728", "#16d5d5", "#7a8831"]
            proj_colors = (list(base) + [f"C{i}" for i in range(len(base), n_proj)])[:n_proj]

        fig, (ax_band, ax_dos) = plt.subplots(
            1,
            2,
            figsize=figsize,
            sharey=True,
            gridspec_kw={"width_ratios": width_ratios},
        )

        for ax in (ax_band, ax_dos):
            ax.minorticks_on()
            ax.tick_params(
                axis="both",
                which="both",
                length=4,
                width=1.0,
                labelbottom=False,
                labelleft=False,
                labelright=False,
                labeltop=False,
            )

        # --- bands
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
            ax_band.set_xlim(float(np.min(k_a)), float(np.max(k_a)))

        ax_band.set_ylim(e_min, e_max)
        ax_band.axhline(fermi_energy_ev, lw=0.8, ls="--", color="k")

        if band_legend:
            handles = [
                Line2D([0], [0], color="orange", lw=0.9, ls="-", label="Spin-up (α)"), 
                Line2D([0], [0], color="k", lw=0.8, ls="--", label="Spin-down (β)"), 
            ]
            ax_band.legend(handles=handles, loc="upper right", frameon=True, framealpha=0.95, fontsize=12)

        # --- DOS
        ax_dos.set_xlim(dos_xmin, dos_xmax)
        ax_dos.set_ylim(e_min, e_max)
        ax_dos.axhline(fermi_energy_ev, lw=0.8, ls="--", color="k")

        for i in range(n_proj):
            ya = cols_a[:, i]
            yb = cols_b[:, i]

            ax_dos.fill_betweenx(energies_ev, 0, ya, color=proj_colors[i], alpha=fill_alpha, linewidth=0)
            ax_dos.plot(ya, energies_ev, color=proj_colors[i], lw=1.0)
            ax_dos.plot(yb, energies_ev, color=proj_colors[i], lw=1.0, ls="--")

        fig.tight_layout()
        return fig, (ax_band, ax_dos)
