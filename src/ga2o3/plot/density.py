from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import (
    LogNorm,
    TwoSlopeNorm,
    SymLogNorm,
    ListedColormap,
    to_rgba,
)

from ga2o3.io.f25 import read_f25_maps
from ga2o3.util import ensure_parent_dir, mpl_rc


def default_transform(arr: np.ndarray) -> np.ndarray:
    """Match the orientation used in your reference figures."""
    arr = np.flip(arr, axis=1)
    arr = np.rot90(arr, k=3)
    arr = np.flip(arr, axis=0)
    arr = np.flip(arr, axis=1)
    return arr


def make_soft_plasma(bg_hex: str = "#131315", low_frac: float = 0.7, gamma: float = 0.9) -> ListedColormap:
    """Start from 'plasma' but replace the darkest end with a soft background."""
    base = plt.colormaps["plasma"](np.linspace(0, 1, 256))
    bg = np.array(to_rgba(bg_hex))
    n = int(256 * low_frac)
    ramp = np.linspace(0.0, 1.0, n) ** gamma
    base[:n, :] = (1 - ramp)[:, None] * bg + ramp[:, None] * base[:n, :]
    return ListedColormap(base)


def plot_charge_density_f25(
    file_path: str | Path,
    out_path: Optional[str | Path] = None,
    *,
    apply_transform: bool = True,
    transform_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    charge_percentiles: Tuple[float, float] = (0.01, 99.2),
    charge_limits: Optional[Tuple[float, float]] = None,
    n_levels_charge: int = 320,
    dpi: int = 900,
    figsize: Tuple[float, float] = (3.5, 3.5),
    show_colorbar: bool = False,
    hide_colorbar_ticks: bool = True,
    cmap: Optional[ListedColormap] = None,
    rc_params: Optional[dict] = None,
):
    """Plot the first 2D map in a CASTEP .f25 file as a log-scaled contour plot."""
    rc = {"axes.linewidth": 0.8}
    if rc_params:
        rc.update(rc_params)

    with mpl_rc(rc):
        maps = read_f25_maps(file_path, max_maps=1)
        if not maps:
            raise ValueError("No MAPN blocks found in the .f25 file.")

        total = maps[0]
        if apply_transform:
            tf = transform_fn or default_transform
            total = tf(total)

        pos = total[total > 0]
        if pos.size == 0:
            raise ValueError("Charge map contains no positive values (required for LogNorm).")

        if charge_limits is None:
            p_lo, p_hi = charge_percentiles
            vmin = float(np.percentile(pos, p_lo))
            vmax = float(np.percentile(pos, p_hi))
        else:
            vmin, vmax = map(float, charge_limits)

        min_pos = float(np.min(pos))
        if vmin <= 0:
            vmin = max(min_pos, 1e-12 * vmax)
        if vmin >= vmax:
            vmin = max(min_pos, 1e-3 * vmax)

        total_clip = total.copy()
        total_clip[total_clip < vmin] = vmin
        total_clip[total_clip > vmax] = vmax

        if cmap is None:
            cmap = make_soft_plasma()

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        levels = np.geomspace(vmin, vmax, n_levels_charge)

        cf = ax.contourf(
            total_clip,
            levels=levels,
            cmap=cmap,
            norm=LogNorm(vmin=vmin, vmax=vmax),
            antialiased=True,
        )

        if show_colorbar:
            cbar = fig.colorbar(cf, ax=ax, pad=0.02, shrink=0.92)
            if hide_colorbar_ticks:
                cbar.set_ticks([])
                cbar.set_label("")
                cbar.ax.tick_params(length=0)
            cbar.outline.set_linewidth(0.6)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal", adjustable="box")
        ax.grid(False)
        fig.tight_layout()

        if out_path is not None:
            out_path = Path(out_path)
            ensure_parent_dir(out_path)
            fig.savefig(out_path, dpi=dpi, bbox_inches="tight", transparent=True)

        return fig, ax


def plot_charge_spin_overlay(
    file_path: str | Path,
    out_path: Optional[str | Path] = None,
    *,
    apply_transform: bool = True,
    transform_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    charge_percentiles: Tuple[float, float] = (0.01, 99.2),
    spin_percentile: float = 99.4,
    spin_norm_mode: str = "quantile",  # "quantile" | "minmax" | "symlog"
    symlog_linthresh_frac: float = 0.02,
    spin_alpha: float = 0.85,
    dpi: int = 900,
    figsize: Tuple[float, float] = (3.5, 3.5),
    charge_cmap: str = "plasma",
    spin_cmap: str = "seismic",
    n_levels_charge: int = 320,
    n_levels_spin: int = 320,
    show_colorbar: bool = False,
    rc_params: Optional[dict] = None,
):
    """Overlay charge (log) and spin (diverging) maps from the first two MAPN blocks."""
    rc = {"axes.linewidth": 0.8}
    if rc_params:
        rc.update(rc_params)

    with mpl_rc(rc):
        maps = read_f25_maps(file_path, max_maps=2)
        if len(maps) < 2:
            raise ValueError("Expected at least two MAPN blocks (charge + spin) in the .f25 file.")

        charge, spin = maps[0], maps[1]
        if apply_transform:
            tf = transform_fn or default_transform
            charge = tf(charge)
            spin = tf(spin)

        # --- charge scaling (LogNorm, percentile)
        pos = charge[charge > 0]
        if pos.size == 0:
            raise ValueError("Charge map contains no positive values (required for LogNorm).")
        vmin = float(np.percentile(pos, charge_percentiles[0]))
        vmax = float(np.percentile(pos, charge_percentiles[1]))
        vmin = max(vmin, float(np.min(pos)))

        charge_clip = np.clip(charge, vmin, vmax)

        # --- spin scaling
        spin_flat = spin[np.isfinite(spin)]
        if spin_flat.size == 0:
            raise ValueError("Spin map contains no finite values.")

        if spin_norm_mode == "minmax":
            smax = float(np.max(np.abs(spin_flat)))
            norm_spin = TwoSlopeNorm(vmin=-smax, vcenter=0.0, vmax=smax)
        elif spin_norm_mode == "symlog":
            smax = float(np.max(np.abs(spin_flat)))
            lin = max(1e-12, symlog_linthresh_frac * smax)
            norm_spin = SymLogNorm(linthresh=lin, vmin=-smax, vmax=smax, base=10)
        else:
            # quantile
            smax = float(np.percentile(np.abs(spin_flat), spin_percentile))
            smax = max(smax, 1e-12)
            norm_spin = TwoSlopeNorm(vmin=-smax, vcenter=0.0, vmax=smax)

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        levels_charge = np.geomspace(vmin, vmax, n_levels_charge)
        ax.contourf(
            charge_clip,
            levels=levels_charge,
            cmap=make_soft_plasma() if charge_cmap == "plasma" else charge_cmap,
            norm=LogNorm(vmin=vmin, vmax=vmax),
            antialiased=True,
            zorder=1,
        )

        # spin overlay
        # use linear-ish spacing for spin contours
        smax_plot = float(getattr(norm_spin, "vmax", np.max(np.abs(spin_flat))))
        levels_spin = np.linspace(-smax_plot, smax_plot, n_levels_spin)
        ax.contourf(
            spin,
            levels=levels_spin,
            cmap=spin_cmap,
            norm=norm_spin,
            alpha=spin_alpha,
            antialiased=True,
            zorder=2,
        )

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal", adjustable="box")
        ax.grid(False)
        fig.tight_layout()

        if out_path is not None:
            out_path = Path(out_path)
            ensure_parent_dir(out_path)
            fig.savefig(out_path, dpi=dpi, bbox_inches="tight", transparent=True)

        return fig, ax
