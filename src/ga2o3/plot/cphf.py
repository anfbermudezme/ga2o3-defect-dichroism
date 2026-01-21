from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Literal, Optional, Sequence

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import LogFormatterMathtext, LogLocator, NullFormatter, ScalarFormatter

from ga2o3.constants import HC_EVN_M
from ga2o3.dataset import Ga2O3Dataset
from ga2o3.io.cphf import add_absorption_columns, apply_energy_correction, parse_eps_components
from ga2o3.util import ensure_parent_dir, mpl_rc


CorrectionMode = Literal["none", "scissor", "ratio"]
YScale = Literal["linear", "log"]


@dataclass(frozen=True)
class CaseStyle:
    color: str
    marker: str
    label: Optional[str] = None


DEFAULT_STYLES: Dict[str, CaseStyle] = {
    "pristine": CaseStyle(color="k", marker="o", label="Pristine"),
    "tetra": CaseStyle(color="darkblue", marker="^", label=r"VGa$_I$"),
    "octa": CaseStyle(color="hotpink", marker="s", label=r"VGa$_{II}$"),
}


_DEFAULT_RC = {
    "font.family": "serif",
    "axes.unicode_minus": True,
    "axes.linewidth": 1.5,
}


def _case_key(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "")


def _load_case_df(
    ds: Ga2O3Dataset,
    case: str,
    *,
    pattern: str = "*.out",
    wl_min: Optional[float] = 465.0,
    wl_min_apply_to: Literal["orig", "corrected"] = "orig",
    round_wavelength: bool = True,
    correction: CorrectionMode = "ratio",
    scissor_shift_ev: float = 0.0,
    ratio_factor: float = 1.0,
    extra_shift_ev: float = 0.0,
) -> "pd.DataFrame":
    import pandas as pd

    files = ds.cphf_outputs(case, pattern=pattern)
    if not files:
        return pd.DataFrame()

    # oldest -> newest (newest wins when de-duplicating wavelength points)
    files.sort(key=lambda p: p.stat().st_mtime)

    by_wl_orig: Dict[float, dict] = {}

    for f in files:
        df_f = parse_eps_components(f, round_wavelength=round_wavelength)
        if df_f.empty:
            continue

        df_f = add_absorption_columns(df_f)
        df_f = apply_energy_correction(
            df_f,
            mode=correction,
            scissor_shift_ev=scissor_shift_ev,
            ratio_factor=ratio_factor,
            extra_shift_ev=extra_shift_ev,
        )

        for _, row in df_f.iterrows():
            wl_orig = float(row.get("wavelength_orig_nm", row["wavelength_nm"]))
            wl_corr = float(row["wavelength_nm"])

            if wl_min is not None:
                if wl_min_apply_to == "orig" and wl_orig < wl_min:
                    continue
                if wl_min_apply_to == "corrected" and wl_corr < wl_min:
                    continue

            by_wl_orig[wl_orig] = dict(row)

    if not by_wl_orig:
        return pd.DataFrame()

    df = pd.DataFrame([by_wl_orig[k] for k in sorted(by_wl_orig.keys())])
    df["case"] = case
    return df.sort_values("wavelength_nm").reset_index(drop=True)


def _compute_global_ylims(case_dfs: Dict[str, "pd.DataFrame"], cols: Sequence[str], y_scale: YScale):
    import pandas as pd

    vals = []
    for df in case_dfs.values():
        if df is None or df.empty:
            continue
        for col in cols:
            if col not in df.columns:
                continue
            a = df[col].to_numpy(dtype=float)
            a = a[np.isfinite(a)]
            if y_scale == "log":
                a = a[a > 0]
            if a.size:
                vals.append(a)
    if not vals:
        return None

    allv = np.concatenate(vals)
    vmin = float(np.min(allv))
    vmax = float(np.max(allv))
    if y_scale == "log":
        ymin = 10 ** math.floor(math.log10(vmin))
        ymax = 10 ** math.ceil(math.log10(vmax))
        return (ymin, ymax)
    return (0.0, 1.05 * vmax)


def _style_axis(ax: plt.Axes, *, y_scale: YScale, y_lims=None):
    ax.set_box_aspect(1)

    ax.grid(True, which="major", axis="both", linestyle="-", linewidth=0.4, alpha=0.15, color="black", zorder=0)

    for s in ax.spines.values():
        s.set_color("black")
        s.set_linewidth(1.5)

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="both", which="both", width=1.2, length=8)

    if y_scale == "log":
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
        ax.yaxis.set_minor_formatter(NullFormatter())
    else:
        sci = ScalarFormatter(useMathText=True)
        sci.set_scientific(True)
        sci.set_powerlimits((0, 0))
        ax.yaxis.set_major_formatter(sci)
        ax.set_ylim(bottom=0)

    if y_lims is not None:
        ax.set_ylim(*y_lims)


def _plot_core_shell(ax: plt.Axes, x, y, *, color: str, marker: str, y_scale: YScale, marker_size: float = 8.0):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if y_scale == "log":
        m &= (y > 0)
    x = x[m]
    y = y[m]
    if x.size == 0:
        return

    line_width = 1.5
    edge_width = 2.0
    gap_scale = 0.75
    core_scale = 0.55

    ax.plot(x, y, color=color, linewidth=line_width, linestyle="-", zorder=2)

    bg = ax.get_facecolor()

    ax.plot(x, y, linestyle="None", marker=marker, markersize=marker_size,
            markerfacecolor=color, markeredgecolor=color, markeredgewidth=edge_width, zorder=3)
    ax.plot(x, y, linestyle="None", marker=marker, markersize=marker_size * gap_scale,
            markerfacecolor=bg, markeredgecolor=bg, markeredgewidth=0.0, zorder=4)
    ax.plot(x, y, linestyle="None", marker=marker, markersize=marker_size * core_scale,
            markerfacecolor=color, markeredgecolor=color, markeredgewidth=0.0, zorder=5)


def _legend_core_shell(ax: plt.Axes, plot_order: Sequence[str], styles: Dict[str, CaseStyle]):
    bg = ax.get_facecolor()

    handles = []
    labels = []

    marker_size = 8.0
    edge_width = 2.0
    gap_scale = 0.75
    core_scale = 0.55
    line_width = 1.5

    for case in plot_order:
        st = styles.get(_case_key(case), CaseStyle(color="k", marker="o", label=case))
        label = st.label or case

        h_line = Line2D([0], [0], color=st.color, lw=line_width)
        h_shell = Line2D([0], [0], linestyle="None", marker=st.marker, markersize=marker_size,
                         markerfacecolor=st.color, markeredgecolor=st.color, markeredgewidth=edge_width)
        h_gap = Line2D([0], [0], linestyle="None", marker=st.marker, markersize=marker_size * gap_scale,
                       markerfacecolor=bg, markeredgecolor=bg, markeredgewidth=0.0)
        h_core = Line2D([0], [0], linestyle="None", marker=st.marker, markersize=marker_size * core_scale,
                        markerfacecolor=st.color, markeredgecolor=st.color, markeredgewidth=0.0)

        handles.append((h_line, h_shell, h_gap, h_core))
        labels.append(label)

    ax.legend(handles, labels, handler_map={tuple: HandlerTuple(ndivide=None)},
              frameon=True, fancybox=True, framealpha=1.0, edgecolor="black",
              fontsize=14, handlelength=2.0, handletextpad=0.6, loc="best")


def plot_absorption_overlays(
    ds: Ga2O3Dataset,
    *,
    cases: Sequence[str] = ("Pristine", "TETRA", "OCTA"),
    wl_min: Optional[float] = 465.0,
    wl_min_apply_to: Literal["orig", "corrected"] = "orig",
    correction: CorrectionMode = "ratio",
    eg_pbe: float = 2.62,
    eg_target: float = 4.43,
    ratio_factor: Optional[float] = None,
    scissor_shift_ev: Optional[float] = None,
    extra_shift_ev: float = 0.0,
    y_scale: YScale = "log",
    use_global_ylims: bool = True,
    add_legend: bool = True,
    figsize: float = 6.0,
    dpi: int = 600,
    styles: Optional[Dict[str, CaseStyle]] = None,
    rc_params: Optional[dict] = None,
) -> Dict[str, plt.Figure]:
    """Create absorption overlays (iso, xx, yy, zz) from CPHF outputs.

    Returns
    -------
    dict with keys: "iso", "xx", "yy", "zz" mapping to matplotlib Figure objects.
    """
    import pandas as pd

    styles2 = dict(DEFAULT_STYLES)
    if styles:
        styles2.update({ _case_key(k): v for k, v in styles.items() })

    if ratio_factor is None:
        ratio_factor = eg_target / eg_pbe
    if scissor_shift_ev is None:
        scissor_shift_ev = eg_target - eg_pbe

    case_dfs: Dict[str, pd.DataFrame] = {}
    for case in cases:
        ckey = _case_key(case)
        df = _load_case_df(
            ds,
            case,
            wl_min=wl_min,
            wl_min_apply_to=wl_min_apply_to,
            correction=correction,
            ratio_factor=float(ratio_factor),
            scissor_shift_ev=float(scissor_shift_ev),
            extra_shift_ev=float(extra_shift_ev),
        )
        case_dfs[case] = df

    ycols = ["alpha_iso_cm^-1", "alpha_xx_cm^-1", "alpha_yy_cm^-1", "alpha_zz_cm^-1"]
    y_lims = _compute_global_ylims(case_dfs, ycols, y_scale) if use_global_ylims else None

    rc = dict(_DEFAULT_RC)
    if rc_params:
        rc.update(rc_params)

    figs: Dict[str, plt.Figure] = {}

    with mpl_rc(rc):
        # isotropic
        fig, ax = plt.subplots(figsize=(figsize, figsize), dpi=dpi)
        for case in cases:
            df = case_dfs.get(case)
            if df is None or df.empty:
                continue
            st = styles2.get(_case_key(case), CaseStyle(color="k", marker="o", label=case))
            _plot_core_shell(ax, df["wavelength_nm"].values, df["alpha_iso_cm^-1"].values,
                             color=st.color, marker=st.marker, y_scale=y_scale)
        _style_axis(ax, y_scale=y_scale, y_lims=y_lims)
        if add_legend:
            _legend_core_shell(ax, cases, styles2)
        fig.tight_layout()
        figs["iso"] = fig

        # directional
        for comp in ("xx", "yy", "zz"):
            fig, ax = plt.subplots(figsize=(figsize, figsize), dpi=dpi)
            col = f"alpha_{comp}_cm^-1"
            for case in cases:
                df = case_dfs.get(case)
                if df is None or df.empty:
                    continue
                st = styles2.get(_case_key(case), CaseStyle(color="k", marker="o", label=case))
                _plot_core_shell(ax, df["wavelength_nm"].values, df[col].values,
                                 color=st.color, marker=st.marker, y_scale=y_scale)
            _style_axis(ax, y_scale=y_scale, y_lims=y_lims)
            if add_legend and comp == "zz":
                _legend_core_shell(ax, cases, styles2)
            fig.tight_layout()
            figs[comp] = fig

    return figs


def save_absorption_figs(figs: Dict[str, plt.Figure], outdir: str | Path, *, fmt: str = "png", dpi: int = 600):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for key, fig in figs.items():
        path = outdir / f"absorption_{key}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
