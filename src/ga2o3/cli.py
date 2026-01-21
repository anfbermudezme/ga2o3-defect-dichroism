from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Sequence

import matplotlib

# Safe default for headless environments (clusters, CI).
matplotlib.use("Agg")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="ga2o3", description="Ga2O3 plotting utilities")
    sub = p.add_subparsers(dest="cmd", required=True)

    # --- bands ---------------------------------------------------------------
    pb = sub.add_parser("bands", help="Plot CRYSTAL band structure (alpha/beta)")
    pb.add_argument("--root", type=Path, default=None, help="Data root (defaults to $GA2O3_DATA or cwd)")
    pb.add_argument("--case", required=True, help="Case name (e.g. pristine, tetra, octa, O_1)")
    pb.add_argument("--out", required=True, type=Path, help="Output image path")
    pb.add_argument("--emin", type=float, default=-1.0)
    pb.add_argument("--emax", type=float, default=5.0)
    pb.add_argument("--shift", type=float, default=0.0, help="Energy shift to subtract (eV)")
    pb.add_argument("--legend", action="store_true", help="Show alpha/beta legend")

    # --- dos -----------------------------------------------------------------
    pd = sub.add_parser("dos", help="Plot spin-resolved DOS from CRYSTAL .DOSS files")
    pd.add_argument("--root", type=Path, default=None, help="Data root (defaults to $GA2O3_DATA or cwd)")
    pd.add_argument("--case", required=True, help="Case name (e.g. pristine, tetra, octa)")
    pd.add_argument("--out", required=True, type=Path, help="Output image path")
    pd.add_argument("--emin", type=float, default=-0.5)
    pd.add_argument("--emax", type=float, default=0.0)
    pd.add_argument("--xmin", type=float, default=-1000.0)
    pd.add_argument("--xmax", type=float, default=1000.0)
    pd.add_argument("--keep-last", action="store_true", help="Keep last DOS column (do not drop)")
    pd.add_argument("--alpha", type=Path, default=None, help="Explicit ALPHA .DOSS file")
    pd.add_argument("--beta", type=Path, default=None, help="Explicit BETA .DOSS file")

    # --- cphf ----------------------------------------------------------------
    pc = sub.add_parser("cphf", help="Plot absorption overlays from discrete CPHF .out files")
    pc.add_argument("--root", type=Path, default=None, help="Data root (defaults to $GA2O3_DATA or cwd)")
    pc.add_argument("--outdir", required=True, type=Path, help="Directory to write absorption_{iso,xx,yy,zz}.png")
    pc.add_argument("--cases", nargs="+", default=["Pristine", "TETRA", "OCTA"], help="Cases to plot")
    pc.add_argument("--wl-min", type=float, default=465.0, help="Minimum wavelength (nm) filter")
    pc.add_argument("--yscale", choices=["linear", "log"], default="log")
    pc.add_argument("--correction", choices=["none", "scissor", "ratio"], default="ratio")
    pc.add_argument("--eg-pbe", type=float, default=2.62)
    pc.add_argument("--eg-target", type=float, default=4.43)
    pc.add_argument("--ratio-factor", type=float, default=None, help="Override Eg_target/Eg_pbe")
    pc.add_argument("--scissor-shift", type=float, default=None, help="Override Eg_target-Eg_pbe (eV)")
    pc.add_argument("--extra-shift", type=float, default=0.0, help="Extra energy shift for ratio mode (eV)")
    pc.add_argument("--dpi", type=int, default=600)
    pc.add_argument("--figsize", type=float, default=6.0)
    pc.add_argument("--no-legend", action="store_true")

    # --- density -------------------------------------------------------------
    pe = sub.add_parser("density", help="Plot charge density from a CASTEP .f25 file")
    pe.add_argument("file", type=Path, help=".f25 file")
    pe.add_argument("--out", required=True, type=Path, help="Output image path")
    pe.add_argument("--pmin", type=float, default=0.01)
    pe.add_argument("--pmax", type=float, default=99.2)
    pe.add_argument("--dpi", type=int, default=900)

    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    # Lazy imports so the CLI starts quickly.
    from ga2o3.dataset import Ga2O3Dataset

    root = args.root if hasattr(args, "root") else None
    if root is None:
        root = Path(os.getenv("GA2O3_DATA", ".")).expanduser().resolve()
    ds = Ga2O3Dataset(root=root)

    if args.cmd == "bands":
        from ga2o3.plot.bands import plot_bands_only

        alpha, beta = ds.band_files(args.case)
        fig, _ = plot_bands_only(
            alpha,
            beta,
            e_min=args.emin,
            e_max=args.emax,
            energy_shift=args.shift,
            band_legend=args.legend,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, dpi=600, bbox_inches="tight")
        return 0

    if args.cmd == "dos":
        from ga2o3.plot.dos import plot_spin_resolved_dos

        alpha = args.alpha
        beta = args.beta
        if alpha is None or beta is None:
            alpha2, beta2 = ds.dos_spin_files(args.case)
            alpha = alpha or alpha2
            beta = beta or beta2

        fig, _ = plot_spin_resolved_dos(
            alpha,
            beta,
            e_range=(args.emin, args.emax),
            dos_range=(args.xmin, args.xmax),
            drop_last_projection=not args.keep_last,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, dpi=600, bbox_inches="tight")
        return 0

    if args.cmd == "cphf":
        from ga2o3.plot.cphf import plot_absorption_overlays, save_absorption_figs

        figs = plot_absorption_overlays(
            ds,
            cases=args.cases,
            wl_min=args.wl_min,
            y_scale=args.yscale,
            correction=args.correction,
            eg_pbe=args.eg_pbe,
            eg_target=args.eg_target,
            ratio_factor=args.ratio_factor,
            scissor_shift_ev=args.scissor_shift,
            extra_shift_ev=args.extra_shift,
            dpi=args.dpi,
            figsize=args.figsize,
            add_legend=not args.no_legend,
        )
        save_absorption_figs(figs, args.outdir, fmt="png", dpi=args.dpi)
        return 0

    if args.cmd == "density":
        from ga2o3.plot.density import plot_charge_density_f25

        fig, _ = plot_charge_density_f25(
            args.file,
            out_path=args.out,
            charge_percentiles=(args.pmin, args.pmax),
            dpi=args.dpi,
        )
        return 0

    raise SystemExit(f"Unknown command: {args.cmd}")
