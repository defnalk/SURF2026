#!/usr/bin/env python
"""
analyze_hydration.py
====================

Skeleton analysis script for the SURF 2026 ion-adsorption project.

Given a LAMMPS data file (topology) and a DCD trajectory (coordinates),
compute:
    1. Radial distribution functions (RDFs) for
         - Na+ .. O(water)
         - Ca2+ .. O(water)
         - Na+ .. Cl-
         - Ca2+ .. Cl-
    2. First-shell coordination numbers by integrating g(r) r^2 dr
       out to the first minimum of each RDF.

The script is intentionally modular: each quantity is produced by a
small function so that additional analyses (density profiles, residence
times, cluster statistics, ...) can be added as the summer progresses.

Usage
-----
    python analyze_hydration.py \
        --data   ../../simulations/monovalent/system.data \
        --traj   ../../simulations/monovalent/prod_298K.dcd \
        --outdir ../results/monovalent_298K

Author: Defne Nihal Ertugrul -- Fong Lab, Caltech SURF 2026
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rdf as mda_rdf

# np.trapz was deprecated in NumPy 2.0 in favour of np.trapezoid; fall
# back to the old name for environments still on NumPy 1.x.
_trapezoid = getattr(np, "trapezoid", np.trapz)


# ---------------------------------------------------------------------
# Configuration (edit atom-type names to match your Moltemplate output)
# ---------------------------------------------------------------------
#
# NOTE: LAMMPS data files label atoms by numeric type, not element.
# MDAnalysis will expose them as Universe.atoms.types == '1', '2', ...
# Update the mapping below once you know your type IDs.
ATOM_TYPES = {
    "OW":  "1",   # SPC/E water oxygen
    "HW":  "2",   # SPC/E water hydrogen
    "Na":  "3",
    "Ca":  "4",
    "Cl":  "5",
}

RDF_PAIRS = [
    ("Na", "OW"),
    ("Ca", "OW"),
    ("Na", "Cl"),
    ("Ca", "Cl"),
]

RDF_RANGE = (0.0, 10.0)   # Angstrom
RDF_NBINS = 200


# ---------------------------------------------------------------------
# Core routines
# ---------------------------------------------------------------------
def load_universe(data_file: Path, traj_file: Path) -> mda.Universe:
    """Load a LAMMPS data file + DCD trajectory into an MDAnalysis Universe."""
    logging.info("Loading topology: %s", data_file)
    logging.info("Loading trajectory: %s", traj_file)
    u = mda.Universe(
        str(data_file),
        str(traj_file),
        topology_format="DATA",
        format="DCD",
    )
    logging.info(
        "Loaded %d atoms, %d frames (dt = %.3f ps)",
        u.atoms.n_atoms, u.trajectory.n_frames, u.trajectory.dt,
    )
    return u


def select(u: mda.Universe, species: str) -> mda.AtomGroup:
    """Return the AtomGroup for a species name in ATOM_TYPES."""
    t = ATOM_TYPES[species]
    ag = u.select_atoms(f"type {t}")
    if ag.n_atoms == 0:
        raise ValueError(f"No atoms of species '{species}' (type {t}) found.")
    return ag


def compute_rdf(
    u: mda.Universe,
    species_a: str,
    species_b: str,
    rmin: float = RDF_RANGE[0],
    rmax: float = RDF_RANGE[1],
    nbins: int = RDF_NBINS,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (r, g(r)) for the A-B pair."""
    a = select(u, species_a)
    b = select(u, species_b)
    logging.info("RDF %s-%s: %d x %d atoms", species_a, species_b, a.n_atoms, b.n_atoms)

    r = mda_rdf.InterRDF(a, b, nbins=nbins, range=(rmin, rmax), exclusion_block=None)
    r.run(start=start, stop=stop, step=step, verbose=True)
    return r.results.bins, r.results.rdf


def coordination_number(
    r: np.ndarray,
    g: np.ndarray,
    rho_b: float,
    r_cut: float | None = None,
) -> tuple[float, float]:
    """
    Integrate  n(r) = 4*pi*rho_b * int_0^{r_cut} g(r') r'^2 dr'.

    If ``r_cut`` is None, use the first minimum of g(r) after its first peak.

    Returns (r_cut, n_coord).
    """
    if r_cut is None:
        # Find first maximum, then the first minimum beyond it.
        # Simple-minded: use diff-sign changes. Good enough for a skeleton.
        dg = np.diff(g)
        maxima = np.where((dg[:-1] > 0) & (dg[1:] <= 0))[0] + 1
        if maxima.size == 0:
            raise RuntimeError("No peak found in g(r); cannot auto-select r_cut.")
        first_max = maxima[0]
        minima = np.where((dg[first_max:-1] < 0) & (dg[first_max + 1:] >= 0))[0]
        if minima.size == 0:
            raise RuntimeError("No minimum found after first peak of g(r).")
        r_cut = r[first_max + minima[0] + 1]

    mask = r <= r_cut
    integrand = g[mask] * r[mask] ** 2
    n_coord = 4.0 * np.pi * rho_b * _trapezoid(integrand, r[mask])
    return float(r_cut), float(n_coord)


def number_density(u: mda.Universe, species: str) -> float:
    """Bulk number density of a species, averaged over the trajectory box."""
    ag = select(u, species)
    volumes = np.array([ts.volume for ts in u.trajectory])
    return ag.n_atoms / volumes.mean()


# ---------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data",   required=True, type=Path, help="LAMMPS .data file")
    p.add_argument("--traj",   required=True, type=Path, help="LAMMPS .dcd trajectory")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory")
    p.add_argument("--start",  type=int, default=None, help="First frame (inclusive)")
    p.add_argument("--stop",   type=int, default=None, help="Last frame (exclusive)")
    p.add_argument("--step",   type=int, default=1,    help="Frame stride")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
    )

    args.outdir.mkdir(parents=True, exist_ok=True)

    u = load_universe(args.data, args.traj)

    summary_rows = []
    for a, b in RDF_PAIRS:
        try:
            r, g = compute_rdf(u, a, b, start=args.start, stop=args.stop, step=args.step)
        except ValueError as exc:
            logging.warning("Skipping RDF %s-%s: %s", a, b, exc)
            continue

        out_csv = args.outdir / f"rdf_{a}_{b}.csv"
        np.savetxt(out_csv, np.column_stack([r, g]), header="r[A] g(r)", comments="")
        logging.info("Wrote %s", out_csv)

        rho_b = number_density(u, b)
        try:
            r_cut, n_coord = coordination_number(r, g, rho_b)
        except RuntimeError as exc:
            logging.warning("Coordination number for %s-%s failed: %s", a, b, exc)
            r_cut, n_coord = float("nan"), float("nan")

        logging.info("  %s-%s : r_cut = %.3f A, N_coord = %.3f", a, b, r_cut, n_coord)
        summary_rows.append((a, b, r_cut, n_coord))

    summary_csv = args.outdir / "coordination_summary.csv"
    with summary_csv.open("w") as f:
        f.write("species_a,species_b,r_cut_A,n_coord\n")
        for a, b, rc, nc in summary_rows:
            f.write(f"{a},{b},{rc:.4f},{nc:.4f}\n")
    logging.info("Wrote %s", summary_csv)


if __name__ == "__main__":
    main()
