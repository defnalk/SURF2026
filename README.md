# SURF 2026  Enthalpic vs. Entropic Drivers of Ion Adsorption at Charged Surfaces

**Caltech Summer Undergraduate Research Fellowship (SURF) 2026**
**Host Lab:** Fong Lab, Caltech
**Student:** Defne Nihal Ertugrul

---

## 1. Scientific Objective

This project uses classical molecular dynamics (MD) to disentangle the
**enthalpic (ΔH)** and **entropic (ΔS)** contributions to the adsorption free
energy (ΔG = ΔH − TΔS) of monovalent (Na⁺) and divalent (Ca²⁺) cations at a
charged planar electrode in explicit water.

Three primary systems are studied:

| System       | Electrolyte       | Purpose                                        |
|--------------|-------------------|------------------------------------------------|
| Monovalent   | NaCl only         | Baseline for Na⁺ adsorption thermodynamics     |
| Divalent     | CaCl₂ only        | Baseline for Ca²⁺ adsorption thermodynamics    |
| Mixed        | NaCl + CaCl₂      | Competitive adsorption between Na⁺ and Ca²⁺    |

Each system is simulated at a **baseline temperature (T₀ ≈ 298 K)** and an
**elevated temperature (T₀ + 40 K)** so that the temperature dependence of the
adsorption free energy can be used to decompose ΔG into ΔH and ΔS via a
finite difference van 't Hoff analysis.

---

## 2. Repository Organization

```
SURF2026/
├── README.md
├── .gitignore
├── environment.yml              # Conda environment for Python analysis
│
├── setup/                       # System construction (Weeks 1–2)
│   ├── packmol/                 # Packmol input files (.inp)
│   ├── moltemplate/             # Moltemplate .lt topology files
│   └── forcefields/             # SPC/E water, Joung–Cheatham ions, etc.
│
├── simulations/                 # LAMMPS input decks (Weeks 3–6)
│   ├── monovalent/              # NaCl-only systems
│   ├── divalent/                # CaCl2-only systems
│   └── mixed/                   # NaCl + CaCl2 competitive systems
│
├── hpc_scripts/                 # Purdue Anvil Slurm job scripts
│   └── submit_anvil.slurm
│
├── analysis/                    # Post-processing (Weeks 6–9)
│   ├── scripts/
│   │   └── analyze_hydration.py # RDF + coordination-number skeleton
│   ├── notebooks/               # Jupyter notebooks for exploratory analysis
│   └── results/                 # Plots, CSVs (small files only)
│
├── visualization/               # Weeks 7–9
│   ├── vmd/                     # .tcl scripts: PBC wrap, representations
│   └── ovito/                   # OVITO Python modifiers
│
├── data_management/             # Automated Anvil ↔ macOS backups
│   ├── rclone_backup.sh
│   └── com.defne.surf.backup.plist   # macOS launchd agent
│
└── docs/                        # Weekly reports, figures, final poster
```

The directory layout mirrors the **10 week SURF workflow**:

| Weeks | Stage                          | Directories                              |
|-------|--------------------------------|------------------------------------------|
| 1/2   | Literature + system building   | `setup/`, `docs/`                        |
| 3/4   | Equilibration + validation     | `simulations/`, `hpc_scripts/`           |
| 5/6   | Production runs (2 T, 3 sys)   | `simulations/`, `hpc_scripts/`           |
| 6/8   | Structural analysis            | `analysis/`, `visualization/`            |
| 8/9   | Thermodynamic decomposition    | `analysis/notebooks/`                    |
| 10    | Poster + final report          | `docs/`                                  |

---

## 3. Installation & Environments

### 3.1 Python analysis environment

```bash
conda env create -f environment.yml
conda activate surf2026
```

### 3.2 LAMMPS on Purdue Anvil

```bash
module load gcc/11 openmpi/4 lammps/2Aug2023
```

See `hpc_scripts/submit_anvil.slurm` for a ready to edit Slurm template.

### 3.3 System building tools (local)

- [Packmol](https://m3g.github.io/packmol/)
- [Moltemplate](https://www.moltemplate.org/)
- [VMD](https://www.ks.uiuc.edu/Research/vmd/)
- [OVITO](https://www.ovito.org/)

---

## 4. Simulation Workflow

1. **Build** the simulation box with Packmol (water + electrodes + ions)
   and assign force field parameters with Moltemplate.
2. **Minimize** and **equilibrate** in LAMMPS: NVT → NPT → NVT with the
   charged walls held rigid.
3. **Production** NVT runs at T₀ and T₀ + 40 K for each of the three
   electrolyte compositions.
4. **Sync** trajectories from Anvil to local storage nightly via
   `data_management/rclone_backup.sh` (scheduled by launchd).
5. **Analyze** with `analysis/scripts/analyze_hydration.py` using
   MDAnalysis: ion density profiles, RDFs, hydration shell coordination
   numbers, and ion clustering statistics.
6. **Decompose** ΔG(T) → ΔH, ΔS via two temperature finite difference.
7. **Visualize** in VMD (custom `.tcl` PBC wrap scripts) and OVITO.

---

## 5. Reproducibility Notes

- LAMMPS input decks are version controlled; trajectory outputs are not
  (see `.gitignore`). Raw trajectories live on Anvil scratch and are
  mirrored via rclone.
- Random seeds for `velocity create` and `fix langevin` are recorded in
  each input deck.
- Module and compiler versions are logged at the top of every Slurm job
  via `module list`.

---

## 6. Acknowledgments

Supported by the Caltech SURF program and the Fong Lab. Computations
performed on **Purdue Anvil** via ACCESS.
