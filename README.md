# GlycoDYNAMICS

**A JSON-driven workflow for preparing, running, and analyzing glycan REST-RECT molecular dynamics simulations.**

GlycoDYNAMICS connects glycan metadata and CHARMM-GUI simulation systems to automated GROMACS minimization/equilibration, REST-RECT (Replica Exchange with Solute Tempering and Replica Exchange with Collective-variable Tempering) setup, HPC production, and postprocessing. It can also export selected results to a separate [GlycoLIBRARY](https://github.com/IsabellGrothaus/GlycoLIBRARY) database.

The main entry point is [`Workflow.ipynb`](Workflow.ipynb); its implementation is in [`scripts/`](scripts). Each glycan has its own `<glycan_id>.json`, and the workflow discovers glycan folders recursively under `data/`. There is **no central `config.yaml` or `glycans.yaml`**.

## Workflow at a glance

```text
Glycan structure (IUPAC and/or CHARMM-GUI PDB)
              |
              v
  data/.../<glycan_id>/
  + CHARMM-GUI GROMACS files
              |
              v
  create_glycan_json()          glycan metadata and graph
              |
              v
  init_glycan_simulation()      simulation parameters stored in JSON
              |
              v
  generate_glycan_simulation()
       |                 |
       v                 v
  01_min_equ/           02_REST-RECT/
  min/equilibrate       PLUMED, REST scaling, TOPO0...N, job.sh
                              |
                              v
                       Run job.sh on HPC
                              |
                              v
                       Download outputs
                              |
                              v
                       analyze_glycan()
                              |
                              v
                          03_lib/
                  COLVAR, torsions, FES, trajectory
                              |
                              v
                  copy_glycan_library() (optional)
```

## Requirements and dependencies

### Python

Use **Python 3.10+** (the code uses modern type annotations and `str.removeprefix`). The scripts use the Python standard library as well as:

| Package | Used for |
|---|---|
| `numpy` | Numerical processing and free-energy profiles |
| `pandas` | Processing PLUMED COLVAR files |
| `MDAnalysis` | Trajectory processing, alignment and PBC handling |
| `matplotlib` | Glycan graph plotting |
| `networkx` | Glycan graph construction and plotting |
| `jupyterlab` / `ipykernel` | Running the example notebook |

For example:

```bash
git clone https://github.com/IsabellGrothaus/GlycoDYNAMICS.git
cd GlycoDYNAMICS
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy pandas MDAnalysis matplotlib networkx jupyterlab ipykernel
jupyter lab Workflow.ipynb
```

These Python packages do **not** install GROMACS or PLUMED. A working molecular-simulation environment is required separately.

### Molecular simulation software

- **CHARMM-GUI**: prepare the solvated glycan/GROMACS system, including its topology, coordinates, and index file. The workflow expects the extracted `gromacs/` folder in a directory whose name contains `charmm-gui`.
- **GROMACS**: `gmx` and `gmx_mpi` commands are used for energy minimization, equilibration, topology preparation, trajectory conversion, and production.
- **PLUMED with REST-RECT support**: the generated PLUMED input uses the `METAD ... RECT=` option. Make sure the chosen PLUMED build supports this functionality and is compatible with GROMACS; the analysis also calls `plumed sum_hills`.
- **Bash and an HPC scheduler (SLURM for the provided template)**: `templates/generate_HREMD.sh` prepares replica input, and `templates/job.sh` submits/runs the multi-replica production simulation.
- **Optional Perl**: needed if using `helper/demux.pl` for replica-exchange demultiplexing.

The included `templates/job.sh` is configured for a particular A100-based cluster and currently references **GROMACS 2022.5 / PLUMED 2.8.3** and **12 MPI ranks/replicas**. These are example environment settings, not portable defaults for all systems. Inspect and adapt the module name, GPU request, wall time, number of MPI ranks, and resource allocation before submission.

**Important:** In the current `scripts/setup.py`, `generate_glycan_simulation()` and `analyze_glycan()` default to loading modules named `chem/GROMACS` and `chem/PLUMED` through `bash -lc`. Adjust the modules in the code or supply the correct site-specific module names. Passing `None` does not disable module loading in the current implementation.

## Repository structure

```text
GlycoDYNAMICS/
├── Workflow.ipynb             # notebook interface for the full workflow
├── scripts/
│   ├── generate_metadata.py   # IUPAC/PDB -> per-glycan metadata JSON and graph
│   ├── setup.py               # min/equilibration, REST-RECT setup and analysis
│   └── generate_library.py    # export results and update library index
├── templates/                # .mdp files, generate_HREMD.sh, job.sh, ...
├── helper/                   # optional helper utilities (e.g. demux.pl)
├── test/notebooks/           # development/test notebooks
└── data/                     # local glycan data root (create if absent)
    └── N-glycans/
        └── High-mannose/
            └── Man5/
                ├── Man5.json
                ├── Man5.png
                ├── <downloaded-charmm-gui-folder>/
                │   └── gromacs/
                ├── 01_min_equ/
                ├── 02_REST-RECT/
                ├── 03_lib/
                └── rest_rect_setup.log
```

The `data/` tree can also include `N-glycans/Complex`, `N-glycans/Hybrid`, `O-glycans`, polymers, or other organizational folders. Category paths are flexible; the final glycan directory name must match the JSON filename and its `glycan_id` value. Glycan IDs must be **unique across the entire data root**.

## Adding a new glycan to `data/`

### 1. Create a directory and obtain the CHARMM-GUI files

For a glycan named `Man5`:

```bash
mkdir -p data/N-glycans/High-mannose/Man5
```

Prepare the molecular system in [CHARMM-GUI](https://www.charmm-gui.org/) and extract its GROMACS output into a subdirectory of that glycan folder:

```text
data/N-glycans/High-mannose/Man5/
└── my_charmm-gui_download/
    └── gromacs/
        ├── step3_input.gro
        ├── topol.top
        ├── index.ndx
        └── toppar/...
```

The code searches for exactly one `*charmm-gui*/gromacs/` directory per glycan. Ensure the CHARMM-GUI package contains all necessary topology, restraint and coordinate files. Verify residue names and atom numbering, particularly for modified sugars and glycopeptides.

### 2. Generate the metadata JSON

From the repository root, in `Workflow.ipynb` (or another Python session):

```python
from pathlib import Path
import sys
sys.path.insert(0, str(Path("scripts").resolve()))

from generate_metadata import create_glycan_json, plot_graph
from setup import (
    init_glycan_simulation,
    generate_glycan_simulation,
    status,
    analyze_glycan,
)

ROOT = Path("data")
glycan_id = "Man5"

# Infer the glycan information from the CHARMM-GUI PDB when possible.
create_glycan_json(ROOT, glycan_id)
plot_graph(ROOT, glycan_id)
```

Alternatively, supply an explicit IUPAC-condensed description and/or `pdb_path`:

```python
create_glycan_json(
    ROOT,
    glycan_id,
    iupac="<your valid IUPAC-condensed glycan string>",
    pdb_path="<path to CHARMM-GUI structure.pdb>",
)
```

The generated JSON contains the glycan identity and structural metadata. **Check the resulting JSON and glycan graph manually** before running a simulation, especially the `monosaccharides`, `pdb_residues`, `pucker_ring_atoms`, `graph_edges`, `forcefield_residue_names` and glycosidic linkage assignments. The metadata generator does not guarantee that candidate force-field residue names match every CHARMM-GUI topology. Existing JSON files are skipped unless `overwrite=True`; overwriting discards changes stored in the file, so use this option carefully.

### 3. Add simulation settings to the same JSON

```python
init_glycan_simulation(
    folder_path=ROOT,
    glycan_id=glycan_id,
    n_replica=12,
    min_temperature=310.15,
    max_temperature=800.0,
    structure_file="step5_production.gro",
    mdp_file="production.mdp",
    template_files=["*.mdp", "generate_HREMD.sh", "job.sh"],
    forcefield_file="toppar/forcefield.itp",
)
```

This **updates the `simulation` object of `Man5.json`** without deleting the structural metadata. It does not generate a separate YAML file.

Relevant `simulation` keys:

| Key | Purpose |
|---|---|
| `n_replica` | Number of Hamiltonian replicas (default 12) |
| `min_temperature`, `max_temperature` | Temperature range used for replica scaling, in K |
| `structure_file` | Equilibrated structure for REST-RECT setup |
| `mdp_file` | Production `.mdp` used to generate replica inputs |
| `template_files` | Files/patterns copied from `templates/` |
| `forcefield_file` | Force-field file relative to `01_min_equ/` |
| `rect_factors` (optional) | Per-replica RECT scaling factors; must match `n_replica` |
| `metad` (optional) | Metadynamics settings (`pace`, `tau`, `sigma`) |

The default RECT factor list contains **12 entries**; if changing `n_replica`, specify an appropriately sized `rect_factors` list in the glycan JSON and adapt `job.sh` accordingly.

### 4. Generate minimization, equilibration and REST-RECT inputs

```python
status(glycan_id, ROOT)

rest_dir = generate_glycan_simulation(
    data_dir=ROOT,
    glycan_id=glycan_id,
    template_dir="templates",
    # Supply installed command names and site-specific modules as required:
    gmx="gmx",
    gmx_mpi="gmx_mpi",
    gromacs_module="chem/GROMACS",
)
print(rest_dir)
status(glycan_id, ROOT)
```

This function performs the following steps:

1. Finds the uniquely named glycan directory and validates its JSON metadata.
2. Copies CHARMM-GUI GROMACS inputs into `01_min_equ/` without overwriting existing files.
3. Runs GROMACS minimization (`step4.0_minimization`), equilibration (`step4.1_equilibration`) and preparation (`step5_production`), skipping completed steps.
4. Creates `reference.pdb` and generates `plumed_allangles.dat` from the glycan's connectivity and ring atoms. It defines glycosidic `phi`/`psi`, relevant `omega` torsions, ring puckering and applicable glycopeptide torsions.
5. Processes the GROMACS topology, identifies glycan solute atoms using JSON-derived CHARMM residue names, and runs `generate_HREMD.sh` to construct replica directories.
6. Organizes the generated input files under `02_REST-RECT/`, including `TOPO0`, `TOPO1`, ..., `TOPO11` for a 12-replica setup.
7. Appends messages and command output to `rest_rect_setup.log`.

**Sialic acid caution:** the setup detects sialic acid residues from JSON metadata and has a specific, validated Neu5Ac force-field correction. Other sialic-acid chemistries are not automatically supported by that correction and may raise an error. The code also has special torsion definitions for Neu5Ac `alpha2-8` linkages; verify other linkages before production.

**Resuming:** repeated setup calls skip completed minimization/equilibration and completed REST setup. If partially generated `TOPO*` folders are present, the code intentionally stops rather than deleting or overwriting them. Inspect the folders and log before continuing.

### 5. Run REST-RECT production on an HPC cluster

Once `generate_glycan_simulation()` completes, copy the prepared `02_REST-RECT/` directory to the cluster, preserving its contents. For example:

```bash
# Run on your workstation; replace host and path with your actual destination.
rsync -av data/N-glycans/High-mannose/Man5/02_REST-RECT/ \
    user@cluster:/path/to/Man5/02_REST-RECT/
```

On the cluster, inspect the copied `job.sh` and adjust the resource request, GROMACS/PLUMED module, replica count and GPU/CPU settings. Then submit:

```bash
cd /path/to/Man5/02_REST-RECT
sbatch job.sh
```

The supplied script runs multi-directory GROMACS/PLUMED replica exchange with `-multidir`, `-replex 400` and `-hrex`. It is **not executed by the Python setup function**. Check the cluster's module system and ensure the production PLUMED library supports the generated `RECT` syntax. Use the checkpoint restart command in `job.sh` if continuing an interrupted production run, after confirming matching checkpoints and replica directories.

### 6. Download the completed production outputs

Copy the `02_REST-RECT/` simulation results back to the same local glycan directory, preserving `TOPO*` subdirectories:

```bash
rsync -av user@cluster:/path/to/Man5/02_REST-RECT/ \
    data/N-glycans/High-mannose/Man5/02_REST-RECT/
```

At minimum, `analyze_glycan()` checks for the following nonempty files in `02_REST-RECT/TOPO0/`:

```text
md.log
COLVAR.0
traj_comp.xtc
```

Further postprocessing requires the corresponding reference structure, topology, `HILLS_*` files and any other inputs used by GROMACS or PLUMED. Check that the full run completed and that the other replicas have their expected outputs; the current readiness check primarily inspects `TOPO0`.

### 7. Analyze the simulation

```python
from setup import analyze_glycan, status

analysis_dir = analyze_glycan(
    data_dir=ROOT,
    glycan_id=glycan_id,
    gmx="gmx",
    gromacs_module="chem/GROMACS",
    plumed_module="chem/PLUMED",
)
print(analysis_dir)
status(glycan_id, ROOT)
```

Postprocessing collects exchange statistics, optionally runs `helper/demux.pl`, converts and aligns the trajectory with GROMACS/MDAnalysis, extracts COLVAR torsion data, and reconstructs 1D free-energy profiles using `plumed sum_hills`. Outputs are stored in `03_lib/`, for example:

```text
03_lib/
├── angles.dat
├── omega_angles.dat
├── colvar.dat
├── pucker.dat                         # when the theta COLVAR is available
├── fes_phi1_2_0.dat                  # example FES filename
├── FES_phi1_2_0_1D/                 # intermediate FES profiles
├── conformations_rest-rect_aligned.xtc
└── reference_rest-rect.pdb
```

The function skips selected existing processed files, and the log records completed or partially completed analysis. Note that `plumed sum_hills` and the FES processing settings are defined in `scripts/setup.py`; review them before applying the workflow to systems requiring different collective-variable ranges or binning.

## Exporting results to GlycoLIBRARY (optional)

After analysis, copy glycan metadata, plots, and selected analysis files into an external library. The destination should **not** be inside `data/`:

```python
from generate_library import (
    copy_glycan_library,
    generate_glycan_json,
    display_glycan_library,
)

DATABASE = Path("/path/to/GlycoLIBRARY")
copy_glycan_library(library_path=DATABASE, data_path=ROOT)
generate_glycan_json(DATABASE)  # rebuild index (also done during copy)
display_glycan_library(DATABASE)
```

The export retains the category hierarchy under `DATABASE/library/`, copies eligible files from `03_lib/` and builds/updates `DATABASE/glycan.json`. It can register glycan metadata before analysis is complete, with explicit missing-file warnings. Simulation source directories are not moved or deleted.

## Status and troubleshooting

Run:

```python
status("Man5", Path("data"))
```

The result reports five stages: `charmm_gui`, `min_equ`, `rest_rect_setup`, `rest_rect_production`, and `analysis`.

| Symptom | What to check |
|---|---|
| Glycan not found or duplicate ID | Exactly one directory named `<glycan_id>` exists below `data/` |
| Missing JSON / validation error | Filename and `glycan_id` match; required fields are present; graph node IDs match residue/pucker dictionaries |
| CHARMM-GUI folder not found | An extracted directory matching `*charmm-gui*/gromacs/` is present under the glycan folder |
| Missing atom or PLUMED atom mapping | PDB residue numbering, atom names, graph edges and ring definitions agree with `reference.pdb` |
| Module / `gmx` / `gmx_mpi` not found | Change the site-specific module names and executable paths; verify GROMACS is installed |
| `generate_HREMD.sh` or `job.sh` missing | Confirm required files are present under `templates/` |
| Existing partial `TOPO*` setup | Inspect incomplete files and `rest_rect_setup.log`; setup intentionally refuses destructive replacement |
| Unsupported sialic-acid correction | The current automated force-field patch is validated for Neu5Ac only |
| `analyze_glycan()` reports incomplete production | Download `TOPO0/md.log`, `COLVAR.0`, `traj_comp.xtc` and remaining analysis inputs |
| Missing FES output | Check `TOPO0/HILLS_*`, PLUMED availability and `03_lib/FES_*/` |

## Implementation notes

- **Per-glycan metadata:** structural and simulation information live together in `<glycan_id>.json`. There is no separate configuration YAML.
- **Glycan naming:** glycan folder names must be unique below `data/`, irrespective of category.
- **Metadata verification:** automatic reconstruction and candidate residue names require scientific validation against CHARMM-GUI topology and the intended glycan structure.
- **Reproducibility:** retain the original CHARMM-GUI files, validated JSON, template files, program/module versions, `job.sh`, and `rest_rect_setup.log` for each system.
- **HPC boundary:** input generation and analysis are local operations, while production is explicitly submitted by the user on the HPC system.

## Related resources

- [GROMACS documentation](https://manual.gromacs.org/)
- [PLUMED documentation](https://www.plumed.org/doc)
- [CHARMM-GUI](https://www.charmm-gui.org/)
- [MDAnalysis](https://www.mdanalysis.org/)
- [GlycoLIBRARY](https://github.com/IsabellGrothaus/GlycoLIBRARY)

---

**Note:** This README describes the repository's current JSON-based `scripts/setup.py` workflow rather than the outdated YAML-based description in the earlier README. Check cluster-specific settings and scientific parameter choices before running new production simulations.
