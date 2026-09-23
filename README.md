# REST-RECT setup workflow

This is a streamlined version of the original `M_Automated_REST-RECT_setup.ipynb`.

## Structure

```text
REST-RECT/
├── config/
│   └── glycans.yaml
├── notebooks/
│   └── REST-RECT_setup.ipynb
├── scripts/
│   └── rest_rect_setup.py
└── templates/
    ├── *.mdp
    ├── generate_HREMD.sh
    └── job.sh
```

The simulation directories (`M/Man5`, `M/Man6`, ...) can remain outside this project directory if the paths in `glycans.yaml` point to them.

## Important compatibility point

The original notebook imported project-specific functions from `setup_RESTRECT.py`:

- `write_plumed`
- `define_solute`
- `check_and_delete_files_and_directories`
- `remove_lines_and_insert`

The new workflow deliberately does **not** reimplement these scientific helper functions because their exact PLUMED/REST2 behavior was not contained in the uploaded notebook. Keep your existing `setup_RESTRECT.py` available on `PYTHONPATH` (for example, next to `scripts/rest_rect_setup.py`, or in the working directory).

The new code handles configuration, validation, command execution, logging, metadata and file organisation around those existing functions.

## Templates

Copy your existing `.mdp` files, `generate_HREMD.sh`, and `job.sh` into `templates/`.

The notebook will copy the templates into the selected simulation directory before executing the workflow.

## First test

Use one glycan, e.g. `Man9`, and run the preview cell first. Then execute. Once that has been verified, change `selected_glycans` to a list of systems.

## Cluster

The original notebook used:

```bash
ml chem/GROMACS
gmx ...
gmx_mpi ...
```

Set `USE_MODULE_LOAD = True` in the notebook for that environment. On a machine where GROMACS is already available on `PATH`, use `False`.


## Update: additional glycans from C_Automated_REST-RECT_setup.ipynb

The configuration now also contains the 23 glycans found in the C notebook:
A1F, A2, A2_26NE5AC, A2F, A2F2, A2G1, FA2G1, Fuc1_Neu4,
Gal1_Fuc1, Gal2_Fuc1_Neu1, Gal2_Fuc2_Neu1, Gal3_Fuc1,
Gal3_Fuc1_Neu1, Gal4_Fuc1, Gal6_Fuc1_Neu3, Man3_Fuc1,
Man3_Gal1, Man3_Gal1_Fuc1, NA2, NA2F, NGA2, NGA2F, NGA3F.

A new `add_glycan_to_yaml()` helper is available in `scripts/rest_rect_setup.py`.
It validates the definition before atomically writing it to `config/glycans.yaml`.
