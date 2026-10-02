"""Reusable setup engine for the REST-RECT glycan workflow.

The scientific operations intentionally follow the original notebooks:
1. copy the required template files
2. run minimization, equilibration and production preprocessing/runs
3. create reference.pdb
4. create plumed_allangles.dat using the existing setup_RESTRECT implementation
5. create processed.top
6. define the REST2 solute
7. run generate_HREMD.sh
8. organise the resulting files

Some functions are working for the atom names of the CHARMM36m FF. Might need to be adapted for the AMBER FF family
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml


class SetupError(RuntimeError):
    """Raised when a REST-RECT setup cannot be completed safely."""


def load_config(config_file: str | Path) -> dict[str, Any]:
    """Load the YAML configuration."""
    config_file = Path(config_file).expanduser().resolve()
    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    with config_file.open() as fh:
        config = yaml.safe_load(fh)
    if not isinstance(config, dict) or "defaults" not in config or "glycans" not in config:
        raise SetupError("YAML must contain 'defaults' and 'glycans'.")
    return config


def merge_config(config: Mapping[str, Any], glycan_name: str) -> dict[str, Any]:
    """Merge global defaults with glycan-specific settings."""
    glycans = config.get("glycans", {})
    if glycan_name not in glycans:
        raise KeyError(f"Unknown glycan '{glycan_name}'. Available: {list(glycans)}")
    params = deepcopy(dict(config["defaults"]))
    params.update(deepcopy(dict(glycans[glycan_name])))
    return params


def validate_glycan_config(name: str, params: Mapping[str, Any]) -> list[str]:
    """Return validation errors without touching the filesystem."""
    errors: list[str] = []
    required = ["path", "n_replica", "min_temperature", "max_temperature",
                "structure_file", "mdp_file", "peptide", "label", "ltype",
                "carbon_t", "carbon_p", "monosaccharide_types"]
    for key in required:
        if key not in params:
            errors.append(f"{name}: missing '{key}'")

    if errors:
        return errors

    n_links = len(params["label"])
    for key in ("ltype", "carbon_t"):
        if len(params[key]) != n_links:
            errors.append(f"{name}: len({key})={len(params[key])} but len(label)={n_links}")

    if len(params["carbon_p"]) < max(max(x) for x in params["label"]):
        errors.append(
            f"{name}: carbon_p has {len(params['carbon_p'])} entries but the highest residue label is "
            f"{max(max(x) for x in params['label'])}"
        )

    if int(params["n_replica"]) < 2:
        errors.append(f"{name}: n_replica must be >= 2")
    if float(params["min_temperature"]) <= 0 or float(params["max_temperature"]) <= 0:
        errors.append(f"{name}: temperatures must be positive")
    if float(params["max_temperature"]) < float(params["min_temperature"]):
        errors.append(f"{name}: max_temperature must be >= min_temperature")

    return errors


def validate_config(config: Mapping[str, Any]) -> dict[str, list[str]]:
    """Validate every configured glycan and return only systems with errors."""
    problems = {}
    for name in config["glycans"]:
        errors = validate_glycan_config(name, merge_config(config, name))
        if errors:
            problems[name] = errors
    return problems



def add_glycan_to_yaml(
    config_file: str | Path,
    name: str,
    path: str | Path,
    label: list[list[int]],
    ltype: list[int],
    carbon_t: list[int],
    carbon_p: list[int],
    *,
    monosaccharide_types: list[str] | str | None = None,
    peptide: str | None = None,
    n_replica: int | None = None,
    min_temperature: float | None = None,
    max_temperature: float | None = None,
    structure_file: str | None = None,
    mdp_file: str | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Add one glycan definition to an existing YAML configuration.

    Only parameters that differ from the global defaults need to be supplied.
    The function validates the new glycan before writing the YAML file.

    Example
    -------
    add_glycan_to_yaml(
        "../config/glycans.yaml",
        name="A2G2S2",
        path="../C/A2G2S2/",
        label=[[1, 2], [2, 3], [3, 4]],
        ltype=[4, 4, 6],
        carbon_t=[1, 1, 1],
        carbon_p=[1, 1, 1, 1],
    )
    """
    config_path = Path(config_file).expanduser().resolve()
    config = load_config(config_path)

    if name in config["glycans"] and not overwrite:
        raise SetupError(
            f"Glycan '{name}' already exists. Use overwrite=True to replace it."
        )

    if isinstance(monosaccharide_types, str):
        monosaccharide_types = [
            x.strip() for x in monosaccharide_types.split(",") if x.strip()
        ]

    glycan: dict[str, Any] = {
        "path": str(path),
        "label": label,
        "ltype": ltype,
        "carbon_t": carbon_t,
        "carbon_p": carbon_p,
    }

    optional = {
        "monosaccharide_types": monosaccharide_types,
        "peptide": peptide,
        "n_replica": n_replica,
        "min_temperature": min_temperature,
        "max_temperature": max_temperature,
        "structure_file": structure_file,
        "mdp_file": mdp_file,
    }
    glycan.update({key: value for key, value in optional.items() if value is not None})

    # Validate against defaults before modifying the file.
    candidate = deepcopy(config)
    candidate["glycans"] = deepcopy(config["glycans"])
    candidate["glycans"][name] = glycan

    params = merge_config(candidate, name)
    errors = validate_glycan_config(name, params)
    if errors:
        raise SetupError("Invalid glycan configuration:\n" + "\n".join(errors))

    config["glycans"][name] = glycan

    # Atomic write: write a temporary file first, then replace the original.
    temp_path = config_path.with_suffix(config_path.suffix + ".tmp")
    try:
        with temp_path.open("w") as fh:
            fh.write(
                "# REST-RECT glycan configuration\n"
                "#\n"
                "# Paths are relative to the repository root.\n"
                "# Keep residue/linkage order identical to the GlyCONFORMER labelling scheme.\n\n"
            )
            yaml.safe_dump(
                config,
                fh,
                sort_keys=False,
                allow_unicode=True,
                default_flow_style=False,
            )
        temp_path.replace(config_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()

    print(f"Added glycan '{name}' to {config_path}")
    return glycan


def available_glycans(config: Mapping[str, Any]) -> list[str]:
    return list(config["glycans"].keys())


def overview_dataframe(config: Mapping[str, Any]):
    """Return a compact pandas overview for notebook display."""
    import pandas as pd

    rows = []
    for name in available_glycans(config):
        p = merge_config(config, name)
        rows.append({
            "glycan": name,
            "path": p["path"],
            "replicas": p["n_replica"],
            "T_min [K]": p["min_temperature"],
            "T_max [K]": p["max_temperature"],
            "peptide": p["peptide"],
            "residues": len(p["carbon_p"]),
            "linkages": len(p["label"]),
            "monosaccharide types": ", ".join(p["monosaccharide_types"]),
        })
    return pd.DataFrame(rows)


def _resolve_repo_root(config_file: Path, repo_root: str | Path | None) -> Path:
    if repo_root is not None:
        return Path(repo_root).expanduser().resolve()
    # config/<file>.yaml -> repository root
    return config_file.resolve().parent.parent


def _resolve_simulation_path(repo_root: Path, configured_path: str | Path) -> Path:
    p = Path(configured_path).expanduser()
    if p.is_absolute():
        return p
    # The paths in the original notebook were relative to the notebook/repo root.
    return (repo_root / p).resolve()


#def _copy_template_files(template_dir: Path, destination: Path, patterns: Iterable[str]) -> list[Path]:
#    copied = []
#    for pattern in patterns:
#        matches = list(template_dir.glob(pattern))
#        if not matches:
#            raise SetupError(f"No template matched '{pattern}' in {template_dir}")
#        for source in matches:
#            if source.is_file():
#                target = destination / source.name
#                shutil.copy2(source, target)
#                copied.append(target)
#    return copied


def run_command(
    command: list[str],
    cwd: Path,
    *,
    dry_run: bool = False,
    env: Mapping[str, str] | None = None,
    log_file: Path | None = None,
) -> None:
    """Run a command with a visible command line and proper error propagation."""
    print(f"$ {' '.join(command)}")
    if dry_run:
        return

    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)

    result = subprocess.run(
        command,
        cwd=cwd,
        env=merged_env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with log_file.open("a") as fh:
            fh.write(f"$ {' '.join(command)}\n")
            fh.write(result.stdout or "")
            fh.write("\n")
    print(result.stdout or "", end="")
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command, result.stdout)

def preview_setup(name: str, config: Mapping[str, Any], *, repo_root: str | Path | None = None) -> None:
    """Print exactly what would be done without executing anything."""
    p = merge_config(config, name)
    errors = validate_glycan_config(name, p)
    if errors:
        raise SetupError("\n".join(errors))

    config_root = Path(repo_root).resolve() if repo_root else Path.cwd().resolve()
    sim_path = _resolve_simulation_path(config_root, p["path"])
    print("=" * 68)
    print(f"REST-RECT setup: {name}")
    print("=" * 68)
    print(f"Simulation directory : {sim_path}")
    print(f"Replicas             : {p['n_replica']}")
    print(f"Temperature          : {p['min_temperature']} -> {p['max_temperature']} K")
    print(f"Structure            : {p['structure_file']}")
    print(f"REST-RECT MDP        : {p['mdp_file']}")
    print(f"Peptide              : {p['peptide']}")
    print(f"Residues             : {len(p['carbon_p'])}")
    print(f"Linkages             : {len(p['label'])}")
    print(f"Monosaccharides      : {', '.join(p['monosaccharide_types'])}")
    print("Linkages:")
    print("  " + ", ".join(f"{a}-{b} (type {t})" for (a, b), t in zip(p["label"], p["ltype"])))
    print("=" * 68)

def glygen(
    name: str,
    config: Mapping[str, Any],
    repo_root: str | Path | None = None,
    template_dir: str | Path | None = None,
    use_module_load: bool = False,
    dry_run: bool = False,
    force: bool = False,
) -> Path:
    """Run the complete setup for one configured glycan.

    The glycan directory must already exist and must contain the
    characteristic CHARMM-GUI files:

        step3_input.*
        topol.top
        index.ndx

    The function also reports independently whether minimisation/
    equilibration, REST-RECT production, and analysis have already
    been performed.

    ``use_module_load=True`` prepends the commands with the shell
    module-loading environment used in the original notebook.

    ``dry_run=True`` performs all filesystem/configuration checks but
    does not modify files or execute commands.
    """
    # GROMACS preprocessing and MD stages.
    if use_module_load:
    # Environment Modules are shell functionality. Use bash -lc only when
    # explicitly requested, so normal Python execution remains portable.
        prefix = "ml chem/GROMACS && "
    else:
        prefix = " "

    config_file_root = (
        Path(repo_root).expanduser().resolve()
        if repo_root
        else Path.cwd().resolve()
    )

    p = merge_config(config, name)
    errors = validate_glycan_config(name, p)
    if errors:
        raise SetupError("\n".join(errors))

    sim_path = _resolve_simulation_path(
        config_file_root,
        p["path"],
    )

    print(f"\nGlycan: {name}")
    print(f"Simulation directory: {sim_path}")

    template_dir = (
        Path(template_dir).expanduser().resolve()
        if template_dir
        else config_file_root / "templates"
    )

    if not template_dir.exists():
        raise SetupError(
            f"Template directory does not exist: {template_dir}"
        )

    if not sim_path.is_dir():
        raise SetupError(
            f"Simulation directory for glycan '{name}' does not exist:\n"
            f"    {sim_path}\n\n"
            f"Please check the 'path' entry in glycans.yaml."
        )
    else:
        print("✓ Glycan directory exists")

    charmm_gui_dirs = list(sim_path.glob("charmm-gui*"))
    charmm_gui_dir = charmm_gui_dirs[0]
    gromacs_dir = charmm_gui_dir / "gromacs"
    min_equ_dir = sim_path / "01_min_equ"
    rest_rect_dir = sim_path / "02_REST-RECT"
    analysis_dir = sim_path / "03_gl_data"
    
    if charmm_gui_dir.is_dir():
        print("✓ CHARMM-GUI directory detected")
        if min_equ_dir.is_dir():
            print("✓ Minimisation/equilibration already performed")
            if rest_rect_dir.is_dir():
                print("✓ REST-RECT simulation already setup")
                if analysis_dir.is_dir():
                    print("✓ Analysis already performed")
                else:
                    print("✗ Analysis not yet performed")
                    #to be continued!
            else:
                _setup_rest_rect(sim_path, template_dir, p, min_equ_dir, name, dry_run)
        else:
            _perform_min_equ(gromacs_dir, sim_path, template_dir, p, dry_run)
            _setup_rest_rect(sim_path, template_dir, p, min_equ_dir, name, dry_run)
    else:
        print("✗ CHARMM-GUI directory not present")

def _cluster_command(binary: str, *args: str) -> list[str]:
    return ["bash", "-lc", " ".join([binary, *args])]
    
    return [binary, *args]
def _perform_min_equ(gromacs_dir: Path, sim_path: Path, template_dir: Path, p: dict, dry_run: bool = False):
    print("✗ Minimisation/equilibration not yet performed")
    if not gromacs_dir.is_dir():
        raise SetupError(
            f"CHARMM-GUI gromacs directory not found: {gromacs_dir}"
        )
    if not dry_run:
        shutil.copytree(
            gromacs_dir,
            sim_path,
            dirs_exist_ok=True,
        )
    log_file = sim_path / "rest_rect_setup.log"
    
    if not dry_run:
        for filename in [
            "step4.0_minimization.mdp",
            "step4.1_equilibration.mdp",
            "step5_production.mdp",
        ]:
            shutil.copy2(
                template_dir / filename,
                sim_path / filename,
            )
        with log_file.open("a") as fh:
            fh.write(f"\n=== {datetime.now(timezone.utc).isoformat()} | {name} ===\n")
            fh.write(json.dumps(p, indent=2, default=str) + "\n")
    # Optional Neu5Ac force-field edit from the original notebook.
    if p.get("modify_forcefield_for_neu5ac", False) and not dry_run:
        remove_lines_and_insert(sim_path, p.get("forcefield_file", "toppar/forcefield.itp"))
    if dry_run:
        print(f"\n--- {name}: dry run ---")
    run_command(_cluster_command("gmx",
        "grompp", "-f", "step4.0_minimization.mdp", "-o", "step4.0_minimization.tpr",
        "-c", "step3_input.gro", "-r", "step3_input.gro", "-p", "topol.top", "-n", "index.ndx"),
        sim_path, dry_run=dry_run, log_file=log_file)
    run_command(_cluster_command("gmx_mpi", "mdrun", "-deffnm", "step4.0_minimization"),
        sim_path, dry_run=dry_run, log_file=log_file)

    run_command(_cluster_command("gmx",
        "grompp", "-f", "step4.1_equilibration.mdp", "-o", "step4.1_equilibration.tpr",
        "-c", "step4.0_minimization.gro", "-r", "step4.0_minimization.gro",
        "-p", "topol.top", "-n", "index.ndx"),
        sim_path, dry_run=dry_run, log_file=log_file)
    run_command(_cluster_command("gmx_mpi", "mdrun", "-deffnm", "step4.1_equilibration"),
        sim_path, dry_run=dry_run, log_file=log_file)

    run_command(_cluster_command("gmx",
        "grompp", "-f", "step5_production.mdp", "-o", "step5_production.tpr",
        "-c", "step4.1_equilibration.gro", "-r", "step4.1_equilibration.gro",
        "-p", "topol.top", "-n", "index.ndx"),
        sim_path, dry_run=dry_run, log_file=log_file)
    run_command(_cluster_command("gmx_mpi", "mdrun", "-deffnm", "step5_production"),
        sim_path, dry_run=dry_run, log_file=log_file)
    min_equ = sim_path / "01_min_equ"
    if not dry_run:
            min_equ.mkdir(exist_ok=True)
    min_equ_names = [
    "step3_input.gro","step3_input.pdb","step3_input.psf",
    "step4.0_minimization.mdp", "step4.0_minimization.tpr", "step4.0_minimization.gro",
    "step4.0_minimization.edr", "step4.0_minimization.log", "step4.0_minimization.mdp",
    "step4.0_minimization.trr",
    "step4.1_equilibration.mdp", "step4.1_equilibration.tpr", "step4.1_equilibration.gro",
    "step4.1_equilibration.edr", "step4.1_equilibration.log", "step4.1_equilibration.cpt",
    "step4.1_equilibration.xtc",
    "step5_production.mdp", "step5_production.tpr", "step5_production.gro",
    "step5_production.edr", "step5_production.log", "step5_production.cpt", "step5_production.xtc",
    "topol.top", "topol.tpr", "toppar", "index.ndx"
    ]
    for name in min_equ_names:
        src = sim_path / name
        if src.exists():
            if dry_run:
                print(f"MOVE {src.name} -> 01_min_equ/")
            else:
                shutil.move(str(src), str(min_equ / src.name))

def _setup_rest_rect(sim_path, template_dir, p, min_equ_dir, name, dry_run):
    print("✗ REST-RECT simulation not yet setup:")
    log_file = sim_path / "rest_rect_setup.log"
    if not dry_run:
        for filename in [
        "production.mdp",
        "job.sh",
        "partial_tempering_pairs.sh",
        "generate_HREMD.sh"
        ]:
            shutil.copy2(
                template_dir / filename,
                sim_path / filename,
            )
        with log_file.open("a") as fh:
            fh.write(f"\n=== {datetime.now(timezone.utc).isoformat()} | {name} ===\n")
            fh.write(json.dumps(p, indent=2, default=str) + "\n")

    files_to_check = ["plumed_allangles.dat", "reference.pdb", "processed.top", "processed_m.top"]
    if not dry_run:
        check_and_delete_files_and_directories(sim_path, files_to_check)

    structure_file = min_equ_dir / p["structure_file"]
    index_file = min_equ_dir/ "index.ndx"
    topol_file = min_equ_dir/ "topol.top"

    run_command(_cluster_command("gmx", "editconf", "-f", str(structure_file.relative_to(sim_path)),
        "-o", "reference.pdb", "-resnr", "1"),
        sim_path, dry_run=dry_run, log_file=log_file)

    if not dry_run:
        write_plumed(
            "plumed_allangles.dat",
            p["peptide"],
            p["label"],
            p["ltype"],
            p["carbon_t"],
            p["carbon_p"],
            str(sim_path)
        )

    run_command(_cluster_command("gmx", "grompp", "-f", p["mdp_file"],
        "-c", str(structure_file.relative_to(sim_path)), "-pp", "processed.top", "-n", str(index_file.relative_to(sim_path)), "-p", str(topol_file.relative_to(sim_path))),
        sim_path, dry_run=dry_run, log_file=log_file)

    if not dry_run:
        define_solute(",".join(p["monosaccharide_types"]), str(sim_path))

    run_command(["bash", "./generate_HREMD.sh", str(structure_file.relative_to(sim_path)), str(p["n_replica"]),
                 str(p["min_temperature"]), str(p["max_temperature"]), p["mdp_file"], str(index_file.relative_to(sim_path))],
                sim_path, dry_run=dry_run, log_file=log_file)

    _organize_output(sim_path, dry_run=dry_run)

    if not dry_run:
        _write_metadata(sim_path, name, p)

    print(f"\nCompleted: {name}")
    return sim_path

def _organize_output(sim_path: Path, *, dry_run: bool = False) -> None:
    """Organise files exactly once, without destructive wildcard deletes."""
    
    rest_rect = sim_path / "02_REST-RECT"

    if not dry_run:
        rest_rect.mkdir(exist_ok=True)

    rest_names = [
        "generate_HREMD.sh", "index.ndx", "job.sh", "mdout.mdp", "production.mdp",
        "plumed_allangles.dat", "reference.pdb", "processed.top", "processed_m.top",
        "partial_tempering_pairs.sh", "topol.tpr"
    ]
    for name in rest_names:
        src = sim_path / name
        if src.exists():
            if dry_run:
                print(f"MOVE {src.name} -> 02_REST-RECT/")
            else:
                shutil.move(str(src), str(rest_rect / src.name))

    for src in sim_path.glob("TOPO*"):
        if src.is_dir():
            destination = rest_rect / src.name
        if dry_run:
            print(f"MOVE {src.name} -> 02_REST-RECT/")

        else:
            shutil.move(
                str(src),
                str(rest_rect / src.name)
            )

    if not dry_run:
        readme = sim_path / "README"
        if readme.exists():
            readme.unlink()


def _write_metadata(sim_path: Path, name: str, params: Mapping[str, Any]) -> None:
    metadata = {
        "glycan": name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": dict(params),
    }
    with (sim_path / "metadata.yaml").open("w") as fh:
        yaml.safe_dump(metadata, fh, sort_keys=False)


__all__ = [
    "SetupError",
    "load_config",
    "merge_config",
    "validate_glycan_config",
    "validate_config",
    "available_glycans",
    "overview_dataframe",
    "preview_setup",
    "glygen",
    "run_command",
]

def update_glycan_in_yaml(
    config_file: str | Path,
    name: str,
    **updates,
) -> dict[str, Any]:
    """Update one or more parameters of an existing glycan."""

    config_path = Path(config_file).expanduser().resolve()
    config = load_config(config_path)

    if name not in config["glycans"]:
        raise SetupError(f"Glycan '{name}' does not exist in the YAML file.")

    # Update requested parameters
    for key, value in updates.items():
        config["glycans"][name][key] = value

    # Validate after modification
    params = merge_config(config, name)
    errors = validate_glycan_config(name, params)

    if errors:
        raise SetupError(
            "Invalid glycan configuration:\n" +
            "\n".join(errors)
        )

    # Write updated YAML
    with config_path.open("w") as fh:
        fh.write(
            "# REST-RECT glycan configuration\n"
            "#\n"
            "# Paths are relative to the repository root.\n"
            "# Keep residue/linkage order identical to the GlyCONFORMER "
            "labelling scheme.\n\n"
        )

        yaml.safe_dump(
            config,
            fh,
            sort_keys=False,
            allow_unicode=True,
            default_flow_style=False,
        )

    print(f"Updated glycan '{name}' in {config_path}")

    return config["glycans"][name]

def write_plumed(filename, peptide, label, ltype, carbon_t, carbon_p, directory):
    # Check if peptide is either "ASA", "ATA" or "ANA"
    if peptide not in ["ASA", "ANA", "ATA", "none"]:
        raise ValueError("Peptide must be either 'ASA', 'ATA', 'ANA' or 'none'")

    # Check if label and type lists are the same length
    if len(label) != len(ltype):
        raise ValueError("Label and type lists must be the same length")

    # Check if directory exists
    if not os.path.exists(directory):
        raise FileNotFoundError(f"The directory {directory} does not exist")

    # Check if we have write permissions to the directory
    if not os.access(directory, os.W_OK):
        raise PermissionError(f"Do not have write permissions for the directory {directory}")

    try:
        with open(os.path.join(directory, filename), 'w') as f:
            if peptide == "ASA":
                f.write("MOLINFO STRUCTURE=reference.pdb\n#\nphi_gb: TORSION ATOMS=@O5-4,@C1-4,@OG-2,@CB-2\npsi_gb: TORSION ATOMS=@C1-4,@OG-2,@CB-2,@CA-2\n")
                f.write("chi1: TORSION ATOMS=@OG-2,@CB-2,@CA-2,@N-2\n#\n")
                f.write("metaD_phi_gb: METAD ARG=phi_gb TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_phi_gb\n")
                f.write("metaD_psi_gb: METAD ARG=psi_gb TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_psi_gb\n")
                f.write("metaD_chi1: METAD ARG=chi1 TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_chi1\n#\n")
            elif peptide == "ANA":
                f.write("MOLINFO STRUCTURE=reference.pdb\n#\nphi_gb: TORSION ATOMS=@O5-4,@C1-4,@ND2-2,@CG-2\npsi_gb: TORSION ATOMS=@C1-4,@ND2-2,@CG-2,@CB-2\n")
                f.write("chi2: TORSION ATOMS=@ND2-2,@CG-2,@CB-2,@CA-2\nchi1: TORSION ATOMS=@CG-2,@CB-2,@CA-2,@N-2\n#\n")
                f.write("metaD_phi_gb: METAD ARG=phi_gb TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_phi_gb\n")
                f.write("metaD_psi_gb: METAD ARG=psi_gb TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_psi_gb\n")
                f.write("metaD_chi1: METAD ARG=chi1 TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_chi1\n")
                f.write("metaD_chi2: METAD ARG=chi2 TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_chi2\n#\n")
            elif peptide == "ATA":
                f.write("MOLINFO STRUCTURE=reference.pdb\n#\nphi_gb: TORSION ATOMS=@O5-4,@C1-4,@OG1-2,@CB-2\npsi_gb: TORSION ATOMS=@C1-4,@OG1-2,@CB-2,@CA-2\n")
                f.write("chi1: TORSION ATOMS=@OG1-2,@CB-2,@CA-2,@N-2\n#\n")
                f.write("metaD_phi_gb: METAD ARG=phi_gb TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_phi_gb\n")
                f.write("metaD_psi_gb: METAD ARG=psi_gb TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_psi_gb\n")
                f.write("metaD_chi1: METAD ARG=chi1 TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_chi1\n#\n")
            elif peptide == "none":
                f.write("MOLINFO STRUCTURE=reference.pdb\n#\n")


            if peptide in ["ASA", "ANA", "ATA"]:
                for e,t,l in zip(range(len(label)),ltype,carbon_t): 
                    if t == 6 and l == 1:
                        f.write(f"phi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O5-{label[e][1]+3},@C1-{label[e][1]+3},@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3}\n")
                        f.write(f"psi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@C1-{label[e][1]+3},@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3},@C{ltype[e]-1}-{label[e][0]+3}\n")
                        f.write(f"omega{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3},@C{ltype[e]-1}-{label[e][0]+3},@C{ltype[e]-2}-{label[e][0]+3}\n#\n")
                    elif t == 6 and l == 2:
                        f.write(f"phi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O6-{label[e][1]+3},@C2-{label[e][1]+3},@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3}\n")
                        f.write(f"psi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@C2-{label[e][1]+3},@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3},@C{ltype[e]-1}-{label[e][0]+3}\n")
                        f.write(f"omega{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3},@C{ltype[e]-1}-{label[e][0]+3},@C{ltype[e]-2}-{label[e][0]+3}\n#\n")
                    elif t != 6 and l == 2:
                        f.write(f"phi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O6-{label[e][1]+3},@C2-{label[e][1]+3},@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3}\n")
                        f.write(f"psi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@C2-{label[e][1]+3},@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3},@C{ltype[e]-1}-{label[e][0]+3}\n#\n")
                    else:
                        f.write(f"phi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O5-{label[e][1]+3},@C1-{label[e][1]+3},@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3}\n")
                        f.write(f"psi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@C1-{label[e][1]+3},@O{ltype[e]}-{label[e][0]+3},@C{ltype[e]}-{label[e][0]+3},@C{ltype[e]-1}-{label[e][0]+3}\n#\n")
                for i,l in zip(range(1,len(carbon_p)+1),carbon_p):
                    if l == 1:
                        f.write(f"puck{i}: PUCKERING ATOMS=@O5-{i+3},@C1-{i+3},@C2-{i+3},@C3-{i+3},@C4-{i+3},@C5-{i+3}\n")
                    else:
                        f.write(f"puck{i}: PUCKERING ATOMS=@O6-{i+3},@C2-{i+3},@C3-{i+3},@C4-{i+3},@C5-{i+3},@C6-{i+3}\n")

                f.write(f"#\n")
            
            else:
                for e,t,l in zip(range(len(label)),ltype,carbon_t): 
                    if t == 6 and l == 1:
                        f.write(f"phi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O5-{label[e][1]},@C1-{label[e][1]},@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]}\n")
                        f.write(f"psi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@C1-{label[e][1]},@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]},@C{ltype[e]-1}-{label[e][0]}\n")
                        f.write(f"omega{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]},@C{ltype[e]-1}-{label[e][0]},@C{ltype[e]-2}-{label[e][0]}\n#\n")
                    elif t == 6 and l == 2:
                        f.write(f"phi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O6-{label[e][1]},@C2-{label[e][1]},@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]}\n")
                        f.write(f"psi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@C2-{label[e][1]},@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]},@C{ltype[e]-1}-{label[e][0]}\n")
                        f.write(f"omega{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]},@C{ltype[e]-1}-{label[e][0]},@C{ltype[e]-2}-{label[e][0]}\n#\n")
                    elif t != 6 and l == 2:
                        f.write(f"phi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O6-{label[e][1]},@C2-{label[e][1]},@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]}\n")
                        f.write(f"psi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@C2-{label[e][1]},@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]},@C{ltype[e]-1}-{label[e][0]}\n#\n")
                    else:
                        f.write(f"phi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@O5-{label[e][1]},@C1-{label[e][1]},@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]}\n")
                        f.write(f"psi{label[e][0]}_{label[e][1]}: TORSION ATOMS=@C1-{label[e][1]},@O{ltype[e]}-{label[e][0]},@C{ltype[e]}-{label[e][0]},@C{ltype[e]-1}-{label[e][0]}\n#\n")
                for i,l in zip(range(1,len(carbon_p)+1),carbon_p):
                    if l == 1:
                        f.write(f"puck{i}: PUCKERING ATOMS=@O5-{i},@C1-{i},@C2-{i},@C3-{i},@C4-{i},@C5-{i}\n")
                    else:
                        f.write(f"puck{i}: PUCKERING ATOMS=@O6-{i},@C2-{i},@C3-{i},@C4-{i},@C5-{i},@C6-{i}\n")
                
                f.write(f"#\n")

            for e,t in zip(range(len(label)),ltype):
                if t == 6:
                    f.write(f"metaD_phi{label[e][0]}_{label[e][1]}: METAD ARG=phi{label[e][0]}_{label[e][1]} TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_phi{label[e][0]}_{label[e][1]}\n")
                    f.write(f"metaD_psi{label[e][0]}_{label[e][1]}: METAD ARG=psi{label[e][0]}_{label[e][1]} TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_psi{label[e][0]}_{label[e][1]}\n")
                    f.write(f"metaD_omega{label[e][0]}_{label[e][1]}: METAD ARG=omega{label[e][0]}_{label[e][1]} TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_omega{label[e][0]}_{label[e][1]}\n#\n")
                else:
                    f.write(f"metaD_phi{label[e][0]}_{label[e][1]}: METAD ARG=phi{label[e][0]}_{label[e][1]} TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_phi{label[e][0]}_{label[e][1]}\n")
                    f.write(f"metaD_psi{label[e][0]}_{label[e][1]}: METAD ARG=psi{label[e][0]}_{label[e][1]} TAU=4.0 PACE=500 SIGMA=0.35 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP=310.15 RECT=1,1.2,1.46,1.82,2.3,2.94,3.78,4.89,6.34,8.23,10.7,14 FILE=HILLS_psi{label[e][0]}_{label[e][1]}\n#\n")
            f.write(f"PRINT ARG=* STRIDE=5000 FILE=COLVAR\n#\nPRINT ARG=*.theta,*.phi STRIDE=5000 FILE=COLVAR_theta\n#\nPRINT ARG=*.bias STRIDE=5000 FILE=COLVAR_metaD")

    except IOError:
        print(f"Error writing to file {filename}")

def define_solute(target_words_input, path):

    target_words = [word.strip() for word in target_words_input.split(",")]

    try:
        with open(os.path.join(path, "processed.top"), "r") as infile:
            lines = infile.readlines()

        # Process each line and print the modified output
        with open(os.path.join(path, "processed_m.top"), "w") as outfile:
            for line in lines:
                # Split the line into columns
                columns = line.split()
                # Check if the fourth column matches any of the target words
                if len(columns) >= 4 and columns[3] in target_words:
                # Append an underscore to the second column
                    columns[1] += "_"
                # Join the modified columns back into a line
                modified_line = " ".join(columns)                
                outfile.write(modified_line + "\n")

        print(f"Lines saved")
    except FileNotFoundError:
        print("File not found. Please check the file path.")

def check_and_delete_files_and_directories(directory_path, files_to_check):
    
    # Check and delete the files
    for file_name in files_to_check:
        file_path = os.path.join(directory_path, file_name)
        
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                print(f"{file_name} found and deleted.")
            except OSError as e:
                print(f"Error deleting {file_name}: {e}")
        else:
            print(f"{file_name} does not exist in the directory.")

    # Check and delete directories starting with "TOPO"
    for item in os.listdir(directory_path):
        if item.startswith("TOPO") and os.path.isdir(os.path.join(directory_path, item)):
            dir_path = os.path.join(directory_path, item)
            try:
                shutil.rmtree(dir_path)  # Use shutil.rmtree to delete directories and their contents
                print(f"Directory {item} found and deleted.")
            except OSError as e:
                print(f"Error deleting directory {item}: {e}")

def remove_lines_and_insert(path, filename):
    # Read the entire file into memory
    with open(os.path.join(path, filename), 'r') as file:
        lines = file.readlines()

    # Search for the lines to be removed
    lines_to_remove = [
    " CC3261  CC3161  CC3161   NC2D1     9  1.800000e+02  8.266329e+00      1\n",
    " CC3261  CC3161  CC3161   NC2D1     9  0.000000e+00  7.073889e+00      2\n",
	" CC3261  CC3161  CC3161   NC2D1     9  0.000000e+00  4.167682e+00      3\n",
	" CC3261  CC3161  CC3161   NC2D1     9  0.000000e+00  1.317123e+00      6\n",
	"  NC2D1  CC3161  CC3163   CC312     9  1.800000e+02  8.903970e+00      1\n",
	"  NC2D1  CC3161  CC3163   CC312     9  0.000000e+00  5.750490e+00      2\n",
	"  NC2D1  CC3161  CC3163   CC312     9  0.000000e+00  2.450987e+00      3\n",
	"  NC2D1  CC3161  CC3163   CC312     9  0.000000e+00  3.709116e+00      6\n",
	"  NC2D1  CC3161  CC3163  OC3C61     9  0.000000e+00  7.568856e+00      1\n",
	"  NC2D1  CC3161  CC3163  OC3C61     9  1.800000e+02  2.875663e+00      2\n",
	"  NC2D1  CC3161  CC3163  OC3C61     9  1.800000e+02  5.568904e+00      3\n",
	"  NC2D1  CC3161  CC3163  OC3C61     9  1.800000e+02  2.367307e+00      6\n"
    ]

    # Remove the specified lines
    updated_lines = [line for line in lines if line not in lines_to_remove]

    # Insert the new lines after the specified comment line
    comment_line = "; i	j	k	l	func	phi0	Kphi	mult\n"
    new_lines = [
        "  NC2D1  CC3161  CC3161  CC3261     9     0.000000     0.836800     3\n",
        "  CC312  CC3163  CC3161   NC2D1     9     0.000000     0.836800     3\n",
        " OC3C61  CC3163  CC3161   NC2D1     9     0.000000     0.836800     3\n"
    ]

    # Find the index of the comment line
    comment_line_index = updated_lines.index(comment_line)

    # Insert the new lines after the comment line
    updated_lines[comment_line_index + 1:comment_line_index + 1] = new_lines

    # Write the updated content back to the file
    with open(os.path.join(path, filename), 'w') as file:
        file.writelines(updated_lines)
    
    # Print messages indicating the changes
    print("Lines removed:")
    for line in lines_to_remove:
        if line in lines:
            print(line.strip())

    print("\nLines added:")
    for line in new_lines:
        print(line.strip())
