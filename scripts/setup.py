"""GlycoGEN: per-glycan JSON driven, resumable CHARMM-GUI -> REST-RECT workflow.

Public API
----------
    status(glycan_id, data_dir)
    init_glycan_simulation()
    generate_glycan_simulation(glycan_id, data_dir, template_dir, ...)
    analyze_glycan(glycan_id, data_dir, ...)

All chemical information is read from the unique <glycan_id>/<glycan_id>.json
found recursively below data_dir.
No config.yaml, metadata.yaml, label/ltype/carbon lists, or dry_run.

Important: the user must run the generated REST-RECT jobs on the HPC separately.
The analysis is only attempted after trajectory outputs have been downloaded.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path



class SetupError(RuntimeError):
    pass


SIALIC_NAMES = {"Neu", "Neu5Ac", "Neu5Gc", "Kdn", "Kdo", "Sia", "Pse", "Leg", "Aci", "4eLeg"}
# CHARMM-specific correction in the original setup.py is validated for Neu5Ac,
# not for arbitrary nonulosonic acids.
NEU5AC_NAMES = {"Neu5Ac"}
REST_RECT_FACTORS = [1, 1.2, 1.46, 1.82, 2.3, 2.94, 3.78, 4.89, 6.34, 8.23, 10.7, 14]
MIN_EQU_STEPS = (
    ("step4.0_minimization", "step3_input.gro", "step3_input.gro"),
    ("step4.1_equilibration", "step4.0_minimization.gro", "step4.0_minimization.gro"),
    ("step5_production", "step4.1_equilibration.gro", "step4.1_equilibration.gro"),
)

def _log(path: Path, message: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(f"[{datetime.now(timezone.utc).isoformat()}] {message}\n")
    print(message)


def _cmd(args, cwd: Path, log: Path, *, module: str | None = None, input_text: str | None = None):
    args = [str(a) for a in args]
    if module:
        # bash -l makes environment modules available on sites that initialize
        # them in login shells; module needs to exist at the target site.
        actual = ["bash", "-lc", f"module load {shlex.quote(module)} && {shlex.join(args)}"]
    else:
        actual = args
    _log(log, f"RUN: {shlex.join(actual)}")
    result = subprocess.run(actual, cwd=cwd, input=input_text, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, check=False)
    with log.open("a", encoding="utf-8") as f:
        f.write(result.stdout or "")
        f.write("\n")
    if result.stdout:
        print(result.stdout[-2500:])
    if result.returncode:
        raise SetupError(f"Command failed with exit code {result.returncode}: {shlex.join(actual)}. See {log}")


def find_glycan_directory(root: str | Path, glycan_id: str) -> Path:
    """Find the unique glycan directory anywhere below the library root.

    Accepts a direct glycan directory as ``root`` as well as the library root.
    Directories are matched by their exact name, not by a substring; ambiguous
    matches raise an error instead of silently selecting the wrong glycan.
    """
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Glycan library directory not found: {root}")
    if not glycan_id or glycan_id in (".", "..") or Path(glycan_id).name != glycan_id:
        raise ValueError(f"Invalid glycan_id: {glycan_id!r}")

    matches = []
    if root.name == glycan_id:
        matches.append(root)
    # rglob searches all library branches, e.g. N-glycans/High-mannose/Man5.
    matches.extend(p for p in root.rglob(glycan_id) if p.is_dir())
    unique = sorted(set(matches))
    if not unique:
        raise FileNotFoundError(f"Glycan directory {glycan_id!r} not found under {root}")
    if len(unique) > 1:
        raise SetupError(
            f"Multiple directories named {glycan_id!r} found under {root}: "
            + ", ".join(str(path) for path in unique)
        )
    return unique[0]


def _load(glycan_id: str, data_dir: str | Path):
    folder = find_glycan_directory(data_dir, glycan_id)
    source = folder / f"{glycan_id}.json"
    if not source.is_file():
        raise FileNotFoundError(f"Missing metadata: {source}")
    data = json.loads(source.read_text(encoding="utf-8"))
    if data.get("glycan_id") != glycan_id:
        raise SetupError(f"glycan_id does not match filename: {source}")
    for key in ("monosaccharides", "pdb_residues", "pucker_ring_atoms", "graph_edges", "simulation"):
        if key not in data:
            raise SetupError(f"{source}: missing {key!r}")
    ids = set(data["monosaccharides"])
    if ids != set(data["pdb_residues"]) or ids != set(data["pucker_ring_atoms"]):
        raise SetupError("monosaccharides, pdb_residues and pucker_ring_atoms keys must match")
    sim = data["simulation"]
    for key in ("n_replica", "min_temperature", "max_temperature", "structure_file",
                "mdp_file", "template_files", "forcefield_file"):
        if key not in sim:
            raise SetupError(f"{source}: missing simulation.{key}")
    if int(sim["n_replica"]) < 2 or float(sim["min_temperature"]) <= 0 or float(sim["max_temperature"]) <= float(sim["min_temperature"]):
        raise SetupError(f"{source}: invalid replica number or temperature range")
    for e in data["graph_edges"]:
        p, c = str(e["parent"]), str(e["child"])
        if p not in ids or c not in ids or not re.fullmatch(r"[abAB][12]-\d+", str(e["linkage"])):
            raise SetupError(f"Invalid graph edge: {e}")
        donor = int(str(e["linkage"]).split("-")[0][1:])
        sugar = data["monosaccharides"][c]
        if (sugar in SIALIC_NAMES) != (donor == 2):
            raise SetupError(f"Donor C{donor} inconsistent with child monosaccharide {sugar}: {e}")
    return folder, data


def _nonempty(p: Path) -> bool:
    return p.is_file() and p.stat().st_size > 0


def _charmm_gromacs(folder: Path) -> Path | None:
    hits = sorted(p / "gromacs" for p in folder.glob("*charmm-gui*") if (p / "gromacs").is_dir())
    if len(hits) > 1:
        raise SetupError(f"Multiple CHARMM-GUI GROMACS directories: {hits}")
    return hits[0] if hits else None


def _min_done(folder: Path) -> bool:
    directory = folder / "01_min_equ"
    return all(_nonempty(directory / f"{step}.gro") and _nonempty(directory / f"{step}.log")
               for step, _, _ in MIN_EQU_STEPS)


def _plumed_missing_definitions(path: Path, data: dict) -> list[str]:
    """Check expected glycosidic torsions and their METAD biases, not just file existence."""
    if not _nonempty(path):
        return ["PLUMED input missing"]
    _, torsions = _torsion_definitions(data)
    peptide_lines, peptide_torsions = _peptide_torsions(data)
    text = path.read_text(encoding="utf-8")
    definitions = set(re.findall(r"^([A-Za-z][\w]*):\s*TORSION\s", text, re.M))
    biases = set(re.findall(r"^metaD_([A-Za-z][\w]*):\s*METAD\s", text, re.M))
    expected = set(torsions + peptide_torsions)
    return ([f"TORSION {name}" for name in sorted(expected - definitions)]
            + [f"METAD {name}" for name in sorted(expected - biases)])


def _validate_plumed_atoms(plumed: Path, pdb: Path):
    """Ensure every @ATOM-RESID used by PLUMED occurs exactly once in reference PDB."""
    if not pdb.is_file():
        raise FileNotFoundError(pdb)
    from collections import Counter
    counts = Counter()
    for line in pdb.open(encoding="utf-8", errors="replace"):
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                counts[(line[12:16].strip(), int(line[22:26]))] += 1
            except ValueError:
                continue
    required = set((atom, int(resid)) for atom, resid in re.findall(
        r"@([A-Za-z0-9]+)-(\d+)", plumed.read_text(encoding="utf-8")))
    invalid = [(atom, resid, counts[(atom, resid)]) for atom, resid in sorted(required)
               if counts[(atom, resid)] != 1]
    if invalid:
        raise SetupError(f"PLUMED atom references missing or ambiguous in {pdb}: {invalid[:20]}")


def _rest_ready(folder: Path, data: dict) -> bool:
    rest = folder / "02_REST-RECT"
    n = int(data["simulation"]["n_replica"])
    return (all(_nonempty(rest / x) for x in ("plumed_allangles.dat", "reference.pdb", "processed_m.top", "job.sh"))
            and not _plumed_missing_definitions(rest / "plumed_allangles.dat", data)
            and all((rest / f"TOPO{i}" / "topol.tpr").is_file() for i in range(n)))


def _production_done(folder: Path, data: dict) -> bool:
    rest = folder / "02_REST-RECT"
    n = int(data["simulation"]["n_replica"])
    return (_nonempty(rest / "TOPO0" / "md.log")
            and _nonempty(rest / "TOPO0" / "COLVAR.0")
            and _nonempty(rest / "TOPO0" / "traj_comp.xtc"))


def _analysis_done(folder: Path, data: dict) -> bool:
    out = folder / "03_lib"
    hills = list((folder / "02_REST-RECT" / "TOPO0").glob("HILLS_*"))
    expected = [out / f"fes_{h.name.removeprefix('HILLS_').replace('.', '_')}.dat" for h in hills]
    return (bool(hills) and all(map(_nonempty, expected))
            and all(_nonempty(out / f) for f in ("colvar.dat", "angles.dat", "omega_angles.dat",
                                                 "conformations_rest-rect_aligned.xtc", "reference_rest-rect.pdb")))


def status(glycan_id: str, data_dir: str | Path = "data") -> dict[str, bool]:
    folder, data = _load(glycan_id, data_dir)
    result = {
        "charmm_gui": _charmm_gromacs(folder) is not None,
        "min_equ": _min_done(folder),
        "rest_rect_setup": _rest_ready(folder, data),
        "rest_rect_production": _production_done(folder, data),
        "analysis": _analysis_done(folder, data),
    }
    print(f"\n{glycan_id}: " + ", ".join(f"{k}={'✓' if v else '—'}" for k, v in result.items()))
    return result


def _copy_templates(source: Path, dest: Path, patterns: list[str], *, required: tuple[str, ...] = ()):
    copied = set()
    for pattern in patterns:
        hits = [p for p in source.glob(pattern) if p.is_file()]
        if not hits and not any(ch in pattern for ch in "*?[]"):
            raise FileNotFoundError(source / pattern)
        for item in hits:
            shutil.copy2(item, dest / item.name)
            copied.add(item.name)
    for filename in required:
        if not (dest / filename).is_file():
            raise FileNotFoundError(f"Missing required template: {filename} (searched {source})")
    return copied


def _apply_neu5ac_correction(path: Path, log: Path):
    """Apply original CHARMM36 Neu5Ac dihedral replacement once, and fail if source layout changes."""
    if not path.is_file():
        raise FileNotFoundError(path)
    text = path.read_text()
    # Original correction was scoped to two specific atom-type patterns.
    inserts = [
        "  NC2D1  CC3161  CC3161  CC3261     9     0.000000     0.836800     3",
        "  CC312  CC3163  CC3161   NC2D1     9     0.000000     0.836800     3",
        " OC3C61  CC3163  CC3161   NC2D1     9     0.000000     0.836800     3",
    ]
    if all(x in text for x in inserts):
        _log(log, "Neu5Ac force-field modification already present")
        return
    # Avoid assuming the original CHARMM FF parameter ordering if version differs.
    old_atoms = {("CC3261", "CC3161", "CC3161", "NC2D1"),
                 ("NC2D1", "CC3161", "CC3163", "CC312"),
                 ("NC2D1", "CC3161", "CC3163", "OC3C61")}
    lines = text.splitlines(keepends=True)
    match_lines = [i for i, line in enumerate(lines)
                   if (tokens := line.split()) and tuple(tokens[:4]) in old_atoms and len(tokens) >= 8]
    if not match_lines or "; i" not in text:
        raise SetupError(f"Cannot safely apply Neu5Ac correction to {path}; FF layout differs from original")
    # Original code removed precisely 12 Fourier terms.
    if len(match_lines) != 12:
        raise SetupError(f"Expected 12 old Neu5Ac terms in {path}, found {len(match_lines)}")
    filtered = [line for i, line in enumerate(lines) if i not in set(match_lines)]
    idx = next((i for i, line in enumerate(filtered) if re.match(r"\s*;\s*i\s+j\s+k\s+l\s+func", line)), None)
    if idx is None:
        raise SetupError("Could not locate dihedrals section header for modification")
    filtered[idx + 1:idx + 1] = [line + "\n" for line in inserts]
    backup = path.with_suffix(path.suffix + ".original")
    if not backup.exists():
        shutil.copy2(path, backup)
    path.write_text("".join(filtered))
    _log(log, f"Applied Neu5Ac CHARMM correction to {path} (backup: {backup})")


def _minimize(folder: Path, gui: Path, templates: Path, data: dict, log: Path,
              gmx: str, gmx_mpi: str, module: str | None):
    target = folder / "01_min_equ"
    target.mkdir(exist_ok=True)
    # Populate only absent CHARMM-GUI files; never overwrite an interrupted run.
    if gui.resolve() != target.resolve():
        for source in gui.rglob("*"):
            if source.is_file():
                destination = target / source.relative_to(gui)
                if not destination.exists():
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, destination)
    for name in ("step4.0_minimization.mdp", "step4.1_equilibration.mdp", "step5_production.mdp"):
        destination = target / name
        if not destination.exists():
            source = templates / name
            if not source.is_file():
                raise FileNotFoundError(source)
            shutil.copy2(source, destination)

    present_sia = set(data["monosaccharides"].values()) & SIALIC_NAMES
    if present_sia:
        if present_sia - NEU5AC_NAMES:
            raise SetupError(f"Sialic acid(s) {present_sia - NEU5AC_NAMES} detected, but only Neu5Ac FF correction is validated")
        ff_path = target / data["simulation"]["forcefield_file"]
        _apply_neu5ac_correction(ff_path, log)

    for step, coordinate, restraint in MIN_EQU_STEPS:
        if _nonempty(target / f"{step}.gro") and _nonempty(target / f"{step}.log"):
            _log(log, f"Skipping completed {step}")
            continue
        if not (target / coordinate).is_file():
            raise SetupError(f"Input for {step} missing: {target / coordinate}")
        _cmd([gmx, "grompp", "-f", f"{step}.mdp", "-o", f"{step}.tpr",
              "-c", coordinate, "-r", restraint, "-p", "topol.top", "-n", "index.ndx"],
             target, log, module=module)
        _cmd([gmx_mpi, "mdrun", "-deffnm", step], target, log, module=module)
        if not _nonempty(target / f"{step}.gro"):
            raise SetupError(f"{step} ended without a nonempty .gro output")


def _torsion_definitions(data: dict):
    """Build glycosidic torsions and ring puckers from the glycan JSON.

    Neu5Ac alpha2-8 edges additionally use omega7, omega8 and omega9,
    following the supplied plumed_polySia3.dat atom definitions. Ordinary
    1-6 linkages retain their single omega torsion.
    """
    node_resids = {str(k): int(v) for k, v in data["pdb_residues"].items()}
    lines = []
    names = []
    for edge in data["graph_edges"]:
        parent, child = str(edge["parent"]), str(edge["child"])
        match = re.fullmatch(r"([abAB])([12])-(\d+)", str(edge["linkage"]))
        if match is None:
            raise SetupError(f"Unsupported glycosidic linkage: {edge['linkage']}")
        _, donor, acceptor = match.groups()
        d, a = int(donor), int(acceptor)
        r_parent, r_child = node_resids[parent], node_resids[child]
        ring_oxygen = data["pucker_ring_atoms"][child][0]
        suffix = f"{parent}_{child}"

        phi, psi = f"phi{suffix}", f"psi{suffix}"
        lines.extend([
            f"{phi}: TORSION ATOMS=@{ring_oxygen}-{r_child},@C{d}-{r_child},@O{a}-{r_parent},@C{a}-{r_parent}",
            f"{psi}: TORSION ATOMS=@C{d}-{r_child},@O{a}-{r_parent},@C{a}-{r_parent},@C{a-1}-{r_parent}",
        ])
        names.extend((phi, psi))

        if a == 8 and d == 2:
            # These atom definitions have been supplied and validated by the
            # user for Neu5Ac alpha2-8, not for other 2-8 sugar chemistries.
            if (data["monosaccharides"][parent] != "Neu5Ac"
                    or data["monosaccharides"][child] != "Neu5Ac"):
                raise SetupError(
                    f"No validated omega7/8/9 atom definitions for {edge['linkage']} "
                    f"between {data['monosaccharides'][child]} and "
                    f"{data['monosaccharides'][parent]} (edge {suffix})"
                )
            omega_atoms = {
                "omega7": ("C8", "C7", "C6", "O6"),
                "omega8": ("O8", "C8", "C7", "C6"),
                "omega9": ("O8", "C8", "C9", "O9"),
            }
            for prefix, atoms in omega_atoms.items():
                name = f"{prefix}_{suffix}"
                lines.append(f"{name}: TORSION ATOMS=" + ",".join(
                    f"@{atom}-{r_parent}" for atom in atoms
                ))
                names.append(name)
        elif a == 6:
            name = f"omega{suffix}"
            lines.append(
                f"{name}: TORSION ATOMS=@O6-{r_parent},@C6-{r_parent},"
                f"@C5-{r_parent},@C4-{r_parent}"
            )
            names.append(name)

    for node, atoms in data["pucker_ring_atoms"].items():
        if len(atoms) != 6:
            raise SetupError(f"Ring of node {node} has {len(atoms)} atoms, expected 6")
        residue = node_resids[str(node)]
        lines.append(f"puck{node}: PUCKERING ATOMS=" + ",".join(
            f"@{atom}-{residue}" for atom in atoms
        ))
    return lines, names


def _peptide_torsions(data: dict):
    peptide = data.get("peptide")
    if peptide in (None, "", "none"):
        return [], []
    if len(peptide) != 3 or peptide[1] not in "NST":
        raise SetupError(f"Cannot define glycopeptide torsions automatically for {peptide!r}; expected tripeptide with central N/S/T")
    root = min(data["pdb_residues"], key=lambda n: int(data["pdb_residues"][n]))
    root_id = int(data["pdb_residues"][root])
    sugar_donor = "C2" if data["monosaccharides"][root] in SIALIC_NAMES else "C1"
    oxygen = data["pucker_ring_atoms"][root][0]
    # A tripeptide in the CHARMM-GUI PDB precedes glycan residues, as in supplied examples.
    # A validated topology-linked protein residue identifier should be added for other organizations.
    pep_res = root_id - 2
    if pep_res < 1:
        raise SetupError("Cannot infer peptide attachment residue index from JSON PDB ordering")
    if peptide[1] == "N":
        attachment, next_atom, next_next = "ND2", "CG", "CB"
        spec = [
            f"phi_gb: TORSION ATOMS=@{oxygen}-{root_id},@{sugar_donor}-{root_id},@{attachment}-{pep_res},@{next_atom}-{pep_res}",
            f"psi_gb: TORSION ATOMS=@{sugar_donor}-{root_id},@{attachment}-{pep_res},@{next_atom}-{pep_res},@{next_next}-{pep_res}",
            f"chi2: TORSION ATOMS=@ND2-{pep_res},@CG-{pep_res},@CB-{pep_res},@CA-{pep_res}",
            f"chi1: TORSION ATOMS=@CG-{pep_res},@CB-{pep_res},@CA-{pep_res},@N-{pep_res}",
        ]
        return spec, ["phi_gb", "psi_gb", "chi1", "chi2"]
    attachment = "OG" if peptide[1] == "S" else "OG1"
    return [
        f"phi_gb: TORSION ATOMS=@{oxygen}-{root_id},@{sugar_donor}-{root_id},@{attachment}-{pep_res},@CB-{pep_res}",
        f"psi_gb: TORSION ATOMS=@{sugar_donor}-{root_id},@{attachment}-{pep_res},@CB-{pep_res},@CA-{pep_res}",
        f"chi1: TORSION ATOMS=@{attachment}-{pep_res},@CB-{pep_res},@CA-{pep_res},@N-{pep_res}",
    ], ["phi_gb", "psi_gb", "chi1"]


def _write_plumed(folder: Path, data: dict):
    pucker_and_glycan, biased = _torsion_definitions(data)
    peptide_def, peptide_bias = _peptide_torsions(data)
    sim = data["simulation"]
    factors = sim.get("rect_factors", REST_RECT_FACTORS)
    if len(factors) != int(sim["n_replica"]):
        raise SetupError("simulation.rect_factors must contain one RECT factor per replica; default list has 12")
    temp = float(sim["min_temperature"])
    # NOTE: HILLS / kinetics reflect these choices; adjustable only via JSON.
    metad = sim.get("metad", {})
    pace, tau, sigma = int(metad.get("pace", 500)), float(metad.get("tau", 4.0)), float(metad.get("sigma", 0.35))
    lines = ["MOLINFO STRUCTURE=reference.pdb", "#"] + peptide_def + pucker_and_glycan + ["#"]
    for name in [*peptide_bias, *biased]:
        lines.append(f"metaD_{name}: METAD ARG={name} TAU={tau} PACE={pace} SIGMA={sigma} "
                     f"GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200 TEMP={temp} "
                     f"RECT={','.join(str(f) for f in factors)} FILE=HILLS_{name}")
    lines += [f"PRINT ARG={','.join([*peptide_bias, *biased])} STRIDE=5000 FILE=COLVAR",
              "PRINT ARG=*.theta,*.phi STRIDE=5000 FILE=COLVAR_theta",
              "PRINT ARG=*.bias STRIDE=5000 FILE=COLVAR_metaD"]
    (folder / "plumed_allangles.dat").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _solute_names(data: dict):
    """Use actual CHARMM/PDB residue names in JSON, not separate YAML monosaccharide_types."""
    names = set()
    for node, sugar in data["monosaccharides"].items():
        field = data.get("forcefield_residue_names", {}).get(node, {}).get("charmm", {})
        if isinstance(field, dict):
            names.update(v for v in field.values() if v)
        elif field:
            names.add(field)
        # User's PDB residue names can be truncated; include them for diagnostics,
        # but matching processed.top is based on actual column-four residue names.
    return names


def _define_solute(processed: Path, output: Path, data: dict):
    target = _solute_names(data)
    if not target:
        raise SetupError("JSON has no CHARMM residue names for solute selection")
    lines = processed.read_text().splitlines(keepends=True)
    out = []
    found = set()
    in_atoms = False
    for line in lines:
        stripped = line.strip()
        if re.match(r"^\s*\[\s*atoms\s*\]", line):
            in_atoms = True
        elif re.match(r"^\s*\[", line):
            in_atoms = False
        if in_atoms and stripped and not stripped.startswith(";"):
            base, sep, comment = line.partition(";")
            parts = base.split()
            if len(parts) >= 4 and parts[3] in target:
                parts[1] = parts[1].rstrip("_") + "_"
                found.add(parts[3])
                line = " ".join(parts) + (" ;" + comment if sep else "\n")
        out.append(line)
    if not found:
        raise SetupError(f"No JSON-derived CHARMM residues found in [ atoms ] of {processed}; names={sorted(target)}")
    output.write_text("".join(out))
    return found


def _organize_rest(folder: Path, rest: Path):
    rest.mkdir(exist_ok=True)
    for name in ("generate_HREMD.sh", "partial_tempering_pairs.sh", "job.sh", "mdout.mdp",
                 "production.mdp", "plumed_allangles.dat", "reference.pdb", "processed.top",
                 "processed_m.top", "topol.tpr"):
        src = folder / name
        if src.is_file():
            destination = rest / name
            if destination.exists():
                raise SetupError(f"Refusing to overwrite existing REST setup file: {destination}")
            src.replace(destination)
    for topo in folder.glob("TOPO[0-9]*"):
        if topo.is_dir():
            dst = rest / topo.name
            if dst.exists():
                raise SetupError(f"Refusing to overwrite existing replica directory: {dst}")
            shutil.move(str(topo), str(dst))


def _setup_rest(folder: Path, data: dict, templates: Path, log: Path,
                gmx: str, module: str | None):
    min_equ = folder / "01_min_equ"
    rest = folder / "02_REST-RECT"
    rest.mkdir(exist_ok=True)
    sim = data["simulation"]
    # Build in the simulation root because generate_HREMD.sh expects relative paths.
    _copy_templates(templates, folder, list(sim["template_files"]),
                    required=("generate_HREMD.sh", "job.sh"))
    for extra in ("partial_tempering_pairs.sh",):
        src = templates / extra
        if src.exists():
            shutil.copy2(src, folder / extra)
    structure = min_equ / sim["structure_file"]
    index_file = min_equ / "index.ndx"
    topol_file = min_equ / "topol.top"
    for input_file in (structure, index_file, topol_file):
        if not _nonempty(input_file):
            raise FileNotFoundError(input_file)
    _cmd([gmx, "editconf", "-f", str(structure), "-o", "reference.pdb", "-resnr", "1"],
         folder, log, module=module)
    _write_plumed(folder, data)
    _validate_plumed_atoms(folder / "plumed_allangles.dat", folder / "reference.pdb")
    _cmd([gmx, "grompp", "-f", str(folder / sim["mdp_file"]), "-c", str(structure),
          "-pp", "processed.top", "-p", str(topol_file), "-n", str(index_file)],
         folder, log, module=module)
    found = _define_solute(folder / "processed.top", folder / "processed_m.top", data)
    _log(log, f"Solute residue names in processed topology: {sorted(found)}")
    _cmd(["bash", "./generate_HREMD.sh", str(structure), str(sim["n_replica"]),
          str(sim["min_temperature"]), str(sim["max_temperature"]), sim["mdp_file"], str(index_file)],
         folder, log)
    _organize_rest(folder, rest)
    if not _rest_ready(folder, data):
        raise SetupError(f"REST-RECT setup incomplete in {rest}; inspect {log}")

def init_glycan_simulation(
    folder_path: str | Path,
    glycan_id: str,
    n_replica: int = 12,
    min_temperature: float = 310.15,
    max_temperature: float = 800.0,
    structure_file: str = "step5_production.gro",
    mdp_file: str = "production.mdp",
    template_files: list[str] | None = None,
    forcefield_file: str = "toppar/forcefield.itp",
) -> Path:
    """
    Add or update simulation configuration in an existing glycan JSON file.

    The configuration is stored under the 'simulation' key.
    All other existing JSON entries are preserved.

    Parameters
    ----------
    glycan_id : str
        Glycan identifier, e.g. 'Man5'.
    folder_path : str | Path
        Library root or the glycan directory containing <glycan_id>.json.
    n_replica : int
        Number of REST-RECT replicas.
    min_temperature, max_temperature : float
        Minimum and maximum replica temperatures in K.
    structure_file, mdp_file : str
        Simulation input filenames.
    template_files : list[str] | None
        Template files to include.
    forcefield_file : str
        Relative path to the force-field file.

    Returns
    -------
    Path
        Path to the updated JSON file.
    """

    json_path = find_glycan_directory(folder_path, glycan_id) / f"{glycan_id}.json"

    if not json_path.is_file():
        raise FileNotFoundError(
            f"Glycan JSON file not found: {json_path}"
        )

    if template_files is None:
        template_files = [
            "*.mdp",
            "generate_HREMD.sh",
            "job.sh",
        ]

    # Read existing JSON
    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    # Simulation configuration
    simulation = {
        "n_replica": n_replica,
        "min_temperature": min_temperature,
        "max_temperature": max_temperature,
        "structure_file": structure_file,
        "mdp_file": mdp_file,
        "template_files": template_files,
        "forcefield_file": forcefield_file,
    }

    # Preserve existing simulation metadata not listed above
    data.setdefault("simulation", {}).update(simulation)

    # Write updated JSON
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)
        f.write("\n")

    print(f"Updated simulation configuration: {json_path}")

    return json_path

def generate_glycan_simulation(data_dir: str, glycan_id: str, template_dir: str | Path = "templates",
           *, gmx: str = "gmx", gmx_mpi: str = "gmx_mpi",
           gromacs_module: str | None = None) -> Path:
    
    """Prepare min/equ and REST-RECT for the glycan found recursively under data_dir."""
    if gromacs_module is None:
        gromacs_module = "chem/GROMACS"
    
    folder, data = _load(glycan_id, data_dir)
    templates = Path(template_dir).expanduser().resolve()
    if not templates.is_dir():
        raise FileNotFoundError(templates)
    log = folder / "rest_rect_setup.log"
    _log(log, f"START setup: {glycan_id}")
    _log(log, f"Metadata JSON: {folder / (glycan_id + '.json')}")
    states = status(glycan_id, data_dir)
    if not states["min_equ"]:
        gui = _charmm_gromacs(folder)
        if gui is None and not (folder / "01_min_equ" / "step3_input.gro").exists():
            raise SetupError(f"Cannot find CHARMM-GUI GROMACS folder for {glycan_id}")
        _minimize(folder, gui or (folder / "01_min_equ"), templates, data, log,
                  gmx, gmx_mpi, gromacs_module)
    if not _min_done(folder):
        raise SetupError("Minimization/equilibration did not finish")
    if not _rest_ready(folder, data):
        # Avoid destroying replicas already prepared by a previous interrupted setup.
        replicas = list((folder / "02_REST-RECT").glob("TOPO*"))
        if replicas:
            missing = _plumed_missing_definitions(folder / "02_REST-RECT" / "plumed_allangles.dat", data)
            if missing:
                raise SetupError(
                    "Existing REST replica directories use an incomplete PLUMED input "
                    f"({len(missing)} missing entries; examples: {missing[:6]}). "
                    "For simulations NOT yet started, regenerate plumed_allangles.dat "
                    "using _write_plumed(folder, data), then replace the existing file "
                    "in 02_REST-RECT. Already-running or completed simulations cannot "
                    "be retroactively fixed; a new production run is required. "
                    "No replica directories were deleted."
                )
            raise SetupError("Partial REST setup with existing TOPO replica directories detected. "
                             "Inspect them and restore/finish setup manually; no automatic deletion.")
        _setup_rest(folder, data, templates, log, gmx, gromacs_module)
    else:
        _log(log, "REST-RECT setup already complete; skipping")
    _log(log, "✓ Setup ready. Transfer 02_REST-RECT to HPC, execute job.sh, and download outputs before analyze_glycan().")
    return folder / "02_REST-RECT"


def _plumed_table(path: Path):
    import pandas as pd
    with path.open(encoding="utf-8") as f:
        first = next((line for line in f if line.strip().startswith("#! FIELDS")), None)
    if first:
        fields = first.split()[2:]
        return pd.read_csv(path, sep=r"\s+", comment="#", names=fields)
    return pd.read_csv(path, sep=r"\s+", comment="#")


def _colvar_analysis(topo: Path, out: Path):
    colvar = _plumed_table(topo / "COLVAR.0")
    if "time" not in colvar.columns:
        raise SetupError("COLVAR.0 missing time column")
    torsions = [str(c) for c in colvar.columns if re.fullmatch(r"(?:phi|psi)\d+_\d+|omega(?:7|8|9)?_?\d+_\d+", str(c))]
    if not torsions:
        raise SetupError("No named glycosidic dihedrals found in COLVAR.0")
    colvar[["time", *torsions]].to_csv(out / "colvar.dat", sep=" ", index=False)
    (out / "angles.dat").write_text(" ".join(torsions) + "\n")
    (out / "omega_angles.dat").write_text(" ".join(c for c in torsions if c.startswith("omega")) + "\n")
    pucker = topo / "COLVAR_theta.0"
    if pucker.is_file():
        df = _plumed_table(pucker)
        cols = [c for c in df if "theta" in str(c) or "phi" in str(c)]
        df[["time", *cols]].to_csv(out / "pucker.dat", sep=" ", index=False)


def _align_trajectory(topo: Path, out: Path, data: dict, log: Path,
                      gmx: str, module: str | None):
    # Requires GROMACS output plus MDAnalysis; avoid assumptions about peptide residues.
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
    import MDAnalysis.transformations as trans
    pdb = topo / "reference_rest-rect.pdb"
    traj = topo / "conformations_rest-rect.xtc"
    target = out / "conformations_rest-rect_aligned.xtc"
    if not _nonempty(pdb):
        _cmd([gmx, "editconf", "-f", "reference.pdb", "-o", str(pdb)],
             topo, log, module=module)
    if not _nonempty(traj):
        _cmd([gmx, "trjconv", "-f", "traj_comp.xtc", "-s", "topol.tpr",
              "-o", str(traj), "-skip", "5"], topo, log, module=module, input_text="0\n")
    u = mda.Universe(str(pdb), str(traj), guess_bonds=True)
    selection = "protein" if data.get("peptide") else "all"
    ref_atoms = u.select_atoms(selection)
    if len(ref_atoms) < 3:
        selection = "all"
    try:
        u.trajectory.add_transformations(trans.unwrap(u.atoms),
                                         trans.center_in_box(u.select_atoms(selection), wrap=False))
    except (ValueError, AttributeError) as exc:
        _log(log, f"Skipping PBC unwrap: {exc}")
    align.AlignTraj(u, u, select=selection, in_memory=True).run()
    with mda.Writer(str(target), n_atoms=u.atoms.n_atoms) as writer:
        for _ in u.trajectory:
            writer.write(u.atoms)
    shutil.copy2(pdb, out / "reference_rest-rect.pdb")


def _fes(topo: Path, out: Path, log: Path, *, plumed_module: str | None, temperature: float):
    import numpy as np
    hills_files = sorted(topo.glob("HILLS_*"))
    if not hills_files:
        raise SetupError(f"No HILLS_* files found in {topo}")
    kbt = 0.008314462618 * temperature
    for hills in hills_files:
        suffix = hills.name.removeprefix("HILLS_").replace(".", "_")
        target = out / f"fes_{suffix}.dat"
        if _nonempty(target):
            _log(log, f"Skipping existing FES {target.name}")
            continue
        fes_dir = out / f"FES_{suffix}_1D"
        fes_dir.mkdir(parents=True, exist_ok=True)
        # Generate directly into a dedicated directory, not TOPO0.
        _cmd(["plumed", "sum_hills", "--hills", str(hills.resolve()), "--stride", "1000",
              "--kt", str(kbt), "--bin", "99", "--min", "0", "--max", "0.5"],
             fes_dir, log, module=plumed_module)
        profiles = fes_dir / "profiles" / "MIN2ZERO"
        profiles.mkdir(parents=True, exist_ok=True)
        fes_files = sorted(fes_dir.glob("fes_*.dat"),
                           key=lambda p: int(p.stem.split("_")[-1])
                           if p.stem.split("_")[-1].isdigit() else -1)
        if not fes_files:
            raise SetupError(f"PLUMED sum_hills returned no FES for {hills}")
        # Shift every FES using global min of last time-window FES (original convention).
        def read_profile(path):
            return np.loadtxt(path, comments="#")
        last = read_profile(fes_files[-1])
        offset = float(np.min(last[:, 1]))
        for i, f in enumerate(fes_files):
            arr = read_profile(f)
            arr[:, 1] -= offset
            np.savetxt(profiles / f"fes_{i:05d}.dat", arr, fmt="%.9f")
        shutil.copy2(profiles / f"fes_{len(fes_files)-1:05d}.dat", target)


def _exchange_stats(folder: Path, topo: Path, log: Path):
    mdlog = topo / "md.log"
    lines = mdlog.read_text(errors="replace").splitlines()
    index = next((i for i, line in enumerate(lines) if "Replica exchange statistics" in line), None)
    if index is not None:
        text = "\n".join(lines[index:index + 20])
        # Repeated calls should not append the same block.
        if text not in log.read_text(errors="replace"):
            _log(log, text)
    demux = Path("helper/demux.pl")
    if demux.is_file():
        _cmd([str(demux.resolve()), str(mdlog.resolve())], folder, log)
        for name in ("replica_index.xvg", "replica_temp.xvg"):
            src = folder / name
            if src.is_file():
                src.replace(folder / "02_REST-RECT" / name)


def analyze_glycan(data_dir: str, glycan_id: str, *,
                   gmx: str = "gmx", gromacs_module: str | None = None,
                   plumed_module: str | None = None) -> Path:
    if gromacs_module is None:
        gromacs_module = "chem/GROMACS"
    if plumed_module is None:
        plumed_module = "chem/PLUMED"
    
    """Analyze completed, downloaded REST-RECT output, preserving completed artifacts."""
    folder, data = _load(glycan_id, data_dir)
    log = folder / "rest_rect_setup.log"
    if not _production_done(folder, data):
        raise SetupError(f"Production outputs for {glycan_id} are incomplete. Download HPC results first.")
    out = folder / "03_lib"
    out.mkdir(exist_ok=True)
    topo = folder / "02_REST-RECT" / "TOPO0"
    if _analysis_done(folder, data):
        _log(log, f"Analysis already complete: {glycan_id}")
        return out
    _log(log, f"Starting analysis: {glycan_id}")
    _exchange_stats(folder, topo, log)
    if not (_nonempty(out / "reference_rest-rect.pdb") and _nonempty(out / "conformations_rest-rect_aligned.xtc")):
        _align_trajectory(topo, out, data, log, gmx, gromacs_module)
    if not _nonempty(out / "colvar.dat"):
        _colvar_analysis(topo, out)
    _fes(topo, out, log, plumed_module=plumed_module,
         temperature=float(data["simulation"]["min_temperature"]))
    _log(log, f"Analysis {'complete' if _analysis_done(folder, data) else 'partially complete'}: {glycan_id}")
    return out


__all__ = ["SetupError", "find_glycan_directory", "status", "init_glycan_simulation", "generate_glycan_simulation", "analyze_glycan"]
