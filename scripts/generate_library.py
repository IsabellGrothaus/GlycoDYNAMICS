"""Copy selected GlycoGEN results to GlycoLIBRARY and index the library.

Layout
------
DATABASE/
    glycan.json                 # automatically updated library index
    library/
        N-glycans/Complex/A2/
            A2.json, A2.png, fes_*.dat, colvar.dat, ...

Source results are expected in data/.../<glycan_id>/03_lib/.
Only selected files are copied, all directly into the destination glycan folder.
The glycan directory hierarchy relative to ``data_path`` is retained.
"""

from __future__ import annotations

import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# Files sourced exclusively from each glycan's 03_lib directory.
RESULT_NAMES = frozenset({
    "colvar.dat", "angles.dat", "pucker.dat", "omega_angles.dat",
})
RESULT_PATTERNS = ("fes_*.dat", "*.xtc", "*.pdb")


def _library_dir(database_path: str | Path) -> Path:
    """Accept the DATABASE root or its ``library`` folder."""
    root = Path(database_path).expanduser().resolve()
    return root if root.name == "library" else root / "library"


def _selected_results(folder: Path) -> list[Path]:
    """Select only files immediately inside 03_lib (no intermediate folders)."""
    found = {p for name in RESULT_NAMES for p in [folder / name] if p.is_file()}
    for pattern in RESULT_PATTERNS:
        found.update(p for p in folder.glob(pattern) if p.is_file())
    return sorted(found, key=lambda p: p.name)


def _atomic_write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    try:
        temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


TORSION_RE = re.compile(r"^(phi|psi|omega(?:7|8|9)?)(?:_)?(\d+)_(\d+)$", re.I)


def _colvar_columns(path: Path) -> tuple[list[str], int | None, list[str]]:
    """Read PLUMED #! FIELDS or plain column headers without loading trajectories."""
    warnings: list[str] = []
    if not path.is_file():
        return [], None, ["colvar.dat missing"]
    with path.open(encoding="utf-8", errors="replace") as handle:
        first = next((line.strip() for line in handle if line.strip()), "")
    if not first:
        return [], 0, ["colvar.dat is empty"]
    if first.startswith("#! FIELDS"):
        columns = first.split()[2:]
        # PLUMED data rows are left untouched; count with streaming iteration.
        with path.open(encoding="utf-8", errors="replace") as handle:
            count = sum(bool(line.strip()) and not line.lstrip().startswith("#") for line in handle)
    elif first.startswith("#"):
        return [], None, ["COLVAR header not recognized; expected '#! FIELDS' or a plain header"]
    else:
        columns = first.split()
        with path.open(encoding="utf-8", errors="replace") as handle:
            count = sum(bool(line.strip()) and not line.lstrip().startswith("#") for line in handle) - 1
    if len(set(columns)) != len(columns):
        warnings.append("Duplicate column names in colvar.dat")
    if "time" not in columns:
        warnings.append("colvar.dat has no time column")
    return columns, count, warnings


def register_glycan(glycan_directory: str | Path) -> dict[str, Any]:
    """Register actual library files into glycan_id.json, without changing graph IDs.

    A torsion named phi1_2 maps to fes_phi1_2_0.dat (likewise for omega).

    No canonicalization or simulation is performed. Existing graph ordering is used.
    Missing files produce explicit warnings, never fabricated successful validation.
    """
    folder = Path(glycan_directory).expanduser().resolve()
    glycan_id = folder.name
    json_path = folder / f"{glycan_id}.json"
    if not json_path.is_file():
        raise FileNotFoundError(json_path)
    data = json.loads(json_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Metadata must be a JSON object: {json_path}")
    edges = data.get("graph_edges", [])
    if not isinstance(edges, list):
        raise ValueError(f"graph_edges must be a list in {json_path}")
    edge_map = {}
    warnings = []
    for edge in edges:
        parent, child = str(edge["parent"]), str(edge["child"])
        edge_id = f"{parent}_{child}"
        if "id" in edge and str(edge["id"]) != edge_id:
            warnings.append(f"Edge ID mismatch: {edge['id']} != {edge_id}")
        if edge_id in edge_map:
            warnings.append(f"Duplicate edge: {edge_id}")
        edge_map[edge_id] = edge

    columns, frame_count, parse_warnings = _colvar_columns(folder / "colvar.dat")
    warnings.extend(parse_warnings)
    actual = {p.name for p in folder.iterdir() if p.is_file()}
    torsions = []
    for name in columns:
        match = TORSION_RE.fullmatch(name)
        if not match:
            continue
        torsion_type, parent, child = match.groups()
        edge_id = f"{parent}_{child}"
        if edge_id not in edge_map:
            warnings.append(f"COLVAR torsion {name} does not match any graph edge")
            continue
        # Production FES files carry a replica/index suffix (_0). Preserve it.
        fes = f"fes_{name}_0.dat"
        if fes not in actual:
            warnings.append(f"Missing FES for {name}: {fes}")
        torsions.append({
            "id": len(torsions), "name": name, "type": torsion_type.lower(),
            "edge_id": edge_id, "colvar_column": name,
            "fes_file": fes if fes in actual else None,
        })

    torsion_names = {item["name"] for item in torsions}
    for edge_id, edge in edge_map.items():
        for prefix in ("phi", "psi"):
            required = f"{prefix}{edge_id}"
            if required not in torsion_names and columns:
                warnings.append(f"Missing COLVAR torsion: {required}")
        linkage = str(edge.get("linkage", ""))
        if linkage.endswith("-6") and not any(
            t["edge_id"] == edge_id and t["type"] == "omega" for t in torsions
        ) and columns:
            warnings.append(f"Missing omega torsion for 1-6 linkage {edge_id}")
        if linkage.lower().endswith("2-8") and columns:
            for omega in ("omega7", "omega8", "omega9"):
                if not any(t["edge_id"] == edge_id and t["type"] == omega for t in torsions):
                    warnings.append(f"Missing {omega} for 2-8 linkage {edge_id}")

    expected_torsions = None
    if columns:
        expected_torsions = 2 * len(edge_map) + sum(
            1 if str(e.get("linkage", "")).endswith("-6") else
            3 if str(e.get("linkage", "")).lower().endswith("2-8") else 0
            for e in edges
        )
    matching = (bool(columns) and len(torsions) == expected_torsions
                and not any("Missing COLVAR torsion" in w or "Missing omega" in w or "does not match" in w
                            for w in warnings))
    all_fes = bool(torsions) and all(t["fes_file"] for t in torsions)
    xtcs = sorted(n for n in actual if n.endswith(".xtc"))
    pdbs = sorted(n for n in actual if n.endswith(".pdb"))
    if not xtcs:
        warnings.append("No XTC trajectory found")
    if not pdbs:
        warnings.append("No PDB reference found")
    if not (folder / "pucker.dat").is_file():
        warnings.append("Missing pucker.dat")
    if not (folder / "angles.dat").is_file():
        warnings.append("Missing angles.dat")
    if not (folder / "omega_angles.dat").is_file():
        warnings.append("Missing omega_angles.dat")
    if not columns and not edge_map:
        status = "metadata_only"
    elif not columns:
        status = "awaiting_analysis"
    elif (matching and all_fes and xtcs and pdbs and not parse_warnings
          and all(n in actual for n in ("angles.dat", "omega_angles.dat", "pucker.dat"))):
        status = "complete"
    else:
        status = "incomplete"
    data["analysis"] = {
        "schema_version": 1,
        "status": status,
        "files": {
            "colvar": "colvar.dat" if "colvar.dat" in actual else None,
            "angles": "angles.dat" if "angles.dat" in actual else None,
            "omega_angles": "omega_angles.dat" if "omega_angles.dat" in actual else None,
            "pucker": "pucker.dat" if "pucker.dat" in actual else None,
            "trajectories": xtcs,
            "references": pdbs,
        },
        "colvar_rows": frame_count,
        "dihedrals": torsions,
        "validation": {
            "colvar_columns_match": matching if columns else None,
            "fes_mapping_complete": all_fes if columns else None,
            "trajectory_frames_match": None,
            "conformer_states_match": None,
        },
        "warnings": sorted(set(warnings)),
    }
    _atomic_write_json(json_path, data)
    return data["analysis"]


def generate_glycan_json(database_path: str | Path) -> Path:
    """Rebuild DATABASE/glycan.json from JSON files physically present in library.

    Records glycan IDs, IUPAC sequences, relative paths, category, linkage count,
    reconstruction status, and completeness indicators. A file that cannot be
    parsed is listed under ``warnings`` rather than being silently ignored.

    This indexes the *current* library, so existing items remain in the index
    even if they weren't part of the latest copy invocation. No files are deleted.
    """
    library_dir = _library_dir(database_path)
    library_dir.mkdir(parents=True, exist_ok=True)
    entries: list[dict[str, Any]] = []
    issues: list[str] = []
    seen: dict[str, str] = {}

    for path in sorted(library_dir.rglob("*.json")):
        # Only <glycan_id>/<glycan_id>.json identifies a glycan; ignore auxiliary JSON.
        if path.parent.name != path.stem:
            continue
        relative = path.parent.relative_to(library_dir)
        glycan_id = path.stem
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(data, dict):
                raise ValueError("metadata JSON must be an object")
        except (OSError, ValueError, UnicodeError) as exc:
            issues.append(f"{relative.as_posix()}: cannot read {path.name}: {exc}")
            continue

        if glycan_id in seen:
            issues.append(f"Duplicate glycan ID {glycan_id!r}: {seen[glycan_id]} and {relative.as_posix()}")
        seen[glycan_id] = relative.as_posix()

        if data.get("glycan_id") not in (None, glycan_id):
            issues.append(f"{relative.as_posix()}: JSON glycan_id mismatches directory name")

        files = sorted(p.name for p in path.parent.iterdir() if p.is_file())
        reconstruction = data.get("reconstruction") or {}
        if not isinstance(reconstruction, dict):
            reconstruction = {}
        required = [f"{glycan_id}.json", f"{glycan_id}.png", *sorted(RESULT_NAMES)]
        # FES / trajectory / PDB are variable-size collections, so report their presence.
        missing = [name for name in required if name not in files]
        if not any(name.startswith("fes_") and name.endswith(".dat") for name in files):
            missing.append("fes_*.dat")
        if not any(name.endswith(".xtc") for name in files):
            missing.append("*.xtc")
        if not any(name.endswith(".pdb") for name in files):
            missing.append("*.pdb")

        edges = data.get("graph_edges")
        entry = {
            "glycan_id": glycan_id,
            "iupac": data.get("iupac"),
            "path": relative.as_posix(),
            "category": list(relative.parts[:-1]),
            "monosaccharide_count": len(data["monosaccharides"]) if isinstance(data.get("monosaccharides"), dict) else None,
            "linkage_count": len(edges) if isinstance(edges, list) else None,
            "reconstruction_status": reconstruction.get("status"),
            "user_checked": reconstruction.get("user_checked"),
            "files": files,
            "missing_files": missing,
            "analysis_status": (data.get("analysis") or {}).get("status"),
            "dihedral_count": len((data.get("analysis") or {}).get("dihedrals", [])),
            "analysis_warnings": (data.get("analysis") or {}).get("warnings", []),
        }
        if not entry["iupac"]:
            issues.append(f"{relative.as_posix()}: missing IUPAC")
        entries.append(entry)

    # Ensure deterministic ordering and let consumers access IDs/IUPAC directly.
    entries.sort(key=lambda item: (item["path"].casefold(), item["glycan_id"]))
    index = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "glycan_count": len(entries),
        "glycan_ids": [entry["glycan_id"] for entry in entries],
        "iupacs": {entry["glycan_id"]: entry["iupac"] for entry in entries},
        "glycans": entries,
        "warnings": issues,
    }
    index_path = library_dir.parent / "glycan.json"
    _atomic_write_json(index_path, index)
    print(f"Updated index: {index_path} ({len(entries)} glycans; {len(issues)} warnings)")
    return index_path


def display_glycan_library(database_path: str | Path, *, show_files: bool = False) -> dict[str, Any]:
    """Read and display DATABASE/glycan.json; return the parsed index."""
    index_path = _library_dir(database_path).parent / "glycan.json"
    if not index_path.is_file():
        raise FileNotFoundError(f"Library index not found: {index_path}. Run generate_glycan_json() first.")
    index = json.loads(index_path.read_text(encoding="utf-8"))
    print(f"GlycoLIBRARY: {index.get('glycan_count', 0)} glycans")
    print(f"Index: {index_path}")
    print(f"Updated: {index.get('generated_at_utc', 'unknown')}")
    for entry in index.get("glycans", []):
        print(f"\n{entry['glycan_id']}  [{entry['path']}]")
        print(f"  IUPAC: {entry.get('iupac') or 'not available'}")
        print(f"  Reconstruction: {entry.get('reconstruction_status') or 'not recorded'}"
              f"; user checked: {entry.get('user_checked')}")
        print(f"  Analysis: {entry.get('analysis_status') or 'not registered'}; "
              f"dihedrals: {entry.get('dihedral_count', 0)}")
        if entry.get("missing_files"):
            print(f"  Missing: {', '.join(entry['missing_files'])}")
        if show_files:
            print(f"  Files: {', '.join(entry.get('files', []))}")
    for warning in index.get("warnings", []):
        print(f"WARNING: {warning}")
    return index


def copy_glycan_library(
    library_path: str | Path,
    data_path: str | Path = "data",
    overwrite: bool = True,
) -> None:
    """Copy metadata and available 03_lib/03_gl_data results; refresh glycan.json.

    ``library_path`` is DATABASE (or DATABASE/library).
    Preserves folder hierarchy, e.g. data/N-glycans/Complex/A2/03_lib
    -> DATABASE/library/N-glycans/Complex/A2/<selected files>.

    Glycan folders are discovered using <glycan_id>.json. Matching results
    are copied from the immediate 03_lib or legacy 03_gl_data directory, if any.
    Never moves or deletes source files, and does not delete old library files.
    """
    source_root = Path(data_path).expanduser().resolve()
    library_dir = _library_dir(library_path)
    if not source_root.is_dir():
        raise FileNotFoundError(f"Data directory not found: {source_root}")
    if library_dir == source_root or source_root in library_dir.parents:
        raise ValueError("Destination library must be outside the source data directory")

    library_dir.mkdir(parents=True, exist_ok=True)
    # Discover glycan folders by their own metadata, rather than requiring analysis.
    # This allows inclusion of glycans before their simulations are finished.
    metadata_files = sorted(
        p for p in source_root.rglob("*.json")
        if p.is_file() and p.stem == p.parent.name
        and not any(part in {"01_min_equ", "02_REST-RECT", "03_lib", "03_gl_data"}
                    or "charmm-gui" in part.lower() for part in p.relative_to(source_root).parts[:-1])
    )
    copied = skipped = 0
    glycan_count = 0
    print(f"Source data directory: {source_root}")
    print(f"Destination library: {library_dir}")
    print(f"Glycan metadata files discovered: {len(metadata_files)}")
    try:
        if not metadata_files:
            print("WARNING: No <glycan_id>/<glycan_id>.json files found under the source directory.")
            print("Check that data_path points to the root of your GlycoGEN data library.")
        for json_file in metadata_files:
            glycan_folder = json_file.parent
            glycan_id = glycan_folder.name
            relative = glycan_folder.relative_to(source_root)
            destination = library_dir / relative
            sources = [json_file]
            png = glycan_folder / f"{glycan_id}.png"
            if png.is_file():
                sources.append(png)
            else:
                print(f"WARNING: missing {png}")

            # Prefer current results; also recognize the older directory name.
            results = glycan_folder / "03_lib"
            legacy = glycan_folder / "03_gl_data"
            if results.is_dir():
                sources.extend(_selected_results(results))
            elif legacy.is_dir():
                print(f"INFO: Using legacy results directory: {legacy}")
                sources.extend(_selected_results(legacy))
            else:
                print(f"INFO: No analysis results yet for {glycan_id}; copying metadata only")

            destination.mkdir(parents=True, exist_ok=True)
            unique_sources = sorted(set(sources), key=lambda p: p.name)
            for source in unique_sources:
                target = destination / source.name
                if target.exists() and not overwrite:
                    skipped += 1
                    continue
                shutil.copy2(source, target)
                copied += 1
            registration = register_glycan(destination)
            print(f"  Registration: {registration['status']}; {len(registration['dihedrals'])} torsions")
            glycan_count += 1
            print(f"✓ {relative.as_posix()}: {len(unique_sources)} eligible files")
    finally:
        # Re-register all existing destinations, including those not copied this run.
        for metadata in sorted(library_dir.rglob("*.json")):
            if metadata.parent.name == metadata.stem:
                try:
                    register_glycan(metadata.parent)
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    print(f"WARNING: could not register {metadata.parent}: {exc}")
        generate_glycan_json(library_dir)
    print(f"Finished: {copied} copied, {skipped} skipped, {glycan_count} glycan folders processed")


__all__ = ["copy_glycan_library", "generate_glycan_json", "display_glycan_library", "register_glycan"]
