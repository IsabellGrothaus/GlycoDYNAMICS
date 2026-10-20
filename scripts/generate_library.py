# do not copy FES files!!! Adapt this in script!

from pathlib import Path
import shutil

def copy_glycan_library(
    library_path: str | Path,
    data_path: str | Path = "data",
    overwrite: bool = True,
) -> None:
    """
    Recursively find all '03_gl_data' directories under data_path and copy
    their contents to the GlycoLIBRARY/library directory, preserving the
    original directory hierarchy but removing '03_gl_data' from the path.

    Example
    -------
    Source:
        data/N-glycans/Complex/A2/03_gl_data/fes_phi1_2.dat

    Destination:
        GlycoLIBRARY/library/N-glycans/Complex/A2/fes_phi1_2.dat

    Parameters
    ----------
    library_path : str | Path
        Absolute path to GlycoLIBRARY (the 'library' subdirectory is appended).
    data_path : str | Path
        Source data directory, relative to the notebook or absolute.
    overwrite : bool
        Whether to overwrite existing destination files.
    """

    data_path = Path(data_path).resolve()
    library_dir = Path(library_path).expanduser().resolve() / "library"

    if not data_path.is_dir():
        raise FileNotFoundError(f"Data directory not found: {data_path}")

    # Find all 03_gl_data directories recursively
    gl_dirs = sorted(
        path for path in data_path.rglob("03_gl_data")
        if path.is_dir()
    )

    if not gl_dirs:
        print("No '03_gl_data' directories found.")
        return

    copied = 0
    skipped = 0

    for gl_dir in gl_dirs:

        # Preserve hierarchy relative to data/, excluding 03_gl_data
        relative_path = gl_dir.parent.relative_to(data_path)
        destination_dir = library_dir / relative_path
        destination_dir.mkdir(parents=True, exist_ok=True)

        # Copy files, including any nested subdirectories
        for source_file in gl_dir.rglob("*"):
            if not source_file.is_file():
                continue

            destination_file = destination_dir / source_file.relative_to(gl_dir)
            destination_file.parent.mkdir(parents=True, exist_ok=True)

            if destination_file.exists() and not overwrite:
                skipped += 1
                continue

            shutil.copy2(source_file, destination_file)
            copied += 1

        print(f"✓ {relative_path}: copied to {destination_dir}")

    print(
        f"\nFinished: {copied} files copied, "
        f"{skipped} skipped, "
        f"{len(gl_dirs)} glycan directories processed."
    )