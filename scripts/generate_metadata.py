"""IUPAC-condensed glycan metadata + optional PDB verification.

Supports bracket branches and common monosaccharides; rejects unsupported chemistry.
Force-field residue codes are names/candidates, NOT a guarantee of a valid topology.
"""
from __future__ import annotations
import json
import re
from collections import defaultdict
from pathlib import Path
import warnings
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.patches import (
    Circle,
    Rectangle,
    RegularPolygon,
    Polygon,
    FancyArrowPatch,
)

RING_PYR = ('O5','C1','C2','C3','C4','C5')
RING_SIA = ('O6','C2','C3','C4','C5','C6')
RING_FUR = ('O4','C1','C2','C3','C4')

# Registry: CHARMM-GUI CSML identifiers (full names), GLYCAM06 naming family
# reference: https://glycam.org/docs/othertoolsservice/2016/06/09/3d-snfg-list-of-residue-names/
# and https://www.charmm-gui.org/?doc=archive&lib=csml
# GLYCAM family letter is not an unambiguous residue code without acceptor positions.
REGISTRY = {
 'GlcNAc': dict(charmm={'a':'AGLCNA','b':'BGLCNA'}, glycam='Y', pdb=['NAG','NDG','SGN','4YS'], ring=RING_PYR),
 'GalNAc': dict(charmm={'a':'AGALNA','b':'BGALNA'}, glycam='V', pdb=['NGA'], ring=RING_PYR),
 'Man': dict(charmm={'a':'AMAN','b':'BMAN'}, glycam='M', pdb=['MAN','BMA'], ring=RING_PYR),
 'Gal': dict(charmm={'a':'AGAL','b':'BGAL'}, glycam='L', pdb=['GAL','GLA'], ring=RING_PYR),
 'Glc': dict(charmm={'a':'AGLC','b':'BGLC'}, glycam='G', pdb=['GLC','BGC','MAL'], ring=RING_PYR),
 'Fuc': dict(charmm={'a':'AFUC','b':'BFUC'}, glycam='F', pdb=['FUC','FUL'], ring=RING_PYR, default_L=True),
 'Xyl': dict(charmm={'a':'AXYL','b':'BXYL'}, glycam='X', pdb=['XYL','XYS','XYP'], ring=RING_PYR),
 'GlcA': dict(charmm={'a':'AGLCA','b':'BGLCA'}, glycam='Z', pdb=['GCU'], ring=RING_PYR),
 'GalA': dict(charmm={'a':'AGALA','b':'BGALA'}, glycam='O', pdb=['ADA'], ring=RING_PYR),
 'IdoA': dict(charmm={'a':'AIDOA','b':'BIDOA'}, glycam='U', pdb=['IDS'], ring=RING_PYR, default_L=True),
 'Rha': dict(charmm={'a':'ARHM','b':'BRHM'}, glycam='H', pdb=['RAM'], ring=RING_PYR, default_L=True),
 'Ara': dict(charmm={'a':'AARB','b':'BARB'}, glycam=None, pdb=['ARA','AHR'], ring=RING_PYR),
 'ManNAc': dict(charmm={}, glycam='W', pdb=[], ring=RING_PYR),
 'GlcN': dict(charmm={'b':'BGLCN'}, glycam=None, pdb=['GCS'], ring=RING_PYR),
 'GalN': dict(charmm={'a':'AGALN','b':'BGALN'}, glycam=None, pdb=[], ring=RING_PYR),
 'FucNAc': dict(charmm={'a':'AFUCNA','b':'BFUCNA'}, glycam=None, pdb=[], ring=RING_PYR, default_L=True),
 'Neu5Ac': dict(charmm={'a':'ANE5AC','b':'BNE5AC'}, glycam='S', pdb=['SIA'], ring=RING_SIA, donor=2),
 'Neu5Gc': dict(charmm={}, glycam=None, pdb=[], ring=RING_SIA, donor=2),
 'Neu': dict(charmm={}, glycam=None, pdb=['NEU'], ring=RING_SIA, donor=2),
 'Kdn': dict(charmm={'b':'BKDN'}, glycam=None, pdb=['KDN'], ring=RING_SIA, donor=2),
 'Kdo': dict(charmm={}, glycam=None, pdb=['KDO'], ring=RING_SIA, donor=2),
 'Fru': dict(charmm={'a':'AFRU','b':'BFRU'}, glycam=None, pdb=['FRU'], ring=None, donor=2),
}
# Additional GLYCAM mapping must be checked against the specific version of its prep files.
# Position code: occupied hydroxyl sites on the *acceptor* residue (up to 4,6 etc).
POSITION_CODES = {(): '0',(2,):'2',(3,):'3',(4,):'4',(6,):'6',
 (2,3):'Z',(2,4):'Y',(2,6):'X',(3,4):'W',(3,6):'V',(4,6):'U',
 (2,3,4):'T',(2,3,6):'S',(2,4,6):'R',(3,4,6):'Q',(2,3,4,6):'P'}
# Restricted because other sites and derivatives require additional residue-library checks.

NAME_RE = '|'.join(sorted(map(re.escape, REGISTRY), key=len, reverse=True))
TOKEN_RE = re.compile(r'\s*(\[|\]|(?:' + NAME_RE + r')(?:\(([ab?\u03b1\u03b2])([12])-(\d+|\?)\))?)')


def _parse_iupac(iupac):
    tokens, pos = [], 0
    while pos < len(iupac):
        if not iupac[pos:].strip(): break
        match = TOKEN_RE.match(iupac, pos)
        if not match:
            raise ValueError(f'Unsupported monosaccharide/linkage around {iupac[pos:pos+40]!r}')
        whole = match.group(1)
        if whole in ('[',']'): tokens.append(whole)
        else:
            name = whole.split('(')[0]
            tokens.append((name, (match.group(2), int(match.group(3)), match.group(4)) if match.group(2) else None))
        pos = match.end()
    nodes, edges, pending, stack = {}, [], [], []
    for tok in tokens:
        if tok == '[':
            stack.append(pending); pending = []
        elif tok == ']':
            if not stack or len(pending) != 1: raise ValueError('Invalid nested branch')
            pending = stack.pop() + pending
        else:
            name, incoming = tok
            nid = len(nodes)
            nodes[nid] = name
            for child, linkage in pending:
                if linkage is None: raise ValueError('Non-reducing residue missing linkage')
                a, donor, acceptor = linkage
                required_donor = REGISTRY[nodes[child]].get('donor', 1)
                if donor != required_donor:
                    raise ValueError(f'{nodes[child]} uses donor C{required_donor}, not C{donor}')
                edges.append((nid, child, f'{a.lower().replace("α","a").replace("β","b")}{donor}-{acceptor}'))
            pending = [(nid, incoming)]
    if stack or len(pending) != 1 or pending[0][1] is not None:
        raise ValueError('Unbalanced branches or unsupported reducing-end linkage')
    root = pending[0][0]
    kids = defaultdict(list)
    for p,c,link in edges: kids[p].append((c,link))
    ordering, seen = [], set()
    def visit(n):
        if n in seen: raise ValueError('Cyclic structure')
        seen.add(n); ordering.append(n)
        for c,l in sorted(kids[n], key=lambda entry: -(int(entry[1].split('-')[1]) if entry[1].split('-')[1].isdigit() else -1)):
            visit(c)
    visit(root)
    if len(ordering) != len(nodes): raise ValueError('Disconnected glycan graph')
    ids = {n:str(i) for i,n in enumerate(ordering,1)}
    monos = {ids[n]:nodes[n] for n in ordering}
    links = dict(sorted(((f'{ids[p]}_{ids[c]}',l) for p,c,l in edges),key=lambda pair: tuple(map(int,pair[0].split('_')))))
    anomers = {ids[c]:l[0] for p,c,l in edges}
    return monos, links, anomers


def _pdb_residues(pdb_path):
    result, lookup = [], {}
    with Path(pdb_path).open(encoding='utf-8', errors='replace') as f:
        for line in f:
            if not line.startswith(('ATOM  ','HETATM')): continue
            segment = line[72:76].strip() or line[21:22].strip()
            key = (segment,line[22:26].strip(),line[26:27].strip())
            if key not in lookup:
                lookup[key] = len(result)
                result.append(dict(segment=segment, resSeq=line[22:26].strip(),
                                   resName=line[17:21].strip(),atoms=set(),
                                   topology_residue_number=len(result)+1))
            result[lookup[key]]['atoms'].add(line[12:16].strip())
    return result


# One-letter protein residue names (common PDB and CHARMM variants).
AMINO_ACID_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'HSD': 'H', 'HSE': 'H', 'HSP': 'H', 'HID': 'H', 'HIE': 'H',
    'HIP': 'H', 'CYX': 'C', 'CYM': 'C', 'ASH': 'D', 'GLH': 'E',
    'LYN': 'K',
}


def _pdb_records(pdb_path):
    """PDB residues with atom positions, in file order; excludes waters only when filtered."""
    residues, seen = [], {}
    with Path(pdb_path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            segment = line[72:76].strip() or line[21:22].strip()
            key = (segment, line[21:22].strip(), line[22:26].strip(), line[26:27].strip())
            if key not in seen:
                seen[key] = len(residues)
                residues.append({"segment": segment, "resSeq": line[22:26].strip(),
                                 "resName": line[17:21].strip(), "atoms": set(),
                                 "coords": {}, "topology_residue_number": len(residues)+1})
            res = residues[seen[key]]
            atom = line[12:16].strip()
            res["atoms"].add(atom)
            if atom not in res["coords"]:
                try:
                    res["coords"][atom] = tuple(float(line[k:k+8]) for k in (30,38,46))
                except ValueError:
                    pass
    return residues


def _pdb_glycan_edges(pdb_residues, bond_range=(1.20, 1.80)):
    """Find potential glycosidic C1/C2 -> O2/O3/O4/O6/etc contacts.

    Return (parent_index, child_index, donor_carbon, acceptor_position, distance).
    Geometrical evidence only; not a replacement for a bond/topology record.
    """
    from math import dist
    edges = []
    for parent, acceptor in enumerate(pdb_residues):
        for child, donor in enumerate(pdb_residues):
            if parent == child:
                continue
            for donor_pos in (1, 2):
                carbon = donor["coords"].get(f"C{donor_pos}")
                if carbon is None:
                    continue
                for atom, position in acceptor["coords"].items():
                    match = re.fullmatch(r"O(\d+)", atom)
                    if not match:
                        continue
                    # Exclude ring oxygen and carbonyl/other non-acceptor oxygens.
                    acceptor_pos = int(match.group(1))
                    if (acceptor_pos == 5 and "O5" in acceptor["coords"] and
                            "C5" in acceptor["coords"] and "C1" in acceptor["coords"]):
                        continue
                    if (acceptor_pos == 6 and "O6" in acceptor["coords"] and
                            "C6" in acceptor["coords"] and "C2" in acceptor["coords"] and
                            "C1" in acceptor["coords"] and "O5" not in acceptor["coords"]):
                        continue
                    distance = dist(carbon, position)
                    if bond_range[0] <= distance <= bond_range[1]:
                        edges.append((parent, child, donor_pos, acceptor_pos, round(distance, 3)))
    return edges


def _map_iupac_to_pdb(monos, links, pdb_residues, *, bond_range=(1.20, 1.80)):
    """Find unique labeled graph isomorphism; output mapping PDB-index -> IUPAC ID.

    Reject missing, extra, or ambiguous geometry rather than guessing mappings.
    """
    pdb_edges = _pdb_glycan_edges(pdb_residues, bond_range)
    if len(pdb_edges) != len(monos) - 1:
        raise ValueError(f"PDB has {len(pdb_edges)} candidate inter-residue bonds, "
                         f"IUPAC expects {len(monos)-1}. Use explicit topology or inspect distances. "
                         f"Detected: {pdb_edges}")
    expected = {}
    for pair, label in links.items():
        a, b = pair.split("_")
        match = re.fullmatch(r"[ab](\d+)-(\d+)", label)
        if match is None:
            raise ValueError(f"Cannot geometrically map non-numeric linkage {label!r}")
        expected[(a,b)] = tuple(map(int, match.groups()))
    observed = {(p,c): (don,acc) for p,c,don,acc,_ in pdb_edges}
    if len(observed) != len(pdb_edges):
        raise ValueError("Ambiguous PDB bond candidates between same residue pair")
    nodes = list(monos)
    possibilities = {i: [node for node in nodes if _pdb_identity_match(monos[node],res)]
                     for i,res in enumerate(pdb_residues)}
    if any(not x for x in possibilities.values()):
        raise ValueError(f"Some PDB residues cannot match IUPAC identities: {possibilities}")
    solutions = []
    def search(assigned, used):
        if len(solutions) > 1:
            return
        if len(assigned) == len(pdb_residues):
            solutions.append(dict(assigned))
            return
        i = min((x for x in possibilities if x not in assigned),
                key=lambda x: len([j for j in possibilities[x] if j not in used]))
        for node in possibilities[i]:
            if node in used:
                continue
            if any(observed.get((i,j)) != expected.get((node, other)) or
                   observed.get((j,i)) != expected.get((other, node))
                   for j,other in assigned.items()):
                continue
            assigned[i] = node
            search(assigned, used | {node})
            del assigned[i]
    search({}, set())
    if len(solutions) != 1:
        raise ValueError(f"Expected unique IUPAC/PDB graph mapping; found "
                         f"{'no' if not solutions else 'multiple'} solutions. "
                         "Supply authoritative bonds/topology or correct PDB/IUPAC.")
    return solutions[0], pdb_edges



def _is_carbohydrate_residue(res):
    """Recognize glycan-like residues using sugar atom signatures and names.

    Excludes peptide, solvent and ions, and tolerates PDB-truncated CHARMM IDs.
    The graph and IUPAC are subsequently checked to avoid false positives.
    """
    name = res['resName'].upper()
    if name in AMINO_ACID_CODES or name in {
        'HOH', 'WAT', 'TIP3', 'TIP', 'SOL', 'SOD', 'CLA', 'NA', 'CL',
        'POT', 'K', 'MG', 'CAL', 'CA', 'ZN'
    }:
        return False
    atoms = res['atoms']
    if not {'C1', 'C2', 'C3', 'C4'}.issubset(atoms):
        return False
    # 5- and 6-member sugar rings, including Neu5Ac C2/O6 ring
    return 'O5' in atoms or 'O6' in atoms or 'O4' in atoms


def detect_glycan_residues(pdb_path, expected_count=None, glycan_segment=None,
                           bond_range=(1.20, 1.80)):
    """Auto-select a glycan connected component in PDB residue order.

    If ``glycan_segment`` is provided it acts as an explicit constraint.
    Otherwise candidate carbohydrate residues are found independently of segID.
    Raises if none or multiple components meet the expected residue count.
    """
    records = _pdb_records(pdb_path)
    candidates = [r for r in records if _is_carbohydrate_residue(r)
                  and (glycan_segment is None or r['segment'] == glycan_segment)]
    if not candidates:
        raise ValueError('No carbohydrate residues identified in PDB')
    contacts = _pdb_glycan_edges(candidates, bond_range=bond_range)
    neighbors = {i: set() for i in range(len(candidates))}
    for p, c, *_ in contacts:
        neighbors[p].add(c)
        neighbors[c].add(p)
    components, visited = [], set()
    for i in range(len(candidates)):
        if i in visited:
            continue
        stack, group = [i], []
        visited.add(i)
        while stack:
            j = stack.pop()
            group.append(j)
            for k in neighbors[j]:
                if k not in visited:
                    visited.add(k)
                    stack.append(k)
        components.append([candidates[j] for j in sorted(group)])
    matching = [comp for comp in components
                if expected_count is None or len(comp) == expected_count]
    if len(matching) != 1:
        details = [(len(c), [(r['segment'], r['resSeq']) for r in c[:3]])
                   for c in components]
        raise ValueError(
            f'Expected one glycan component with {expected_count or "any"} residues; '
            f'found {len(matching)} matching components. Components: {details}. '
            'Supply glycan_segment or explicit connectivity if necessary.'
        )
    return matching[0]


def detect_linked_tripeptide(pdb_path, glycan_segment=None,
                              reducing_pdb_resseq=None, bond_range=(1.2,1.9)):
    """Detect a 3-residue Asn/Ser/Thr-linked glycopeptide via geometric bond.

    Returns its one-letter sequence (e.g. ANA, GSG, GTG) or None.
    """
    from math import dist
    residues = _pdb_records(pdb_path)
    sugar_candidates = [r for r in residues if _is_carbohydrate_residue(r)
                        and (glycan_segment is None or r['segment'] == glycan_segment)
                        and (reducing_pdb_resseq is None or
                             r['resSeq'] == str(reducing_pdb_resseq))]
    if reducing_pdb_resseq is None:
        sugar_candidates = sugar_candidates[:1]
    peptide_groups = defaultdict(list)
    for r in residues:
        if r['resName'].upper() in AMINO_ACID_CODES and (glycan_segment is None or r['segment'] != glycan_segment):
            peptide_groups[r['segment']].append(r)
    matches = []
    for group in peptide_groups.values():
        if len(group) != 3:
            continue
        for aa in group:
            name = aa['resName'].upper()
            protein_atom = {'ASN':'ND2', 'SER':'OG', 'THR':'OG1'}.get(name)
            if not protein_atom or protein_atom not in aa['coords']:
                continue
            for sugar in sugar_candidates:
                # N-glycans typically use GlcNAc; O-glycans can start with
                # GalNAc/GlcNAc/Man/Fuc/Xyl and other reducing-end sugars.
                for donor in ('C1','C2'):
                    if donor not in sugar['coords']:
                        continue
                    d = dist(aa['coords'][protein_atom], sugar['coords'][donor])
                    if bond_range[0] <= d <= bond_range[1]:
                        matches.append(''.join(AMINO_ACID_CODES[x['resName'].upper()] for x in group))
    if len(matches) > 1:
        raise ValueError(f'Ambiguous glycan-linked tripeptide: {matches}')
    return matches[0] if matches else None


def _pdb_identity_match(name, res):
    """A truncated 4-char CHARMM name needs atom-based distinction."""
    spec = REGISTRY[name]
    code = res['resName'].upper()
    full = {v.upper() for v in spec['charmm'].values()}
    if code not in set(spec['pdb']) | full | {x[:4] for x in full}:
        return False
    atoms = res['atoms']
    # For CHARMM-GUI truncated four-letter IDs BGLC/AGLC and BGAL/AGAL,
    # require N-acetyl signature for NAc; reject it for ordinary sugar.
    if name in ('GlcNAc','GalNAc','Glc','Gal'):
        has_acetamido = ('N2' in atoms or {'N', 'CT', 'C', 'O'}.issubset(atoms))
        if name.endswith('NAc') and not has_acetamido: return False
        if name in ('Glc','Gal') and has_acetamido: return False
    return True


def _forcefield_metadata(monos, links, anomers):
    accepted_positions = defaultdict(set)
    for key, linkage in links.items():
        parent, child = key.split('_')
        p = linkage.split('-')[1]
        if p.isdigit(): accepted_positions[parent].add(int(p))
    out = {}
    for node, name in monos.items():
        spec = REGISTRY[name]
        anomer = anomers.get(node)
        charmm = dict(spec['charmm'])
        # GLYCAM codes depend on occupied hydroxyls and alpha/beta configuration.
        glycam_code = None
        positions = tuple(sorted(accepted_positions[node]))
        if anomer in ('a','b') and spec.get('glycam') and positions in POSITION_CODES and (name != 'Neu5Ac' or positions == ()):
            glycam_code = POSITION_CODES[positions] + spec['glycam'] + anomer.upper()
            # Candidate codes must be validated against specific GLYCAM release;
            # GLYCAM anomer case must additionally reflect enantiomeric D/L form;
            # canonical L-Fuc/L-Rha/L-IdoA have lowercase middle-letter codes.
            if spec.get('default_L'):
                # Enantiomer must be specified explicitly for a validated exact code.
                glycam_code = None
        out[node] = dict(charmm=charmm,
                         glycam06_candidate=glycam_code,
                         glycam06_family=spec.get('glycam'),
                         pdb_ccd_names=spec['pdb'],
                         outgoing_anomer=anomer)
    return out


def find_glycan_directory(ROOT, glycan_id):
    """Locate a unique glycan directory anywhere under ROOT (or ROOT itself).

    Match directory names exactly to avoid confusing Man5 with e.g. Man50.
    A duplicate glycan ID is an error rather than an arbitrary selection.
    """
    root = Path(ROOT).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Glycan library root not found: {root}")
    matches = ([root] if root.name == str(glycan_id) else [])
    matches.extend(p for p in root.rglob(str(glycan_id)) if p.is_dir() and p.name == str(glycan_id))
    matches = sorted(set(matches))
    if not matches:
        raise FileNotFoundError(f"No directory named {glycan_id!r} found under {root}")
    if len(matches) > 1:
        raise ValueError(f"Multiple directories named {glycan_id!r} found under {root}: "
                         + ", ".join(map(str, matches)))
    return matches[0]


def find_glycan_pdb(glycan_directory):
    """Prefer step3_input.pdb in the glycan folder, then CHARMM-GUI/gromacs."""
    folder = Path(glycan_directory)
    direct = folder / "step3_input.pdb"
    if direct.is_file():
        return direct
    candidates = sorted(
        child / "gromacs" / "step3_input.pdb"
        for child in folder.iterdir()
        if child.is_dir() and "charmm-gui" in child.name.lower()
        and (child / "gromacs" / "step3_input.pdb").is_file()
    )
    if not candidates:
        raise FileNotFoundError(
            f"No step3_input.pdb in {folder} or in a *charmm-gui*/gromacs subfolder"
        )
    if len(candidates) > 1:
        raise ValueError(f"Multiple CHARMM-GUI PDBs found; pass pdb_path explicitly: {candidates}")
    return candidates[0]


def create_graph_nodes(record):
    """Build Man5-style node records without the obsolete yaml_position."""
    modifications = record.get("terminal_modifications", [])
    # Terminal modifications may be given as node IDs, or dictionaries with node/id.
    modified = set()
    for item in modifications:
        if isinstance(item, dict):
            node = item.get("node", item.get("id"))
            if node is not None:
                modified.add(str(node))
        elif isinstance(item, (str, int)):
            modified.add(str(item))
    return [
        {
            "id": str(node),
            "mono": sugar,
            "terminal_modification": str(node) in modified,
            "pdb_resid": record.get("pdb_residues", {}).get(str(node)),
            "pucker_ring_atoms": record["pucker_ring_atoms"][str(node)],
        }
        for node, sugar in record["monosaccharides"].items()
    ]


def _find_charmm_topology_files(pdb_path):
    """Find the GROMACS topology files belonging to the selected PDB."""
    pdb_path = Path(pdb_path).resolve()
    if pdb_path.parent.name == "gromacs":
        gromacs = pdb_path.parent
    else:
        candidates = sorted(
            d / "gromacs" for d in pdb_path.parent.iterdir()
            if d.is_dir() and "charmm-gui" in d.name.lower()
            and (d / "gromacs").is_dir()
        )
        if len(candidates) > 1:
            raise ValueError(f"Multiple CHARMM-GUI GROMACS directories for {pdb_path}: {candidates}")
        if not candidates:
            return []
        gromacs = candidates[0]
    files = sorted((gromacs / "toppar").rglob("*.itp")) if (gromacs / "toppar").is_dir() else []
    files += sorted(gromacs.glob("*.itp"))
    files += sorted(gromacs.glob("*.top"))
    return list(dict.fromkeys(files))


def _topology_residue_signatures(files):
    """Read GROMACS [ atoms ] atom names grouped by residue and molecule.

    The residue-name column is not restricted to four characters, unlike PDB.
    Parsing only the [ atoms ] section avoids accidentally reading atomtypes.
    """
    from collections import defaultdict
    entries = []
    for path in files:
        section = None
        molecule = None
        residues = defaultdict(set)
        for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw.partition(";")[0].strip()
            if not line or line.startswith("#"):
                continue
            hit = re.match(r"^\[\s*([^\]]+)\s*\]$", line)
            if hit:
                section = hit.group(1).strip().lower()
                if section == "moleculetype":
                    if residues:
                        entries.extend((name, frozenset(atoms), path) for (number, name), atoms in residues.items())
                        residues.clear()
                    molecule = None
                continue
            if section == "moleculetype" and molecule is None:
                molecule = line.split()[0]
            elif section == "atoms" and molecule is not None:
                cols = line.split()
                if len(cols) >= 5 and cols[0].isdigit() and cols[2].lstrip("-").isdigit():
                    # nr type resnr residue atom cgnr charge mass
                    residues[(cols[2], cols[3].upper())].add(cols[4].upper())
        if residues:
            entries.extend((name, frozenset(atoms), path) for (number, name), atoms in residues.items())
    return entries


def _topology_identity(residue, signatures):
    """Resolve a residue through full GROMACS residue names AND its atom signature.

    Require all PDB heavy atoms to match a topology residue signature. Missing
    hydrogen atoms are tolerated, but different heavy-atom compositions are not.
    Return (name, anomer, evidence), or None when topology is inconclusive.
    """
    code = residue["resName"].upper()
    pdb_atoms = {a.upper() for a in residue["atoms"] if not a.upper().startswith(("H", "D"))}
    compatible = set()
    evidence = set()
    for resname, atom_names, source in signatures:
        if code != resname and code != resname[:4]:
            continue
        heavy = {a for a in atom_names if not a.startswith(("H", "D"))}
        if not pdb_atoms.issubset(heavy):
            continue
        for name, spec in REGISTRY.items():
            for anomer, full in spec["charmm"].items():
                if resname == full.upper():
                    compatible.add((name, anomer))
                    evidence.add(f"{source.name}:{resname}")
    if len(compatible) == 1:
        name, anomer = next(iter(compatible))
        return name, anomer, sorted(evidence)
    return None


def _identify_pdb_monosaccharide(residue):
    """Resolve a CHARMM/PDB carbohydrate residue; never choose ambiguous identities."""
    matches = [name for name in REGISTRY if _pdb_identity_match(name, residue)]
    if len(matches) != 1:
        raise ValueError(
            f"Cannot uniquely identify residue {residue['resName']} "
            f"{residue['segment']}:{residue['resSeq']}: candidates={matches}. "
            "Check chemical modifications and residue names manually."
        )
    return matches[0]


def _pdb_anomer(residue, name):
    """Only use anomer labels unambiguously encoded by CHARMM residue names."""
    code = residue['resName'].upper()
    charmm = REGISTRY[name]['charmm']
    matches = {anomer for anomer, full in charmm.items()
               if code in (full.upper(), full[:4].upper())}
    return next(iter(matches)) if len(matches) == 1 else None


def _serialize_glycan_iupac(monos, links, root):
    """Serialize a reducing-end rooted tree using parser-compatible bracket syntax."""
    children = defaultdict(list)
    for key, linkage in links.items():
        parent, child = key.split('_')
        children[parent].append((child, linkage))

    visited = set()

    def render(node):
        if node in visited:
            raise ValueError('Cycle detected in reconstructed glycan')
        visited.add(node)
        branches = sorted(children[node], key=lambda x: (-int(x[1].split('-')[1]), int(x[0])))
        if not branches:
            return monos[node]
        fragments = [render(child) + f'({label})' for child, label in branches]
        return fragments[0] + ''.join(f'[{frag}]' for frag in fragments[1:]) + monos[node]

    result = render(root)
    if len(visited) != len(monos):
        raise ValueError('Reconstructed glycan contains disconnected nodes')
    return result


def reconstruct_iupac_from_pdb(pdb_path, *, glycan_segment=None,
                               bond_range=(1.20, 1.80)):
    """Reconstruct connectivity, anomers and a provisional IUPAC from CHARMM PDB.

    Status meanings:
      verified: reserved for externally/topologically verified reconstruction;
      partially_verified: PDB geometry plus unambiguous CHARMM residue identities;
      requires_manual_inspection: unknown anomers or possible modifications.

    No 'verified' status is assigned from coordinates alone.
    Returns (iupac, details) where details includes glycan records and warnings.
    """
    pdb_path = Path(pdb_path).expanduser().resolve()
    residues = detect_glycan_residues(pdb_path, glycan_segment=glycan_segment,
                                      bond_range=bond_range)
    topology_files = _find_charmm_topology_files(pdb_path)
    topology_signatures = _topology_residue_signatures(topology_files)
    resolved = [_topology_identity(r, topology_signatures) for r in residues]
    names = [item[0] if item else _identify_pdb_monosaccharide(r)
             for r, item in zip(residues, resolved)]
    contacts = _pdb_glycan_edges(residues, bond_range=bond_range)
    if len(contacts) != len(residues) - 1:
        raise ValueError(f'Expected {len(residues)-1} bonds, found {len(contacts)}: {contacts}')

    monos = {str(i + 1): name for i, name in enumerate(names)}
    anomers = {}
    links = {}
    warnings_list = []
    if topology_signatures:
        warnings_list.append("Residue identities matched to topology [ atoms ] records by atom names; PDB bonds still geometry-inferred")
    else:
        warnings_list.append("No usable CHARMM-GUI topology [ atoms ] records were found")
    child_nodes = set()
    for p, c, donor, acceptor, distance in contacts:
        parent, child = str(p + 1), str(c + 1)
        if child in child_nodes:
            raise ValueError(f'Residue {child} has multiple parent bonds')
        child_nodes.add(child)
        required_donor = REGISTRY[names[c]].get('donor', 1)
        if donor != required_donor:
            raise ValueError(f'Unexpected donor C{donor} for {names[c]} (expected C{required_donor})')
        anomer = resolved[c][1] if resolved[c] else _pdb_anomer(residues[c], names[c])
        if anomer is None:
            warnings_list.append(f'Unknown anomer for residue {child} ({residues[c]["resName"]})')
        else:
            anomers[child] = anomer
        links[f'{parent}_{child}'] = f'{anomer or "?"}{donor}-{acceptor}'

    roots = set(monos) - child_nodes
    if len(roots) != 1:
        raise ValueError(f'Expected one reducing-end root, got {sorted(roots)}')
    root = roots.pop()
    # The geometry identifies a potential C--O bond, not its force-field definition.
    warnings_list.append('Inter-residue bonds inferred from PDB distances, not verified against topology')
    # Do not silently interpret modified sugars as their unmodified registry entries.
    all_residues = _pdb_records(pdb_path)
    suspected = [f"{r['resName']}:{r['resSeq']}" for r in all_residues
                 if (r['resName'].upper() in {'SO4', 'SUL', 'SULF', 'PO4', 'PHO'}
                     or any(atom.upper().startswith(('SUL', 'SO4', 'PO4'))
                            for atom in r['atoms']))]
    if suspected:
        warnings_list.append('Possible sulfate/phosphate modifications: ' + ', '.join(suspected))
    for i, (r, name) in enumerate(zip(residues, names), 1):
        ring = REGISTRY[name]['ring']
        if ring is None or not set(ring).issubset(r['atoms']):
            raise ValueError(f'Residue {i} has no validated ring atoms for {name}')
    iupac = _serialize_glycan_iupac(monos, links, root)
    # Do not claim verified when stereochemistry or modifications are unresolved.
    status = ('requires_manual_inspection' if any(
        x.startswith(('Unknown anomer', 'Possible sulfate')) for x in warnings_list)
              else 'partially_verified')
    details = {
        'status': status,
        'warnings': warnings_list,
        'source_pdb': str(pdb_path),
        'monosaccharides': monos,
        'linkages': links,
        'anomers': anomers,
        'pdb_records': residues,
        'pdb_contacts': contacts,
        'root_node': root,
        'topology_sources': [str(f) for f in topology_files],
        'topology_residues_resolved': sum(x is not None for x in resolved),
    }
    return iupac, details


def create_glycan_json(ROOT, glycan_id, iupac=None, pdb_path=None, output_path=None,
                       *, glycan_segment=None, overwrite=False,
                       pdb_node_order=None):
    """Write glycan metadata with PDB-order IDs and IUPAC-derived labeled edges.

    With a PDB, atom-coordinate bond candidates provide the mapping.
    Never infer mapping simply from equal numbers of residues or PDB order.
    """
    glycan_dir = find_glycan_directory(ROOT, glycan_id)
    if pdb_path is None:
        pdb_path = find_glycan_pdb(glycan_dir)
    else:
        pdb_path = Path(pdb_path).expanduser().resolve()
        if not pdb_path.is_file():
            raise FileNotFoundError(pdb_path)
    if output_path is None:
        output_path = glycan_dir / f"{glycan_id}.json"

    # Auto reconstruction bypasses the IUPAC->PDB isomorphism mapper.
    reconstructed = None
    if iupac is None:
        iupac, reconstructed = reconstruct_iupac_from_pdb(
            pdb_path, glycan_segment=glycan_segment)
        monos = reconstructed['monosaccharides']
        links = reconstructed['linkages']
        anomers = reconstructed['anomers']
    else:
        monos, links, anomers = _parse_iupac(iupac)
    out = Path(output_path or f'{glycan_id}.json')
    if out.exists() and not overwrite:
        print(f"File already exists: {out}. Skipping.")
        return out
    record = dict(glycan_id=glycan_id, iupac=iupac,
                  terminal_modifications=[], peptide=None)
    if pdb_path is not None:
        pdb = detect_glycan_residues(pdb_path, expected_count=len(monos), glycan_segment=glycan_segment)
        if len(pdb) != len(monos):
            raise ValueError(f'Expected {len(monos)} glycan residues, got {len(pdb)}')
        if reconstructed is not None:
            index_to_old = {i: str(i+1) for i in range(len(pdb))}
            contacts = reconstructed['pdb_contacts']
        else:
            index_to_old, contacts = _map_iupac_to_pdb(monos, links, pdb)
        if pdb_node_order is not None:
            requested = [str(x) for x in pdb_node_order]
            if requested != [index_to_old[i] for i in range(len(pdb))]:
                raise ValueError('Provided pdb_node_order conflicts with topology-based mapping')
        remap = {index_to_old[i]: str(i+1) for i in range(len(pdb))}
        record['monosaccharides'] = {remap[index_to_old[i]]: monos[index_to_old[i]]
                                     for i in range(len(pdb))}
        record['linkages'] = {f'{remap[a]}_{remap[b]}': label
                              for key,label in links.items()
                              for a,b in [key.split('_')]}
        new_anomers = {remap[k]:v for k,v in anomers.items()}
        record['pdb_residues'] = {str(i+1):r['topology_residue_number']
                                   for i,r in enumerate(pdb)}
        record['pdb_original_resseq'] = {str(i+1):r['resSeq'] for i,r in enumerate(pdb)}
        record['pdb_residue_names'] = {str(i+1):r['resName'] for i,r in enumerate(pdb)}
        record['pdb_segments'] = {str(i+1):r['segment'] for i,r in enumerate(pdb)}
        for i,r in enumerate(pdb):
            node = str(i+1)
            ring = REGISTRY[record['monosaccharides'][node]]['ring']
            if ring is not None and not set(ring).issubset(r['atoms']):
                raise ValueError(f'PDB node {node}: missing ring atoms {set(ring)-r["atoms"]}')
        record['pucker_ring_atoms'] = {node:list(REGISTRY[name]['ring']) if REGISTRY[name]['ring'] else None
                                        for node,name in record['monosaccharides'].items()}
        record['pdb_mapping_method'] = ('direct_PDB_graph_reconstruction' if reconstructed is not None
                                        else 'unique_graph_isomorphism_with_interresidue_distance_bonds')
        record['pdb_glycosidic_bonds'] = [
            dict(parent=str(p+1), child=str(c+1), donor_carbon=f'C{don}',
                 acceptor_oxygen=f'O{acc}', distance_angstrom=d)
            for p,c,don,acc,d in contacts]
        root_old = reconstructed['root_node'] if reconstructed is not None else next(iter(monos))
        reducing_res = pdb[int(remap[root_old])-1]
        record['peptide'] = detect_linked_tripeptide(
            pdb_path, glycan_segment=reducing_res['segment'],
            reducing_pdb_resseq=reducing_res['resSeq'])
    else:
        record['monosaccharides'] = monos
        record['linkages'] = links
        new_anomers = anomers
        record['pdb_residues'] = {}
        record['pucker_ring_atoms'] = {node:list(REGISTRY[name]['ring']) if REGISTRY[name]['ring'] else None
                                        for node,name in monos.items()}
        record['pdb_mapping_method'] = None
    if reconstructed is not None:
        record['reconstruction'] = {
            'status': reconstructed['status'],
            'warnings': reconstructed['warnings'],
            'source_pdb': reconstructed['source_pdb'],
            'user_checked': False,
            'topology_sources': reconstructed.get('topology_sources', []),
            'topology_residues_resolved': reconstructed.get('topology_residues_resolved', 0),
        }
        print(f"{glycan_id}: IUPAC reconstruction {reconstructed['status']}: {iupac}")
        for warning in reconstructed['warnings']:
            print(f'  - {warning}')
    else:
        record['reconstruction'] = {
            'status': 'provided_by_user',
            'warnings': [],
            'source_pdb': str(pdb_path),
            'user_checked': True,
        }
    record['forcefield_residue_names'] = _forcefield_metadata(
        record['monosaccharides'],record['linkages'],new_anomers)
    # Graph edges use stable PDB-order node identifiers when a PDB is supplied.
    record['graph_edges'] = [
        {'id':key,'parent':a,'child':b,'linkage':l}
        for key,l in record['linkages'].items()
        for a,b in [key.split('_')]
    ]
    record['graph_nodes'] = create_graph_nodes(record)
    children = {edge['child'] for edge in record['graph_edges']}
    roots = set(record['monosaccharides']) - children
    if len(roots) != 1:
        raise ValueError(f"Expected exactly one reducing-end root, found {roots}")
    record['root_node'] = roots.pop()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return out


#def update_glycan_json(path, **fields):
#    path = Path(path)
#    record = json.loads(path.read_text(encoding='utf-8'))
#    record.update(fields)
#    path.write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')


#def update_glycan_section(path, section, updates):
#    path = Path(path)
#    record = json.loads(path.read_text(encoding='utf-8'))
#    if section not in record: record[section] = {}
#    if not isinstance(record[section],dict): raise TypeError(f'{section} is not a JSON object')
#    record[section].update(updates)
#    path.write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')


def plot_graph(
    ROOT,
    glycan_id,
    save_path=None,
    figsize=(12, 8),
    node_size=950,
    font_size=20,
    show_labels=True,
    show=True,
):
    """
    Plot a glycan graph using SNFG (Symbol Nomenclature for Glycans).

    Nodes:
        SNFG shapes and colors based on monosaccharide identity.
        Node IDs correspond to PDB residue order.

    Edges:
        Directed parent -> child.
        Labeled with glycosidic linkages (e.g. a1-3, b1-4, a2-8).

    Parameters
    ----------
    ROOT : str or Path
        Path to the glycan library.
    glycan_id : str
        Name of the glycan.
    save_path : str or Path, optional
        Save figure as PNG, PDF, SVG, etc.
    figsize : tuple
        Figure dimensions.
    node_size : int
        Size of SNFG symbols.
    font_size : int
        Font size for labels.
    show_labels : bool
        Whether to display node IDs and monosaccharide names.
    show : bool
        Whether to display the figure.

    Returns
    -------
    fig, ax
        Matplotlib figure and axes.
    """

    # ========================================================
    # SNFG nomenclature
    # Reference:
    # https://www.ncbi.nlm.nih.gov/glycans/snfg.html
    # ========================================================

    COLORS = {
        "white": "#FFFFFF",
        "blue": "#0072BC",
        "green": "#00A651",
        "yellow": "#FFD400",
        "orange": "#F47920",
        "pink": "#F69EA1",
        "purple": "#A54399",
        "lightblue": "#8FCCE9",
        "brown": "#A17A4D",
        "red": "#ED1C24",
    }

    # SNFG shape and color definitions
    SNFG = {}

    def register(names, shape, color):
        for name in names:
            SNFG[name] = (shape, COLORS[color])

    # Hexoses: circles
    register(["Glc"], "circle", "blue")
    register(["Man"], "circle", "green")
    register(["Gal"], "circle", "yellow")
    register(["Gul"], "circle", "orange")
    register(["Alt"], "circle", "pink")
    register(["All"], "circle", "purple")
    register(["Tal"], "circle", "lightblue")
    register(["Ido"], "circle", "brown")

    # N-acetylhexosamines: squares
    register(["GlcNAc"], "square", "blue")
    register(["ManNAc"], "square", "green")
    register(["GalNAc"], "square", "yellow")
    register(["GulNAc"], "square", "orange")
    register(["AltNAc"], "square", "pink")
    register(["AllNAc"], "square", "purple")
    register(["TalNAc"], "square", "lightblue")
    register(["IdoNAc"], "square", "brown")

    # Hexosamines: crossed squares
    register(["GlcN"], "crossed_square", "blue")
    register(["ManN"], "crossed_square", "green")
    register(["GalN"], "crossed_square", "yellow")
    register(["GulN"], "crossed_square", "orange")
    register(["AltN"], "crossed_square", "pink")
    register(["AllN"], "crossed_square", "purple")
    register(["TalN"], "crossed_square", "lightblue")
    register(["IdoN"], "crossed_square", "brown")

    # Uronic acids: divided diamonds
    register(["GlcA"], "divided_diamond", "blue")
    register(["ManA"], "divided_diamond", "green")
    register(["GalA"], "divided_diamond", "yellow")
    register(["GulA"], "divided_diamond", "orange")
    register(["AltA"], "divided_diamond", "pink")
    register(["AllA"], "divided_diamond", "purple")
    register(["TalA"], "divided_diamond", "lightblue")
    register(["IdoA"], "divided_diamond", "brown")

    # Deoxyhexoses: triangles
    register(["Qui"], "triangle", "blue")
    register(["Rha"], "triangle", "green")
    register(["6dGul"], "triangle", "orange")
    register(["6dAlt"], "triangle", "pink")
    register(["6dTal"], "triangle", "lightblue")
    register(["Fuc"], "triangle", "red")

    # N-acetyl deoxyhexoses: divided triangles
    register(["QuiNAc"], "divided_triangle", "blue")
    register(["RhaNAc"], "divided_triangle", "green")
    register(["6dAltNAc"], "divided_triangle", "pink")
    register(["6dTalNAc"], "divided_triangle", "lightblue")
    register(["FucNAc"], "divided_triangle", "red")

    # Pentoses: stars
    register(["Ara"], "star", "green")
    register(["Lyx"], "star", "yellow")
    register(["Xyl"], "star", "orange")
    register(["Rib"], "star", "pink")

    # Sialic acids and related nonulosonic acids: diamonds
    register(["Kdn"], "diamond", "green")
    register(["Neu5Ac"], "diamond", "purple")
    register(["Neu5Gc"], "diamond", "lightblue")
    register(["Neu"], "diamond", "brown")
    register(["Sia"], "diamond", "red")

    # Other ketoses / acidic sugars
    register(["Fru"], "pentagon", "green")
    register(["Kdo"], "hexagon", "green")
    register(["Pse"], "flat_diamond", "green")
    register(["Leg"], "flat_diamond", "yellow")
    register(["Aci"], "flat_diamond", "pink")
    register(["4eLeg"], "flat_diamond", "lightblue")

    # ========================================================
    # Load JSON
    # ========================================================

    glycan_dir = find_glycan_directory(ROOT, glycan_id)
    json_path = glycan_dir / f"{glycan_id}.json"

    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    monosaccharides = data["monosaccharides"]
    edges = data["graph_edges"]

    # ========================================================
    # Construct graph
    # ========================================================

    G = nx.DiGraph()

    for node_id, sugar in monosaccharides.items():
        G.add_node(str(node_id), sugar=sugar)

    for edge in edges:
        parent = str(edge["parent"])
        child = str(edge["child"])

        G.add_edge(
            parent,
            child,
            linkage=edge["linkage"],
        )

    if not nx.is_directed_acyclic_graph(G):
        raise ValueError(
            "The glycan graph contains a cycle. "
            "Expected a directed acyclic glycan graph."
        )

    # ========================================================
    # Tree layout, reducing end on the left
    # ========================================================

    roots = [n for n in G.nodes if G.in_degree(n) == 0]

    if len(roots) != 1:
        raise ValueError(
            f"Expected exactly one glycan root, found {roots}"
        )

    root = roots[0]

    if not nx.is_arborescence(G):
        raise ValueError(
            "Expected a connected, rooted glycan tree."
        )

    pos = {}
    leaf_counter = [0]

    def layout(node, depth=0):
        children = list(G.successors(node))

        # Higher acceptor position first
        def acceptor_position(child):
            linkage = G.edges[node, child]["linkage"]

            try:
                return int(str(linkage).split("-")[-1])
            except (ValueError, IndexError):
                return -1

        children.sort(
            key=lambda c: (
                -acceptor_position(c),
                int(c) if c.isdigit() else c,
            )
        )

        if not children:
            y = leaf_counter[0]
            leaf_counter[0] += 1
        else:
            child_y = [
                layout(child, depth + 1)
                for child in children
            ]
            y = sum(child_y) / len(child_y)

        pos[node] = (depth * 2.2, -y * 1.6)

        return y

    layout(root)

    # ========================================================
    # Figure
    # ========================================================

    fig, ax = plt.subplots(figsize=figsize)

    ax.set_aspect("equal")
    ax.axis("off")

    symbol_radius = 0.33 * np.sqrt(node_size / 950)

    # ========================================================
    # Draw edges first
    # ========================================================

    for parent, child, attributes in G.edges(data=True):

        x1, y1 = pos[parent]
        x2, y2 = pos[child]

        dx = x2 - x1
        dy = y2 - y1

        length = np.hypot(dx, dy)

        ux = dx / length
        uy = dy / length

        start = (
            x1 + ux * symbol_radius,
            y1 + uy * symbol_radius,
        )

        end = (
            x2 - ux * symbol_radius,
            y2 - uy * symbol_radius,
        )

        arrow = FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.5,
            color="black",
            zorder=1,
        )

        ax.add_patch(arrow)

        # Edge linkage label
        mx = (start[0] + end[0]) / 2
        my = (start[1] + end[1]) / 2

        ax.text(
            mx,
            my + 0.17,
            attributes["linkage"],
            ha="center",
            va="center",
            fontsize=font_size,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.9,
                "pad": 1.5,
            },
            zorder=5,
        )

    # ========================================================
    # SNFG symbol drawing
    # ========================================================

    def draw_symbol(x, y, shape, color, radius):

        common = {
            "facecolor": color,
            "edgecolor": "black",
            "linewidth": 1.3,
            "zorder": 3,
        }

        if shape == "circle":
            patch = Circle((x, y), radius, **common)

        elif shape in {"square", "crossed_square"}:
            patch = Rectangle(
                (x - radius, y - radius),
                2 * radius,
                2 * radius,
                **common,
            )

        elif shape in {"triangle", "divided_triangle"}:
            patch = RegularPolygon(
                (x, y),
                numVertices=3,
                radius=radius * 1.35,
                orientation=np.pi / 2,
                **common,
            )

        elif shape in {"diamond", "divided_diamond"}:
            patch = Polygon(
                [
                    (x, y + radius * 1.3),
                    (x + radius * 1.3, y),
                    (x, y - radius * 1.3),
                    (x - radius * 1.3, y),
                ],
                closed=True,
                **common,
            )

        elif shape == "flat_diamond":
            patch = Polygon(
                [
                    (x - radius, y),
                    (x, y + radius * 0.65),
                    (x + radius, y),
                    (x, y - radius * 0.65),
                ],
                closed=True,
                **common,
            )

        elif shape == "star":
            vertices = []
            for i in range(10):
                angle = np.pi / 2 + i * np.pi / 5
                r = radius * (1.25 if i % 2 == 0 else 0.55)
                vertices.append(
                    (x + r * np.cos(angle), y + r * np.sin(angle))
                )

            patch = Polygon(vertices, closed=True, **common)

        elif shape in {"pentagon", "hexagon"}:
            patch = RegularPolygon(
                (x, y),
                numVertices=5 if shape == "pentagon" else 6,
                radius=radius * 1.2,
                orientation=np.pi / 2,
                **common,
            )

        else:
            # Unknown residue: white hexagon with abbreviation
            patch = RegularPolygon(
                (x, y),
                numVertices=6,
                radius=radius * 1.2,
                **common,
            )

        ax.add_patch(patch)

        # Crossed square for hexosamines
        if shape == "crossed_square":
            ax.plot(
                [x - radius, x + radius],
                [y - radius, y + radius],
                color="black",
                linewidth=1.2,
                zorder=4,
            )
            ax.plot(
                [x - radius, x + radius],
                [y + radius, y - radius],
                color="black",
                linewidth=1.2,
                zorder=4,
            )

        # Divided SNFG shapes
        if shape in {"divided_diamond", "divided_triangle"}:
            ax.plot(
                [x, x],
                [y - radius * 0.8, y + radius * 0.8],
                color="black",
                linewidth=1.2,
                zorder=4,
            )

    # ========================================================
    # Draw nodes
    # ========================================================

    for node, attributes in G.nodes(data=True):

        sugar = attributes["sugar"]
        x, y = pos[node]

        if sugar not in SNFG:
            warnings.warn(
                f"No SNFG definition for '{sugar}'. "
                "Using a white hexagon.",
                stacklevel=2,
            )

        shape, color = SNFG.get(
            sugar,
            ("hexagon", COLORS["white"]),
        )

        draw_symbol(
            x,
            y,
            shape,
            color,
            symbol_radius,
        )

        if sugar not in SNFG:
            ax.text(
                x,
                y,
                sugar,
                fontsize=font_size * 0.75,
                ha="center",
                va="center",
                zorder=6,
            )

        if show_labels:
            ax.text(
                x,
                y - symbol_radius * 1.7,
                f"{node}: {sugar}",
                fontsize=font_size,
                ha="center",
                va="top",
                zorder=6,
            )

    # ========================================================
    # Final formatting
    # ========================================================

    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]

    ax.set_xlim(min(xs) - 1.0, max(xs) + 1.0)
    ax.set_ylim(min(ys) - 1.3, max(ys) + 1.0)

    fig.tight_layout()

    if save_path is None:
        save_path = glycan_dir / f"{glycan_id}.png"

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight",
        )

    if show:
        plt.show()

    return fig, ax