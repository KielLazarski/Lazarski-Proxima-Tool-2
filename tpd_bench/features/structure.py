"""Minimal structure parsing and per-residue geometry, pure Python.

The degron scorer needs secondary structure, relative solvent accessibility
and pLDDT for every residue. DSSP supplies the first two and is not available
without third-party libraries, so both are approximated from CA geometry.

These are PROXIES. They are good enough to rank a proteome and not good enough
to make a structural claim about one protein. Anything the sweep nominates
should have its degron re-checked with real DSSP before anyone acts on it.

Secondary structure from CA(i)-CA(i+3) and CA(i)-CA(i+4) distances:
  alpha helix  CA(i)-CA(i+3) about 5.0-5.5 A, CA(i)-CA(i+4) about 6.2 A
  beta strand  extended, CA(i)-CA(i+3) about 9.5-11 A
This catches regular secondary structure reliably and is weakest exactly where
DSSP is also weakest, at the ends of elements.

Solvent accessibility from neighbour count: residues are buried in proportion
to how many other residues sit within 10 A. A standard proxy, monotonic with
true RSA, and it does not need surface-area integration.

pLDDT comes straight from the B-factor column, which is where AlphaFold puts it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .degron import Residue

# three-letter to one-letter, including the residues that actually appear
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

NEIGHBOUR_RADIUS = 10.0
# Neighbour counts at which a residue reads as fully buried / fully exposed.
MIN_HELIX_RUN = 4
BURIED_AT = 24
EXPOSED_AT = 8


@dataclass
class Atom:
    resseq: int
    resname: str
    name: str
    x: float
    y: float
    z: float
    bfactor: float


def parse_pdb(text: str) -> list[Atom]:
    """Column-based PDB parse of CA atoms only. Ignores everything else."""
    out = []
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        name = line[12:16].strip()
        if name != "CA":
            continue
        resname = line[17:20].strip()
        if resname not in AA3:
            continue
        try:
            out.append(
                Atom(
                    resseq=int(line[22:26]),
                    resname=resname,
                    name=name,
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    bfactor=float(line[60:66] or 0.0),
                )
            )
        except ValueError:
            continue
    return out


def _dist(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def secondary_structure(cas: list[Atom]) -> list[str]:
    """H / E / C per residue, from CA spacing.

    A helix must SUSTAIN its signal. Tight turns -- including the beta-hairpin
    turns that carry degrons -- produce short CA(i)-CA(i+3) distances that
    match the helix criterion for one or two residues, so a per-residue test
    labels every hairpin turn as helix and the degron scorer then finds
    nothing. Requiring a run of consecutive helical residues separates a real
    helix from a turn that merely looks like the start of one.
    """
    n = len(cas)
    ss = ["C"] * n

    helical = [False] * n
    extended = [False] * n
    for i in range(n):
        d3 = _dist(cas[i], cas[i + 3]) if i + 3 < n else None
        d4 = _dist(cas[i], cas[i + 4]) if i + 4 < n else None
        if d3 is not None and 4.5 <= d3 <= 6.3 and d4 is not None and 5.2 <= d4 <= 7.4:
            helical[i] = True
        if d3 is not None and d3 >= 9.0:
            extended[i] = True

    # a helix needs MIN_HELIX_RUN consecutive helical starts
    i = 0
    while i < n:
        if helical[i]:
            j = i
            while j < n and helical[j]:
                j += 1
            if (j - i) >= MIN_HELIX_RUN:
                for k in range(i, min(j + 3, n)):
                    ss[k] = "H"
            i = j
        else:
            i += 1

    for i in range(n):
        if extended[i]:
            for k in range(i, min(i + 3, n)):
                if ss[k] == "C":
                    ss[k] = "E"
    return ss


def rsa_proxy(cas: list[Atom]) -> list[float]:
    """0 (buried) to 1 (exposed), from neighbour density.

    O(n^2) over CA atoms. Fine to a few thousand residues; a grid would be
    needed above that, and human proteins rarely get there.
    """
    n = len(cas)
    out = []
    r2 = NEIGHBOUR_RADIUS ** 2
    for i in range(n):
        a = cas[i]
        count = 0
        for j in range(n):
            if i == j:
                continue
            b = cas[j]
            if (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2 <= r2:
                count += 1
        frac = (BURIED_AT - count) / (BURIED_AT - EXPOSED_AT)
        out.append(max(0.0, min(1.0, frac)))
    return out


def to_residues(pdb_text: str) -> list[Residue]:
    """PDB text to the Residue list the degron scorer consumes."""
    cas = parse_pdb(pdb_text)
    if not cas:
        return []
    ss = secondary_structure(cas)
    rsa = rsa_proxy(cas)
    return [
        Residue(
            index=a.resseq,
            aa=AA3[a.resname],
            ss=ss[i],
            rsa=rsa[i],
            plddt=a.bfactor,  # AlphaFold stores pLDDT here
        )
        for i, a in enumerate(cas)
    ]
