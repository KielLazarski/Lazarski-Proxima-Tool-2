"""Tests for geometry-derived secondary structure and burial.

The hairpin test exists because of a real bug: tight turns produce short
CA(i)-CA(i+3) distances that match the helix criterion, so a per-residue test
labelled every beta-hairpin turn as helix. The degron scorer then found
nothing anywhere, silently, and a whole sweep would have returned zeros with
no error. Requiring a sustained helical run fixes it, and this test keeps it
fixed.
"""

import math

from .degron import score_protein
from .structure import parse_pdb, rsa_proxy, secondary_structure, to_residues


def at(i, rn, x, y, z, b=95.0):
    return f"ATOM  {i:5d}  CA  {rn} A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{b:6.2f}           C"


def ideal_helix(n=20):
    out = []
    for k in range(n):
        a = math.radians(100 * k)
        out.append(at(k + 1, "ALA", 2.3 * math.cos(a), 2.3 * math.sin(a), 1.5 * k))
    return "\n".join(out)


def beta_hairpin(turn=(("SER", 1.2, 2.2), ("GLY", 2.6, 3.4), ("GLY", 4.0, 3.4), ("SER", 5.0, 2.2))):
    lines, i, N = [], 1, 7
    for k in range(N):
        lines.append(at(i, "VAL", 0.0, 0.0, 3.35 * k)); i += 1
    z0 = 3.35 * (N - 1)
    for rn, dx, dz in turn:
        lines.append(at(i, rn, dx, 0.0, z0 + dz)); i += 1
    for k in range(N):
        lines.append(at(i, "LYS" if k == 2 else "ILE", 4.8, 0.0, z0 - 3.35 * k)); i += 1
    return "\n".join(lines)


def run():
    checks = []

    ss = "".join(secondary_structure(parse_pdb(ideal_helix())))
    checks.append(("ideal helix reads as helix", ss.count("H") / len(ss) > 0.8, ss))

    ss = "".join(secondary_structure(parse_pdb(beta_hairpin())))
    checks.append(("hairpin strands read as strand", ss.count("E") >= 10, ss))
    checks.append(("hairpin TURN is coil, not helix", "H" not in ss, ss))

    s_hp, hits = score_protein(to_residues(beta_hairpin()))
    checks.append(("hairpin yields a degron", s_hp > 0.5, f"{s_hp}"))
    checks.append((
        "G-loop detected in the turn",
        any("G-loop" in f for h in hits for f in h.features),
        str(hits[0].features) if hits else "no hits",
    ))

    s_hx, _ = score_protein(to_residues(ideal_helix()))
    checks.append(("all-helix protein yields no degron", s_hx == 0.0, f"{s_hx}"))

    noG = beta_hairpin(turn=(("SER", 1.2, 2.2), ("ALA", 2.6, 3.4), ("ALA", 4.0, 3.4), ("SER", 5.0, 2.2)))
    s_noG, _ = score_protein(to_residues(noG))
    checks.append(("G-loop turn outranks non-G turn", s_hp > s_noG, f"{s_hp} vs {s_noG}"))

    rsa = rsa_proxy(parse_pdb(ideal_helix()))
    checks.append(("lone helix reads exposed", sum(rsa) / len(rsa) > 0.6, f"{sum(rsa)/len(rsa):.2f}"))

    ok = 0
    for name, passed, detail in checks:
        print(f"  [{'PASS' if passed else 'FAIL'}] {name:36s} {detail}")
        ok += passed
    print(f"\n{ok}/{len(checks)} passed")
    return 0 if ok == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(run())
