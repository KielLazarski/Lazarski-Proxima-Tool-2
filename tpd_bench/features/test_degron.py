"""Tests for the degron scorer, against cases with known answers.

Synthetic structures, not real AlphaFold models -- the sandbox cannot reach
AlphaFold DB. These test the SCORING LOGIC: that a canonical hairpin outranks
a non-canonical loop, that a buried degron scores zero, that a missing lysine
is penalised, and that low model confidence propagates. Replace with real
models before trusting any number on a real protein.
"""

from .degron import Residue, score_protein


def build(spec, plddt=95.0, rsa=0.5):
    """spec: list of (aa, ss) tuples."""
    return [Residue(i + 1, aa, ss, rsa, plddt) for i, (aa, ss) in enumerate(spec)]


def hairpin(turn="SGSG", lys=True, plddt=95.0, rsa=0.5):
    spec = [("V", "E")] * 4 + [("W", "E")]
    spec += [(c, "C") for c in turn]
    spec += [("V", "E")] * 4
    spec += [("K", "C")] if lys else [("A", "C")]
    return build(spec, plddt, rsa)


def run():
    checks = []

    s_g, _ = score_protein(hairpin("SGSG"))
    s_nog, _ = score_protein(hairpin("SASA"))
    checks.append(("G-loop scores above non-G turn", s_g > s_nog, f"{s_g} vs {s_nog}"))

    s_lys, _ = score_protein(hairpin("SGSG", lys=True))
    s_nolys, _ = score_protein(hairpin("SGSG", lys=False))
    checks.append(("missing lysine penalised", s_lys > s_nolys, f"{s_lys} vs {s_nolys}"))

    s_hi, _ = score_protein(hairpin("SGSG", plddt=95))
    s_lo, _ = score_protein(hairpin("SGSG", plddt=55))
    checks.append(("low pLDDT lowers score", s_hi > s_lo, f"{s_hi} vs {s_lo}"))

    s_exp, _ = score_protein(hairpin("SGSG", rsa=0.6))
    s_bur, _ = score_protein(hairpin("SGSG", rsa=0.02))
    checks.append(("buried degron scores zero", s_bur == 0.0 and s_exp > 0, f"{s_exp} vs {s_bur}"))

    # non-canonical loop should be found, but rank below a real hairpin
    loop = build([("A", "E")] * 2 + [("S", "C")] * 7 + [("K", "C")] + [("A", "E")] * 2)
    s_loop, hits = score_protein(loop)
    checks.append(("non-canonical loop detected", s_loop > 0, f"{s_loop}"))
    checks.append(("hairpin outranks non-canonical", s_g > s_loop, f"{s_g} vs {s_loop}"))

    flat = build([("A", "H")] * 30)
    s_flat, _ = score_protein(flat)
    checks.append(("all-helix protein scores zero", s_flat == 0.0, f"{s_flat}"))

    ok = 0
    for name, passed, detail in checks:
        print(f"  [{'PASS' if passed else 'FAIL'}] {name:38s} {detail}")
        ok += passed
    print(f"\n{ok}/{len(checks)} passed")
    return 0 if ok == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(run())
