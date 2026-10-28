"""Run the sweep. Just: py RUN_SWEEP.py

Edit LIMIT below if you want a different size.
Ctrl+C is safe at any point -- work is cached and resumes next run.
"""

import sys

LIMIT = 500      # proteins to scan. None = whole proteome (hours).
GENES = None     # or e.g. "STAT6,IRF5,PNPLA3" to test specific genes

if __name__ == "__main__":
    try:
        from tpd_bench.features import sweep
    except Exception as exc:
        print("Could not load the bench:", exc)
        print("Are you in the folder that contains the 'tpd_bench' directory?")
        sys.exit(1)

    genes = [g.strip() for g in GENES.split(",")] if GENES else None
    print(f"Running sweep (limit={LIMIT}). Ctrl+C is safe.\n")
    code = sweep.run(LIMIT, None, genes)
    print("\nDone. Send me sweep_results.csv")
    sys.exit(code)
