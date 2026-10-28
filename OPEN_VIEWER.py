"""Build the viewer and open it. Just: py OPEN_VIEWER.py"""

import pathlib
import sys
import webbrowser

if __name__ == "__main__":
    try:
        from tpd_bench.viewer import build
    except Exception as exc:
        print("Could not load the bench:", exc)
        print("Are you in the folder that contains the 'tpd_bench' directory?")
        sys.exit(1)
    p = build().resolve()
    print(f"Built {p.name}")
    print("(This is the same file GitHub Pages serves. If you changed anything")
    print(" the page shows, commit it.)")
    try:
        webbrowser.open(p.as_uri())
        print("Opening it in your browser.")
    except Exception:
        print(f"Open this file in a browser: {p}")
    print("\nIt opens on the bundled demo sweep, so there is real output to look at")
    print("straight away. To browse your own run instead, drop sweep_results.csv")
    print("onto the Nominations tab.")
