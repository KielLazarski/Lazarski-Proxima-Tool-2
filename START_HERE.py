"""Run this file. It does everything and writes one report.

    python3 START_HERE.py

Then send me report.txt.
"""

import io
import platform
import sys
import traceback
from contextlib import redirect_stdout

OUT = "report.txt"


def section(title, fn, log):
    log.write("\n" + "=" * 70 + "\n")
    log.write(title + "\n")
    log.write("=" * 70 + "\n")
    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            code = fn()
        log.write(buf.getvalue())
        # These return an exit code rather than raising, so a non-zero result
        # is a failure even though nothing blew up.
        return code in (0, None)
    except Exception:
        log.write(buf.getvalue())
        log.write("\n--- this step hit an error ---\n")
        log.write(traceback.format_exc())
        return False


def main():
    log = io.StringIO()
    log.write(f"python {platform.python_version()} on {platform.system()}\n")

    if sys.version_info < (3, 9):
        print("You need Python 3.9 or newer. Yours is", platform.python_version())
        print("That's the only problem -- install a newer Python and run this again.")
        return 2

    steps = []
    try:
        from tpd_bench import ablate, calibrate
        from tpd_bench.corpus import build as corpus_build
        from tpd_bench.features import doctor, test_degron, test_structure

        steps = [
            ("1. Does the bench score correctly?", calibrate.run),
            ("2. Which scoring terms matter?", ablate.run),
            ("3. Scaffolding corpus", corpus_build.run),
            ("4. Degron scorer tests", test_degron.run),
            ("4b. Structure parsing tests", test_structure.run),
            ("5. Can we reach the data sources?", doctor.run),
        ]
    except Exception:
        log.write("\nCould not load the bench:\n" + traceback.format_exc())
        log.write("\nMost likely you are not in the folder that contains 'tpd_bench'.\n")
        with open(OUT, "w", encoding="utf-8") as fh:
            fh.write(log.getvalue())
        print("Something went wrong loading the code. Details are in", OUT)
        print("Send me that file.")
        return 1

    print("Running 5 checks. This takes a minute.\n")
    results = []
    for title, fn in steps:
        print("  " + title)
        results.append(section(title, fn, log))

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(log.getvalue())

    offline_ok = all(results[:-1])
    online_ok = results[-1]

    print("\n" + "-" * 60)
    if offline_ok:
        print("The bench itself works.")
    else:
        print("Something is wrong with the bench itself.")
    if online_ok:
        print("Data sources are reachable. Good to go.")
    else:
        print("Data sources are NOT reachable from this machine.")
        print("That is the thing to fix, and it may not be fixable from here.")
    print("-" * 60)
    print(f"\nWritten to {OUT} -- send me that file.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
