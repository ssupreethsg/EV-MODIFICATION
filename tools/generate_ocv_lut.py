"""Regenerates the OCV lookup tables in src/ev_modify.c from the Chen2020 OCP
functions and the stoichiometry windows in the header.

    python tools/generate_ocv_lut.py           # rewrite the tables in place
    python tools/generate_ocv_lut.py --check   # exit 1 if the tables are out of date
"""

import argparse
import re
import sys

import ev_model as m

BEGIN = "/* BEGIN GENERATED OCV TABLES (tools/generate_ocv_lut.py) */"
END = "/* END GENERATED OCV TABLES */"


def _table(name, values):
    rows = []
    for start in range(0, len(values), 8):
        rows.append("    " + ", ".join(f"{v:.6f}f" for v in values[start : start + 8]))
    return f"static const float {name}[EV_LUT_POINTS] = {{\n" + ",\n".join(rows) + "\n};"


def generated_block():
    graphite = m.chen2020_graphite_ocp(m.anode_stoich(m.LUT_SOC))
    nmc = m.chen2020_nmc_ocp(m.cathode_stoich(m.LUT_SOC))
    return "\n".join(
        [
            BEGIN,
            f"#define EV_LUT_POINTS {m.LUT_POINTS}U",
            _table("LUT_GRAPHITE_OCP", graphite),
            _table("LUT_NMC_OCP", nmc),
            END,
        ]
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="exit 1 if the tables are out of date")
    args = parser.parse_args()
    with open(m.C_SOURCE, encoding="utf-8") as handle:
        source = handle.read()
    pattern = re.compile(re.escape(BEGIN) + r".*?" + re.escape(END), re.DOTALL)
    if not pattern.search(source):
        print(f"generated-table markers not found in {m.C_SOURCE}", file=sys.stderr)
        return 2
    updated = pattern.sub(lambda _: generated_block(), source)
    if args.check:
        if updated != source:
            print(
                "OCV tables are out of date: run python tools/generate_ocv_lut.py",
                file=sys.stderr,
            )
            return 1
        print("OCV tables are up to date")
        return 0
    with open(m.C_SOURCE, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(updated)
    print(f"Wrote {m.LUT_POINTS}-point OCV tables to {m.C_SOURCE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
