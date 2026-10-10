"""Checks the C OCV model.

1. The stoichiometry windows conserve lithium for the Chen2020 electrode capacities.
2. Full-cell OCV is 2.5 V at SOC 0 and 4.2 V at SOC 1 (the voltage limits defining SOC).
3. The tables in the C source equal the Chen2020 OCP functions (tools/generate_ocv_lut.py output).
4. Linear interpolation error on a dense grid is below its bound.
5. ev_cli ocv reproduces the Python interpolation and its slope (requires the built ev_cli).
6. If PyBaMM is installed, the windows match pybamm.lithium_ion.get_min_max_stoichiometries.

Usage: python tests/test_ocv.py [--cli PATH]
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools"),
)
import ev_model as m  # noqa: E402

# Chen2020 electrode loading (electrode height x width cancels in the ratio)
NEG_CAPACITY = 0.75 * 85.2e-6 * 33133.0  # active fraction * thickness * c_max
POS_CAPACITY = 0.665 * 75.6e-6 * 63104.0

FAILURES = []


def check(condition, message):
    print(("PASS: " if condition else "FAIL: ") + message)
    if not condition:
        FAILURES.append(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cli", help="path to ev_cli (default: $EV_CLI or ./build)")
    args = parser.parse_args()

    neg_li = NEG_CAPACITY * (m.ANODE_STOICH_SOC100 - m.ANODE_STOICH_SOC0)
    pos_li = POS_CAPACITY * (m.CATHODE_STOICH_SOC0 - m.CATHODE_STOICH_SOC100)
    check(
        abs(neg_li / pos_li - 1.0) < 1e-3,
        f"lithium conserved across the SOC window (ratio {neg_li / pos_li:.5f})",
    )

    ocv0 = float(m.full_cell_ocv_V(0.0))
    ocv1 = float(m.full_cell_ocv_V(1.0))
    check(abs(ocv0 - 2.5) < 1e-3, f"OCV(SOC=0) = {ocv0:.4f} V (2.5 V)")
    check(
        abs(ocv1 - m.MAX_CELL_VOLTAGE_V) < 1e-3,
        f"OCV(SOC=1) = {ocv1:.4f} V ({m.MAX_CELL_VOLTAGE_V} V)",
    )

    graphite, nmc = m.luts()
    table_error = max(
        np.max(np.abs(m.chen2020_graphite_ocp(m.anode_stoich(m.LUT_SOC)) - graphite)),
        np.max(np.abs(m.chen2020_nmc_ocp(m.cathode_stoich(m.LUT_SOC)) - nmc)),
    )
    check(
        table_error < 1e-6,
        f"C tables equal the Chen2020 OCP functions (max {table_error:.1e} V)",
    )

    fine = np.linspace(0.0, 1.0, 20001)
    interp_error = np.abs(m.full_cell_ocv_V(fine) - m.exact_full_cell_ocv_V(fine))
    check(
        np.max(interp_error) < 3e-3,
        f"interpolation error {np.max(interp_error) * 1000:.2f} mV < 3 mV everywhere",
    )
    check(
        np.max(interp_error[fine >= 0.1]) < 3e-4,
        f"interpolation error {np.max(interp_error[fine >= 0.1]) * 1000:.3f} mV < 0.3 mV above 10 % SOC",
    )

    try:
        cli = m.find_ev_cli(args.cli)
    except FileNotFoundError as error:
        print(f"FAIL: {error}")
        return 1
    check(cli is not None, "ev_cli found (build it, pass --cli or set EV_CLI)")
    if cli is not None:
        socs = np.linspace(0.0, 1.0, 1001)
        rows = np.array(m.run_cli(cli, ["ocv"], [f"{s:.9f}" for s in socs]), dtype=float)
        cli_error = float(np.max(np.abs(rows[:, 0] - m.full_cell_ocv_V(socs))))
        check(cli_error < 1e-5, f"ev_cli ocv matches Python (max {cli_error:.1e} V)")
        # dOCV/dSOC (the EKF Jacobian) is the slope of the table segment; test between nodes
        mid = (np.arange(m.LUT_POINTS - 1) + 0.5) / (m.LUT_POINTS - 1)
        slope_rows = np.array(m.run_cli(cli, ["ocv"], [f"{s:.9f}" for s in mid]), dtype=float)
        expected = np.diff(nmc - graphite) * (m.LUT_POINTS - 1)
        slope_error = float(np.max(np.abs(slope_rows[:, 1] - expected)))
        check(
            slope_error < 1e-3,
            f"ev_cli dOCV/dSOC matches the table slope (max {slope_error:.1e} V)",
        )

    try:
        import pybamm
    except ImportError:
        print("SKIP: PyBaMM not installed, window cross-check skipped")
    else:
        x0, x100, y100, y0 = pybamm.lithium_ion.get_min_max_stoichiometries(pybamm.ParameterValues("Chen2020"))
        window_error = max(
            abs(x0 - m.ANODE_STOICH_SOC0),
            abs(x100 - m.ANODE_STOICH_SOC100),
            abs(y0 - m.CATHODE_STOICH_SOC0),
            abs(y100 - m.CATHODE_STOICH_SOC100),
        )
        check(
            window_error < 1e-5,
            f"windows match PyBaMM (max difference {window_error:.1e})",
        )

    print("OCV TEST: " + ("PASS" if not FAILURES else f"FAIL ({len(FAILURES)})"))
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    sys.exit(main())
