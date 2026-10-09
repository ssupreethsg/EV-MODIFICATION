"""Verifies the C OCV lookup tables against the PyBaMM Chen2020 OCP functions,
and (if built) the C ocv_proxy against the same interpolation in Python."""
import subprocess
import sys

import numpy as np

import ev_model as m

LUT_TOLERANCE_V = 1.0e-4
PROXY_TOLERANCE_V = 1.0e-4


def main():
    sto_an = m.ANODE_STOICH_SOC0 + m.LUT_SOC * (m.ANODE_STOICH_SOC100 - m.ANODE_STOICH_SOC0)
    sto_ca = m.CATHODE_STOICH_SOC0 + m.LUT_SOC * (m.CATHODE_STOICH_SOC100 - m.CATHODE_STOICH_SOC0)
    graphite_err = np.max(np.abs(m.chen2020_graphite_ocp(sto_an) - m.LUT_GRAPHITE))
    nmc_err = np.max(np.abs(m.chen2020_nmc_ocp(sto_ca) - m.LUT_NMC))
    print(f"Graphite LUT max error: {graphite_err:.2e} V")
    print(f"NMC LUT max error:      {nmc_err:.2e} V")
    ok = graphite_err < LUT_TOLERANCE_V and nmc_err < LUT_TOLERANCE_V

    for soc in (0.0, 0.5, 0.8, 1.0):
        print(f"SOC {soc:.1f}: {float(m.full_cell_ocv_V(soc)):.3f} V")

    proxy = m.find_ocv_proxy()
    if proxy is None:
        print("ocv_proxy not built; skipping C proxy check")
    else:
        socs = np.linspace(0.0, 1.0, 257)
        out = subprocess.run([proxy], input="\n".join(f"{s:.9f}" for s in socs) + "\n",
                             text=True, capture_output=True, check=True).stdout.split()
        proxy_err = np.max(np.abs(np.array(out, dtype=float) - m.full_cell_ocv_V(socs)))
        print(f"C proxy vs Python max error: {proxy_err:.2e} V")
        ok = ok and len(out) == len(socs) and proxy_err < PROXY_TOLERANCE_V

    print("OCV TEST: PASS" if ok else "OCV TEST: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
