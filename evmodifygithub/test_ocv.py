import numpy as np

def graphite(sto):
    return (
        1.9793 * np.exp(-39.3631 * sto)
        + 0.2482
        - 0.0909 * np.tanh(29.8538 * (sto - 0.1234))
        - 0.04478 * np.tanh(14.9159 * (sto - 0.2769))
        - 0.0205 * np.tanh(30.4444 * (sto - 0.6103))
    )

def nmc(sto):
    return (
        -0.8090 * sto
        + 4.4875
        - 0.0428 * np.tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * np.tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * np.tanh(15.9308 * (sto - 0.3120))
    )

x_an_0 = 0.0279
x_an_100 = 0.9014
x_ca_0 = 0.9287
x_ca_100 = 0.2479

for soc in [0.0, 0.5, 0.8, 1.0]:
    sto_an = x_an_0 + soc * (x_an_100 - x_an_0)
    sto_ca = x_ca_0 + soc * (x_ca_100 - x_ca_0)
    ocv = nmc(sto_ca) - graphite(sto_an)
    print(f"SOC {soc:.1f}: {ocv:.3f} V")
