/* Reads one SOC (0..1) per line from stdin, writes the C model's full-cell OCV [V]. */
#include <stdio.h>
#include "spme_v7_apex_solver.h"

int main(void) {
    double soc;
    while (scanf("%lf", &soc) == 1) {
        float ocv = 0.0f;
        ev_full_cell_ocv_V_and_dv_dsoc((float)soc, &ocv, NULL);
        printf("%.6f\n", (double)ocv);
    }
    return 0;
}
