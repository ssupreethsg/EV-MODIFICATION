#include <stdio.h>
#include <stdlib.h>
#include "spme_v7_apex_solver.h"

int main(void) {
    float soc;
    while (scanf("%f", &soc) == 1) {
        float ocv;
        ev_full_cell_ocv_V_and_dv_dsoc(soc, &ocv, NULL);
        printf("%f\n", ocv);
        fflush(stdout);
    }
    return 0;
}
