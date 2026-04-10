/* 
    Main file for testing the C implementation of Fast Iterative Filtering (FIF)
    
    written based on  https://github.com/Acicone/FIF/blob/master/Example_v8.m

    Authors: Igor Bertello, Emanuele Papini
    Affiliation(s): IAPS - INAF, University of L'Aquila (Italy)
    
    Revised: Chuck Coleman
    Affiliation: Timely Analytics, LLC

*/

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include "Fif.h"

#define lung 1000
#define maxIMF 10 

static void free_fif_result(Fif_t *head)
{
    if (head == NULL) return;

    Fif_t *node = head->next;
    free(head->dati);
    head->dati = NULL;
    head->next = NULL;

    while (node != NULL) {
        Fif_t *next = node->next;
        free(node->dati);
        free(node);
        node = next;
    }
}

int main(int varn, char *varc[]) {
    (void)varn;
    (void)varc;

    Fif_t IMFs;
    Fif_t *IMFappo;
    double f[lung];
    double t, dt;
    double pi;
    int numIMF = maxIMF;
    char fname[32];
    int exit_code = 0;
    FILE *fil = NULL;

    memset(&IMFs, 0, sizeof(IMFs));

    dt = 0.001; 
    pi = acos(-1.0);

    // create the signal to analyse
    srand((unsigned)time(NULL));
    for (int i = 0; i < lung; i++) {
        t = i * dt;
        f[i] = (2 * pow(t - 0.5, 2) + 0.2) * sin(20 * pi * t + 0.2 * cos(40 * pi * t)) +
               4 * pow(t - 0.5, 2) + 1;
        // f[i] = 5*sin(i*dt*pi/500) + 5*sin(i*pi/200);
    }

    snprintf(fname, sizeof(fname), "signal.dat");
    fil = fopen(fname, "w");
    if (fil == NULL) {
        perror("fopen signal.dat");
        return 1;
    }
    for (int i = 0; i < lung; i++) {
        fprintf(fil, "%f\n", f[i]);
    }
    fclose(fil);
    fil = NULL;


    // RUN FIF
    IMFs = FIF_v2_1(f, lung, &numIMF);

    /* Summary prints suppressed for cleaner Valgrind output on macOS. */

    IMFappo = &IMFs;

    for (int j = 0; j < numIMF && IMFappo != NULL; j++) {
        snprintf(fname, sizeof(fname), "IMF%d.dat", j + 1);
        fil = fopen(fname, "w");
        if (fil == NULL) {
            perror("fopen IMF output");
            exit_code = 1;
            goto cleanup;
        }

        if (IMFappo->dati != NULL) {
            for (int i = 0; i < lung; i++) {
                fprintf(fil, "%f\n", IMFappo->dati[i]);
            }
        }

        fclose(fil);
        fil = NULL;
        IMFappo = IMFappo->next;
    }

    // Write residual only if there is an additional node
    if (IMFappo != NULL && IMFappo->dati != NULL) {
        snprintf(fname, sizeof(fname), "residual.dat");
        fil = fopen(fname, "w");
        if (fil == NULL) {
            perror("fopen residual.dat");
            exit_code = 1;
            goto cleanup;
        }

        for (int i = 0; i < lung; i++) {
            fprintf(fil, "%f\n", IMFappo->dati[i]);
        }
        fclose(fil);
        fil = NULL;
    }

cleanup:
    if (fil != NULL) {
        fclose(fil);
    }
    free_fif_result(&IMFs);
    fftwl_cleanup();
    return exit_code;
}
