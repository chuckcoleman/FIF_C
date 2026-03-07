/* 
    
    Collection of FFT functions needed by Fif.c
    
    This file is part of the C implementation of Fast Iterative Filtering (FIF)
    

    Authors: Igor Bertello, Emanuele Papini
    Affiliation(s): IAPS - INAF
    
    Revised: Chuck Coleman
    Affiliation: Timely Analytics, LLC

    Dependencies: FFTW3

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "FFT.h"

/*
 * Forward FFT of a real signal x of length N.
 * Returns an fftwl_complex array of length N allocated with fftwl_malloc.
 * Caller must free with fftwl_free().
 */
fftwl_complex *fft_dir(double *x, int N)
{
    if (x == NULL || N <= 0) {
        return NULL;
    }

    fftwl_complex *in  = (fftwl_complex *)fftwl_malloc(sizeof(fftwl_complex) * N);
    fftwl_complex *out = (fftwl_complex *)fftwl_malloc(sizeof(fftwl_complex) * N);

    if (in == NULL || out == NULL) {
        if (in)  fftwl_free(in);
        if (out) fftwl_free(out);
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        in[i][0] = (long double)x[i];
        in[i][1] = 0.0L;
    }

    fftwl_plan p = fftwl_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    if (p == NULL) {
        fftwl_free(in);
        fftwl_free(out);
        return NULL;
    }

    fftwl_execute(p);
    fftwl_destroy_plan(p);
    fftwl_free(in);

    return out;
}

/*
 * Inverse FFT of a complex signal x of length N.
 * Returns a real array of length N allocated with malloc.
 * Caller must free with free().
 */
double *fft_inv(fftwl_complex *x, int N)
{
    if (x == NULL || N <= 0) {
        return NULL;
    }

    fftwl_complex *in  = (fftwl_complex *)fftwl_malloc(sizeof(fftwl_complex) * N);
    fftwl_complex *out = (fftwl_complex *)fftwl_malloc(sizeof(fftwl_complex) * N);

    if (in == NULL || out == NULL) {
        if (in)  fftwl_free(in);
        if (out) fftwl_free(out);
        return NULL;
    }

    memcpy(in, x, sizeof(fftwl_complex) * N);

    fftwl_plan p = fftwl_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (p == NULL) {
        fftwl_free(in);
        fftwl_free(out);
        return NULL;
    }

    fftwl_execute(p);
    fftwl_destroy_plan(p);
    fftwl_free(in);

    double *y = (double *)malloc(sizeof(double) * N);
    if (y == NULL) {
        fftwl_free(out);
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        y[i] = (double)(out[i][0] / (long double)N);
    }

    fftwl_free(out);
    return y;
}

/*
 * Real part of FFT of a real signal x of length N.
 * Returns a real array of length N allocated with malloc.
 * Caller must free with free().
 */
double *realFFT(double *x, int N)
{
    if (x == NULL || N <= 0) {
        return NULL;
    }

    fftwl_complex *fx = fft_dir(x, N);
    if (fx == NULL) {
        return NULL;
    }

    double *y = (double *)malloc(sizeof(double) * N);
    if (y == NULL) {
        fftwl_free(fx);
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        y[i] = (double)fx[i][0];
    }

    fftwl_free(fx);
    return y;
}