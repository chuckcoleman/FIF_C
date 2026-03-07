#ifndef FFT_H
#define FFT_H

#include <fftw3.h>

fftwl_complex *fft_dir(double *x, int N);
double *fft_inv(fftwl_complex *x, int N);
double *realFFT(double *x, int N);

#endif