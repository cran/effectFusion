// Copyright 2016 Nick Polson, James Scott, and Jesse Windle.
// Extracted from the R package "BayesLogit" (Version: 0.6).

////////////////////////////////////////////////////////////////////////////////

#ifndef __LOGITWRAPPER__
#define __LOGITWRAPPER__

extern "C" {

    // RPG

    void rpg_hybrid(double *x, double *h, double *z, int* num);

    // Default Logistic

    void gibbs(double *wp, double *betap,                            // Posterior
               double *yp, double *tXp, double *np,                  // Data
               double *m0p, double *P0p,                             // Prior
               int *N, int *P,                                       // Dim
               int *samp, int *burn);                                // MCMC


    void combine(double *yp, double *tXp, double *np,
                 int *N, int *P);

}

#endif
