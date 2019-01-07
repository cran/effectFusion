// Copyright 2016 Nick Polson, James Scott, and Jesse Windle.
// Extracted from the R package "BayesLogit" (Version: 0.6).

////////////////////////////////////////////////////////////////////////////////

#ifdef USE_R
#include "R.h"
#include "Rmath.h"
#endif

#include "LogitWrapper.h"
#include "Logit.h"
#include "RNG.h"
#include "PolyaGamma.h"
#include "PolyaGammaAlt.h"
#include "PolyaGammaSP.h"
#include <exception>
#include <stdio.h>

using std::pow;
using std::fabs;
using std::sqrt;
using std::log;
using std::exp;

////////////////////////////////////////////////////////////////////////////////
// PolyaGamma //
////////////////////////////////////////////////////////////////////////////////


void rpg_hybrid(double *x, double *h, double *z, int* num)
{
    RNG r;
    PolyaGamma dv;
    PolyaGammaAlt alt;
    PolyaGammaSP sp;

#ifdef USE_R
    GetRNGstate();
#endif

    for(int i=0; i < *num; ++i){
        double b = h[i];
        if (b > 170) {
            double m = dv.pg_m1(b,z[i]);
            double v = dv.pg_m2(b,z[i]) - m*m;
            x[i] = r.norm(m, sqrt(v));
        }
        else if (b > 13) {
            sp.draw(x[i], b, z[i], r);
        }
        else if (b==1 || b==2) {
            x[i] = dv.draw((int)b, z[i], r);
        }
        else if (b > 1) {
            x[i] = alt.draw(b, z[i], r);
        }
        else if (b > 0) {
            x[i] = dv.draw_sum_of_gammas(b, z[i], r);
        }
        else {
            x[i] = 0.0;
        }
    }

#ifdef USE_R
    PutRNGstate();
#endif
}


////////////////////////////////////////////////////////////////////////////////
// POSTERIOR INFERENCE //
////////////////////////////////////////////////////////////////////////////////

// Posterior by Gibbs.
//------------------------------------------------------------------------------
void gibbs(double *wp, double *betap,                            // Posterior
           double *yp, double *tXp, double *np,                  // Data
           double *m0p, double *P0p,                             // Prior
           int *N, int *P,                                       // Dim
           int *samp, int *burn)                                 // MCMC
{

    // Set up data.
    Matrix y (yp,  *N,  1);
    Matrix tX(tXp, *P, *N);
    Matrix n (np,  *N,  1);
    Matrix m0(m0p, *P,  1);
    Matrix P0(P0p, *P, *P);

    // Declare posteriors.
    Matrix w, beta;

    // Random number generator.
    RNG r;

#ifdef USE_R
    GetRNGstate();
#endif

    // Logit Gibbs
    try{
        Logit logit(y, tX, n);
        logit.set_prior(m0, P0);
        // logit.compress();

        // Set the correct dimensions after combining data.
        w.resize(logit.get_N(), 1, *samp);
        beta.resize(logit.get_P(), 1, *samp);
        MatrixFrame w_mf    (wp   , w.rows()   , w.cols()   , w.mats());
        MatrixFrame beta_mf (betap, beta.rows(), beta.cols(), beta.mats());

        // Copy values to test code.  Must be using known value with correct dim.
        // w.copy(w_mf);
        // beta.copy(beta_mf);

        // Run gibbs.
        logit.gibbs(w, beta, *samp, *burn, r);

        // Copy values to return.
        w_mf.copy(w);
        beta_mf.copy(beta);

        // Adjust for combined data.
        *N = w.rows();
    }
    catch (std::exception& e) {
        Rprintf("Error: %s\n", e.what());
        Rprintf("Aborting Gibbs sampler.\n");
    }

#ifdef USE_R
    PutRNGstate();
#endif
} // gibbs



////////////////////////////////////////////////////////////////////////////////

// combine_data
//------------------------------------------------------------------------------
void combine(double *yp, double *tXp, double *np,                  // Data
             int *N, int *P)
{

    // Set up data.
    Matrix y (yp,  *N,  1);
    Matrix tX(tXp, *P, *N);
    Matrix n (np,  *N,  1);

    // Logit Gibbs
    try{
        Logit logit(y, tX, n);
        logit.compress();
        logit.get_data(y, tX, n);

        // Copy.
        MatrixFrame y_mf    (yp ,  y.rows(),  y.cols(),  y.mats());
        MatrixFrame tX_mf   (tXp, tX.rows(), tX.cols(), tX.mats());
        MatrixFrame n_mf    (np ,  n.rows(),  n.cols(),  n.mats());

        y_mf.copy(y);
        tX_mf.copy(tX);
        n_mf.copy(n);

        *N = tX.cols();
    }
    catch (std::exception& e) {
        Rprintf("Error: %s\n", e.what());
        Rprintf("Aborting combine.\n");
    }

} // combine_data

