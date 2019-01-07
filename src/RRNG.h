/// Copyright 2016 Nick Polson, James Scott, and Jesse Windle.
// Extracted from the R package "BayesLogit" (Version: 0.6).

//////////////////////////////////////////////////////////////////////

// YOU MUST ALWAYS CALL GetRNGSeed() and PutRNGSeed() WHEN USING THESE FUNCTIONS!!!

//////////////////////////////////////////////////////////////////////

#ifndef __BASICRNG__
#define __BASICRNG__

#include "R.h"
#include "Rmath.h"
// #include "Matrix.h"

class BasicRNG {

 public:

  // Random variates.
  double unif  ();                             // Uniform
  double expon_mean(double mean);                  // Exponential
  double expon_rate(double rate);                  // Exponential
  double chisq (double df);                    // Chisq
  double norm  (double sd);                    // Normal
  double norm  (double mean , double sd);      // Normal
  double gamma_scale (double shape, double scale); // Gamma_Scale
  double gamma_rate  (double shape, double rate);  // Gamma_Rate
  double igamma(double shape, double scale);   // Inv-Gamma
  double flat  (double a=0  , double b=1  );   // Flat
  double beta  (double a=1.0, double b=1.0);   // Beta

  int bern  (double p);                     // Bernoulli

  // CDF
  static double p_norm (double x, int use_log=0);
  static double p_gamma_rate(double x, double shape, double rate, int use_log=0);

  // Density
  static double d_beta(double x, double a, double b);

  // Utility
  static double Gamma (double x, int use_log=0);

}; // BasicRNG

#endif
