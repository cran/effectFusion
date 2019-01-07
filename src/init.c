#include <R.h>
#include <R_ext/Rdynload.h>

extern void combine(double *yp, double *tXp, double *np, int *N, int *P);
extern void gibbs(double *wp, double *betap, double *yp, double *tXp,
		  double *np, double *m0p, double *P0p, int *N, int *P, 
		  int *samp, int *burn);                                
extern void rpg_hybrid(double *x, double *h, double *z, int* num);

#define CDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CMethodDef cMethods[] = {
    CDEF(combine, 5),
    CDEF(gibbs, 11),
    CDEF(rpg_hybrid, 4),
    {NULL, NULL, 0}
};

void R_init_effectFusion(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}
