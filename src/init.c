#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* Declarations (prototypes) of the functions */

extern void get_EYs_C (double *xs, int *Ns, int *nmodel, 
                double *beta, double *sig2, double *ps, double *gamma, double *epsilon, double *EX,
                /* Sortie */ double *EYs, double *EXs, double *phis);
  
extern void strata_bh_opti_C(double *xnoc, int *Nnoc, double *bhfull, int *L, int *takenone, int *takeall,
                      int *Nc, double *EYc, double *q1, double *q2, double *q3,
                      int *nmodel, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, double *EX, double *EX2,
                      int *findn, int *n, double *CV, double *rhL, double *biaspenalty, int *takealladjust, int *dotests, int *minNh,
                      /* Sortie */ int *NhOK, int *nhOK, /* uniquement utile a l'algo de Kozak*/ 
                      double *phih, double *psih, double *gammah, double *ah, double *U2, double *U, double *V, 
                      /* uniquement utile a l'algo de Sethi*/
                      int *stratumIDnoc, int *Nh, double *EYh, double *VYh, double *TY, double *TAY,  
                      double *nhnonint, int *takeallout, double *nh, double *optinhnonint, double *optinh);

extern void get_RRMSE_C (double *biaspenalty, double *TY, double *TAY, int *Nh, double *VYh, double *nhcalcul, 
                  double *rhL, int *L,int *takenone,
                  /* sortie */ double *RRMSE);

extern void get_n_C (double *nhcalcul, int *L, int *Nc,
              /* sortie */ double *n);

extern void complete_enum_C(int *pbhsol, int *nsol, int *L, double *x1noc, int *N1noc, double *xnoc, int *Nnoc, 
                     int *takenone, int *takeall, int *Nc, double *EYc, double *q1, double *q2, double *q3,
                     int *nmodel, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, 
                     double *EX, double *EX2, int *findn, int *n, double *CV, double *rhL, 
                     double *biaspenalty, int *minNh,
                     /* Sortie */ double *soldetail);

extern void pbh2bhfull_C(int *pbh, int *L, double *x1noc, int *N1noc, double *bhfull);

extern void algo_Kozak_C(double *combin2try, int *ncombin, int *L, double *x1noc, int *N1noc, double *xnoc, int *Nnoc, 
                  int *takenone, int *takeall, int *Nc, double *EYc, double *q1, double *q2, double *q3,
                  int *nmodel, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, 
                  double *EX, double *EX2, int *findn, int *n, double *CV, double *rhL, 
                  double *biaspenalty, int *minNh, int *maxiter, int *idoptinh, int *rep,
                  /* Sortie */ double *rundetail, double *iterdetail, int *nrowiter);

extern void verif_takeall_C(double *nhnonint, int *Nh, int *L, int *takenone, 
                     /* Sortie et entree */ int *takeall,
                     /* Sortie */ int *valid);

extern void get_momentY_C(double *xnoc, int *stratumIDnoc, int *Nnoc, int *Nh, int *L, int *Nc, double *EYc, int *takenone,
                   int *nmodel, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, 
                   double *EX, double *EX2,
                   /* Sortie */ double *EYh, double *VYh, double *phih, double *psih, double *TY, double *TAY);
                       

/* Construction of the registration table */

static const R_CMethodDef CEntries[]  = {
  {"get_EYs_C", (DL_FUNC) &get_EYs_C, 12},
  {"strata_bh_opti_C", (DL_FUNC) &strata_bh_opti_C, 47},
  {"get_RRMSE_C", (DL_FUNC) &get_RRMSE_C, 10},
  {"get_n_C", (DL_FUNC) &get_n_C, 4},
  {"complete_enum_C", (DL_FUNC) &complete_enum_C, 29},
  {"pbh2bhfull_C", (DL_FUNC) &pbh2bhfull_C, 5},
  {"algo_Kozak_C", (DL_FUNC) &algo_Kozak_C, 34},
  {"verif_takeall_C", (DL_FUNC) &verif_takeall_C, 6},
  {"get_momentY_C", (DL_FUNC) &get_momentY_C, 22},
  {NULL, NULL, 0}
};

void attribute_visible R_init_stratification(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
