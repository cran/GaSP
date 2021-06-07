#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "RCconvert.h"
#include "model.h"
#include "kriging.h"
#include "alex.h"

/* For debugging memory leakage */
extern size_t nPointers;

/* Used by optimizer */
extern real LambdaPrior;
extern size_t FitCritNum;
extern real SPVarPropMin;
extern real SPVarPropMax;

/* Used by correlation families */
/* Power-exponential only */
extern real AlphaMax; /* p > 1 */
extern real AlphaMin;
/* Matern only */
extern size_t derivMax; /* Codes infinity! */
extern size_t derivMin;
/* Matern and power-exponential */
extern real ThetaStandMax;
extern real ThetaStandMin;

extern int ErrNum;
boolean DesignJob = NO;

boolean isCorParNull;

int CalcFit(const Matrix *X, const real *y,
            const LinModel *RegMod, const LinModel *SPMod,
            size_t CorFamNum, boolean RanErr,
            real *SPVar, real *ErrVar, size_t Tries, real CritLogLikeDiff,
            real LogLikeTol, size_t ModCompCritNum, Matrix *CorPar, real **Beta, real **Summary)
{
  int ErrReturn;
  KrigingModel KrigMod;
  real NegLogLike;
  real CondNum, CVRootMSE;
  size_t m;

  m = MatNumRows(X);

  *Beta = AllocReal(ModDF(RegMod), NULL);
  *Summary = AllocReal(5, NULL);

  /* Perform a fit for each response. */
  ErrReturn = OK;
  /* Set up kriging model. */
  KrigModAlloc(m, MatNumCols(X), RegMod, SPMod, CorFamNum, RanErr, &KrigMod);
  KrigModData(m, NULL, X, y, &KrigMod);

  ProgressInit((int)Tries);
  ErrNum = FitBest(&KrigMod, Tries, CritLogLikeDiff, LogLikeTol, ModCompCritNum,
                   *Beta, CorPar,
                   SPVar, ErrVar, &NegLogLike,
                   &CVRootMSE, &CondNum, m, RanErr);
  if (ErrNum == OK)
  {
    Summary[0][0] = -NegLogLike;
    Summary[0][1] = CondNum;
    Summary[0][2] = CVRootMSE;
    Summary[0][3] = SPVar[0];
    Summary[0][4] = ErrVar[0];
  }
  KrigModFree(&KrigMod);

  if (ErrNum != OK)
    ErrReturn = ErrNum;

  return ErrReturn;
}

int FitBest(KrigingModel *KrigMod, size_t Tries, real CritLogLikeDiff,
            real LogLikeTol, size_t ModCompCritNum, real *Beta,
            Matrix *CorPar, real *SPVar, real *ErrVar, real *NegLogLike,
            real *CVRootMSE, real *CondNum, size_t nCasesXY, boolean RanErr)
{
  boolean Better;
  int ErrNum, ErrThisTry;
  Matrix RegCorPar;
  real CondNumTry, CVRootMSETry, MaxErr;
  real NegLogLikeTry;
  real *YHatCV;
  size_t IndexMaxErr;
  unsigned nEvalsTry;

  YHatCV = AllocReal(nCasesXY, NULL);
  ErrNum = !OK;
  *CondNum = NA_REAL;
  *CVRootMSE = REAL_MAX;
  *NegLogLike = REAL_MAX;
  for (size_t j = 0; j < Tries; j++)
  {
    MLEStart(KrigMod, &RegCorPar);

    if (j == 0)
    {
      /* First try: If SPModMat contains correlation   */
      /* parameters, then use them as starting values. */
      if (!isCorParNull)
      {
        MatCopy(CorPar, KrigCorPar(KrigMod));
      }

      if (RanErr && *SPVar != NA_REAL && *ErrVar != NA_REAL)
      {
        KrigMod->SPVarProp = *SPVar / (*SPVar + *ErrVar);

        /* User-supplied variances must satisfy nugget condition */
        if (KrigMod->SPVarProp > SPVarPropMax)
          KrigMod->SPVarProp = SPVarPropMax;
      }
    }

    ErrThisTry = MLEFit(&RegCorPar, KrigMod, LogLikeTol,
                        CritLogLikeDiff, j + 1, &NegLogLikeTry,
                        &CondNumTry, &nEvalsTry);
    MatFree(&RegCorPar);

    Better = FALSE;
    if (ErrThisTry == OK &&
        (ErrThisTry = CalcCV(KrigMod, YHatCV, NULL)) == OK)
    {
      CVRootMSETry = RootMSE(nCasesXY, YHatCV, KrigY(KrigMod),
                             &MaxErr, &IndexMaxErr);
      switch (ModCompCritNum)
      {
      case MOD_COMP_CRIT_CV:
        if (CVRootMSETry < *CVRootMSE)
          Better = TRUE;
        break;

      case MOD_COMP_CRIT_LIKE:
        if (NegLogLikeTry < *NegLogLike)
          Better = TRUE;
        break;

      default:
        CodeBug(ILLEGAL_COND_TXT);
      }
    }

    if (ErrThisTry == OK && Better)
    {
      /* One good try is sufficient. */
      ErrNum = OK;

      /* Best parameters so far. */

      *CVRootMSE = CVRootMSETry;
      VecCopy(KrigMod->Beta,
              ModDF(KrigRegMod(KrigMod)), Beta);

      MatCopy(KrigCorPar(KrigMod), CorPar);

      *SPVar = KrigMod->SigmaSq * KrigMod->SPVarProp;
      *ErrVar = KrigMod->SigmaSq * (1.0 - KrigMod->SPVarProp);
      *NegLogLike = NegLogLikeTry;
      *CondNum = CondNumTry;
    }
    tick(1);
  }

  AllocFree(YHatCV);

  return ErrNum;
}

SEXP fit(SEXP x_R, SEXP y_R, SEXP reg_mod, SEXP sp_mod,
         SEXP corFamNum, SEXP ranErr,
         SEXP corpar, SEXP spVar, SEXP errVar, SEXP nugget,
         SEXP tries_R, SEXP seed_R, SEXP fit_objective,
         SEXP theta_standardized_min, SEXP theta_standardized_max, SEXP alpha_min, SEXP alpha_max,
         SEXP derivatives_min, SEXP derivatives_max,
         SEXP logLikeTol, SEXP critLogLikeDiff,
         SEXP lambda_prior, SEXP mod_comp_num)
{
  size_t CorFamNum = asInteger(corFamNum);
  boolean RanErr = asLogical(ranErr);
  real SPVar = asReal(spVar);
  real ErrVar = asReal(errVar);
  SPVarPropMax = 1.0 - asReal(nugget);
  size_t Tries = asInteger(tries_R);
  int Seed = asInteger(seed_R);
  RandInit(Seed, Seed, Seed);
  FitCritNum = asInteger(fit_objective);
  ThetaStandMax = asReal(theta_standardized_max);
  ThetaStandMin = asReal(theta_standardized_min);
  AlphaMax = asReal(alpha_max);
  AlphaMin = asReal(alpha_min);
  derivMax = asInteger(derivatives_max);
  derivMin = asInteger(derivatives_min);
  real LogLikeTol = asReal(logLikeTol);
  real CritLogLikeDiff = asReal(critLogLikeDiff);
  LambdaPrior = asReal(lambda_prior);
  size_t ModCompCritNum = asInteger(mod_comp_num);
  matrix X, CorPar;
  real *y, *Beta;
  string *xName;
  LinModel RegMod;
  LinModel SPMod;
  string *RegMod_Term, *SPMod_Term;
  real *Summary;

  MatrixDFAlloc(&X, x_R);
  if (Rf_length(corpar) == 1)
  {
    isCorParNull = 1;
    MatAlloc(MatNumCols(&X), 2, RECT, &CorPar);
  }
  else
  {
    isCorParNull = 0;
    MatrixDFAlloc(&CorPar, corpar);
  }
  if (SPVar < 0)
    SPVar = NA_REAL;
  if (ErrVar < 0)
    ErrVar = NA_REAL;
  RealVecAlloc(&y, y_R);
  RegModDFAlloc(&RegMod_Term, reg_mod);
  RegModDFAlloc(&SPMod_Term, sp_mod);
  GetColName(&xName, x_R);

  /* Set up RegMod and SPMod */
  ErrNum = ModParse1((size_t)Rf_length(VECTOR_ELT(reg_mod, 0)), RegMod_Term, "RegressionModel", &RegMod);
  if (ErrNum == OK)
  {
    ErrNum = ModParse2(MatNumCols(&X), xName, NULL, "RegressionModel", &RegMod);
  }
  if (ErrNum == OK)
  {
    ErrNum = ModParse1((size_t)Rf_length(VECTOR_ELT(sp_mod, 0)), SPMod_Term, "StochasticProcessModel", &SPMod);
  }
  if (ErrNum == OK)
  {
    ErrNum = ModParse2(MatNumCols(&X), xName, NULL, "StochasticProcessModel", &SPMod);
  }
  else
  {
    AllocFree(y);
    AllocFree(RegMod_Term);
    AllocFree(SPMod_Term);
    AllocFree(xName);
    MatFree(&X);
    MatFree(&CorPar);
    ModFree(&RegMod);
    ModFree(&SPMod);
    Rf_error("Regression model and Stochastic Process model setup failed.");
  }

  int result = CalcFit(&X, y, &RegMod, &SPMod, CorFamNum, RanErr,
                       &SPVar, &ErrVar,
                       Tries, CritLogLikeDiff, LogLikeTol, ModCompCritNum, &CorPar, &Beta, &Summary);
  SEXP results = PROTECT(allocVector(VECSXP, 3));
  if (result == OK)
  {

    SEXP summary = RealVecConstructor(&Summary, 5);
    SET_VECTOR_ELT(results, 0, summary);
    SEXP beta = RealVecConstructor(&Beta, ModDF(&RegMod));
    SET_VECTOR_ELT(results, 1, beta);

    SEXP y_rowName = VECTOR_ELT(sp_mod, 0);
    SEXP y_colName = PROTECT(allocVector(STRSXP, 2));
    if (CorFamNum)
    {
      SET_STRING_ELT(y_colName, 0, mkChar("Theta"));
      SET_STRING_ELT(y_colName, 1, mkChar("Derivatives"));
    }
    else
    {
      SET_STRING_ELT(y_colName, 0, mkChar("Theta"));
      SET_STRING_ELT(y_colName, 1, mkChar("Alpha"));
    }
    SEXP CorPardf = MatrixDFConstructor(&CorPar, y_rowName, y_colName);
    SET_VECTOR_ELT(results, 2, CorPardf);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  AllocFree(y);
  AllocFree(RegMod_Term);
  AllocFree(SPMod_Term);
  AllocFree(xName);
  AllocFree(Beta);
  AllocFree(Summary);
  MatFree(&X);
  MatFree(&CorPar);
  ModFree(&RegMod);
  ModFree(&SPMod);
  if (result != OK)
  {
    Rf_error("GaSP Fit failed.");
  }
  return results;
}
