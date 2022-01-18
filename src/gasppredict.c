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

extern int ErrNum;

int CalcPred(const Matrix *X, const real *y,
             const LinModel *RegMod, const LinModel *SPMod,
             size_t CorFamNum, boolean RanErr, const Matrix *CorPar,
             real SPVar, real ErrVar, const Matrix *XPred, boolean GenPredCoefs,
             Matrix *YPred, real **PredCoef)
/* Arguments up to ErrVar are extracted from the R GaSP object.
Memory for outputs YPred and PredCoef has to be allocated by the
calling function.
*/
{
  int ErrReturn;
  KrigingModel KrigMod;
  real *ResTildeTilde, *SE, *yHat;
  size_t m;

  m = MatNumRows(XPred);

  ErrReturn = OK;
  /* Allocate and set up and kriging model. */
  KrigModAlloc(MatNumRows(X), MatNumCols(X), RegMod, SPMod, CorFamNum, RanErr, &KrigMod);
  KrigModData(MatNumRows(X), NULL, X, y, &KrigMod);
  ErrNum = KrigModSetUp(CorPar, SPVar, ErrVar, &KrigMod);

  if (ErrNum == OK)
  {
    yHat = AllocReal(m, NULL);
    SE = AllocReal(m, NULL);

    ErrNum = KrigPredSE(&KrigMod, XPred, yHat, SE);
    /*populate yPred */
    MatAlloc(m, 2, RECT, YPred);
    VecCopy(yHat, m, YPred->Elem[0]);
    VecCopy(SE, m, YPred->Elem[1]);
    AllocFree(yHat);
    AllocFree(SE);
  }

  if (GenPredCoefs && ErrNum == OK)
  {
    ResTildeTilde = AllocReal(MatNumRows(X), NULL);

    /* Compute prediction coefficients. */
    /* This seems to be unstable!       */
    ErrNum = TriBackSolve(KrigChol(&KrigMod),
                          KrigMod.ResTilde, ResTildeTilde);
    if (GenPredCoefs && ErrNum == OK)
    { /*populate pred coefs */
      *PredCoef = AllocReal(MatNumRows(X), NULL);
      VecCopy(ResTildeTilde, MatNumRows(X), *PredCoef);
    }
    AllocFree(ResTildeTilde);
  }

  KrigModFree(&KrigMod);

  if (ErrNum != OK)
    ErrReturn = ErrNum;

  return ErrReturn;
}

SEXP predict(SEXP reg_mod, SEXP sp_mod, SEXP ranErr, SEXP corFamNum,
             SEXP x_R, SEXP y_R, SEXP xPred, SEXP generate_coefficients,
             SEXP spVar, SEXP errVar, SEXP corpar)
{
  boolean RanErr = asLogical(ranErr);
  size_t CorFamNum = asInteger(corFamNum);
  real SPVar = asReal(spVar);
  real ErrVar = asReal(errVar);
  matrix X, XPred, CorPar, YPred;
  real *y, *PredCoef;
  string *xName;
  string *RegMod_Term, *SPMod_Term;
  LinModel RegMod;
  LinModel SPMod;
  boolean GenPredCoefs = asLogical(generate_coefficients);

  MatrixDFAlloc(&X, x_R);
  MatrixDFAlloc(&XPred, xPred);
  MatrixDFAlloc(&CorPar, corpar);
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
    StrFree(&RegMod_Term, (size_t)Rf_length(VECTOR_ELT(reg_mod, 0)));
    StrFree(&SPMod_Term, (size_t)Rf_length(VECTOR_ELT(sp_mod, 0)));
    StrFree(&xName, (size_t)Rf_length(getAttrib(x_R, R_NamesSymbol)));
    MatFree(&CorPar);
    MatFree(&YPred);
    MatFree(&X);
    MatFree(&XPred);
    ModFree(&RegMod);
    ModFree(&SPMod);
    Rf_error("Regression model and Stochastic Process model setup failed.");
  }

  int result = CalcPred(&X, y, &RegMod, &SPMod, CorFamNum, RanErr, &CorPar, SPVar, ErrVar, &XPred, GenPredCoefs, &YPred, &PredCoef);
  SEXP reslist = PROTECT(allocVector(VECSXP, 2));
  if (result == OK)
  {
    SEXP y_rowName = PROTECT(getAttrib(xPred, R_RowNamesSymbol));
    SEXP y_colName = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(y_colName, 0, mkChar("Pred"));
    SET_STRING_ELT(y_colName, 1, mkChar("SE"));
    SEXP yPreddf = MatrixDFConstructor(&YPred, y_rowName, y_colName);

    SET_VECTOR_ELT(reslist, 0, yPreddf);
    UNPROTECT(2);
    if (GenPredCoefs)
    {
      SEXP predcoefvec = RealVecConstructor(&PredCoef, MatNumRows(&X));
      SET_VECTOR_ELT(reslist, 1, predcoefvec);
      AllocFree(PredCoef);
    }
  }

  UNPROTECT(1);
  AllocFree(y);
  StrFree(&RegMod_Term, (size_t)Rf_length(VECTOR_ELT(reg_mod, 0)));
  StrFree(&SPMod_Term, (size_t)Rf_length(VECTOR_ELT(sp_mod, 0)));
  StrFree(&xName, (size_t)Rf_length(getAttrib(x_R, R_NamesSymbol)));
  MatFree(&CorPar);
  MatFree(&YPred);
  MatFree(&X);
  MatFree(&XPred);
  ModFree(&RegMod);
  ModFree(&SPMod);
  if (result != OK)
  {
    Rf_error("GaSP Predict failed.");
  }
  return reslist;
}
