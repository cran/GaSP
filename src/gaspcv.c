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

int ErrNum;

int CVHelper(const Matrix *X, const real *y,
             const LinModel *RegMod, const LinModel *SPMod,
             size_t CorFamNum, boolean RanErr, const Matrix *CorPar,
             real SPVar, real ErrVar,
             Matrix *CV)
{
  int ErrReturn;
  KrigingModel KrigMod;
  real *SE, *YHatCV;
  size_t m;

  m = MatNumRows(X);

  /* Compute cross-validation predictions for each response. */
  ErrReturn = OK;

  /* Allocations. */
  YHatCV = AllocReal(m, NULL);
  SE = AllocReal(m, NULL);

  /* Set up kriging model. */
  KrigModAlloc(m, MatNumCols(X), RegMod, SPMod, CorFamNum, RanErr, &KrigMod);
  KrigModData(m, NULL, X, y, &KrigMod);

  /* SPModMat contains the correlation parameters. */
  ErrNum = KrigModSetUp(CorPar, SPVar, ErrVar, &KrigMod);

  if (ErrNum == OK)
    ErrNum = CalcCV(&KrigMod, YHatCV, SE);

  if (ErrNum == OK)
  {
    /* Put cross validations and standard errors */
    /* in CV.                                    */

    MatAlloc(m, 2, RECT, CV);
    VecCopy(YHatCV, m, CV->Elem[0]);
    VecCopy(SE, m, CV->Elem[1]);
  }

  KrigModFree(&KrigMod);

  if (ErrNum != OK)
    ErrReturn = ErrNum;

  AllocFree(YHatCV);
  AllocFree(SE);

  return ErrReturn;
}

int CalcCV(KrigingModel *KrigMod, real *YHatCV, real *SE)
{
  int ErrNum;
  Matrix C, FTilde;
  Matrix *Chol, *F, *Q, *R;
  real c, s, t;
  real *Col, *Beta, *f, *r, *RBeta, *ResTilde, *Y, *YTilde;
  size_t i, ii, j, k, m, n;

  Y = KrigY(KrigMod);
  F = KrigF(KrigMod);
  Chol = KrigChol(KrigMod);
  Q = KrigQ(KrigMod);
  R = KrigR(KrigMod);

  /* Use workspace in KrigMod. */
  f = KrigMod->fRow;
  r = KrigMod->r;
  RBeta = KrigMod->RBeta;
  Beta = KrigMod->Beta;
  ResTilde = KrigMod->ResTilde;

  n = MatNumRows(F);
  k = MatNumCols(F);

  if (n == 0)
    return OK;
  else if (n == 1)
  {
    YHatCV[0] = NA_REAL;
    if (SE != NULL)
      SE[0] = NA_REAL;
    return OK;
  }

  MatAlloc(n, n, UP_TRIANG, &C);
  MatAlloc(n, k, RECT, &FTilde);
  YTilde = AllocReal(n, NULL);

  MatPutNumRows(Q, n - 1);

  /* Put correlation matrix in C. */
  KrigCorMat(0, NULL, KrigMod);
  MatCopy(Chol, &C);

  /* Overwrite correlation matrix with Cholesky decomposition. */
  if (TriCholesky(Chol, 0, Chol) != OK)
  {
    Error("Ill-conditioned Cholesky factor.\n");
    ErrNum = NUMERIC_ERR;
  }
  else
    ErrNum = OK;

  /* Compute FTilde and YTilde for all n rows. */
  if (ErrNum == OK)
    ErrNum = KrigSolve(Chol, F, Y, &FTilde, YTilde);

  /* Delete case i and predict Y[i]. */
  for (i = n - 1, ii = 0; ii < n && ErrNum == OK; ii++, i--)
  {

    /* Permute adjacent columns of Chol until  */
    /* column i is moved to the last column.   */
    for (j = i; j < n - 1; j++)
    {
      TriPerm(j, j + 1, Chol, &c, &s);

      /* Apply the same rotation to YTilde and FTilde. */
      t = c * YTilde[j] + s * YTilde[j + 1];
      YTilde[j + 1] = -s * YTilde[j] + c * YTilde[j + 1];
      YTilde[j] = t;
      for (m = 0; m < k; m++)
      {
        Col = MatCol(&FTilde, m);
        t = c * Col[j] + s * Col[j + 1];
        Col[j + 1] = -s * Col[j] + c * Col[j + 1];
        Col[j] = t;
      }
    }

    /* Correlations between case i and the other cases.  */
    /* Note that cases after i are now in reverse order. */
    for (j = 0; j < i; j++)
      r[j] = MatElem(&C, j, i);
    for (j = 0; j < n - 1 - i; j++)
      r[i + j] = MatElem(&C, i, n - 1 - j);

    /* Linear model terms for case i. */
    MatRow(F, i, f);

    /* Pretend we have only n - 1 cases. */
    MatPutNumRows(Chol, n - 1);
    MatPutNumCols(Chol, n - 1);
    MatPutNumRows(&FTilde, n - 1);

    /* Gram-Schmidt QR orthogonalization of FTilde. */
    if (QRLS(&FTilde, YTilde, Q, R, RBeta, ResTilde) != OK)
    {
      Error("Cannot perform QR decomposition.\n");
      ErrNum = NUMERIC_ERR;
    }

    else
    {
      /* Leave-one-out beta's can be obtained as follows. */
      /*
               if (TriBackSolve(R, RBeta, Beta) != OK)
                    Error("Cannot compute regression beta's.\n");
               else
               {
                    for (j = 0; j < k; j++)
                         Output(" %e", Beta[j]);
                    Output("\n");
               }
               */

      if (SE != NULL)
      {
        /* Standard error required.             */
        /* KrigMod->SigmaSq is not updated.     */
        /* RAve = 1.0 for epsilon contribution. */
        ErrNum = KrigYHatSE(KrigMod, 1.0, f, r,
                            &YHatCV[i], &SE[i]);
      }
      else
        /* No standard error. */
        ErrNum = KrigYHatSE(KrigMod, 1.0, f, r,
                            &YHatCV[i], NULL);
    }

    /* Restore sizes of Chol and FTilde. */
    MatPutNumRows(Chol, n);
    MatPutNumCols(Chol, n);
    MatPutNumRows(&FTilde, n);
  }

  if (ErrNum != OK)
    for (i = 0; i < n; i++)
      YHatCV[i] = SE[i] = NA_REAL;

  MatPutNumRows(Q, n);

  MatFree(&C);
  MatFree(&FTilde);
  AllocFree(YTilde);

  return ErrNum;
}

SEXP crossvalidate(SEXP reg_mod, SEXP sp_mod, SEXP ranErr, SEXP corFamNum,
                   SEXP x_R, SEXP y_R, SEXP spVar, SEXP errVar, SEXP corpar)
{
  boolean RanErr = asLogical(ranErr);
  size_t CorFamNum = asInteger(corFamNum);
  real SPVar = asReal(spVar);
  real ErrVar = asReal(errVar);
  matrix X, CorPar, CV;
  real *y;
  string *xName;
  string *RegMod_Term, *SPMod_Term;
  LinModel RegMod;
  LinModel SPMod;

  MatrixDFAlloc(&X, x_R);
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
    MatFree(&X);
    MatFree(&CorPar);
    ModFree(&RegMod);
    ModFree(&SPMod);
    Rf_error("Regression model and Stochastic Process model setup failed.");
  }
  int result = CVHelper(&X, y, &RegMod, &SPMod, CorFamNum, RanErr, &CorPar, SPVar, ErrVar, &CV);
  SEXP cvdf;
  if (result == OK)
  {
    SEXP y_rowName = PROTECT(getAttrib(x_R, R_RowNamesSymbol));
    SEXP y_colName = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(y_colName, 0, mkChar("Pred"));
    SET_STRING_ELT(y_colName, 1, mkChar("SE"));
    cvdf = MatrixDFConstructor(&CV, y_rowName, y_colName);

    UNPROTECT(2);
    MatFree(&CV);
  }

  AllocFree(y);
  StrFree(&RegMod_Term, (size_t)Rf_length(VECTOR_ELT(reg_mod, 0)));
  StrFree(&SPMod_Term, (size_t)Rf_length(VECTOR_ELT(sp_mod, 0)));
  StrFree(&xName, (size_t)Rf_length(getAttrib(x_R, R_NamesSymbol)));
  MatFree(&X);
  MatFree(&CorPar);
  ModFree(&RegMod);
  ModFree(&SPMod);
  if (result != OK)
  {
    Rf_error("GaSP Cross Validation failed.");
  }
  return cvdf;
}
