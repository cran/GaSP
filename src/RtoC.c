#include "RCconvert.h"

/*   2023.12.05: newline added at end of file */

SEXP RealVecConstructor(real **r, size_t nRows)
{
  SEXP vec = PROTECT(allocVector(REALSXP, (int)nRows));
  double *pvec = REAL(vec);
  for (size_t j = 0; j < nRows; j++)
    pvec[j] = r[0][j];
  UNPROTECT(1);
  return vec;
}

SEXP RealDFConstructor(real **r, SEXP rowName_R, SEXP colName_R, size_t nRows)
{
  SEXP df = PROTECT(allocVector(VECSXP, 1));
  SEXP col = PROTECT(allocVector(REALSXP, (int)nRows));
  double *pcol = REAL(col);
  for (size_t j = 0; j < nRows; j++)
    pcol[j] = r[0][j];
  SET_VECTOR_ELT(df, 0, col);

  setAttrib(df, R_ClassSymbol, ScalarString(mkChar("data.frame")));
  setAttrib(df, R_RowNamesSymbol, rowName_R);
  setAttrib(df, R_NamesSymbol, colName_R);

  UNPROTECT(2);
  return df;
}

SEXP JointEffDFConstructor(matrix *m)
/* 2024.06.23: size_t nRows replaces int */
/* 2024.06.23 size_t i replaces int i (twice) */
/* 2024.06.23 size_t j replaces int j (twice) */
{
  size_t nRows = MatNumRows(m);
  SEXP df = PROTECT(allocVector(VECSXP, 6));
  SEXP colName = PROTECT(allocVector(STRSXP, 6));
  SEXP rowName = PROTECT(allocVector(STRSXP, nRows));
  for (size_t j = 0; j < nRows; j++)
  {
    SET_STRING_ELT(rowName, j, mkChar(StrFromSize_t(j + 1)));
  }
  SET_STRING_ELT(colName, 0, mkChar("Variable.x_i"));
  SET_STRING_ELT(colName, 1, mkChar("Variable.x_j"));
  SET_STRING_ELT(colName, 2, mkChar("x_i"));
  SET_STRING_ELT(colName, 3, mkChar("x_j"));
  SET_STRING_ELT(colName, 4, mkChar("y"));
  SET_STRING_ELT(colName, 5, mkChar("SE"));
  for (size_t i = 0; i < 2; i++)
  {
    SEXP scol = PROTECT(allocVector(STRSXP, nRows));
    string *scol_c = MatStrCol(m, i);
    for (size_t j = 0; j < nRows; j++)
    {
      SET_STRING_ELT(scol, j, mkChar(scol_c[j]));
    }
    SET_VECTOR_ELT(df, i, scol);
    UNPROTECT(1);
  }
  for (size_t i = 2; i < 6; i++)
  {
    SEXP col = PROTECT(allocVector(REALSXP, nRows));
    double *pcol = REAL(col);
    for (size_t j = 0; j < nRows; j++)
    {
      pcol[j] = MatElem(m, j, i + 1);
    }
    SET_VECTOR_ELT(df, i, col);
    UNPROTECT(1);
  }

  setAttrib(df, R_ClassSymbol, ScalarString(mkChar("data.frame")));
  setAttrib(df, R_NamesSymbol, colName);
  setAttrib(df, R_RowNamesSymbol, rowName);

  UNPROTECT(3);
  return df;
}

SEXP MainEffDFConstructor(matrix *m)
/* 2024.06.23: size_t nRows, nCols replaces int */
/* 2024.06.23 size_t i replaces int i  */
/* 2024.06.23 size_t j replaces int j (twice) */
{
  size_t nCols = MatNumCols(m);
  size_t nRows = MatNumRows(m);
  SEXP df = PROTECT(allocVector(VECSXP, nCols - 1));
  SEXP colName = PROTECT(allocVector(STRSXP, 4));
  SEXP rowName = PROTECT(allocVector(STRSXP, nRows));
  for (size_t j = 0; j < nRows; j++)
  {
    SET_STRING_ELT(rowName, j, mkChar(StrFromSize_t(j + 1)));
  }
  SET_STRING_ELT(colName, 0, mkChar("Variable.x_i"));
  SET_STRING_ELT(colName, 1, mkChar("x_i"));
  SET_STRING_ELT(colName, 2, mkChar("y"));
  SET_STRING_ELT(colName, 3, mkChar("SE"));

  SEXP scol = PROTECT(allocVector(STRSXP, nRows));
  string *scol_c = MatStrCol(m, 0);
  for (size_t j = 0; j < nRows; j++)
  {
    SET_STRING_ELT(scol, j, mkChar(scol_c[j]));
  }
  SET_VECTOR_ELT(df, 0, scol);
  UNPROTECT(1);
  for (size_t i = 1; i < 4; i++)
  {
    SEXP col = PROTECT(allocVector(REALSXP, nRows));
    double *pcol = REAL(col);
    for (size_t j = 0; j < nRows; j++)
    {
      pcol[j] = MatElem(m, j, i + 1);
    }
    SET_VECTOR_ELT(df, i, col);
    UNPROTECT(1);
  }

  setAttrib(df, R_ClassSymbol, ScalarString(mkChar("data.frame")));
  setAttrib(df, R_NamesSymbol, colName);
  setAttrib(df, R_RowNamesSymbol, rowName);

  UNPROTECT(3);
  return df;
}

SEXP ANOVAMatrixDFConstructor(matrix *m)
/* 2024.06.23: size_t nRows, nCols replaces int */
/* 2024.06.23: size_t j replaces int (twice) */
{
  size_t nCols = MatNumCols(m);
  size_t nRows = MatNumRows(m);
  SEXP df = PROTECT(allocVector(VECSXP, nCols));
  SEXP colName = PROTECT(allocVector(STRSXP, 1));
  SEXP rowName = PROTECT(allocVector(STRSXP, nRows));
  string *rownames = MatRowNames(m);
  for (size_t j = 0; j < nRows; j++)
  {
    SET_STRING_ELT(rowName, j, mkChar(rownames[j]));
  }

  SET_STRING_ELT(colName, 0, mkChar("y"));
  SEXP col = PROTECT(allocVector(REALSXP, nRows));
  double *pcol = REAL(col);
  for (size_t j = 0; j < nRows; j++)
  {
    pcol[j] = MatElem(m, j, 0);
  }
  SET_VECTOR_ELT(df, 0, col);

  setAttrib(df, R_ClassSymbol, ScalarString(mkChar("data.frame")));
  setAttrib(df, R_NamesSymbol, colName);
  setAttrib(df, R_RowNamesSymbol, rowName);

  UNPROTECT(4);
  return df;
}

SEXP MatrixDFConstructor(matrix *m, SEXP rowName_R, SEXP colName_R)
{
  int nCols = (int)Rf_length(colName_R);
  int nRows = (int)Rf_length(rowName_R);
  SEXP df = PROTECT(allocVector(VECSXP, nCols));
  for (int i = 0; i < nCols; ++i)
  {
    SEXP col = PROTECT(allocVector(REALSXP, nRows));
    double *pcol = REAL(col);
    for (int j = 0; j < nRows; j++)
    {
      pcol[j] = MatElem(m, j, i);
    }
    SET_VECTOR_ELT(df, i, col);
  }

  setAttrib(df, R_ClassSymbol, ScalarString(mkChar("data.frame")));
  setAttrib(df, R_RowNamesSymbol, rowName_R);
  setAttrib(df, R_NamesSymbol, colName_R);

  UNPROTECT(1 + nCols);
  return df;
}

int ANOVAPercAlloc(matrix *ANOVAPerc, matrix *PredReg, const string *xName)
{
  int ErrNum = OK;
  size_t nXVars = MatNumRows(PredReg);
  size_t GroupSize, h, i, j, nEffects, nGroups;
  size_t *GroupVarIndex = AllocSize_t(nXVars, NULL);
  string s;
  for (nGroups = 0, j = 0; j < nXVars; j++)
  {
    RegGroupIndices(PredReg, j, GroupVarIndex);
    if (GroupVarIndex[0] == j)
      nGroups++;
  }
  nEffects = nGroups * (nGroups + 1) / 2;
  MatAllocate(nEffects, 0, RECT, REALC, NULL, YES, ANOVAPerc);
  if (nEffects != MatNumRows(ANOVAPerc))
    ErrNum = INCOMPAT_ERR;

  string *RowName = MatRowNames(ANOVAPerc);

  for (i = 0, j = 0; j < nXVars && ErrNum == OK; j++)
  {
    GroupSize = RegGroupIndices(PredReg, j, GroupVarIndex);
    if (GroupVarIndex[0] != j)
      continue;

    if (GroupSize == 1)
      s = StrDup(xName[j]);
    else
      s = StrPaste(2, GROUP, StrFromSize_t(RegCandGroup(PredReg, j)));

    if (RowName[i] == NULL)
      MatPutRowName(ANOVAPerc, i, s);
    else if (stricmp(RowName[i], s) != 0)
      ErrNum = INCOMPAT_ERR;

    AllocFree(s);
    i++;
  }

  for (h = nGroups, i = 0; i < nGroups - 1; i++)
  {
    for (j = i + 1; j < nGroups && ErrNum == OK; j++, h++)
    {
      s = StrPaste(3, RowName[i], ":", RowName[j]);
      if (RowName[h] == NULL)
        MatPutRowName(ANOVAPerc, h, s);
      else if (stricmp(RowName[h], s) != 0)
        ErrNum = INCOMPAT_ERR;
      AllocFree(s);
    }
  }
  AllocFree(GroupVarIndex);
  return ErrNum;
}

void ColNameCopy(string **s, SEXP colName)
{
  size_t nRows = (size_t)Rf_length(colName);

  for (size_t i = 0; i < nRows; i++)
  {
    if (s[0][i] != NULL)
      AllocFree(s[0][i]);

    s[0][i] = StrDup((string)CHAR(STRING_ELT(colName, i)));
  }
}

void StrFree(string **s, size_t n)
{
  for (size_t i = 0; i < n; i++)
  {
    if (s[0][i] != NULL)
      AllocFree(s[0][i]);
  }
  AllocFree(s[0]);
}

void XDescripAlloc(matrix *m, SEXP df, const string *xName)
/* 2024.06.23 size_t i replaces int i (twice) */
{
  string *colNames;
  GetColName(&colNames, df);
  size_t nCols = (size_t)Rf_length(df);
  size_t nRows = (size_t)Rf_length(VECTOR_ELT(df, 0));
  int *Coltypes = AllocInt(nCols, NULL);
  Coltypes[0] = STRING;
  Coltypes[1] = REALC;
  Coltypes[2] = REALC;
  if (nCols > 3)
  {
    for (size_t i = 3; i < nCols; i++)
    {
      if (stricmp(colNames[i], "Support") == 0)
      {
        Coltypes[i] = STRING;
      }
      else if (stricmp(colNames[i], "NumberLevels") == 0)
      {
        Coltypes[i] = SIZE_T;
      }
      else if (stricmp(colNames[i], "Distribution") == 0)
      {
        Coltypes[i] = STRING;
      }
    }
  }
  MatAllocate(nRows, nCols, RECT, MIXED, Coltypes, YES, m);
  MatPutColName(m, 0, VARIABLE);
  string min_pred = StrPaste(2, MIN, ".Pred");
  string max_pred = StrPaste(2, MAX, ".Pred");
  MatPutColName(m, 1, min_pred);
  MatPutColName(m, 2, max_pred);
  AllocFree(min_pred);
  AllocFree(max_pred);
  VecStrCopy(xName, nRows, MatStrCol(m, 0));
  SEXP v = VECTOR_ELT(df, 1);
  VecCopy(REAL(v), nRows, MatCol(m, 1));
  v = VECTOR_ELT(df, 2);
  VecCopy(REAL(v), nRows, MatCol(m, 2));
  if (nCols > 3)
  {
    for (size_t i = 3; i < nCols; i++)
    {
      if (stricmp(colNames[i], "Support") == 0)
      {
        MatPutColName(m, i, SUPPORT);
        v = VECTOR_ELT(df, i);
        string *s = MatStrCol(m, i);
        ColNameCopy(&s, v);
      }
      else if (stricmp(colNames[i], "NumberLevels") == 0)
      {
        MatPutColName(m, i, NUM_LEVELS);
        v = VECTOR_ELT(df, i);
        int *pcol = INTEGER(v);
        size_t *col = MatSize_tCol(m, i);
        for (size_t i = 0; i < nRows; i++)
          col[i] = (size_t)pcol[i];
      }
      else if (stricmp(colNames[i], "Distribution") == 0)
      {
        MatPutColName(m, i, DISTRIBUTION);
        v = VECTOR_ELT(df, i);
        string *s = MatStrCol(m, i);
        ColNameCopy(&s, v);
      }
    }
  }
  StrFree(&colNames, (size_t)Rf_length(getAttrib(df, R_NamesSymbol)));
  AllocFree(Coltypes);
}

void MatrixDFCopy(matrix *m, SEXP df)
{
  size_t nCols = MatNumCols(m);
  size_t nRows = MatNumRows(m);
  for (size_t i = 0; i < nCols; i++)
  {
    SEXP v = VECTOR_ELT(df, i);
    VecCopy(REAL(v), nRows, m->Elem[i]);
  }
}

void MatrixDFAlloc(matrix *m, SEXP df)
{
  size_t nCols = (size_t)Rf_length(df);
  size_t nRows = (size_t)Rf_length(VECTOR_ELT(df, 0));

  MatAlloc(nRows, nCols, RECT, m);
  MatrixDFCopy(m, df);
}

void RealVecAlloc(real **r, SEXP v)
{
  size_t n = (size_t)Rf_length(v);
  *r = AllocReal(n, NULL);
  VecCopy(REAL(v), n, *r);
}

void RealDFAlloc(real **r, SEXP df)
{
  SEXP v = VECTOR_ELT(df, 0);
  size_t nRows = (size_t)Rf_length(v);
  *r = AllocReal(nRows, NULL);
  VecCopy(REAL(v), nRows, *r);
}

void RegModDFAlloc(string **s, SEXP df)
{
  SEXP v = VECTOR_ELT(df, 0);
  size_t nRows = (size_t)Rf_length(v);
  *s = AllocStr(nRows, NULL);
  for (size_t i = 0; i < nRows; i++)
  {
    if (s[0][i] != NULL)
      AllocFree(s[0][i]);

    s[0][i] = StrDup((string)CHAR(STRING_ELT(v, i)));
  }
}

void GetColName(string **s, SEXP df)
{
  SEXP colName = getAttrib(df, R_NamesSymbol);
  size_t n = (size_t)Rf_length(colName);
  *s = AllocStr(n, NULL);
  ColNameCopy(s, colName);
}

double tickCount = 0;
double totalTasks;
double tickSize;
void ProgressInit(double taskCount)
{
  R_FlushConsole();
  Rprintf("\nProgress: [--------------------------------------------------]");
  tickSize = 50 / taskCount;
  totalTasks = taskCount;
  tickCount = 0;
}

void tick(double times)
/* 2024.06.23: int replaced by size_t three times */
/* 2024.06.23: cast applied to *product* of doubles */
{
  tickCount += times;
  if (tickCount == totalTasks)
  {
    Rprintf("\rProgress: [==================================================]");
    R_FlushConsole();
    Rprintf("\n\n");
    ErrorMatOut();
    R_FlushConsole();
  }
  else
  {
    Rprintf("\rProgress: [");
    /* int tick_f = (int)tickCount * tickSize; */
    /* for (int j = 0; j < tick_f; j++) */
    size_t tick_f = (size_t)(tickCount * tickSize);
    for (size_t j = 0; j < tick_f; j++)
    {
      Rprintf("=");
    }
    R_FlushConsole();
  }
}

