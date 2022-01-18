#include <R.h>
#include <Rinternals.h>

#include "implem.h"
#include "define.h"
#include "matrix.h"
#include "lib.h"

void MatrixDFAlloc(matrix *m, SEXP df);

int ANOVAPercAlloc(matrix *ANOVAPerc, matrix *PredReg, const string *xName);

void StrFree(string **s, size_t n);

void XDescripAlloc(matrix *m, SEXP df, const string *xName);

void RealDFAlloc(real **r, SEXP df);

void RealVecAlloc(real **r, SEXP v);

void RegModDFAlloc(string **s, SEXP df);

void GetColName(string **s, SEXP df);

SEXP JointEffDFConstructor(matrix *m);
  
SEXP MainEffDFConstructor(matrix *m);

SEXP ANOVAMatrixDFConstructor(matrix *m);

SEXP MatrixDFConstructor(matrix *m, SEXP rowName_R, SEXP colName_R);

SEXP RealDFConstructor(real **r, SEXP rowName_R, SEXP colName_R, size_t nRows);

SEXP RealVecConstructor(real **r, size_t nRows);

void ProgressInit(double taskCount);

void tick(double times);
