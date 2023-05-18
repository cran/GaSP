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

Matrix XDescrip;
Matrix ANOVAPerc;
Matrix MainEff;
Matrix JointEff;
Matrix PredReg;

extern int ErrNum;
extern int ErrorSeverityLevel;

int CalcVisualize(const Matrix *X, const real *y,
                  const LinModel *RegMod, const LinModel *SPMod,
                  size_t CorFamNum, boolean RanErr,
                  real *SPVar, real *ErrVar, Matrix *CorPar, real *MainPerc, real *InterPerc, real **Summary)
{
     int ErrReturn;
     KrigingModel KrigMod;
     Matrix GroupVarIndex;
     real *ANOVATot, *Average, *Perc;
     real *SEAve;
     size_t i, m, nNotNA;
     size_t *GroupSize;
     *Summary = AllocReal(3, NULL);

     /* Determine group structure (GroupSize and  */
     /* GroupVarIndex allocated in RegGroupings). */
     RegGroupings(&PredReg, &GroupSize, &GroupVarIndex);

     // /* Add columns to the YDescrip matrix. */
     // ANOVATot Average SEAve
     ANOVATot = AllocReal(1, NULL);
     Average = AllocReal(1, NULL);
     SEAve = AllocReal(1, NULL);
     int *Coltypes = AllocInt(5, NULL);
     Coltypes[0] = STRING;
     Coltypes[1] = STRING;
     Coltypes[2] = REALC;
     Coltypes[3] = REALC;
     Coltypes[4] = REALC;
     MatAllocate(0, 5, RECT, MIXED, Coltypes, YES, &MainEff);
     int *ColtypesI = AllocInt(7, NULL);
     ColtypesI[0] = STRING;
     ColtypesI[1] = STRING;
     ColtypesI[2] = STRING;
     ColtypesI[3] = REALC;
     ColtypesI[4] = REALC;
     ColtypesI[5] = REALC;
     ColtypesI[6] = REALC;
     MatAllocate(0, 7, RECT, MIXED, ColtypesI, YES, &JointEff);
     /* Compute effects for each response. */
     ErrReturn = OK;

     m = MatNumRows(X);

     /* Set up kriging model. */
     KrigModAlloc(m, MatNumCols(X), RegMod, SPMod, CorFamNum, RanErr, &KrigMod);
     KrigModData(m, NULL, X, y, &KrigMod);

     /* SPModMat contains the correlation parameters. */
     ErrNum = KrigModSetUp(CorPar, *SPVar, *ErrVar, &KrigMod);

     if (ErrNum == OK)
     {
          Perc = MatColAdd("y", &ANOVAPerc);
          ErrNum = CompEffects(&KrigMod, "y", &PredReg,
                               GroupSize, &GroupVarIndex, *MainPerc, *InterPerc,
                               Perc, Average, SEAve);
     }

     if (ErrNum != OK)
          ErrReturn = ErrNum;

     *ANOVATot = 0.0;
     for (nNotNA = 0, i = 0; i < MatNumRows(&ANOVAPerc); i++)
          if (Perc[i] != NA_REAL)
          {
               nNotNA++;
               *ANOVATot += Perc[i];
          }
     if (nNotNA == 0)
          *ANOVATot = NA_REAL;
     if (ErrNum == OK)
     {
          Summary[0][0] = *ANOVATot;
          Summary[0][1] = *Average;
          Summary[0][2] = *SEAve;
     }

     AllocFree(GroupSize);
     AllocFree(ANOVATot);
     AllocFree(Average);
     AllocFree(SEAve);
     AllocFree(Coltypes);
     AllocFree(ColtypesI);
     MatFree(&GroupVarIndex);
     KrigModFree(&KrigMod);
     if (ErrNum != OK)
          ErrReturn = ErrNum;

     return ErrReturn;
}

SEXP visualize(SEXP x_R, SEXP y_R, SEXP reg_mod, SEXP sp_mod,
               SEXP corFamNum, SEXP ranErr,
               SEXP spVar, SEXP errVar, SEXP corpar, SEXP mainEffectPct, SEXP interactionEffectPct,
               SEXP x_descrip)
{
     size_t CorFamNum = asInteger(corFamNum);
     boolean RanErr = asLogical(ranErr);
     real SPVar = asReal(spVar);
     real ErrVar = asReal(errVar);
     real MainPerc = asReal(mainEffectPct);
     real InterPerc = asReal(interactionEffectPct);
     matrix X, CorPar;
     real *y;
     string *xName;
     LinModel RegMod;
     LinModel SPMod;
     string *RegMod_Term, *SPMod_Term;
     real *Summary;

     MatrixDFAlloc(&X, x_R);
     MatrixDFAlloc(&CorPar, corpar);
     RealVecAlloc(&y, y_R);
     RegModDFAlloc(&RegMod_Term, reg_mod);
     RegModDFAlloc(&SPMod_Term, sp_mod);
     GetColName(&xName, x_R);
     XDescripAlloc(&XDescrip, x_descrip, xName);

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
     if (ErrNum == OK)
     {
          ErrNum = RegExtract(&XDescrip, X_DESCRIP, "." PRED, &PredReg);
     }
     if (ErrNum == OK)
     {
          ErrNum = ANOVAPercAlloc(&ANOVAPerc, &PredReg, xName);
     }
     if (ErrNum != OK)
     {
          AllocFree(y);
          StrFree(&RegMod_Term, (size_t)Rf_length(VECTOR_ELT(reg_mod, 0)));
          StrFree(&SPMod_Term, (size_t)Rf_length(VECTOR_ELT(sp_mod, 0)));
          StrFree(&xName, (size_t)Rf_length(getAttrib(x_R, R_NamesSymbol)));
          MatFree(&X);
          MatFree(&XDescrip);
          MatFree(&PredReg);
          MatFree(&ANOVAPerc);
          MatFree(&CorPar);
          ModFree(&RegMod);
          ModFree(&SPMod);
          Rf_error("GaSP Visualize c setup failed.");
     }

     int result = CalcVisualize(&X, y, &RegMod, &SPMod, CorFamNum, RanErr,
                                &SPVar, &ErrVar, &CorPar, &MainPerc, &InterPerc, &Summary);
     SEXP results = PROTECT(allocVector(VECSXP, 4));
     if (result == OK)
     {
          SEXP ANOVAdf = ANOVAMatrixDFConstructor(&ANOVAPerc);
          SET_VECTOR_ELT(results, 0, ANOVAdf);
          SEXP MainEffdf = MainEffDFConstructor(&MainEff);
          SET_VECTOR_ELT(results, 1, MainEffdf);
          SEXP JointEffdf = JointEffDFConstructor(&JointEff);
          SET_VECTOR_ELT(results, 2, JointEffdf);
          SEXP summary = RealVecConstructor(&Summary, 3);
          SET_VECTOR_ELT(results, 3, summary);
     }
     AllocFree(y);
     StrFree(&RegMod_Term, (size_t)Rf_length(VECTOR_ELT(reg_mod, 0)));
     StrFree(&SPMod_Term, (size_t)Rf_length(VECTOR_ELT(sp_mod, 0)));
     StrFree(&xName, (size_t)Rf_length(getAttrib(x_R, R_NamesSymbol)));
     AllocFree(Summary);
     MatFree(&X);
     MatFree(&CorPar);
     MatFree(&XDescrip);
     MatFree(&PredReg);
     MatFree(&ANOVAPerc);
     MatFree(&MainEff);
     MatFree(&JointEff);
     ModFree(&RegMod);
     ModFree(&SPMod);
     if (result != OK)
     {
          Rf_error("GaSP Visualise failed.");
     }
     UNPROTECT(1);
     return results;
}

int CompEffects(KrigingModel *KrigMod, const string yName,
                const Matrix *PredReg, const size_t *GroupSize,
                const Matrix *GroupVarIndex, real MainPerc, real InterPerc,
                real *Perc, real *Average, real *SEAve)
{
     boolean *ActiveGroup;
     int ErrNum, SevSave;
     Matrix IndexSP;
     real RAve, SSTot, VarEff, VarEffRow;
     real *Eff = NULL, *fAve, *rAve;
     real *SE = NULL;
     size_t c, kSP, i, i1, i2, j, jj, j1, j2, m, m1, m2;
     size_t n, nGroups, x1Index, x2Index;
     size_t IndexGroup[2];
     size_t *nSPTerms, *xIndex;

     kSP = ModDF(KrigSPMod(KrigMod));
     n = MatNumRows(KrigChol(KrigMod));
     nGroups = MatNumCols(GroupVarIndex);

     /* Workspace in KrigMod. */
     fAve = KrigMod->fRow;
     rAve = KrigMod->r;

     /* Allocations. */
     ActiveGroup = (boolean *)AllocGeneric(nGroups,
                                           sizeof(boolean), NULL);
     nSPTerms = AllocSize_t(nGroups, NULL);
     MatAllocate(kSP, nGroups, RECT, SIZE_T, NULL, NO, &IndexSP);

     for (j = 0; j < nGroups; j++)
          nSPTerms[j] = KrigSPActiveTerms(KrigMod, GroupSize[j],
                                          MatSize_tCol(GroupVarIndex, j),
                                          MatSize_tCol(&IndexSP, j));

     // /* Average predictor w.r.t. all x variables. */
     AvePred(KrigMod, PredReg, 0, NULL, GroupSize, GroupVarIndex,
             nSPTerms, &IndexSP, fAve, rAve, &RAve);

     if ((ErrNum = KrigYHatSE(KrigMod, RAve, fAve, rAve, Average,
                              SEAve)) == OK)
     {
          /* Subtract average prediction from data,   */
          /* so that predictions will have average 0. */
          VecAddScalar(-(*Average), n, KrigY(KrigMod));

          /* Get new decompositions and total sum of squares */
          /* for all variables.                              */
          KrigCorMat(0, NULL, KrigMod);
          if ((ErrNum = KrigDecompose(KrigMod)) == OK)
               ErrNum = CompSSTot(KrigMod, PredReg, GroupSize,
                                  GroupVarIndex, nSPTerms, &IndexSP, &SSTot);
     }

     if (ErrNum == OK)
     {
          /* Determine active groups. */
          for (j = 0; j < nGroups; j++)
          {
               xIndex = MatSize_tCol(GroupVarIndex, j);
               for (jj = 0; jj < GroupSize[j]; jj++)
                    if (KrigIsXActive(KrigMod, xIndex[jj]))
                         /* Active variable found. */
                         break;
               ActiveGroup[j] = (jj < GroupSize[j]);
          }
     }

     if (SSTot < sqrt(EPSILON))
     {
          SSTot = 0.0;
          SevSave = ErrorSeverityLevel;
          ErrorSeverityLevel = SEV_WARNING;
          Rf_error("No variation in predictor");
          ErrorSeverityLevel = SevSave;
     }

     /* Main effects and contributions. */
     for (j = 0; j < nGroups && ErrNum == OK; j++)
     {
          // Rprintf("Variable: %s  Main effect: %d \n", yName, j + 1);

          if (!ActiveGroup[j] && MainPerc > 0.0)
          {
               Perc[j] = 0.0;
               continue;
          }

          x1Index = MatSize_tElem(GroupVarIndex, 0, j);
          m = RegNumLevels(PredReg, x1Index);

          Eff = AllocReal(m, Eff);
          SE = AllocReal(m, SE);

          /* Compute main effect of group j. */
          AnyEffect(KrigMod, PredReg, 1, &j, GroupSize,
                    GroupVarIndex, nSPTerms, &IndexSP, Eff, SE);

          if (SSTot == 0.0)
               Perc[j] = NA_REAL;
          else
          {
               for (VarEff = 0.0, i = 0; i < m; i++)
                    VarEff += RegLevelWt(PredReg, x1Index, i) * Eff[i] * Eff[i];

               Perc[j] = VarEff / SSTot * 100.0;
               if (Perc[j] < sqrt(EPSILON))
                    Perc[j] = 0.0;
          }

          /* Previously required GroupSize[j] == 1. */
          if (Perc[j] != NA_REAL && Perc[j] >= MainPerc)
          {
               /* Generate plotting coordinates. */

               /* Add average prediction back in. */
               VecAddScalar(*Average, m, Eff);

               /* Append effect to MAIN_EFF. */
               AppendEffect(yName, 1, &j, PredReg, GroupSize,
                            GroupVarIndex, Eff, SE, &MainEff);
          }
     }

     /* Joint effects and interaction contributions. */
     for (c = nGroups, j1 = 0; j1 < nGroups - 1; j1++)
     {
          IndexGroup[0] = j1;
          x1Index = MatSize_tElem(GroupVarIndex, 0, j1);
          m1 = RegNumLevels(PredReg, x1Index);

          for (j2 = j1 + 1; j2 < nGroups; j2++, c++)
          {
               // Rprintf("Variable: %s  Joint effect: %d\n", yName, c - nGroups + 1);

               IndexGroup[1] = j2;

               if ((!ActiveGroup[j1] || !ActiveGroup[j2]) &&
                   InterPerc > 0.0)
               {
                    Perc[c] = 0.0;
                    continue;
               }

               x2Index = MatSize_tElem(GroupVarIndex, 0, j2);
               m2 = RegNumLevels(PredReg, x2Index);

               Eff = AllocReal(m1 * m2, Eff);
               SE = AllocReal(m1 * m2, SE);

               /* Compute joint effect of variables j1 and j2. */
               AnyEffect(KrigMod, PredReg, 2, IndexGroup, GroupSize,
                         GroupVarIndex, nSPTerms, &IndexSP, Eff, SE);

               /* Contribution of the *interaction* effect. */
               if (SSTot == 0.0)
                    Perc[c] = NA_REAL;
               else
               {
                    for (VarEff = 0.0, i1 = 0; i1 < m1; i1++)
                    {
                         for (VarEffRow = 0.0, i2 = 0; i2 < m2; i2++)
                              VarEffRow += RegLevelWt(PredReg,
                                                      x2Index, i2) *
                                           Eff[i1 * m2 + i2] * Eff[i1 * m2 + i2];
                         VarEff += RegLevelWt(PredReg, x1Index, i1) * VarEffRow;
                    }

                    Perc[c] = VarEff / SSTot * 100.0 - ((Perc[j1] != NA_REAL) ? Perc[j1] : 0.0) - ((Perc[j2] != NA_REAL) ? Perc[j2] : 0.0);
                    if (Perc[c] < sqrt(EPSILON))
                         Perc[c] = 0.0;
               }

               /* Previously required:                         */
               /* GroupSize[j1] == 1 && and GroupSize[j2] == 1 */
               if (Perc[c] != NA_REAL && Perc[c] >= InterPerc)
               {
                    /* Generate plotting coordinates. */

                    /* Add average back in. */
                    VecAddScalar(*Average, m1 * m2, Eff);

                    /* Append joint effect to JOINT_EFF. */
                    AppendEffect(yName, 2, IndexGroup, PredReg,
                                 GroupSize, GroupVarIndex, Eff, SE,
                                 &JointEff);
               }
          }
     }

     AllocFree(ActiveGroup);
     AllocFree(nSPTerms);
     MatFree(&IndexSP);
     AllocFree(Eff);
     AllocFree(SE);

     return ErrNum;
}

void AvePred(KrigingModel *KrigMod, const Matrix *PredReg,
             size_t nGroups, const size_t *IndexGroup,
             const size_t *GroupSize, const Matrix *GroupVarIndex,
             const size_t *nSPTerms, const Matrix *IndexSP,
             real *fAve, real *rAve, real *RAve)
{
     const LinModel *RegMod, *SPMod;
     Matrix         GAve;
     real           wRw, wRwj, SPVarPropSave;
     real           *f, *fj, *g, *r, *rj, *Rj = NULL, *Wt = NULL, *xRow;
     size_t         i, j, kReg, kSP, m, n;
     size_t         *IndexSPCol, *xIndex;

     n = MatNumRows(KrigChol(KrigMod));
     RegMod = KrigRegMod(KrigMod);
     SPMod = KrigSPMod(KrigMod);
     kReg = ModDF(RegMod);
     kSP = ModDF(SPMod);

     /* Allocations. */
     f = AllocReal(kReg, NULL);
     fj = AllocReal(kReg, NULL);
     g = AllocReal(kSP, NULL);
     r = AllocReal(n, NULL);
     rj = AllocReal(n, NULL);
     MatInit(RECT, REALC, NO, &GAve);

     /* Workspace in KrigMod. */
     xRow = KrigMod->xRow;

     VecInit(1.0, kReg, fAve);
     VecInit(KrigMod->SPVarProp, n, rAve);

     wRw = KrigMod->SPVarProp;

     /* SPVarProp is applied only once. */
     SPVarPropSave = KrigMod->SPVarProp;
     KrigMod->SPVarProp = 1.0;

     for (j = 0; j < MatNumCols(GroupVarIndex); j++)
     {
          if (VecSize_tIndex(j, nGroups, IndexGroup) != INDEX_ERR)
               continue;

          xIndex = MatSize_tCol(GroupVarIndex, j);
          IndexSPCol = MatSize_tCol(IndexSP, j);

          m = RegNumLevels(PredReg, xIndex[0]);
          if (m <= 0)
               Rf_error("AvePred failed.");
          Rj = AllocReal(m, Rj);
          Wt = AllocReal(m, Wt);
          MatReAlloc(m, kSP, &GAve);

          VecInit(0.0, kReg, fj);
          VecInit(0.0, n, rj);

          for (wRwj = 0.0, i = 0; i < m; i++)
          {
               fgrGroup(KrigMod, PredReg, GroupSize[j], xIndex, i,
                        nSPTerms[j], IndexSPCol, xRow, f, g, r);

               Wt[i] = RegLevelWt(PredReg, xIndex[0], i);

               VecAddVec(Wt[i], f, kReg, fj);
               VecAddVec(Wt[i], r, n, rj);

               MatRowPut(g, i, &GAve);
               KrigCorVec(g, &GAve, i, nSPTerms[j], IndexSPCol, YES,
                          KrigMod, Rj);
               wRwj += Wt[i] * (Wt[i] + 2.0 * DotProd(Wt, Rj, i));
          }

          VecMultVec(fj, kReg, fAve);
          VecMultVec(rj, n, rAve);

          wRw *= wRwj;
     }

     *RAve = wRw;

     KrigMod->SPVarProp = SPVarPropSave;

     AllocFree(f);
     AllocFree(fj);
     AllocFree(g);
     AllocFree(r);
     AllocFree(Rj);
     AllocFree(rj);
     AllocFree(Wt);

     MatFree(&GAve);
}

int CompSSTot(KrigingModel *KrigMod, const Matrix *PredReg,
              const size_t *GroupSize, const Matrix *GroupVarIndex,
              const size_t *nSPTerms, const Matrix *IndexSP,
              real *SSTot)
{
     int ErrNum;
     Matrix frfr, frfrj;
     real a;
     real *eVal, *v;
     size_t j, k, n;

     ErrNum = OK;

     n = MatNumRows(KrigChol(KrigMod));
     k = ModDF(KrigRegMod(KrigMod));

     /* Allocations:                                            */
     /* frfr is RECT because it is overwritten by eigenvectors. */
     MatAlloc(k + n, k + n, RECT, &frfr);
     MatAlloc(k + n, k + n, SYM, &frfrj);

     /* Workspace in KrigMod. */
     eVal = KrigMod->fr;

     /* frfr must be symmetric for frfrAve. */
     MatPutShape(&frfr, SYM);
     frfrAve(KrigMod, PredReg, GroupSize, GroupVarIndex, nSPTerms,
             IndexSP, &frfrj, &frfr);
     MatPutShape(&frfr, RECT);

     /* Eigen decomposition, overwriting frfr with eigenvectors. */
     if ((ErrNum = MatEig(YES, &frfr, eVal, &frfr)) != OK)
          Error("Eigen decomposition of averaging moment matrix failed.");

     for (*SSTot = 0.0, j = 0; j < k + n && ErrNum == OK; j++)
     {

          if (eVal[j] < EPSILON * eVal[0])
               break;

          /* Eigenvector j is row j of V'. */
          v = MatCol(&frfr, j);

          /* Overwrite first k elements of v with solution. */
          if ((ErrNum = TriForSolve(KrigR(KrigMod), v, 0, v)) != OK)
               Error("Ill-conditioned expanded-design matrix.\n");

          /* Overwrite next n elements of v with solution. */
          else if ((ErrNum = TriForSolve(KrigChol(KrigMod),
                                         v + k, 0, v + k)) != OK)
               Error("Ill-conditioned correlation matrix.\n");

          else
          {
               a = VecDotProd(k, v, KrigMod->RBeta) + VecDotProd(n, v + k, KrigMod->ResTilde);
               *SSTot += eVal[j] * a * a;
          }
     }

     MatFree(&frfr);
     MatFree(&frfrj);

     return ErrNum;
}

void AnyEffect(KrigingModel *KrigMod, const Matrix *PredReg,
               size_t nGroups, const size_t *IndexGroup,
               const size_t *GroupSize, const Matrix *GroupVarIndex,
               const size_t *nSPTerms, const Matrix *IndexSP,
               real *Eff, real *SE)
{
     real RAve, SPVarPropSave;
     real *fAve, *f, *fj, *g, *r, *rj, *rAve, *xRow;
     size_t i, j, jj, k, n;
     size_t *Level, *nLevels, *xIndex;

     n = MatNumRows(KrigChol(KrigMod));
     k = ModDF(KrigRegMod(KrigMod));

     /* Workspace in KrigMod. */
     fAve = KrigMod->fRow;
     g = KrigMod->gRow;
     rAve = KrigMod->r;
     xRow = KrigMod->xRow;

     /* Allocations. */
     f = AllocReal(k, NULL);
     fj = AllocReal(k, NULL);
     r = AllocReal(n, NULL);
     rj = AllocReal(n, NULL);
     Level = AllocSize_t(nGroups, NULL);
     nLevels = AllocSize_t(nGroups, NULL);

     AvePred(KrigMod, PredReg, nGroups, IndexGroup, GroupSize,
             GroupVarIndex, nSPTerms, IndexSP, fAve, rAve, &RAve);

     /* SPVarProp was applied in AvePred. */
     SPVarPropSave = KrigMod->SPVarProp;
     KrigMod->SPVarProp = 1.0;

     /* Adjust averages for variables in the effect. */

     /* Set up initial group levels. */
     for (j = 0; j < nGroups; j++)
     {
          Level[j] = 0;
          nLevels[j] = RegNumLevels(PredReg,
                                    MatSize_tElem(GroupVarIndex, 0, IndexGroup[j]));
     }

     /* For each level combination. */
     i = 0;
     do
     {
          VecCopy(fAve, k, f);
          VecCopy(rAve, n, r);

          /* Adjust f and r for combination i */
          /* of variables in effect.          */
          for (jj = 0; jj < nGroups; jj++)
          {
               j = IndexGroup[jj];

               xIndex = MatSize_tCol(GroupVarIndex, j);

               fgrGroup(KrigMod, PredReg, GroupSize[j], xIndex,
                        Level[jj], nSPTerms[j],
                        MatSize_tCol(IndexSP, j), xRow, fj, g, rj);

               VecMultVec(fj, k, f);
               VecMultVec(rj, n, r);
          }

          /* Need to transform r with T? */

          KrigYHatSE(KrigMod, RAve, f, r, &Eff[i], &SE[i]);

          i++;
     } while (LevelLex(nGroups, nLevels, Level) != ALL_DONE);

     KrigMod->SPVarProp = SPVarPropSave;

     AllocFree(f);
     AllocFree(fj);
     AllocFree(r);
     AllocFree(rj);
     AllocFree(Level);
     AllocFree(nLevels);

     return;
}

void AppendEffect(const string yName, size_t DegreeEff,
                  const size_t *IndexGroup, const Matrix *PredReg,
                  const size_t *GroupSize, const Matrix *GroupVarIndex,
                  real *Eff, real *SE, Matrix *EffMat)
{
     real l;
     size_t i, j, jj, nRowsNew, nRowsOld, xIndex;
     size_t *Level, *nLevels;
     string s;

     nRowsOld = MatNumRows(EffMat);

     /* Allocations. */
     Level = AllocSize_t(DegreeEff, NULL);
     nLevels = AllocSize_t(DegreeEff, NULL);

     /* Set up initial variable levels. */
     for (nRowsNew = 1, j = 0; j < DegreeEff; j++)
     {
          Level[j] = 0;
          nLevels[j] = RegNumLevels(PredReg,
                                    MatSize_tElem(GroupVarIndex, 0, IndexGroup[j]));
          nRowsNew *= nLevels[j];
     }

     MatReAlloc(nRowsOld + nRowsNew, MatNumCols(EffMat), EffMat);

     /* Put x-variable names in first new row of EffMat. */
     /* Names will be copied for further rows.           */
     for (jj = 0; jj < DegreeEff; jj++)
     {
          j = IndexGroup[jj];
          xIndex = MatSize_tElem(GroupVarIndex, 0, j);
          if (GroupSize[j] == 1)
               MatPutStrElem(EffMat, nRowsOld, jj, RegVar(PredReg, xIndex));
          else
          {
               s = StrPaste(2, GROUP, StrFromSize_t(RegCandGroup(PredReg, xIndex)));
               MatPutStrElem(EffMat, nRowsOld, jj, s);
               AllocFree(s);
          }
     }

     i = 0;
     do
     {
          for (jj = 0; jj < DegreeEff; jj++)
          {
               if (i > 0)
                    MatPutStrElem(EffMat, nRowsOld + i, jj,
                                  MatStrElem(EffMat, nRowsOld + i - 1, jj));

               j = IndexGroup[jj];
               xIndex = MatSize_tElem(GroupVarIndex, 0, j);

               l = (GroupSize[j] == 1) ? RegLevel(PredReg, xIndex, Level[jj]) : (real)(Level[jj] + 1);
               MatPutElem(EffMat, nRowsOld + i, DegreeEff + jj + 1,
                          l);
          }

          MatPutStrElem(EffMat, nRowsOld + i, DegreeEff, yName);
          MatPutElem(EffMat, nRowsOld + i, 2 * DegreeEff + 1, Eff[i]);
          MatPutElem(EffMat, nRowsOld + i, 2 * DegreeEff + 2, SE[i]);
          i++;

     } while (LevelLex(DegreeEff, nLevels, Level) != ALL_DONE);

     AllocFree(Level);
     AllocFree(nLevels);

     return;
}
