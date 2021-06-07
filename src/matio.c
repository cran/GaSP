/*****************************************************************/
/*   MATRIX INPUT-OUTPUT ROUTINES                                */
/*                                                               */
/*   Copyright (c) William J. Welch 1991--95.                    */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "implem.h"
#include "matrix.h"
#include "lib.h"

#define BETWEEN_SPACES   2

#define ROW_ALLOC   100  /* Was 32. */

/*******************************+++*******************************/
void MatWriteBlock(Matrix *M, boolean CaseLabels)
/*****************************************************************/
/*   Purpose:  Write a matrix to a file in blocks of columns     */
/*             that fit on a single output line.                 */
/*                                                               */
/*   Comment:  Only RECT, Labelled matrices can be written.      */
/*                                                               */
/*   Version:  1995 October 25                                   */
/*****************************************************************/
{
     boolean   RightAdj;
     char      *Conversion;
     size_t    CaseWidth, FirstCol, i, j, LastCol;
     size_t    LineWidth, NumCols, NumRows;
     size_t    *ColWidth;
     int       *Precision;
     string    *ColName, s;

     NumRows = MatNumRows(M);
     NumCols = MatNumCols(M);
     ColName = MatColNames(M);

     /* Output the matrix name. */
     Output("%s", (MatText(M) != NULL) ? MatText(M) :
               "Unnamed matrix.\n\n");

     /* Get column widths and precisions. */
     ColWidth   = AllocSize_t(NumCols, NULL);
     Conversion = AllocChar(NumCols, NULL);
     Precision = AllocInt(NumCols, NULL);
     for (j = 0; j < NumCols; j++)
          ColWidth[j] = MatColWidth(M, j, &Precision[j],
                    &Conversion[j]);

     /* Column width for case labels. */
     if (CaseLabels)
          CaseWidth = MatCaseWidth(M, &RightAdj);
     else
          CaseWidth = 0;

     /* Output in blocks of columns that will fit on a line. */
     FirstCol = 0;
     while (FirstCol < NumCols)
     {
          /* Always output at least one column, whatever. */
          LastCol = FirstCol;
          LineWidth = ((CaseLabels) ? CaseWidth + BETWEEN_SPACES : 0)
                    + ColWidth[FirstCol];

          while (LastCol + 1 < NumCols &&
                    LineWidth + BETWEEN_SPACES
                    + ColWidth[LastCol+1] <= OUTPUT_COLS)
               LineWidth += BETWEEN_SPACES + ColWidth[++LastCol];

          for (j = 0; j < LineWidth; j++)
               Output("-");
          Output("\n");

          if (CaseLabels)
               Output((RightAdj) ? "%*s" : "%-*s",
                         (int) CaseWidth, "Case");
          for (j = FirstCol; j <= LastCol; j++)
          {
               if (CaseLabels || j > FirstCol)
                    Output("%*s", BETWEEN_SPACES, "");
               Output((MatColType(M, j) == STRING)
                         ? "%-*s" : "%*s", (int) ColWidth[j],
                         ColName[j]);
          }
          Output("\n");

          for (j = 0; j < LineWidth; j++)
               Output("-");
          Output("\n\n");

          for (i = 0; i < NumRows; i++)
          {
               if (CaseLabels)
                    Output((RightAdj) ? "%*s" : "%-*s",
                              (int) CaseWidth, MatRowName(M, i));
               for (j = FirstCol; j <= LastCol; j++)
               {
                    if (CaseLabels || j > FirstCol)
                         Output( "%*s", BETWEEN_SPACES, "");
                    s = MatElemToStr(M, i, j, Precision[j],
                              Conversion[j]);
                    Output((MatColType(M, j) == STRING)
                              ? "%-*s" : "%*s", (int) ColWidth[j], s);
               }
               Output("\n");
          }

          if ( (FirstCol = LastCol + 1) < NumCols)
               Output("\n");
     }

     AllocFree(ColWidth);
     AllocFree(Conversion);
     AllocFree(Precision);

     return;
}

/*******************************+++*******************************/
size_t MatColWidth(const Matrix *M, size_t j, int *Precision,
               char *Conversion)
/*****************************************************************/
/*   Purpose:    Return the column width necessary for outputting  */
/*               column j of M.                                    */
/*               If the column is of type real, then, on return,   */
/*               *Precision will be the minimum precision that     */
/*               does not lose accuracy, and *Conversion will be   */
/*               'e' or 'f'.                                       */
/*                                                                 */
/*   2020.08.22: NULL replaced by '\0'                             */
/*****************************************************************/
{
     size_t    ExponLen, i, Width;
     string    DecPoint, Expon, s, StrEnd;

     Width = 0;

     if (MatColType(M, j) == REALC)
     {
          *Conversion = 'g';
          *Precision = 0;

          /* Are there any e conversions? */
          for (i = 0; i < MatNumRows(M); i++)
               if (strchr(StrFromReal(MatElem(M, i, j), "",
                         PRECISION, 'g'), 'e') != NULL)
               {
                    *Conversion = 'e';
                    break;
               }

          for (i = 0; i < MatNumRows(M); i++)
          {
               /* Conversion is either e or g (not f) to maintain */
               /* the minimum number of significant digits.       */
               s = StrFromReal(MatElem(M, i, j), "#",
                         PRECISION, *Conversion);

               if (stricmp(s, NOT_AVAIL) == 0)
                    Width = max(strlen(s), Width);
               else
               {
                    ExponLen = 0;
                    if ( (Expon = strchr(s, 'e')) != NULL)
                    {
                         ExponLen = strlen(Expon);

                         /* Delete exponent part of string. */
                         *Expon = '\0';
                    }

                    /* Find decimal point.  This assumes that the */
                    /* real has been converted to a string with a */
                    /* # flag, so that '.' is  always present.    */
                    DecPoint = strchr(s, '.');
                    CodeCheck(DecPoint != NULL);

                    /* Delete trailing zeros to get precision. */
                    StrEnd = s + strlen(s) - 1;
                    while (StrEnd > DecPoint && *StrEnd == '0')
                         *StrEnd-- = '\0';
                    *Precision = max(StrEnd - DecPoint, *Precision);

                    /* Delete digits after decimal point;        */
                    /* if none, delete the decimal point itself. */
                    if (StrEnd == DecPoint)
                         *DecPoint = '\0';
                    else
                         *(DecPoint + 1) = '\0';

                    /* Everything but precision. */
                    Width = max(strlen(s) + ExponLen, Width);
               }
          }

          Width += (size_t) *Precision;
          if (*Conversion == 'g')
               *Conversion = 'f';
     }
     else
          for (i = 0; i < MatNumRows(M); i++)
               Width = max(strlen(MatElemToStr(M, i, j, 'x', -1)),
                         Width);

     Width = max(Width, strlen(MatColName(M, j)));

     return Width;
}

/*******************************+++*******************************/
size_t MatCaseWidth(const Matrix *M, boolean *RightAdj)
/*****************************************************************/
/*   Purpose:  Return the column width necessary for outputting  */
/*             the case column of M.                             */
/*                                                               */
/*   Version:  1995 October 25                                   */
/*****************************************************************/
{
     size_t    CaseNumber, i, Width;
     string    s;

     *RightAdj = TRUE;
     for (Width = 0, i = 0; i < MatNumRows(M); i++)
     {
          s = MatRowName(M, i);

          Width = max(strlen(s), Width);

          if (StrToSize_t(s, &CaseNumber) != OK)
               /* Non-numeric case labels. */
               *RightAdj = FALSE;
     }

     Width = max(Width, strlen("Case"));

     return Width;
}

/*******************************+++*******************************/
string MatElemToStr(const Matrix *M, size_t i, size_t j,
          int Precision, char Conversion)
/*****************************************************************/
/*   Purpose:  Return a string representing element i, j of M.   */
/*             Precision and Conversion are only used if the     */
/*             element is real.                                  */
/*                                                               */
/*   96.01.22: IllegalType replaced by CodeBug.                  */
/*                                                               */
/*   Version:  1996 January 22                                   */
/*****************************************************************/
{
     string    s;

     switch (MatColType(M, j))
     {
          case INTEGERC:
               s = StrFromInt(MatIntElem(M, i, j));
               break;

          case REALC:
               s = StrFromReal(MatElem(M, i, j), "", Precision,
                         Conversion);
               break;

          case SIZE_T:
               s = StrFromSize_t(MatSize_tElem(M, i, j));
               break;

          case STRING:
               s = MatStrElem(M, i, j);
               break;

          default:
               CodeBug("Illegal type");
     }

     return s;
}
