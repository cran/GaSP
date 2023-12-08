/*****************************************************************/
/*   LOW-LEVEL OUTPUT ROUTINES                                   */
/*                                                               */
/*   Copyright (c) William J. Welch 1990--2000.                  */
/*   All rights reserved.                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <R.h>

#include "define.h"
#include "implem.h"
#include "lib.h"

/* 2020.08.20: removed */
/* static FILE *LogFile = NULL; */

static char Buf[MAXTOK + 1];
static Matrix ErrorMat;

boolean ErrorSave = YES;
string ErrorVar = NULL;
int ErrorSeverityLevel = SEV_ERROR;
size_t ErrorTry = 0;

static string SeverityStr[] = SEVERITY_STRS;

/******************************+++********************************/
/*   void      Output(const string Format, ...)                  */
/*                                                               */
/*   Purpose:  Output variable argument list to stdout and to    */
/*             the log file.                                     */
/*             Error(), Fatal(), and Incompatibility() are       */
/*             similar except:                                   */
/*                                                               */
/*                                 precedes the output by        */
/*                                                               */
/*             Error()             "Error/Warning: "             */
/*             Fatal()             "Fatal error: "               */
/*             Incompatibility     "Incompatibility: "           */
/*                                                               */
/*             and Fatal() also calls exit(1).                   */
/*                                                               */
/* 1995.05.02:                                                   */
/* 2000.02.15: Output("\n") added to Fatal().                    */
/* 2023.05.17: vsprintf replaced by vsnprintf                    */
/* 2023.12.05: format argument added to Rprintf and Rf_error     */
/*****************************************************************/

void Output(const string Format, ...)
{
     va_list Args;

     va_start(Args, Format);

     vsnprintf(Buf, MAXTOK + 1, Format, Args);

     Rprintf("%s", Buf);

     va_end(Args);

     return;
}

void Error(const string Format, ...)
{
     va_list Args;

     va_start(Args, Format);

     if (ErrorSave)
          ErrorToMat(SeverityStr[ErrorSeverityLevel], Format, Args);

     va_end(Args);

     return;
}

void Fatal(const string Format, ...)
{
        string message = StrPaste(2, "Fatal error: ", Format);
        Rf_error("%s", message);
}


/******************************+++********************************/
void ErrorToMat(const string Severity, const string Format,
                va_list Args)
/*****************************************************************/
/*   Purpose:  Save error, warning, etc. messages in a matrix.   */
/*                                                               */
/* 1995.03.10:                                                   */
/* 2023.05.17: vsprintf replaced by vsnprintf                    */
/*****************************************************************/
{
     size_t j, nRowsOld, LastTry;
     size_t *Try;
     string LastMess, LastVar, TermPtr;
     string *Message, *Variable;

     if (!MatInitialized(&ErrorMat))
     {
          MatInit(RECT, MIXED, YES, &ErrorMat);
          MatPutText(&ErrorMat, "The following warning/error messages were generated:\n");
     }
     MatPutText(&ErrorMat, "The following warning/error messages were generated:\n");
     Variable = MatStrColFind(&ErrorMat, VARIABLE, NO);
     Try = MatSize_tColFind(&ErrorMat, "Try", NO);
     Message = MatStrColFind(&ErrorMat, "Message", NO);

     nRowsOld = MatNumRows(&ErrorMat);

     LastVar = (Variable != NULL) ? Variable[nRowsOld - 1] : NULL;
     LastTry = (Try != NULL) ? Try[nRowsOld - 1] : 0;
     LastMess = (Message != NULL) ? Message[nRowsOld - 1] : NULL;

     /* Put the message in Buf. */
     vsnprintf(Buf, MAXTOK + 1, Format, Args);

     /* Remove any terminating ".\n". */
     TermPtr = Buf + strlen(Buf) - 2;
     if (stricmp(TermPtr, ".\n") == 0)
          *TermPtr = NULL;

     if (stricmp(ErrorVar, LastVar) == 0 && ErrorTry == LastTry &&
         stricmp(Buf, LastMess) == 0)
          /* Do not repeat same message. */
          return;

     /* Allocate a new row for the new message. */
     MatReAlloc(nRowsOld + 1, MatNumCols(&ErrorMat), &ErrorMat);

     if (ErrorVar != NULL)
     {
          j = MatColumnAdd(VARIABLE, STRING, &ErrorMat);
          MatPutStrElem(&ErrorMat, nRowsOld, j, ErrorVar);
     }

     if (ErrorTry != 0)
     {
          j = MatColumnAdd("Try", SIZE_T, &ErrorMat);
          MatPutSize_tElem(&ErrorMat, nRowsOld, j, ErrorTry);
     }

     j = MatColumnAdd("Severity", STRING, &ErrorMat);
     MatPutStrElem(&ErrorMat, nRowsOld, j, Severity);

     j = MatColumnAdd("Message", STRING, &ErrorMat);
     MatPutStrElem(&ErrorMat, nRowsOld, j, Buf);

     return;
}

/******************************+++********************************/
void ErrorMatOut(void)
/*****************************************************************/
/*   Purpose:  Output matrix of error, warning, etc. messages.   */
/*                                                               */
/*   Version:  1995 May 11                                       */
/*****************************************************************/
{
     if (MatNumRows(&ErrorMat) > 0)
     {
          MatWriteBlock(&ErrorMat, YES);
          Output("\n");
          MatFree(&ErrorMat);
     }

     ErrorVar = NULL;
     ErrorTry = 0;
     ErrorSave = YES;
}
