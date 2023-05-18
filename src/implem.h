/*****************************************************************/
/*   Implementation-dependent definitions                        */
/*                                                               */
/*   Copyright (c) William J. Welch 1991--99.                    */
/*   All rights reserved.                                        */
/*                                                               */
/*   2022.10.10: Define SIZE_T_MAX only if undefined             */
/*****************************************************************/

#ifndef DBL_EPSILON
     #include <float.h>
#endif

#ifndef SIZE_T_MAX
     #define SIZE_T_MAX  UINT_MAX
#endif

#ifndef UINT_MAX
     #include <limits.h>
#endif

#ifdef NULL
     #undef NULL
#endif
#define NULL             0L

typedef double      real;
#define EPSILON     DBL_EPSILON    /* Machine precision. */
#ifndef REAL_MAX
     #define REAL_MAX    DBL_MAX
#endif
#define REAL_MIN    DBL_MIN
#define LN_MAX      (700.0)        /* Approx. ln(REAL_MAX) */
#define LN_MIN      (-700.0)       /* Approx. ln(REAL_MIN) */

#define OUTPUT_COLS  80   /* Line width (in characters) for output. */
#define MAXTOK      256   /* Maximum length for filenames, */
                          /* converted numbers, etc.       */

#define PRECISION     6   /* Default precision for %g, etc. */

