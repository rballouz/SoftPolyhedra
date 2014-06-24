#ifndef FLOATTYPE_INCLUDED
#define FLOATTYPE_INCLUDED

/*
** NOTE: With FPE trapping for Linux/glibc 2.2 and later, problems can
** arise when dealing with machine infinity.  In particular, HUGE_VAL
** (which is actually a macro function) should be used instead of
** DBL_MAX (a fixed constant). and there should be no math ops with
** infinities, like 0*HUGE_VAL.  Unfortunately there does not appear
** to be an equivalent preferred float-precision infinity, so use
** SVID's HUGE macro for now.
*/

#if TRAP_FPE/* || __STDC_VERSION__ >= 199901L*/
#include <math.h>
#undef DBL_MAX
#define DBL_MAX HUGE_VAL
#undef FLT_MAX
#define FLT_MAX HUGE
#else
#include <limits.h>
#include <float.h>
#ifndef FLT_MAX
#define FLT_MAX 3.402823466E+38F
#endif
#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157E+308
#endif
#endif

/*
 * If you change these, see also xdr_FLOAT() in outtype.c
 */
#ifndef SINGLE
#define FLOAT double
#define FLOAT_MAXVAL DBL_MAX
#else
#define FLOAT float
#define FLOAT_MAXVAL FLT_MAX
#endif

#endif
