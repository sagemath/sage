/* sage_setup: distribution = sagemath-ntl
 */
#ifndef SAGE_NTLWRAP_H
#define SAGE_NTLWRAP_H

#ifdef __cplusplus
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>
#include <NTL/HNF.h>
#include <NTL/LLL.h>
using namespace NTL;
#endif

#ifndef __cplusplus
/* Ensure that code passing around pointers to NTL classes can be
 * compiled in pure C mode. */
typedef struct {} ZZ;
typedef struct {} ZZ_p;
typedef struct {} ZZX;
typedef struct {} ZZ_pX;
typedef struct {} zz_p;
typedef struct {} zz_pX;
typedef struct {} ZZ_pEContext;
typedef struct {} ZZ_pE;
typedef struct {} ZZ_pEX;
typedef struct {} GF2X;
typedef struct {} GF2EContext;
typedef struct {} GF2X_c;
typedef struct {} GF2E;
typedef struct {} GF2;
typedef struct {} mat_ZZ;
typedef struct {} mat_GF2;
typedef struct {} mat_GF2E;
#endif

#endif
