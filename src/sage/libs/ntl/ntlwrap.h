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
typedef class {} ZZ;
typedef class {} ZZ_p;
typedef class {} ZZX;
typedef class {} ZZ_pX;
typedef class {} zz_p;
typedef class {} zz_pX;
typedef class {} ZZ_pEContext;
typedef class {} ZZ_pE;
typedef class {} ZZ_pEX;
typedef class {} GF2X;
typedef class {} GF2EContext;
typedef class {} GF2X_c;
typedef class {} GF2E;
typedef class {} GF2;
typedef class {} mat_ZZ;
typedef class {} mat_GF2;
typedef class {} mat_GF2E;
#endif

#endif
