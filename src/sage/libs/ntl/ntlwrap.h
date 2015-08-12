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
struct ZZ;
struct ZZ_p;
struct ZZX;
struct ZZ_pX;
struct zz_p;
struct zz_pX;
struct ZZ_pEContext;
struct ZZ_pE;
struct ZZ_pEX;
struct ZZ_pE;
struct GF2X;
struct GF2EContext;
struct GF2X_c;
struct GF2E;
struct GF2;
typedef struct {} mat_ZZ;
typedef struct {} mat_GF2;
typedef struct {} mat_GF2E;
#endif

#define zz_p_set_from_long( obj1, obj2 )\
        ((obj1) = (obj2))
#define NTL_zz_p_DOUBLE_EQUALS( obj1, obj2 )\
        ((obj1) == (obj2))
#define NTL_zz_pX_DOUBLE_EQUALS( obj1, obj2 )\
        ((obj1) == (obj2))

#endif
