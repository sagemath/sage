/* sage_setup: distribution = sagemath-categories
 */
#ifndef _SAGE_FINITE_RINGS_INTEGER_MOD_LIMITS_H
#define _SAGE_FINITE_RINGS_INTEGER_MOD_LIMITS_H

#include <stdint.h>

#define INTEGER_MOD_INT32_LIMIT 46341          //  = ceil(sqrt(2^31-1))
#define INTEGER_MOD_INT64_LIMIT 2147483647     //  = 2^31-1 for now, should be 3037000500LL = ceil(sqrt(2^63-1))

#endif
