/* sage_setup: distribution = sagemath-objects
 */
/* 3-argument version of Py_VISIT, easier to use from Cython */

#define Py_VISIT3(op, visit, arg)                                       \
    do {                                                                \
        if (op) {                                                       \
            int vret = visit((PyObject *)(op), arg);                    \
            if (vret)                                                   \
                return vret;                                            \
        }                                                               \
    } while (0)

