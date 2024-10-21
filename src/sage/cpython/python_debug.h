/* sage_setup: distribution = sagemath-objects
 */
/*
 * Workaround to handle Python preprocessor macros: Translate defined /
 * undefined into true / false
 */


#include <Python.h>

#ifdef Py_DEBUG
    #define SAGE_Py_DEBUG Py_DEBUG
#else
    #define SAGE_Py_DEBUG 0
#endif

#ifdef Py_REF_DEBUG
    #define SAGE_Py_REF_DEBUG Py_REF_DEBUG
#else
    #define SAGE_Py_REF_DEBUG 0
#endif

#ifdef Py_TRACE_REFS
    #define SAGE_Py_TRACE_REFS Py_TRACE_REFS
    #define if_Py_TRACE_REFS_then_PyObject_INIT(obj, type) PyObject_INIT(obj, type)
#else
    #define SAGE_Py_TRACE_REFS 0
    #define if_Py_TRACE_REFS_then_PyObject_INIT(obj, type)
#endif

#ifdef PYMALLOC_DEBUG
    #define SAGE_PYMALLOC_DEBUG PYMALLOC_DEBUG
#else
    #define SAGE_PYMALLOC_DEBUG 0
#endif

#ifdef WITH_PYMALLOC
    #define SAGE_WITH_PYMALLOC PYMALLOC_DEBUG
#else
    #define SAGE_WITH_PYMALLOC 0
#endif



