/* sage_setup: distribution = sagemath-objects
 */
#include <Python.h>
#include <gmp.h>

static PyTypeObject* Integer = NULL;


/* Replication of the internal structure of a Sage Integer */
struct Integer_Object
{
    PyObject_HEAD;
    void* __pyx_vtab;
    void* _parent;
    mpz_t value;
};


static inline mpz_ptr Integer_AS_MPZ(PyObject* x)
{
    return ((struct Integer_Object*)x)->value;
}
