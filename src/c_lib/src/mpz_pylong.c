/* mpz <-> pylong conversion and "pythonhash" for mpz
 *
 * Author:  Gonzalo Tornaría <tornaria@math.utexas.edu>
 * Date:    March 2006
 * License: GPL v2 or later
 *
 * this is free software: if it breaks, you get to keep all the pieces

AUTHORS:
  -- David Harvey (2007-08-18): added mpz_get_pyintlong function

 */

#include "mpn_pylong.h"
#include "mpz_pylong.h"

/* mpz python hash */
long
mpz_pythonhash (mpz_srcptr z)
{
  long x = mpn_pythonhash(z->_mp_d, abs(z->_mp_size));
  if (z->_mp_size < 0)
    x = -x;
  if (x == -1)
    x = -2;
  return x;
}

/* mpz -> pylong conversion */
PyObject *
mpz_get_pylong(mpz_srcptr z)
{
  py_size_t size = mpn_pylong_size(z->_mp_d, abs(z->_mp_size));
  PyLongObject *l = PyObject_NEW_VAR(PyLongObject, &PyLong_Type, size);

  if (l != NULL)
  {
    mpn_get_pylong(l->ob_digit, size, z->_mp_d, abs(z->_mp_size));
    if (z->_mp_size < 0)
      Py_SIZE(l) = -Py_SIZE(l);
  }

  return (PyObject *) l;
}

/* mpz -> pyint/pylong conversion; if the value fits in a python int, it
returns a python int (optimised for that pathway), otherwise returns
a python long */
PyObject *
mpz_get_pyintlong(mpz_srcptr z)
{
  if (mpz_fits_slong_p(z))
     return PyInt_FromLong(mpz_get_si(z));

  return mpz_get_pylong(z);
}

/* pylong -> mpz conversion */
int
mpz_set_pylong(mpz_ptr z, PyObject * ll)
{
  register PyLongObject * l = (PyLongObject *) ll;
  mp_size_t size;
  int i;

  if (l==NULL || !PyLong_Check(l)) {
    PyErr_BadInternalCall();
    return -1;
  }

  size = mpn_size_from_pylong(l->ob_digit, abs(Py_SIZE(l)));

  if (z->_mp_alloc < size)
    _mpz_realloc (z, size);

  mpn_set_pylong(z->_mp_d, size, l->ob_digit, abs(Py_SIZE(l)));
  z->_mp_size = Py_SIZE(l) < 0 ? -size : size;

  return size;
}

