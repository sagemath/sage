/* sage_setup: distribution = sagemath-objects
 */
/*****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.bray@lri.fr>
#       Copyright (C) 2018 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
*****************************************************************************/

#include <Python.h>
#include <string.h>


static inline PyObject* _cstr_to_str(const char* c, PyObject* encoding, PyObject* errors)
{
    const char* err = NULL;  // Default: strict
    const char* enc = NULL;  // Default: utf-8

    if (errors != Py_None) {
        err = PyUnicode_AsUTF8(errors);
        if (!err) return NULL;
    }

    if (encoding != Py_None) {
        enc = PyUnicode_AsUTF8(encoding);
        if (!enc) return NULL;
    }

    return PyUnicode_Decode(c, strlen(c), enc, err);
}


static inline PyObject* _str_to_bytes(PyObject* s, PyObject* encoding, PyObject* errors)
{
    if (!PyUnicode_Check(s)) {
        PyErr_Format(PyExc_TypeError,
                     "expected str, %s found",
                     Py_TYPE(s)->tp_name);
        return NULL;
    }

    const char* err = NULL;  // Default: strict
    const char* enc = NULL;  // Default: utf-8

    if (errors != Py_None) {
        err = PyUnicode_AsUTF8(errors);
        if (!err) return NULL;
    }

    if (encoding != Py_None) {
        enc = PyUnicode_AsUTF8(encoding);
        if (!enc) return NULL;
    }

    return PyUnicode_AsEncodedString(s, enc, err);
}
