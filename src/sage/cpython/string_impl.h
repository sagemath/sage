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
#if PY_MAJOR_VERSION <= 2
    return PyBytes_FromString(c);
#else
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
#endif
}


static inline PyObject* _str_to_bytes(PyObject* s, PyObject* encoding, PyObject* errors)
{
#if PY_MAJOR_VERSION <= 2
    /* On Python 2, we accept bytes == str as input */
    if (PyBytes_CheckExact(s)) {
        Py_INCREF(s);
        return s;
    }
#endif

    if (!PyUnicode_Check(s)) {
        PyErr_Format(PyExc_TypeError,
#if PY_MAJOR_VERSION >= 3
                     "expected str, %s found",
#else
                     "expected str or unicode, %s found",
#endif
                     Py_TYPE(s)->tp_name);
        return NULL;
    }

    const char* err = NULL;  // Default: strict
    const char* enc = NULL;  // Default: utf-8

#if PY_MAJOR_VERSION <= 2
    if (errors != Py_None) {
        err = PyString_AsString(errors);
        if (!err) return NULL;
    }

    if (encoding != Py_None) {
        enc = PyString_AsString(encoding);
        if (!enc) return NULL;
    }
    else {
        enc = "utf-8";
    }
#else
    if (errors != Py_None) {
        err = PyUnicode_AsUTF8(errors);
        if (!err) return NULL;
    }

    if (encoding != Py_None) {
        enc = PyUnicode_AsUTF8(encoding);
        if (!enc) return NULL;
    }
#endif

    return PyUnicode_AsEncodedString(s, enc, err);
}
