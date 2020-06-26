SAGE_SPKG_CONFIGURE([python3], [
   SAGE_SPKG_DEPCHECK([sqlite libpng bzip2 xz libffi], [
      AS_IF([test $SAGE_PYTHON_VERSION = 2], [
        dnl If we use Python 2 for Sage, we install Python 3 too and do NOT attempt to do
        dnl venv using system python3 over SAGE_LOCAL.
        dnl (In particular, the setuptools and pip install scripts are not prepared for
        dnl handling this situation.)
        sage_spkg_install_python3=yes
      ], [
        dnl Using Python 3 for Sage.  Check if we can do venv with a system python3
        dnl instead of building our own copy.
        check_modules="sqlite3, ctypes, math, hashlib, crypt, readline, socket, zlib, distutils.core"
        AC_CACHE_CHECK([for python3 >= 3.7.3, < 3.8 with modules $check_modules], [ac_cv_path_PYTHON3], [
            AC_MSG_RESULT([])
            AC_PATH_PROGS_FEATURE_CHECK([PYTHON3], [python3.7 python3], [
                AC_MSG_CHECKING([... whether $ac_path_PYTHON3 is good])
                python3_version=`"$ac_path_PYTHON3" --version 2>&1 \
                    | $SED -n -e 's/\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\).*/\1/p'`
                AS_IF([test -n "$python3_version"], [
                    AX_COMPARE_VERSION([$python3_version], [ge], [3.7.3], [
                        AX_COMPARE_VERSION([$python3_version], [lt], [3.8.0], [
                            dnl Because the system python is not used directly but rather in a venv without site-packages,
                            dnl we test whether the module will be available in a venv.
                            dnl Otherwise, some system site-package may be providing this module to the system python.
                            dnl m4_define([conftest_venv], [config-venv]) .... for debugging only
                            rm -rf conftest_venv
                            AS_IF(["$ac_path_PYTHON3" build/bin/sage-venv conftest_venv && conftest_venv/bin/python3 -c "import $check_modules"], [
                                AC_LANG_PUSH([C])
                                AC_LANG_CONFTEST([
                                    AC_LANG_SOURCE([[
#define PY_SSIZE_T_CLEAN
#include <Python.h>
static PyMethodDef SpamMethods[] = {
    {NULL, NULL, 0, NULL}        /* Sentinel */
};
static struct PyModuleDef spammodule = {
    PyModuleDef_HEAD_INIT,
    "spam",   /* name of module */
    NULL,     /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    SpamMethods
};
PyMODINIT_FUNC
PyInit_spam(void)
{
    PyObject *m;

    m = PyModule_Create(&spammodule);
    return m;
}
                                    ]])
                                ])
                                AC_LANG_POP([C])
                                cat > conftest.py <<EOF
from distutils.core import setup
from distutils.extension import Extension
from sys import exit
modules = list((Extension("config_check_distutils", list(("conftest.c",))),))
setup(name="config_check_distutils", ext_modules=modules)
exit(0)
EOF
                                dnl (echo "***ENV***:"; env; echo "***SYSCONFIG***"; conftest_venv/bin/python3 -m sysconfig) >& AS_MESSAGE_LOG_FD
                                echo CC="$CC" CXX="$CXX" conftest_venv/bin/python3 conftest.py --verbose build --build-base=conftest.dir >& AS_MESSAGE_LOG_FD
                                AS_IF([CC="$CC" CXX="$CXX" conftest_venv/bin/python3 conftest.py --verbose build --build-base=conftest.dir >& AS_MESSAGE_LOG_FD 2>&1 ], [
                                    rm -rf conftest.*
                                    AC_LANG_PUSH([C++])
                                    AC_LANG_CONFTEST([
                                        AC_LANG_SOURCE([[
#define PY_SSIZE_T_CLEAN
#include <Python.h>
static PyMethodDef SpamMethods[] = {
    {NULL, NULL, 0, NULL}        /* Sentinel */
};
static struct PyModuleDef spammodule = {
    PyModuleDef_HEAD_INIT,
    "spam",   /* name of module */
    NULL,     /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    SpamMethods
};
PyMODINIT_FUNC
PyInit_spam(void)
{
    PyObject *m;

    m = PyModule_Create(&spammodule);
    return m;
}
// Partial C++11 test, from ax_cxx_compile_stdcxx.m4

  namespace test_noexcept
  {

    int f() { return 0; }
    int g() noexcept { return 0; }

    static_assert(noexcept(f()) == false, "");
    static_assert(noexcept(g()) == true, "");

  }

                                        ]])
                                    ])
                                    AC_LANG_POP([C++])
                                    cat > conftest.py <<EOF
from distutils.core import setup
from distutils.extension import Extension
from sys import exit
modules = list((Extension("config_check_distutils_cxx", list(("conftest.cpp",)),
                          extra_compile_args=list(("-std=c++11",)), language="c++"),))
setup(name="config_check_distutils_cxx", ext_modules=modules)
exit(0)
EOF
                                    AS_IF([CC="$CC" CXX="$CXX" conftest_venv/bin/python3 conftest.py --verbose build --build-base=conftest.dir >& AS_MESSAGE_LOG_FD 2>&1 ], [
                                        ac_cv_path_PYTHON3="$ac_path_PYTHON3"
                                        ac_path_PYTHON3_found=:
                                        AC_MSG_RESULT([yes])
                                        dnl introduction for AC_MSG_RESULT printed by AC_CACHE_CHECK
                                        AC_MSG_CHECKING([for python3 >= 3.7.3, < 3.8 with modules $check_modules])
                                    ], [
                                        AC_MSG_RESULT([no, the version is in the supported range, and the modules can be imported, but distutils cannot build a C++ 11 extension])
                                    ])
                                ], [
                                    AC_MSG_RESULT([no, the version is in the supported range, and the modules can be imported, but distutils cannot build a C extension])
                                ])
                            ], [
                                AC_MSG_RESULT([no, the version is in the supported range but cannot import one of the required modules: $check_modules])
                            ])
                        ], [
                            AC_MSG_RESULT([no, $python3_version is too recent])
                        ])
                    ], [
                        AC_MSG_RESULT([no, $python3_version is too old])
                    ])
                ], [
                    AC_MSG_RESULT([no, "$ac_path_PYTHON3 --version" does not work])
                ])
            ])
        ])
        AS_IF([test -z "$ac_cv_path_PYTHON3"],
              [sage_spkg_install_python3=yes])
      ])
   ])
],, [
    dnl PRE
], [
    dnl POST
    AS_IF([test x$sage_spkg_install_python3 = xno], [PYTHON_FOR_VENV="$ac_cv_path_PYTHON3"])
    AC_SUBST([PYTHON_FOR_VENV])

    dnl These temporary directories are created by the check above
    dnl and need to be cleaned up to prevent the "rm -f conftest*"
    dnl (that a bunch of other checks do) from emitting warnings about
    dnl conftest.dir and conftest_venv being directories.
    rm -rf conftest.dir conftest_venv
])
