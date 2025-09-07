# SAGE_CHECK_PYTHON_FOR_VENV(PYTHON_EXE, MIN_VERSION, LT_VERSION, REQUIRED_MODULES, COMMANDS_IF_GOOD)

AC_DEFUN([SAGE_CHECK_PYTHON_FOR_VENV], [
                m4_pushdef([PYTHON_EXE],       [$1])
                m4_pushdef([MIN_VERSION],      [$2])
                m4_pushdef([LT_VERSION],       [$3])
                m4_pushdef([REQUIRED_MODULES], [$4])
                m4_pushdef([COMMANDS_IF_GOOD], [$5])

                AC_SUBST([SAGE_ARCHFLAGS])

                AC_MSG_CHECKING([... whether ]PYTHON_EXE[ is good])
                python3_version=`"PYTHON_EXE" --version 2>&1 \
                    | $SED -n -e 's/\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\).*/\1/p'`
                AS_IF([test -n "$python3_version"], [
                    AX_COMPARE_VERSION([$python3_version], [ge], MIN_VERSION, [
                        AX_COMPARE_VERSION([$python3_version], [lt], LT_VERSION, [
                            dnl Because the system python is not used directly but rather in a venv without site-packages,
                            dnl we should test whether the module will be available in a venv.
                            dnl Otherwise, some system site-package may be providing this module to the system python.
                            dnl However, on Python >= 3.12, we need setuptools to run our extension compilation tests
                            dnl because distutils has been removed from the standard library.
                            AX_COMPARE_VERSION([$python3_version], [ge], [3.12.0], [
                                conftest_venv_options="--system-site-packages"
                                distutils_core="setuptools"
                                distutils_extension="setuptools.extension"
                            ], [
                                conftest_venv_options=
                                distutils_core="distutils.core"
                                distutils_extension="distutils.extension"
                            ])
                            all_required_modules="]REQUIRED_MODULES[, $distutils_core, $distutils_extension"
                            dnl m4_define([conftest_venv], [config-venv]) .... for debugging only
                            rm -rf conftest_venv
                            AS_IF(["]PYTHON_EXE[" build/bin/sage-venv $conftest_venv_options conftest_venv && conftest_venv/bin/python3 -c "import $all_required_modules" 2>& ]AS_MESSAGE_LOG_FD, [
                                AS_VAR_SET([python3_result], [yes])
                                SAGE_PYTHON_CHECK_DISTUTILS([CC="$CC" CXX="$CXX" conftest_venv/bin/python3], [
                                    SAGE_ARCHFLAGS="unset"
                                ], [
                                    AS_CASE([$host],
                                        [*-*-darwin*], [
                                            dnl #31227: Try if setting ARCHFLAGS to empty fixes it
                                            SAGE_PYTHON_CHECK_DISTUTILS([CC="$CC" CXX="$CXX" ARCHFLAGS="" conftest_venv/bin/python3], [
                                                SAGE_ARCHFLAGS=""
                                            ], [
                                                AS_VAR_SET([python3_result],
                                                           ["no, the version is in the supported range, and the modules can be imported, but $reason (even with ARCHFLAGS set to empty)"])
                                            ])
                                        ], [
                                            AS_VAR_SET([python3_result],
                                                       ["no, the version is in the supported range, and the modules can be imported, but $reason"])
                                        ])
                                ])
                                AS_VAR_IF([python3_result], [yes], [
                                    dnl these commands are expected to call AC_MSG_RESULT
                                    COMMANDS_IF_GOOD
                                ], [
                                    AC_MSG_RESULT([$python3_result])
                                ])
                            ], [
                                AC_MSG_RESULT([no, the version is in the supported range but cannot import one of the required modules: $all_required_modules])
                            ])
                        ], [
                            AC_MSG_RESULT([no, $python3_version is too recent])
                        ])
                    ], [
                        AC_MSG_RESULT([no, $python3_version is too old])
                    ])
                ], [
                    AC_MSG_RESULT([no, "]PYTHON_EXE[ --version" does not work])
                ])


                m4_popdef([PYTHON_EXE])
                m4_popdef([MIN_VERSION])
                m4_popdef([LT_VERSION])
                m4_popdef([REQUIRED_MODULES])
                m4_popdef([COMMANDS_IF_GOOD])

])

dnl distutils test
AC_DEFUN([SAGE_PYTHON_CHECK_DISTUTILS], [
    m4_pushdef([PYTHON_EXE], [$1])
    m4_pushdef([COMMANDS_IF_DISTUTILS_GOOD], [$2])
    m4_pushdef([COMMANDS_IF_DISTUTILS_NOT_GOOD], [$3])
    SAGE_PYTHON_DISTUTILS_C_CONFTEST
    dnl (echo "***ENV***:"; env; echo "***SYSCONFIG***"; conftest_venv/bin/python3 -m sysconfig) >& AS_MESSAGE_LOG_FD
    echo PYTHON_EXE conftest.py --verbose build --build-base=conftest.dir >& AS_MESSAGE_LOG_FD
    AS_IF([PYTHON_EXE conftest.py --verbose build --build-base=conftest.dir >& AS_MESSAGE_LOG_FD 2>&1 ], [
        SAGE_PYTHON_DISTUTILS_CXX_CONFTEST
        echo PYTHON_EXE conftest.py --verbose build --build-base=conftest.dir >& AS_MESSAGE_LOG_FD 2>&1
        AS_IF([PYTHON_EXE conftest.py --verbose build --build-base=conftest.dir >& AS_MESSAGE_LOG_FD 2>&1 ], [
            COMMANDS_IF_DISTUTILS_GOOD], [
            reason="distutils cannot build a C++ 11 extension"
            COMMANDS_IF_DISTUTILS_NOT_GOOD
        ])
    ], [
       reason="distutils cannot build a C extension"
       COMMANDS_IF_DISTUTILS_NOT_GOOD
    ])
    m4_popdef([PYTHON_EXE])
    m4_popdef([COMMANDS_IF_DISTUTILS_GOOD])
    m4_popdef([COMMANDS_IF_DISTUTILS_NOT_GOOD])
])

dnl Write conftest.py and conftest.c
AC_DEFUN([SAGE_PYTHON_DISTUTILS_C_CONFTEST], [
                                rm -rf conftest.*
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

from $distutils_core import setup
from $distutils_extension import Extension
from sys import exit
modules = list((Extension("config_check_distutils", list(("conftest.c",))),))
setup(name="config_check_distutils", ext_modules=modules)
exit(0)
EOF
])

dnl Write conftest.py and conftest.cpp
AC_DEFUN([SAGE_PYTHON_DISTUTILS_CXX_CONFTEST], [
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

from $distutils_core import setup
from $distutils_extension import Extension
from sys import exit
modules = list((Extension("config_check_distutils_cxx", list(("conftest.cpp",)),
                          extra_compile_args=list(("-std=c++11",)), language="c++"),))
setup(name="config_check_distutils_cxx", ext_modules=modules, py_modules=["conftest"])
exit(0)
EOF
])
