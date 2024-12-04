import os
import subprocess
import sys

import pkgconfig

from sage.env import (
    CBLAS_PC_MODULES,
    ECL_CONFIG,
    NTL_INCDIR,
    NTL_LIBDIR,
    OPENMP_CFLAGS,
    OPENMP_CXXFLAGS,
)


def compiler_directives(profile: bool):
    """
    Return a list of Cython directives used for compilation.
    """
    return dict(
        # Do not generate __reduce__ methods
        auto_pickle=False,
        # Do not create __test__ dictionary automatically from docstrings
        autotestdict=False,
        binding=False,
        c_api_binop_methods=True,
        # Do not check for division by 0 (this is about 35% quicker than with check)
        cdivision=True,
        cpow=True,
        # Embed a textual copy of the call signature in the docstring (to support tools like IPython)
        embedsignature=True,
        fast_getattr=True,
        # Use Python 3 (including source code semantics) for module compilation
        language_level="3",
        legacy_implicit_noexcept=True,
        # Enable support for late includes (make declarations in Cython code available to C include files)
        preliminary_late_includes_cy28=True,
        # Add hooks for Python profilers into the compiled C code
        profile=profile,
    )


def compile_time_env_variables():
    """
    Return a list of environmental variables used for compilation.
    """
    return dict(
        PY_PLATFORM=sys.platform,
        # The following two constants are here only for backwards compatibility of user packages
        PY_VERSION_HEX=sys.hexversion,
        PY_MAJOR_VERSION=sys.version_info[0],
    )


def get_cblas_pc_module_name() -> str:
    """
    Return the name of the BLAS libraries to be used.
    """
    cblas_pc_modules = CBLAS_PC_MODULES.split(":")
    return next(blas_lib for blas_lib in cblas_pc_modules if pkgconfig.exists(blas_lib))


default_required_modules = (
    "fflas-ffpack",
    "givaro",
    "gsl",
    "linbox",
    "Singular",
    "libpng",
    "gdlib",
    "m4ri",
    "zlib",
    "cblas",
    "ecl",
)

default_optional_modules = ("lapack",)


def cython_aliases(required_modules=None, optional_modules=None):
    """
    Return the aliases for compiling Cython code. These aliases are
    macros which can occur in ``# distutils`` headers.

    INPUT:

    - ``required_modules`` -- (default: taken from ``default_required_modules``)
      iterable of string values

    - ``optional_modules`` -- (default: taken from ``default_optional_modules``)
      iterable of string values

    EXAMPLES::

        sage: from sage.env import cython_aliases
        sage: cython_aliases()
        {...}
        sage: sorted(cython_aliases().keys())
        ['CBLAS_CFLAGS',
         ...,
         'ZLIB_LIBRARIES']
        sage: cython_aliases(required_modules=('module-that-is-assumed-to-not-exist'))
        Traceback (most recent call last):
        ...
        PackageNotFoundError: ...
        sage: cython_aliases(required_modules=(), optional_modules=('module-that-is-assumed-to-not-exist'))
        {...}

    TESTS:

    We can use ``cython.parallel`` regardless of whether OpenMP is supported.
    This will run in parallel, if OpenMP is supported::

        sage: cython(                                               # optional - sage.misc.cython
        ....: '''
        ....: #distutils: extra_compile_args = OPENMP_CFLAGS
        ....: #distutils: extra_link_args = OPENMP_CFLAGS
        ....: from cython.parallel import prange
        ....:
        ....: cdef int i
        ....: cdef int n = 30
        ....: cdef int sum = 0
        ....:
        ....: for i in prange(n, num_threads=4, nogil=True):
        ....:     sum += i
        ....:
        ....: print(sum)
        ....: ''')
        435
    """
    import itertools

    if required_modules is None:
        required_modules = default_required_modules

    if optional_modules is None:
        optional_modules = default_optional_modules

    aliases = {}

    for lib, required in itertools.chain(
        ((lib, True) for lib in required_modules),
        ((lib, False) for lib in optional_modules),
    ):
        var = lib.upper().replace("-", "") + "_"
        if lib == "cblas":
            lib = get_cblas_pc_module_name()
        if lib == "zlib":
            aliases[var + "CFLAGS"] = ""
            try:
                pc = pkgconfig.parse("zlib")
                libs = pkgconfig.libs(lib)
            except pkgconfig.PackageNotFoundError:
                from collections import defaultdict

                pc = defaultdict(list, {"libraries": ["z"]})
                libs = "-lz"
        elif lib == "ecl":
            try:
                # Determine ecl-specific compiler arguments using the ecl-config script
                ecl_cflags = subprocess.run(
                    [ECL_CONFIG, "--cflags"], check=True, capture_output=True, text=True
                ).stdout.split()
                ecl_libs = subprocess.run(
                    [ECL_CONFIG, "--libs"], check=True, capture_output=True, text=True
                ).stdout.split()
            except subprocess.CalledProcessError:
                if required:
                    raise
                else:
                    continue
            aliases["ECL_CFLAGS"] = list(
                filter(lambda s: not s.startswith("-I"), ecl_cflags)
            )
            aliases["ECL_INCDIR"] = [
                s[2:] for s in filter(lambda s: s.startswith("-I"), ecl_cflags)
            ]
            aliases["ECL_LIBDIR"] = [
                s[2:] for s in filter(lambda s: s.startswith("-L"), ecl_libs)
            ]
            aliases["ECL_LIBRARIES"] = [
                s[2:] for s in filter(lambda s: s.startswith("-l"), ecl_libs)
            ]
            aliases["ECL_LIBEXTRA"] = list(
                filter(lambda s: not s.startswith(("-l", "-L")), ecl_libs)
            )
            continue
        else:
            try:
                aliases[var + "CFLAGS"] = pkgconfig.cflags(lib).split()
                pc = pkgconfig.parse(lib)
                libs = pkgconfig.libs(lib)
            except pkgconfig.PackageNotFoundError:
                if required:
                    raise
                else:
                    continue

        # It may seem that INCDIR is redundant because the -I options are also
        # passed in CFLAGS.  However, "extra_compile_args" are put at the end
        # of the compiler command line.  "include_dirs" go to the front; the
        # include search order matters.
        aliases[var + "INCDIR"] = pc["include_dirs"]
        aliases[var + "LIBDIR"] = pc["library_dirs"]
        aliases[var + "LIBEXTRA"] = list(
            filter(lambda s: not s.startswith(("-l", "-L")), libs.split())
        )
        aliases[var + "LIBRARIES"] = pc["libraries"]

    # uname-specific flags
    UNAME = os.uname()

    def uname_specific(name, value, alternative):
        if name in UNAME[0]:
            return value
        else:
            return alternative

    aliases["LINUX_NOEXECSTACK"] = uname_specific("Linux", ["-Wl,-z,noexecstack"], [])

    # LinBox needs special care because it actually requires C++11 with
    # GNU extensions: -std=c++11 does not work, you need -std=gnu++11
    # (this is true at least with GCC 7.2.0).
    #
    # Further, note that LinBox does not add any C++11 flag in its .pc
    # file (possibly because of confusion between CFLAGS and CXXFLAGS?).
    # This is not a problem in practice since LinBox depends on
    # fflas-ffpack and fflas-ffpack does add such a C++11 flag.
    if "LINBOX_CFLAGS" in aliases:
        aliases["LINBOX_CFLAGS"].append("-std=gnu++11")

    try:
        aliases["M4RI_CFLAGS"].remove("-pedantic")
    except (ValueError, KeyError):
        pass

    # NTL
    aliases["NTL_CFLAGS"] = ["-std=c++11"]
    aliases["NTL_INCDIR"] = [NTL_INCDIR] if NTL_INCDIR else []
    aliases["NTL_LIBDIR"] = [NTL_LIBDIR] if NTL_LIBDIR else []
    aliases["NTL_LIBRARIES"] = ["ntl"]
    aliases["NTL_LIBEXTRA"] = []

    # OpenMP
    aliases["OPENMP_CFLAGS"] = OPENMP_CFLAGS.split()
    aliases["OPENMP_CXXFLAGS"] = OPENMP_CXXFLAGS.split()

    return aliases
