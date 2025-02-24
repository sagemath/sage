import sys


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
        PY_MAJOR_VERSION=sys.version_info[0]
    )
