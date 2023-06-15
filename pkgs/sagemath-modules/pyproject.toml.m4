[build-system]
# Minimum requirements for the build system to execute.
#
# Note we include numpy here to build some modules that cimport numpy,
# but it is not part of the install-requires.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        wheel          \
        sage_setup     \
        pkgconfig      \
        sagemath_environment \
        sagemath_categories \
        cython         \
        gmpy2          \
        numpy          \
        cysignals      \
        memory_allocator   \
                    ')]
build-backend = "setuptools.build_meta"
