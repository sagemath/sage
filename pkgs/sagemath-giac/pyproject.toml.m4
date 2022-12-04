[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        wheel          \
        sage_setup     \
        sagemath_environment \
        sagemath_categories \
        cython         \
        gmpy2          \
        cysignals      \
                    ')]
build-backend = "setuptools.build_meta"
