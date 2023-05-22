[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        sage_setup     \
        sagemath_environment \
        cython         \
        cysignals      \
                    ')]
build-backend = "setuptools.build_meta"
