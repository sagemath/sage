[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        sage_setup     \
        cython         \
        pkgconfig      \
        sagemath_environment \
        sagemath_categories \
                    ')]
build-backend = "setuptools.build_meta"
