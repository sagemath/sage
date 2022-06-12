[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        wheel          \
                    ')]
build-backend = "setuptools.build_meta"
