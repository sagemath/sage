## Pipfile with all dependencies of sagelib and version information as free as possible
## (for developers to generate a dev environment)
[[source]]
name = "pypi"
url = "https://pypi.org/simple"
verify_ssl = true

[dev-packages]
## We do not list packages that are already declared as dependencies (install_requires)
## in pyproject.toml
pycodestyle = "*"
tox = "*"
pytest = "*"
rope = "*"
six = "*"

[packages]
## We do not list packages that are already declared as dependencies (install_requires)
## in pyproject.toml

[packages.e1839a8]
path = "."
