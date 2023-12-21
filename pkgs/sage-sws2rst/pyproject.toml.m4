include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
]
build-backend = "setuptools.build_meta"

[project]
name = "sage-sws2rst"
description = "Sage: Open Source Mathematics Software: SageNB worksheet converter"
license = {text = "GNU General Public License (GPL) v3 or later"}
authors = [{name = "The Sage Developers", email = "sage-support@googlegroups.com"}]
urls = {Homepage = "https://www.sagemath.org"}
dynamic = ["version"]

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
script-files = ["bin/sage-sws2rst"]
include-package-data = false

[tool.setuptools.packages]
find = {namespaces = false}

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
