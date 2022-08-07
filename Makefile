# Main Makefile for Sage.

# The default target ("all") builds Sage and the whole (HTML) documentation.
#
# Target "build" just builds Sage.
#
# See below for targets to build the documentation in other formats,
# to run various types of test suites, and to remove parts of the build etc.

default: all

all: base-toolchain
	$(MAKE) all-start

build: base-toolchain
	$(MAKE) all-build

build-local: base-toolchain
	$(MAKE) all-build-local

build-venv: base-toolchain
	$(MAKE) all-build-venv

start: base-toolchain
	$(MAKE) build-start

sageruntime: base-toolchain
	$(MAKE) all-sageruntime

SAGE_ROOT_LOGS = logs

# The --stop flag below is just a random flag to induce graceful
# breakage with non-GNU versions of make.
# See https://trac.sagemath.org/ticket/24617

# Defer unknown targets to build/make/Makefile
%::
	@if [ -x relocate-once.py ]; then ./relocate-once.py; fi
	$(MAKE) build/make/Makefile --stop
	+build/bin/sage-logger \
		"cd build/make && ./install '$@'" logs/install.log

# CONFIG_FILES lists all files that appear in AC_CONFIG_FILES in configure.ac;
# except for build/make/Makefile-auto, which is unused by the build system
CONFIG_FILES = build/make/Makefile src/bin/sage-env-config build/bin/sage-build-env-config pkgs/sage-conf/_sage_conf/_conf.py

# SPKG_COLLECT_FILES contains all files that influence the SAGE_SPKG_COLLECT macro
SPKG_COLLECT_FILES = build/pkgs/*/type build/pkgs/*/package-version.txt build/pkgs/*/dependencies build/pkgs/*/requirements.txt build/pkgs/*/checksums.ini build/pkgs/*/spkg-install

# If configure was run before, rerun it with the old arguments.
# Otherwise, run configure with argument $PREREQ_OPTIONS.
build/make/Makefile: configure $(SPKG_COLLECT_FILES) $(CONFIG_FILES:%=%.in)
	rm -f config.log
	mkdir -p logs/pkgs
	ln -s logs/pkgs/config.log config.log
	@if [ -x config.status ]; then \
		./config.status --recheck && ./config.status; \
	else \
		echo >&2 '****************************************************************************'; \
		echo >&2 'error: Sage source tree is unconfigured. Please run "./configure" first.'; \
		echo >&2 'note:  Type "./configure --help" to see the available configuration options.'; \
		echo >&2 '****************************************************************************'; \
	        exit 1; \
	fi

# This is used to monitor progress towards Python 3 and prevent
# regressions. Should be removed after the full switch to python3.
#
# As of Sage 9.0: keep the build target for backward compatibility,
# but it just runs "make".
buildbot-python3:
	$(MAKE)

# Preemptively download all source tarballs of normal packages.
download:
	export SAGE_ROOT=$$(pwd) && \
	export PATH=$$SAGE_ROOT/build/bin:$$PATH && \
	sage-package download :all:

dist: build/make/Makefile
	./sage --sdist

pypi-sdists: sage_setup
	./sage --sh build/pkgs/sage_conf/spkg-src
	./sage --sh build/pkgs/sage_sws2rst/spkg-src
	./sage --sh build/pkgs/sage_docbuild/spkg-src
	./sage --sh build/pkgs/sage_setup/spkg-src
	./sage --sh build/pkgs/sagelib/spkg-src
	./sage --sh build/pkgs/sagemath_objects/spkg-src
	./sage --sh build/pkgs/sagemath_categories/spkg-src
	./sage --sh build/pkgs/sagemath_environment/spkg-src
	./sage --sh build/pkgs/sagemath_repl/spkg-src
	@echo "Built sdists are in upstream/"

# Ensuring wheels are present, even for packages that may have been installed
# as editable. Until we have better uninstallation of script packages, we
# just remove the timestamps, which will lead to rebuilds of the packages.
PYPI_WHEEL_PACKAGES = sage_sws2rst sage_setup sagemath_environment sagemath_objects sagemath_repl sagemath_categories
pypi-wheels:
	for a in $(PYPI_WHEEL_PACKAGES); do \
	    rm -f venv/var/lib/sage/installed/$$a-*; \
	done
	for a in $(PYPI_WHEEL_PACKAGES); do \
	    $(MAKE) SAGE_EDITABLE=no $$a; \
	done
	@echo "Built wheels are in venv/var/lib/sage/wheels/"

# sage_docbuild is here, not in PYPI_WHEEL_PACKAGES, because it depends on sagelib
WHEEL_PACKAGES = $(PYPI_WHEEL_PACKAGES) sage_conf sagelib sage_docbuild
wheels:
	for a in $(WHEEL_PACKAGES); do \
	    rm -f venv/var/lib/sage/installed/$$a-*; \
	done
	for a in $(WHEEL_PACKAGES); do \
	    $(MAKE) SAGE_EDITABLE=no $$a; \
	done
	@echo "Built wheels are in venv/var/lib/sage/wheels/"

###############################################################################
# Cleaning up
###############################################################################

SAGE_ROOT = $(CURDIR)
SAGE_SRC = $(SAGE_ROOT)/src

clean:
	@echo "Deleting package build directories..."
	if [ -f "$(SAGE_SRC)"/bin/sage-env-config ]; then \
	    . "$(SAGE_SRC)"/bin/sage-env-config; \
	    if [ -n "$$SAGE_LOCAL" ]; then \
	        rm -rf "$$SAGE_LOCAL/var/tmp/sage/build"; \
	    fi; \
	fi

# "c_lib", ".cython_version", "build" in $(SAGE_SRC) are from old sage versions
# Cleaning .so files (and .c and .cpp files associated with .pyx files) is for editable installs.
# Also cython_debug is for editable installs.
sagelib-clean:
	@echo "Deleting Sage library build artifacts..."
	if [ -d "$(SAGE_SRC)" ]; then \
	    (cd "$(SAGE_SRC)" && \
	     rm -rf c_lib .cython_version cython_debug; \
	     rm -rf build; find . -name '*.pyc' -o -name "*.so" | xargs rm -f; \
	     rm -f $$(find . -name "*.pyx" | sed 's/\(.*\)[.]pyx$$/\1.c \1.cpp/'); \
	     rm -rf sage/ext/interpreters) \
	    && (cd "$(SAGE_ROOT)/build/pkgs/sagelib/src/" && rm -rf build); \
	fi

sage_docbuild-clean:
	(cd "$(SAGE_ROOT)/build/pkgs/sage_docbuild/src" && rm -rf build)

sage_setup-clean:
	(cd "$(SAGE_ROOT)/build/pkgs/sage_setup/src" && rm -rf build)

build-clean: clean doc-clean sagelib-clean sage_docbuild-clean

doc-clean:
	if [ -f "$(SAGE_SRC)"/bin/sage-env-config ]; then \
	    . "$(SAGE_SRC)"/bin/sage-env-config; \
	    if [ -n "$$SAGE_LOCAL" ]; then \
	        rm -rf "$$SAGE_LOCAL/share/doc/sage/inventory"; \
	        rm -rf "$$SAGE_LOCAL/share/doc/sage/doctrees"; \
	    fi; \
	fi; \
	cd "$(SAGE_SRC)/doc" && $(MAKE) clean

# Deleting src/lib is to get rid of src/lib/pkgconfig
# that was forgotten to clean in #29082.
misc-clean:
	@echo "Deleting build artifacts generated by autoconf, automake, make ..."
	rm -rf logs
	rm -rf dist
	rm -rf tmp
	rm -f aclocal.m4 config.log confcache
	rm -rf autom4te.cache
	rm -f build/make/Makefile build/make/Makefile-auto
	rm -rf src/lib

bdist-clean: clean
	$(MAKE) misc-clean
	@echo "Deleting build artifacts generated by configure ..."
	rm -f config.status

distclean: build-clean
	$(MAKE) misc-clean
	@echo "Deleting all remaining output from build system ..."
	rm -rf local
	rm -f src/bin/sage-env-config
	rm -f prefix venv

# Delete all auto-generated files which are distributed as part of the
# source tarball
bootstrap-clean:
	rm -rf config/install-sh config/compile config/config.guess config/config.sub config/missing configure build/make/Makefile-auto.in
	rm -f src/doc/en/installation/*.txt
	rm -rf src/doc/en/reference/spkg/*.rst
	rm -f environment.yml
	rm -f src/environment.yml
	rm -f src/environment-dev.yml
	rm -f environment-optional.yml
	rm -f src/environment-optional.yml
	rm -f src/Pipfile
	rm -f src/pyproject.toml
	rm -f src/requirements.txt
	rm -f src/setup.cfg

# Remove absolutely everything which isn't part of the git repo
maintainer-clean: distclean bootstrap-clean
	rm -rf upstream

# Remove everything that is not necessary to run Sage and pass all its
# doctests.
micro_release:
	$(MAKE) sagelib-clean
	$(MAKE) misc-clean
	@echo "Stripping binaries ..."
	LC_ALL=C find local/lib local/bin local/var/lib/sage/venv-python* -type f -exec strip '{}' ';' 2>&1 | grep -v "File format not recognized" |  grep -v "File truncated" || true
	@echo "Removing sphinx artifacts..."
	rm -rf local/share/doc/sage/doctrees local/share/doc/sage/inventory
	@echo "Removing documentation. Inspection in IPython still works."
	rm -rf local/share/doc local/share/*/doc local/share/*/examples local/share/singular/html
	@echo "Removing unnecessary files & directories - make will not be functional afterwards anymore"
	@# We keep src/sage for some doctests that it expect it to be there and
	@# also because it does not add any weight with rdfind below.
	@# We need src/bin/ for the scripts that invoke Sage
	@# We need sage, the script to start Sage
	@# We need local/, the dependencies and the built Sage library itself.
	@# We keep VERSION.txt.
	@# We keep COPYING.txt so we ship a license with this distribution.
	find . -name . -o -prune ! -name config.status ! -name src ! -name sage ! -name local ! -name VERSION.txt ! -name COPYING.txt ! -name build -exec rm -rf \{\} \;
	cd src && find . -name . -o -prune ! -name sage ! -name bin -exec rm -rf \{\} \;
	if command -v rdfind > /dev/null; then \
		echo "Hardlinking identical files."; \
		rdfind -makeresultsfile false -makehardlinks true .; \
	else \
		echo "rdfind not installed. Not hardlinking identical files."; \
	fi

# Leaves everything that is needed to make the next "make" fast but removes
# all the cheap build artifacts that can be quickly regenerated.
# Trac #30960: We no longer uninstall sagelib.
fast-rebuild-clean: misc-clean
	rm -rf upstream/
	rm -rf build/pkgs/sagelib/src/build/temp.*
	# The .py files in src/build are restored from src/sage without their
	# mtimes changed.
	-find build/pkgs/sagelib/src/build -name '*.py' -exec rm \{\} \;
	# Remove leftovers from ancient branches
	rm -rf src/build

###############################################################################
# Testing
###############################################################################

TEST_LOG = $(SAGE_ROOT_LOGS)/test.log

TEST_FILES = --all

TEST_FLAGS =

# When the documentation is installed, "optional" also includes all tests marked 'sagemath_doc_html',
# see https://trac.sagemath.org/ticket/25345, https://trac.sagemath.org/ticket/26110, and
# https://trac.sagemath.org/ticket/32759
TEST_OPTIONAL = sage,optional

# Keep track of the top-level *test* Makefile target for logging.
TEST_TARGET = $@

TEST = ./sage -t --logfile=$(TEST_LOG) $(TEST_FLAGS) --optional=$(TEST_OPTIONAL) $(TEST_FILES)

test: all
	@echo '### make $(TEST_TARGET): Running $(TEST)' >> $(TEST_LOG)
	$(TEST)

check:
	@$(MAKE) test

testall:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_OPTIONAL="$(TEST_OPTIONAL),external" test

testlong:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) --long" test

testalllong:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) --long" testall

ptest:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) -p" test

ptestall:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_OPTIONAL="$(TEST_OPTIONAL),external" ptest

ptestlong:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) --long" ptest

ptestalllong:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) --long" ptestall

testoptional:
	@echo "'make $@' is deprecated; use 'make test'"
	@$(MAKE) test

testoptionallong:
	@echo "'make $@' is deprecated; use 'make testlong'"
	@$(MAKE) testlong

ptestoptional:
	@echo "'make $@' is deprecated; use 'make ptest'"
	@$(MAKE) ptest

ptestoptionallong:
	@echo "'make $@' is deprecated; use 'make ptestlong'"
	@$(MAKE) ptestlong

# *test*-nodoc targets skip the docbuild and skip all tests that depend on the documentation

test-nodoc: TEST_OPTIONAL := $(TEST_OPTIONAL),!sagemath_doc_html,!sagemath_doc_pdf
test-nodoc: build
	@echo '### make $(TEST_TARGET): Running $(TEST)' >> $(TEST_LOG)
	$(TEST)

check-nodoc:
	@$(MAKE) test-nodoc

testall-nodoc:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_OPTIONAL="$(TEST_OPTIONAL),external" test-nodoc

testlong-nodoc:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) --long" test-nodoc

testalllong-nodoc:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) --long" testall-nodoc

ptest-nodoc:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) -p" test-nodoc

ptestall-nodoc:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_OPTIONAL="$(TEST_OPTIONAL),external" ptest-nodoc

ptestlong-nodoc:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) --long" ptest-nodoc

ptestalllong-nodoc:
	@$(MAKE) TEST_TARGET="$(TEST_TARGET)" TEST_FLAGS="$(TEST_FLAGS) --long" ptestall-nodoc

testoptional-nodoc:
	@echo "'make $@' is deprecated; use 'make test-nodoc'"
	@$(MAKE) test-nodoc

testoptionallong-nodoc:
	@echo "'make $@' is deprecated; use 'make testlong-nodoc'"
	@$(MAKE) testlong-nodoc

ptestoptional-nodoc:
	@echo "'make $@' is deprecated; use 'make ptest-nodoc'"
	@$(MAKE) ptest-nodoc

ptestoptionallong-nodoc:
	@echo "'make $@' is deprecated; use 'make ptestlong-nodoc'"
	@$(MAKE) ptestlong-nodoc

###############################################################################

configure: bootstrap src/doc/bootstrap configure.ac src/bin/sage-version.sh m4/*.m4 build/pkgs/*/spkg-configure.m4 build/pkgs/*/type build/pkgs/*/install-requires.txt build/pkgs/*/package-version.txt build/pkgs/*/distros/*.txt
	./bootstrap -d

install: all
	@echo "******************************************************************"
	@echo "The '$@' target is a no-op; 'make' already does 'make install'"
	@echo "You can change the install prefix from its default"
	@echo "(the subdirectory 'local') by using ./configure --prefix=PREFIX"
	@echo "You can also consider using the binary packaging scripts"
	@echo "from https://github.com/sagemath/binary-pkg"
	@echo "******************************************************************"

# Setting SAGE_PKGCONFIG is only so that make does not exit with
# "This Makefile needs to be invoked by build/make/install".
list:
	@$(MAKE) --silent build/make/Makefile >&2
	@$(MAKE) --silent -f build/make/Makefile SAGE_PKGCONFIG=dummy $@

.PHONY: default build dist install micro_release \
	pypi-sdists pypi-wheels wheels \
	misc-clean bdist-clean distclean bootstrap-clean maintainer-clean \
	test check testoptional testall testlong testoptionallong testallong \
	ptest ptestoptional ptestall ptestlong ptestoptionallong ptestallong \
	buildbot-python3 list \
	doc-clean clean sagelib-clean build-clean
