.PHONY: all sage csage clean

# By default, build csage and sage (in this order)
all:
	$(MAKE) csage
	$(MAKE) sage

sage: sage/libs/pari/auto_gen.pxi sage/ext/interpreters/__init__.py
	python setup.py install

csage:
	cd c_lib && scons -Q install

clean:
	@echo "Deleting Sage library build artifacts..."
	rm -f c_lib/.sconsign.dblite
	find c_lib -name '*.os' | xargs rm -f
	find c_lib -name '*.so' | xargs rm -f
	find c_lib -name '*.dylib' | xargs rm -f
	rm -rf build
	find . -name '*.pyc' | xargs rm -f
	rm -f sage/libs/pari/auto_*
	rm -rf sage/ext/interpreters


# Auto-generated files
sage/libs/pari/auto_gen.pxi: $(SAGE_LOCAL)/share/pari/pari.desc \
        sage/libs/pari/decl.pxi sage_setup/autogen/pari/*.py
	python -c "from sage_setup.autogen.pari import rebuild; rebuild()"

sage/ext/interpreters/__init__.py: sage_setup/autogen/interpreters.py
	python -c "from sage_setup.autogen.interpreters import rebuild; rebuild('sage/ext/interpreters')"
