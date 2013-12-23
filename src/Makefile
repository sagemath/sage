clean:
	@echo "Deleting Sage library build artifacts..."
	rm -f c_lib/.sconsign.dblite
	find c_lib -name '*.os' | xargs rm -f
	find c_lib -name '*.so' | xargs rm -f
	rm -rf build
	find . -name '*.pyc' | xargs rm -f
	find . -name '*.pyx' | sed 's/pyx$$/c/' | xargs rm -f
	find . -name '*.pyx' | sed 's/pyx$$/cpp/' | xargs rm -f
	rm -rf sage/modular/arithgroup/farey_symbol.h
	rm -rf sage/rings/real_mpfi.h
	rm -rf sage/symbolic/pynac.h
