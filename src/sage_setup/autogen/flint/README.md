Autogeneration of flint header files for SageMath
=================================================

The scripts in this folder are responsible for the automatic generation of the
Cython header files in `$SAGE_ROOT/src/sage/libs/flint/`. To make the
autogeneration

1. Obtain a clone of the flint repo, eg `git clone
   https://github.com/flintlib/flint2`

2. Checkout to the appropriate commit, eg `git checkout v2.9.0`

3. Possibly adjust the content of `types.pxd.template` (which will be used to
   generate types.pxd)

4. Set the environment variable `FLINT_GIT_DIR`

5. Run the `flint_autogen.py` script eg `python
   $SAGE_ROOT/src/sage_setup/autogen/flint_autogen.py`. The script writes down
   the headers in the sage source tree `$SAGE_ROOT/src/sage/libs/flint/`.


Additional notes
----------------

- macros in flint documentation are not converted into cython declarations
  (because they lack a signature). The cython signature of flint macros must be
  manually written down in the files contained
  `SAGE_SRC/sage/src/sage_setup/autogen/flint/macros`
  See https://github.com/flintlib/flint/issues/1529.
