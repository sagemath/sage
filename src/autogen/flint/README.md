Autogeneration of flint header files for SageMath
=================================================

The scripts in this folder are responsible for the automatic generation of the
Cython header files in `$SAGE_ROOT/src/sage/libs/flint/`. To make the
autogeneration

1. Obtain a clone of the flint repo, eg `git clone
   https://github.com/flintlib/flint`

2. Checkout to the appropriate commit, eg `git checkout v3.0.1`. The correct
   version can be found in `$SAGE_ROOT/build/pkgs/flint/package-version.txt`

3. Possibly edit the docs at `$FLINT_ROOT/doc/source/*.rst` to match the
   exposed API.txt`. For example, the docs of release `v3.0.1` was incorrect,
   so the commit `3e2c3a3e091106a25ca9c6fba28e02f2cbcd654a` was used instead

4. Possibly adjust the content of `types.pxd.template` (which will be used to
   generate types.pxd)

5. Set the environment variable `FLINT_GIT_DIR`

6. Run the `flint_autogen.py` script e.g. `python
   $SAGE_ROOT/src/autogen/flint_autogen.py`. The script writes down
   the headers in the sage source tree `$SAGE_ROOT/src/sage/libs/flint/`


Additional notes
----------------

- macros in flint documentation are not converted into cython declarations
  (because they lack a signature). The cython signature of flint macros must be
  manually written down in the files contained
  `SAGE_SRC/sage/src/autogen/flint/macros`
  See https://github.com/flintlib/flint/issues/1529.
