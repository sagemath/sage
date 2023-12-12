Autogeneration of flint header files for SageMath
=================================================

The scripts in this folder are responsible for the automatic generation of the Cython header files
in `sage/libs/flint/`. To make the autogeneration

1. Obtain a clone of the flint repo, eg `git clone https://github.com/flintlib/flint2`

2. Checkout to the appropriate commit, eg `git checkout v2.9.0`

3. Possibly adjust the content of `types.pxd.template` (which will be used to generate
   types.pxd)

4. Set the environment variable `FLINT_GIT_DIR`

5. Run the `run.py` script, eg `python autogen.py`

6. Copy all the files generated in `pxd_headers` inside the sage source tree (ie `SAGE_SRC/sage/libs/flint`)
