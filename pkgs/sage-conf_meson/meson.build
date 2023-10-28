conf_data = configuration_data()

conf_data.set('PACKAGE_VERSION', '1.2.3')

conf_data.set('SAGE_MAXIMA', maxima.path())
# Conda's ecl does not have any problems with Maxima, so nothing needs to be set here:
conf_data.set('SAGE_MAXIMA_FAS', '')

# Kenzo cannot yet be provided by the system, so we always use the SAGE_LOCAL path for now.
conf_data.set('SAGE_KENZO_FAS', '\'${prefix}\'/lib/ecl/kenzo.fas')

conf_data.set('SAGE_ARB_LIBRARY', 'arb')

# It can be found, so we don't have to set anything here:
conf_data.set('NTL_INCDIR', '')
conf_data.set('NTL_LIBDIR', '')

conf_data.set('SAGE_ECL_CONFIG', ecl_config.path())

conf_data.set('SAGE_ARCHFLAGS', 'unset')

# not needed when using conda, as we then don't build any pc files
conf_data.set('SAGE_PKG_CONFIG_PATH', '')

if openmp.found()
  conf_data.set('OPENMP_CFLAGS', '-fopenmp')
  conf_data.set('OPENMP_CXXFLAGS', '-fopenmp')
endif

configure_file(input : '../sage-conf_conda/_sage_conf/_conf.py.in',
               output : '_conf.py',
               configuration : conf_data)
