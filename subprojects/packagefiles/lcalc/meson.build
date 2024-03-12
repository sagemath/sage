project('lcalc', 'cpp', version : '2.0.5')

cc = meson.get_compiler('cpp')

# Check for all of the standard headers that we use
required_headers = ['getopt.h', 'limits.h', 'math.h', 'stdlib.h', 'string.h']
foreach header : required_headers
    cc.has_header(header, required : true)
endforeach

# Check for GNU gengetopt
gengetopt = find_program('gengetopt', required : true)

# Check for the mkstemp function
# We invoke mkstemp() through m4 in tests/lib/compare_fp_lists.sh.
# We already check for the necessary stdlib.h above.
cc.has_function('mkstemp', required : true)

# TODO: add options for openmp, pari and precision as in https://gitlab.com/sagemath/lcalc/-/blob/master/configure.ac
# add_project_arguments('openmp args', language : 'cpp')

#subdir('doc')
subdir('src')
#subdir('tests')
