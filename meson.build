project('sage', 'c')
cc = meson.get_compiler('c')

maxima = find_program('maxima', required: true)
ecl_config = find_program('ecl-config', required: true)
openmp = dependency('openmp', required : false)
#ntl = dependency('ntl', required: true) doesn't work, so ask the compiler directly:
ntl = cc.find_library('ntl', required: true)

subdir('pkgs/sage-conf_meson')
