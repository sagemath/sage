doc_sources = [
  'coding_basics.rst',
  'coding_in_cython.rst',
  'coding_in_other.rst',
  'coding_in_python.rst',
  'conf.py',
  'doctesting.rst',
  'downstream.rst',
  'git_advanced.rst',
  'git_background.rst',
  'git_basic.rst',
  'git_setup.rst',
  'github.rst',
  'index.rst',
  'packaging.rst',
  'packaging_sage_library.rst',
  'portability_platform_table.rst',
  'portability_testing.rst',
  'review.rst',
  'sage_manuals.rst',
  'tools.rst',
  'walkthrough.rst',
  'workflows.rst',
  'workspace.rst',
]

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

subdir('static')
subdir('sagenb')
