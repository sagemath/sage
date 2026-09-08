# Generate the "Bibliographic References" ReST document from a
# template.
reference_db = custom_target(
  'reference_db',
  input: [src / 'doc/reference-database.rst', 'index.rst.in'],
  output: 'index.rst',
  command: [
    py,
    src / 'sage_docbuild/generate-references.py',
    '@INPUT@',
    '@OUTPUT@',
  ],
)

doc_src += fs.copyfile('conf.py')
doc_src += reference_db
