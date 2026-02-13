if not get_option('build-docs')
  subdir_done()
endif

warning(
  'Documentation building enabled, generating targets may be slow. To disable this, pass -Dbuild-docs=false.',
)

sphinx_check = py_module.find_installation(required: false, modules: ['sphinx'])
if not sphinx_check.found()
  warning(
    'Python package "sphinx" is not installed. Documentation build not enabled.',
  )
  subdir_done()
endif

doc_src = []
subdir('el')
subdir('ca')
subdir('zh')
subdir('it')
subdir('common')
subdir('hu')
subdir('de')
subdir('en')
subdir('pt')
subdir('ru')
subdir('tr')
subdir('es')
subdir('fr')
subdir('ja')

# TODO: Migrate this completely to meson
doc_bootstrap = custom_target(
  'bootstrap',
  output: ['autogen'],
  command: [files('bootstrap'), meson.current_build_dir()],
  env: {'SAGE_ROOT': root},
  # doc_src is not really a dependency of the bootstrap, but we want to make sure
  # that all the source files are present before running the actual doc build
  # so let's collect all source-related targets here.
  depends: doc_src,
)

references = run_command(
  py,
  [
    src / 'build-docs.py',
    '--no-prune-empty-dirs',
    '--all-documents',
    'reference',
    '--source',
    meson.current_source_dir(),
  ],
  check: true,
).stdout().strip()

reference_inventory = []
reference_html = []
reference_pdf = []
bibliography = []
foreach type : ['inventory', 'html', 'pdf']
  foreach ref : references.splitlines()
    if '/' in ref
      short_ref = ref.split('/')[1]
    else
      short_ref = ref
    endif
    deps = []
    deps += doc_bootstrap
    if short_ref == 'repl'
      deps += doc_src_repl
    endif
    if type == 'html' or type == 'pdf'
      deps += reference_inventory
    endif
    if short_ref != 'references'
      deps += bibliography
    endif
    if short_ref == 'reference_top'
      deps += reference_inventory
      if type == 'html'
        deps += reference_html
      elif type == 'pdf'
        deps += reference_pdf
      endif
    endif
    target = custom_target(
      'doc-' + type + '-reference-' + short_ref,
      output: [type + short_ref],
      command: [
        py,
        src / 'build-docs.py',
        '--no-prune-empty-dirs',
        '--no-pdf-links',
        ref,
        type,
        '-o',
        '@OUTDIR@',
        '--source',
        meson.current_build_dir(),
      ],
      depends: deps,
    )
    if short_ref == 'references'
      bibliography += target
    endif
    if type == 'inventory'
      reference_inventory += target
    elif type == 'html'
      reference_html += target
    elif type == 'pdf'
      reference_pdf += target
    endif
  endforeach
endforeach

other_documents = run_command(
  py,
  [
    src / 'build-docs.py',
    '--no-prune-empty-dirs',
    '--all-documents',
    'all',
    '--source',
    meson.current_source_dir(),
  ],
  check: true,
).stdout().strip()
other_documents_html = []
other_documents_pdf = []
foreach type : ['html', 'pdf']
  foreach doc : other_documents.splitlines()
    short_doc = doc.replace('/', '-')
    target = custom_target(
      'doc-' + type + '-other-' + short_doc,
      output: [type + short_doc],
      command: [
        py,
        src / 'build-docs.py',
        '--no-prune-empty-dirs',
        '--no-pdf-links',
        doc,
        type,
        '-o',
        '@OUTDIR@',
        '--source',
        meson.current_build_dir(),
      ],
      depends: reference_inventory,
    )
    if type == 'html'
      other_documents_html += target
    elif type == 'pdf'
      other_documents_pdf += target
    endif
  endforeach
endforeach

# Custom target for building the complete documentation
alias_target('doc-html', [reference_html, other_documents_html])
alias_target('doc-pdf', [reference_pdf, other_documents_pdf])
