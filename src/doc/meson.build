# Clean target
# clean = custom_target(
#     'clean',
#     output: [],
#     command: [
#         'sh', '-c',
#         '''
#         rm -rf en/reference/*/sage
#         rm -rf en/reference/documentation/sage_docbuild
#         rm -rf en/reference/sage
#         rm -f common/*.pyc
#         '''
#     ]
# )

# Generate sources at build time
# doc_src = custom_target(
#     'doc-src',
#     output: ['en/reference/repl/options.txt'],
#     command: [
#         'sh', '-c',
#         '''
#         mkdir -p en/reference/repl
#         sage -advanced > en/reference/repl/options.txt
#         '''
#     ]
# )



references = run_command(
  py,
  [
    src / 'build-docs.py',
    '--all-documents',
    'reference',
    '--source', meson.current_source_dir()
  ],
  check: true,
).stdout().strip()

reference_inventory = []
reference_html = []
reference_pdf = []
bibliography = []
foreach type : ['inventory', 'html', 'pdf']
    foreach ref : references.splitlines()
        short_ref = ref.split('/')[1]
        deps = []
        if type == 'html' or type == 'pdf'
            deps += reference_inventory
        endif
        if short_ref != 'references'
            deps += bibliography
        endif
        target = custom_target(
            'doc-' + type + '-reference-' + short_ref,
            output: [type + short_ref],
            command: [
                py, 
                src / 'build-docs.py',
                '--no-pdf-links',
                ref,
                type,
                '-o', '@OUTDIR@',
                '--source', meson.current_source_dir()
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
    '--all-documents',
    'all',
    '--source', meson.current_source_dir()
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
                '--no-pdf-links',
                doc,
                type,
                '-o', '@OUTDIR@',
                '--source', meson.current_source_dir()
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
