# Discover at build time so a newly added bibliography needs no reconfiguration.
doc_bibliographies = custom_target(
  'copy-bibliographies',
  output: 'bibliographies.stamp',
  command: [
    py,
    files('../../../tools/copy-bibliographies.py'),
    meson.current_source_dir(),
    meson.current_build_dir(),
    '@OUTPUT@',
  ],
  build_always_stale: true,
)
doc_src += doc_bibliographies
