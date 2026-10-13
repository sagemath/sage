# SPDX-License-Identifier: GPL-2.0-or-later
"""Check that splitting a bibliography preserves its metadata and rendered output."""

import argparse
import hashlib
import subprocess
import tempfile
from pathlib import Path

from _html_compare import main_content, read_html, write_report
from pybtex.database.input.bibtex import Parser


class DatabaseParser(Parser):
    """Reject duplicate keys across files instead of silently skipping them."""

    def __init__(self):
        super().__init__()
        self.seen = set()
        self.data.want_entry = lambda key: True

    def process_entry(self, entry_type, key, fields):
        if key.upper() in self.seen:
            raise ValueError(f'case-insensitive duplicate bibliography key: {key}')
        self.seen.add(key.upper())
        super().process_entry(entry_type, key, fields)


def bibliography(paths):
    files = []
    for path in paths:
        files.extend(sorted(path.glob('*.bib')) if path.is_dir() else [path])
    if not files:
        raise ValueError('no bibliography files found')
    return DatabaseParser().parse_files(files)


def pdf_text(path, gs, directory, name):
    output = directory / name
    subprocess.run([gs, '-q', '-dBATCH', '-dNOPAUSE', '-sDEVICE=txtwrite',
                    f'-sOutputFile={output}', str(path)], check=True, capture_output=True)
    return output.read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before-bib', type=Path, nargs='+', required=True,
                        help='old .bib files or directories containing .bib files')
    parser.add_argument('--after-bib', type=Path, nargs='+', required=True,
                        help='new .bib files or directories containing .bib files')
    parser.add_argument('--before-html', type=Path, required=True)
    parser.add_argument('--after-html', type=Path, required=True)
    parser.add_argument('--before-latex', type=Path)
    parser.add_argument('--after-latex', type=Path)
    parser.add_argument('--before-pdf', type=Path)
    parser.add_argument('--after-pdf', type=Path)
    parser.add_argument('--gs', default='gs', help='Ghostscript executable for PDF text extraction')
    parser.add_argument('--output', type=Path, help='write the full JSON report to this file')
    args = parser.parse_args()
    for kind in ('latex', 'pdf'):
        if bool(getattr(args, f'before_{kind}')) != bool(getattr(args, f'after_{kind}')):
            parser.error(f'--before-{kind} and --after-{kind} must be supplied together')
    try:
        old = bibliography(args.before_bib)
        new = bibliography(args.after_bib)
        old_page = read_html(args.before_html)
        new_page = read_html(args.after_html)
        old_citations = [str(node) for node in old_page.select('div.citation')]
        new_citations = [str(node) for node in new_page.select('div.citation')]
        if not old_citations or not new_citations:
            raise ValueError('both HTML files must contain bibliography citations')
        old_html = args.before_html.read_bytes()
        new_html = args.after_html.read_bytes()
        checks = {
            'entry_metadata_identical': old.entries == new.entries,
            'preamble_identical': old.preamble == new.preamble,
            'citation_html_identical': old_citations == new_citations,
            'main_html_identical': str(main_content(old_page, args.before_html))
            == str(main_content(new_page, args.after_html)),
        }
        report = {
            'summary': {},
            'checks': checks,
            'before_entries': len(old.entries),
            'after_entries': len(new.entries),
            'before_citations': len(old_citations),
            'after_citations': len(new_citations),
            'added_keys': sorted(set(new.entries) - set(old.entries)),
            'missing_keys': sorted(set(old.entries) - set(new.entries)),
            'changed_keys': sorted(key for key in old.entries
                                   if key in new.entries and old.entries[key] != new.entries[key]),
            'entire_html_byte_identical': old_html == new_html,
            'html_sha256_before': hashlib.sha256(old_html).hexdigest(),
            'html_sha256_after': hashlib.sha256(new_html).hexdigest(),
        }
        if args.before_latex is not None:
            checks['latex_byte_identical'] = args.before_latex.read_bytes() == args.after_latex.read_bytes()
        if args.before_pdf is not None:
            with tempfile.TemporaryDirectory(prefix='rst2bib-') as temporary_directory:
                directory = Path(temporary_directory)
                checks['pdf_extracted_text_identical'] = (
                    pdf_text(args.before_pdf, args.gs, directory, 'before.txt')
                    == pdf_text(args.after_pdf, args.gs, directory, 'after.txt'))
            report['pdf_byte_identical'] = args.before_pdf.read_bytes() == args.after_pdf.read_bytes()
        report['summary'] = {'entries': len(new.entries), 'citations': len(new_citations),
                             'ok': all(checks.values())}
        write_report(report, args.output)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    return 0 if report['summary']['ok'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
