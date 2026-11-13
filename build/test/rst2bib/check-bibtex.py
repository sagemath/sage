"""Check field structure and exercise ordinary BibTeX formatting.

This is an offline regression check, not a substitute for checking a reference
against its publication. It deliberately removes Sage's annotation field before
formatting, so preserved prose cannot conceal missing structured metadata.
"""

import argparse
import copy
import csv
import json
import re
from collections import Counter
from pathlib import Path

from pybtex.database import parse_file
from pybtex.plugin import find_plugin


def check(paths, exceptions):
    entries = {}
    errors = []
    for path in paths:
        for key, entry in parse_file(path, 'bibtex').entries.items():
            if key.casefold() in {k.casefold() for k in entries}:
                errors.append(f'{key}: duplicate citation key')
            entries[key] = entry

    authorless = set(exceptions['authorless_web_references'])
    incomplete = exceptions['incomplete_publication_fields']
    formatter = find_plugin('pybtex.style.formatting', 'unsrt')()
    formatted = 0
    observed_incomplete = {}
    dois = {}
    for key, entry in entries.items():
        doi = entry.fields.get('doi', '').strip().casefold()
        if doi:
            if doi in dois:
                errors.append(f'{key}: DOI also used by {dois[doi]}: {doi}')
            dois[doi] = key
        title = entry.fields.get('title', '')
        if not title:
            errors.append(f'{key}: missing title')
        people = entry.persons.get('author') or entry.persons.get('editor')
        if not people and key not in authorless:
            errors.append(f'{key}: missing author/editor without explanation')
        if key in authorless and (people or not entry.fields.get('url')):
            errors.append(f'{key}: stale authorless-web-reference exception')
        for role, names in entry.persons.items():
            for person in names:
                name = str(person)
                if name.startswith('{') and name.endswith('}'):
                    continue  # Explicit corporate author.
                if re.search(r'https?://|\\(?:emph|url|href)|\b(?:Journal|Proceedings|Springer)\b', name):
                    errors.append(f'{key}: publication prose in {role}: {name}')
                if any(re.fullmatch(r'[A-Z]\.', n) for n in person.last_names):
                    errors.append(f'{key}: initial parsed as a family name: {name}')
        for field, value in entry.fields.items():
            if field == 'annotation':
                continue
            if re.search(r'<(?:/?(?:i|b|sub|sup)|mml:)[ >]|&amp;', value):
                errors.append(f'{key}: HTML markup in {field}')
            if value.count('{') != value.count('}'):
                errors.append(f'{key}: unbalanced braces in {field}')
        if re.search(r'\$[^$]*\\[_^]', title):
            errors.append(f'{key}: RST escape inside title mathematics')
        if 'note' in entry.fields and entry.fields['note'] == entry.fields.get('annotation'):
            errors.append(f'{key}: full annotation duplicated as note')

        # Missing publication facts must remain explicit exceptions, not guessed
        # dates or institutions merely to satisfy a formatting style.
        missing = [field for field in incomplete.get(key, {}).get('fields', [])
                   if not entry.fields.get(field)]
        if key in incomplete:
            expected = incomplete[key]['fields']
            if missing != expected:
                errors.append(f'{key}: stale incomplete-fields exception')
            observed_incomplete[key] = missing
        plain = copy.deepcopy(entry)
        plain.fields.pop('annotation', None)
        try:
            formatter.format_entry(key, plain).text.render_as('text')
            formatted += 1
        except Exception as error:
            message = str(error)
            allowed = {f'missing {field} in {key}' for field in missing}
            if message not in allowed:
                errors.append(f'{key}: ordinary formatting failed: {message}')

    for key in (authorless | incomplete.keys()) - entries.keys():
        errors.append(f'{key}: exception refers to absent entry')
    return entries, {
        'entries': len(entries),
        'types': dict(sorted(Counter(e.type for e in entries.values()).items())),
        'with_doi': sum(bool(e.fields.get('doi')) for e in entries.values()),
        'ordinary_formatting_passed': formatted,
        'incomplete_publication_fields': observed_incomplete,
        'errors': errors,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bibfiles', nargs='+', type=Path)
    parser.add_argument('--exceptions', type=Path,
                        default=Path(__file__).with_name('metadata-exceptions.json'))
    parser.add_argument('--output', type=Path)
    parser.add_argument('--tsv', type=Path,
                        help='export author/title/publication fields for manual review')
    args = parser.parse_args()
    entries, report = check(args.bibfiles, json.loads(args.exceptions.read_text()))
    if args.tsv:
        with args.tsv.open('w', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t')
            fields = ['title', 'journal', 'booktitle', 'year', 'doi']
            writer.writerow(['key', 'type', 'author', 'editor', *fields])
            for key, entry in sorted(entries.items()):
                writer.writerow([key, entry.type,
                                 ' and '.join(map(str, entry.persons.get('author', []))),
                                 ' and '.join(map(str, entry.persons.get('editor', []))),
                                 *(entry.fields.get(f, '') for f in fields)])
    result = json.dumps(report, indent=2, ensure_ascii=False) + '\n'
    if args.output:
        args.output.write_text(result)
    print(result, end='')
    return bool(report['errors'])


if __name__ == '__main__':
    raise SystemExit(main())
