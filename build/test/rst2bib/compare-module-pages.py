# SPDX-License-Identifier: GPL-2.0-or-later
"""Compare module-page prose and examples, with explicit migration exclusions."""

import argparse
import re
from difflib import unified_diff
from pathlib import Path

from _html_compare import main_content, normalize, read_html, read_key_map, write_report


def read_page(path, key_map, *, ignore_title=False, ignore_bibliography=False):
    article = main_content(read_html(path), path)
    examples = [node.get_text() for node in article.select('.highlight')]
    citation_count = len(article.select('div.citation'))
    if ignore_bibliography:
        for node in article.select('div.citation'):
            node.decompose()
        for node in article.select('li'):
            if re.fullmatch(r'\[[^]]+\]', node.get_text().strip()):
                node.decompose()
        for node in article.select('p.rubric'):
            if node.get_text() == 'Bibliography':
                node.decompose()
    if ignore_title and article.h1 is not None:
        article.h1.decompose()
    text = normalize(article.get_text(' ', strip=True))
    # Apply each rename once; replacements must not cascade through another key.
    text = re.sub(r'\[([^]]+)\]', lambda match: f'[{key_map.get(match[1], match[1])}]', text)
    return {'text': text, 'examples': examples, 'citations': citation_count}


def page_pairs(before, after, selected):
    if before.is_file() and after.is_file():
        if selected:
            raise ValueError('--page applies to directory inputs only')
        return [(before.name, before, after)], [], []
    if not before.is_dir() or not after.is_dir():
        raise ValueError('before and after must both be files or both be directories')
    if selected:
        paths = set()
        for name in selected:
            path = Path(name)
            if path.is_absolute() or '..' in path.parts:
                raise ValueError(f'--page must name a path inside the input directories: {name}')
            paths.add(path.as_posix())
        old = {name for name in paths if (before / name).is_file()}
        new = {name for name in paths if (after / name).is_file()}
        missing = sorted(paths - new)
        missing_before = sorted(paths - old)
    else:
        old = {path.relative_to(before).as_posix() for path in before.rglob('*.html')}
        new = {path.relative_to(after).as_posix() for path in after.rglob('*.html')}
        missing = sorted(old - new)
        missing_before = sorted(new - old)
    pairs = [(name, before / name, after / name) for name in sorted(old & new)]
    if not pairs and not missing and not missing_before:
        raise ValueError('no HTML pages found to compare')
    return pairs, missing, missing_before


def difference(before, after):
    return '\n'.join(unified_diff(before, after, fromfile='before', tofile='after', lineterm=''))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('before', type=Path, help='old HTML file or directory')
    parser.add_argument('after', type=Path, help='new HTML file or directory')
    parser.add_argument('--page', nargs='+', help='relative page paths to compare within directories')
    parser.add_argument('--key-map', type=Path, help='JSON mapping old citation keys to new keys')
    parser.add_argument('--ignore-title', action='store_true', help='exclude the first h1 heading')
    parser.add_argument('--ignore-bibliography', action='store_true',
                        help='exclude citation blocks, citation-only list items, and bibliography rubrics')
    parser.add_argument('--output', type=Path, help='write the full JSON report to this file')
    args = parser.parse_args()
    try:
        key_map = read_key_map(args.key_map)
        pairs, missing_after, missing_before = page_pairs(args.before, args.after, args.page)
        pages = {}
        for name, before, after in pairs:
            old = read_page(before, key_map, ignore_title=args.ignore_title,
                            ignore_bibliography=args.ignore_bibliography)
            new = read_page(after, {}, ignore_title=args.ignore_title,
                            ignore_bibliography=args.ignore_bibliography)
            prose_equal = old['text'] == new['text']
            examples_equal = old['examples'] == new['examples']
            page = {
                'prose_equal': prose_equal,
                'examples_equal': examples_equal,
                'example_blocks_before': len(old['examples']),
                'example_blocks_after': len(new['examples']),
                'citations_before': old['citations'],
                'citations_after': new['citations'],
            }
            if not prose_equal:
                page['prose_diff'] = difference(old['text'].split(), new['text'].split())
            if not examples_equal:
                page['examples_diff'] = difference('\n\n'.join(old['examples']).splitlines(),
                                                  '\n\n'.join(new['examples']).splitlines())
            pages[name] = page
        changed = [name for name, page in pages.items()
                   if not page['prose_equal'] or not page['examples_equal']]
        report = {
            'summary': {
                'compared': len(pages),
                'changed': len(changed),
                'example_blocks': sum(page['example_blocks_before'] for page in pages.values()),
                'ok': not changed and not missing_before and not missing_after,
            },
            'pages': pages,
            'missing_before': missing_before,
            'missing_after': missing_after,
        }
        write_report(report, args.output)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0 if report['summary']['ok'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
