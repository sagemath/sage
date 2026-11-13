# SPDX-License-Identifier: GPL-2.0-or-later
"""Compare citation text, anchors, links, mathematics, and formatting in HTML."""

import argparse
from collections import Counter
from pathlib import Path

from _html_compare import citations, read_json, read_key_map, write_report


def compare(before, after, key_map, expected, *, ignore_ids=False, allow_added=False):
    targets = {key: key_map.get(key, key) for key in before}
    missing = {key: target for key, target in targets.items() if target not in after}
    added = sorted(set(after) - set(targets.values()))
    differences = {}
    fields = Counter()
    for key, old in before.items():
        if key in missing:
            continue
        new = after[targets[key]]
        delta = {field: {'before': value, 'after': new[field]}
                 for field, value in old.items()
                 if (field != 'id' or not ignore_ids) and value != new[field]}
        if delta:
            differences[key] = delta
            fields.update(delta.keys())
    expected_match = differences == expected
    ok = not missing and (allow_added or not added) and expected_match
    return {
        'summary': {
            'before': len(before),
            'after': len(after),
            'compared': len(before) - len(missing),
            'identical': len(before) - len(missing) - len(differences),
            'different': len(differences),
            'missing': len(missing),
            'added': len(added),
            'fields': dict(fields),
            'expected_differences_match': expected_match,
            'ok': ok,
        },
        'differences': differences,
        'expected_differences': expected,
        'missing': missing,
        'added': added,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('before', type=Path, help='HTML file before migration')
    parser.add_argument('after', type=Path, help='HTML file after migration')
    parser.add_argument('--key-map', type=Path, help='JSON mapping old citation keys to new keys')
    parser.add_argument('--expected-differences', type=Path,
                        help='JSON with the exact accepted field differences (default: none)')
    parser.add_argument('--ignore-ids', action='store_true',
                        help='omit citation-anchor IDs from the comparison')
    parser.add_argument('--allow-added', action='store_true',
                        help='allow extra citations in the after file; missing citations still fail')
    parser.add_argument('--output', type=Path, help='write the full JSON report to this file')
    args = parser.parse_args()
    try:
        report = compare(citations(args.before), citations(args.after), read_key_map(args.key_map),
                         read_json(args.expected_differences),
                         ignore_ids=args.ignore_ids, allow_added=args.allow_added)
        write_report(report, args.output)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0 if report['summary']['ok'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
