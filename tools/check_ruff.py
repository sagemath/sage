#!/usr/bin/env python3
# See README.md for more details

import argparse
import json
import re
import subprocess
import tomllib

# These are not currently relevant for Sage
EXCLUDE_PREFIXES = ['AIR', 'ASYNC', 'DJ', 'FAST', 'PD', 'T20']


def get_ruff_config() -> tuple[list[str], list[str]]:
    # First we read from the ruff config
    with open('pyproject.toml', 'rb') as f:
        toml_config = tomllib.load(f)

    ignored_rules = toml_config['tool']['ruff']['lint']['ignore']
    selected_rules = toml_config['tool']['ruff']['lint']['select']
    for ex in EXCLUDE_PREFIXES:
        if ex not in ignored_rules:
            ignored_rules.append(ex)

    # We also check against the ruff config, because we can't easily check if a rule code belongs to a rule prefix.
    # The ruff config outputs all enabled rules with their full codes.
    ruff_config = subprocess.run(
        ['ruff', 'check', '--show-settings', '--silent'],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout

    enabled_rules = re.compile(r'linter\.rules\.enabled = \[.*?\]', re.DOTALL).search(
        ruff_config
    )
    linter_rules = enabled_rules.group().splitlines()[1:-1]

    for linter_rule in linter_rules:
        name, code = linter_rule.split()
        code = code[1:-2]
        if code not in selected_rules:
            selected_rules.append(code)

    return ignored_rules, selected_rules


def check_prefixes(args):
    ignored_rules, selected_rules = get_ruff_config()
    linters = subprocess.run(
        ['ruff', 'linter'], check=True, text=True, stdout=subprocess.PIPE
    ).stdout
    prefixes = []
    for linter in linters.splitlines():
        codes = linter.split()[0]
        prefixes.extend(codes.split('/'))

    for ep in EXCLUDE_PREFIXES:
        prefixes.remove(ep)

    prefixes_to_check = [
        prefix for prefix in prefixes if prefix not in selected_rules + ignored_rules
    ]

    passed_prefixes = []
    failed_prefixes = []
    for i, prefix in enumerate(prefixes_to_check):
        print(f'Checking {prefix} ({i + 1} / {len(prefixes_to_check)})')
        prefix_check = subprocess.run(
            ['ruff', 'check', '--silent', '--preview', '--select', prefix], check=False
        )
        prefix_passes = prefix_check.returncode == 0
        if prefix_passes:
            passed_prefixes.append(prefix)
        else:
            failed_prefixes.append(prefix)

    print(f'The following {len(failed_prefixes)} rule prefixes have failures:')
    print(failed_prefixes)
    print()
    print(
        f'The following {len(passed_prefixes)} rule prefixes pass on the entire codebase:'
    )
    print(passed_prefixes)


def check_rules(args):
    allow_preview = args.allow_preview
    only_fixable = args.only_fixable

    ignored_rules, selected_rules = get_ruff_config()
    rules = subprocess.run(
        ['ruff', 'rule', '--all', '--output-format', 'json'],
        check=True,
        stdout=subprocess.PIPE,
    )
    rule_data = json.loads(rules.stdout)
    del rules

    codes_to_check = [
        rule['code']
        for rule in rule_data
        if rule['code'] not in (ignored_rules + selected_rules)
        and (not rule['preview'] or allow_preview)
        and (rule['fix_availability'] == 'Always' or not only_fixable)
        and not any(rule['code'].startswith(ep) for ep in EXCLUDE_PREFIXES)
        and 'Removed' not in rule['status']
    ]

    passed_rules = []
    failed_rules = []
    for i, code in enumerate(codes_to_check):
        print(f'Checking {code} ({i + 1} / {len(codes_to_check)})')
        rule_check = subprocess.run(
            ['ruff', 'check', '--silent', '--select', code], check=False
        )
        if rule_check.returncode == 0:
            passed_rules.append(code)
        elif rule_check.returncode == 1:
            failed_rules.append(code)
        else:
            print('Error with rule code', code)
            for rule in rule_data:
                if rule['code'] == code:
                    print(rule)

    print(f'The following {len(failed_rules)} rules have failures:')
    print(failed_rules)
    print()
    print(f'The following {len(passed_rules)} rules pass on the entire codebase:')
    print(passed_rules)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(required=True)

    prefix_parser = subparsers.add_parser(
        'prefix',
        description='Check all ruff rule prefixes, with preview enabled. '
        'Rule prefixes that are already enabled or disabled in pyproject.toml '
        'are skipped, as well as some prefixes that are not relevant to Sage.',
    )
    prefix_parser.set_defaults(func=check_prefixes)

    rule_parser = subparsers.add_parser(
        'rule',
        description='Check all ruff rules. '
        'Rules that are already enabled or disabled in pyproject.toml '
        'are skipped, as well as some rules that are not relevant to Sage.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    rule_parser.add_argument(
        '--only-fixable',
        action='store_true',
        help='Only check rules that are always fixable.',
    )
    rule_parser.add_argument(
        '--allow-preview',
        action='store_true',
        help='Check rules that are in preview mode.',
    )
    rule_parser.set_defaults(func=check_rules)

    args = parser.parse_args()
    args.func(args)
