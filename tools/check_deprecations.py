#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pygithub",
#     "tqdm",
# ]
# ///
# pyright: strict

import argparse
import os
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path

import tqdm
from github.MainClass import Github

# Regex pattern to find deprecation instances
DEPRECATION_PATTERN = re.compile(r'deprecation\((\d+),')


def get_pr_closed_date(github_token: str, pr_number: int) -> datetime:
    g = Github(github_token)
    repo = g.get_repo("sagemath/sage")
    issue = repo.get_issue(number=pr_number)
    return issue.closed_at


def search_deprecations(path: str) -> set[tuple[str, int]]:
    deprecations: set[tuple[str, int]] = set()
    for filepath in Path(path).rglob('*.py*'):
        try:
            with filepath.open('r') as f:
                content = f.read()
                matches = DEPRECATION_PATTERN.findall(content)
                for match in matches:
                    deprecations.add((str(filepath), int(match)))
        except (PermissionError, UnicodeDecodeError):
            pass
    print(f"Found {len(deprecations)} deprecations.")
    return deprecations


def main():
    # Get source directory from command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sourcedir", help="Source directory", nargs="?", default=".", type=Path
    )
    parser.add_argument(
        "--token", help="GitHub API token", default=os.getenv('GITHUB_TOKEN'), type=str
    )
    options = parser.parse_args()

    deprecations = search_deprecations(options.sourcedir)

    one_year_ago = datetime.now(timezone.utc) - timedelta(days=365)
    old_deprecations: set[tuple[str, int, datetime]] = set()
    for filepath, pr_number in tqdm.tqdm(deprecations):
        closed_date = get_pr_closed_date(options.token, pr_number)
        if closed_date and closed_date < one_year_ago:
            old_deprecations.add((filepath, pr_number, closed_date))

    if old_deprecations:
        print("Deprecations over one year ago:")
        for filepath, pr_number, closed_date in old_deprecations:
            print(
                f"File: {filepath}, PR: https://github.com/sagemath/sage/pull/{pr_number}, Closed Date: {closed_date:%Y-%m-%d}"
            )
    else:
        print("No deprecations over one year ago found.")


if __name__ == '__main__':
    main()
