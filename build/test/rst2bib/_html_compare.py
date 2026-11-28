# SPDX-License-Identifier: GPL-2.0-or-later
"""Shared HTML readers for the RST-to-BibTeX comparison commands."""

import copy
import json
import re

from bs4 import BeautifulSoup


def normalize(text):
    return re.sub(r'\s+', ' ', text).strip()


def read_html(path):
    return BeautifulSoup(path.read_text(encoding='utf-8'), 'html.parser')


def main_content(page, path):
    article = page.select_one('article[role="main"]')
    if article is None:
        article = page.select_one('[role="main"]')
    if article is None:
        raise ValueError(f'no main content found in {path}')
    return article


def citations(path):
    result = {}
    for original in read_html(path).select('div.citation'):
        label = original.select_one('.label')
        if label is None or 'id' not in original.attrs:
            raise ValueError(f'citation without a label or anchor in {path}')
        key = label.get_text().strip().strip('[]')
        if key in result:
            raise ValueError(f'duplicate citation {key!r} in {path}')
        node = copy.deepcopy(original)
        for redundant in node.select('.label, .backrefs'):
            redundant.decompose()
        result[key] = {
            'id': original['id'],
            'text': normalize(node.get_text(' ', strip=True)),
            'links': [[a.get_text(), a['href']] for a in node.select('a[href]')],
            'math': [n.get_text() for n in node.select('.math')],
            'emphasis': [normalize(n.get_text(' ', strip=True)) for n in node.select('em')],
            'strong': [normalize(n.get_text(' ', strip=True)) for n in node.select('strong')],
            'literal': [normalize(n.get_text()) for n in node.select('code')],
        }
    if not result:
        raise ValueError(f'no citations found in {path}')
    return result


def read_json(path):
    if path is None:
        return {}
    result = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(result, dict):
        raise ValueError(f'expected a JSON object in {path}')
    return result


def read_key_map(path):
    result = read_json(path)
    if not all(isinstance(key, str) and isinstance(value, str)
               for key, value in result.items()):
        raise ValueError(f'expected citation keys mapped to citation keys in {path}')
    return result


def write_report(report, path):
    text = json.dumps(report, ensure_ascii=False, indent=2) + '\n'
    if path is None:
        print(text, end='')
    else:
        path.write_text(text, encoding='utf-8')
        print(json.dumps(report['summary'], ensure_ascii=False, indent=2))
