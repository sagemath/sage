#!/usr/bin/env python3

"""
Script for generating checksums.ini files from PyPI.
"""

import argparse

import requests

argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('package_name')
    parser.add_argument('version')

    args = parser.parse_args()
    name = args.package_name
    version = args.version

    r = requests.get(f'https://pypi.org/pypi/{name}/{version}/json')
    response = r.json()
    urls = response['urls']

    entries = []
    for url in urls:
        tarball = url['filename']
        sha256 = url['digests']['sha256']
        upstream_url = url['url']

        entries.append(f'tarball={tarball}\nsha256={sha256}\nupstream_url={upstream_url}')
    print('\n\n'.join(entries))
