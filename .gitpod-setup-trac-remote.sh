#!/usr/bin/env bash

# Exit on error
set -e

git remote remove trac 2> /dev/null || true  # might still exists from a previous run/prebuild

# Setup trac as remote
git remote add trac https://github.com/sagemath/sagetrac-mirror.git -t master -t develop
git remote set-url --push trac no-pushing--this-is-a-read-only-archive
