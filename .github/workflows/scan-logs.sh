#!/usr/bin/env bash
# The markup in the output is a GitHub Actions logging command
# https://help.github.com/en/actions/automating-your-workflow-with-github-actions/development-tools-for-github-actions
LOGS=${1-logs}
for a in $(find "$LOGS" -type f -name "*.log"); do
    if tail -100 "$a" 2>/dev/null | grep "^[A-Za-z]*Error" >/dev/null; then
        echo :":"error file=$a:":" ==== ERROR IN LOG FILE $a ====
        echo "::group::$a"
        cat "$a"
        echo "::endgroup::"
    elif tail -20 "$a" 2>/dev/null | grep -E "^(Warning: Error testing|^sage .*doctest.*failed)" >/dev/null; then
        echo :":"warning file=$a:":" ==== TESTSUITE FAILURE IN LOG FILE $a ====
        echo "::group::$a"
        cat "$a"
        echo "::endgroup::"
    fi
done
