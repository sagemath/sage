#!/usr/bin/env bash
for X in standard-python2 minimal; do sed 's/\[standard\]/['$X']/g;s/CI cygwin-standard/CI cygwin-'$X'/g;' ci-cygwin-standard.yml > ci-cygwin-$X.yml; done
