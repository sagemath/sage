#!/bin/bash
set -e
make build
exec "$@"
