#! /bin/sh
set -e

UV_VERSION="${VERSION:-}"
[ "$UV_VERSION" = "latest" ] && UV_VERSION=""
curl -LsSf https://astral.sh/uv/${UV_VERSION}/install.sh | sh
