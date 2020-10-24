#!/bin/sh

# This script gets called from CI to run minimal tests on the sagemath-jupyter
# image.

# Usage: ./test-jupyter.sh IMAGE-NAME [HOST]

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

set -ex
docker run --name sage-jupyter -d "$1" sage-jupyter
echo "Checking that the Jupyter notebook is running…"
sleep 10 # giving the server some time to start
docker logs sage-jupyter
docker run --link sage-jupyter alpine sh -c "apk --update add wget; wget --retry-connrefused --tries=10 --wait=3 http://sage-jupyter:8888"
