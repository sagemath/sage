#!/bin/sh

# This script protects all environment variables that start with "SECRET_".
# It puts them in a temporary file. The name of the variable contains the path
# of that file.  This filename can then safely be used in `cat` even if `set
# -x` has been turned on. Also you can run "export" to understand the
# environment without danger.
# Be careful, however, not to use this like the following:
# docker login $DOCKER_USER $(cat $SECRET_DOCKER_PASS)
# as this would expose the password if `set -x` has been turned on.

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

set -eo pipefail
set +x

function encrypt {
    RET=`mktemp`
    eval " echo \$$1" > "$RET"
    echo $RET
}

for name in `awk 'END { for (name in ENVIRON) { print name; } }' < /dev/null`; do
case "$name" in
  SECRET_*)
    export $name="$(encrypt $name)"
    echo "Protected $name"
    ;;
esac
done

unset encrypt
