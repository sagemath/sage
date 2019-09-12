#!/bin/sh

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

set +e -x

docker info
docker run docker sh -c "
  set -x
  uname -a
  df -h
  cat /proc/cpuinfo
  cat /proc/meminfo
  cat /proc/sys/vm/overcommit_memory
  cat /proc/sys/vm/overcommit_ratio"
