#!/bin/sh

# Source this to set CPUTHREADS (the number of apparent cores) and RAMTHREADS
# (free RAM divided by the maximum amount needed per thread typically)
# From this this script infers reasonable defaults for SAGE_NUM_THREADS and
# MAKEFLAGS.

# We do exactly the same for CPUTHREADS_DOCBUILD, RAMTHREADS_DOCBUILD,
# SAGE_NUM_THREADS_DOCBUILD, MAKEFLAGS_DOCBUILD. As the docbuild needs
# substantially more RAM as of May 2018.

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

if [ -z "$CPUTHREADS" ]; then
    # Determine the number of threads that can run simultaneously on this system
    # (we might not have nproc available.)
    # Note that this value is incorrect for some CI providers (notably CircleCI:
    # https://circleci.com/docs/2.0/configuration-reference/#resource_class) which
    # provision fewer vCPUs than shown in /proc/cpuinfo. So it is probably better
    # to set CPUTHREADS manuall in your CI configuration.
    CPUTHREADS=`docker run docker cat /proc/cpuinfo | grep -E '^processor' | wc -l`
fi
if [ -z "$CPUTHREADS_DOCBUILD" ]; then
    CPUTHREADS_DOCBUILD=$CPUTHREADS
fi

if [ -z "$RAMTHREADS" ]; then
    RAMTHREADS=$(( `docker run docker cat /proc/meminfo | grep MemTotal | awk '{ print $2 }'` / 1048576 ))
    if [ $RAMTHREADS = 0 ];then
        RAMTHREADS=1;
    fi
fi
if [ -z "$RAMTHREADS_DOCBUILD" ]; then
    RAMTHREADS_DOCBUILD=$(( `docker run docker cat /proc/meminfo | grep MemTotal | awk '{ print $2 }'` / 2097152 ))
    if [ $RAMTHREADS_DOCBUILD = 0 ];then
        RAMTHREADS_DOCBUILD=1;
    fi
fi

# On CI machines with their virtual CPUs, it seems to be quite beneficial to
# overcommit on CPU usage. We only need to make sure that we do not exceed RAM
# (as there is no swap.)
if [ $CPUTHREADS -lt $RAMTHREADS ]; then
    export SAGE_NUM_THREADS=$((CPUTHREADS + 1))
else
    export SAGE_NUM_THREADS=$RAMTHREADS
fi
if [ $CPUTHREADS_DOCBUILD -lt $RAMTHREADS_DOCBUILD ]; then
    export SAGE_NUM_THREADS_DOCBUILD=$((CPUTHREADS_DOCBUILD + 1))
else
    export SAGE_NUM_THREADS_DOCBUILD=$RAMTHREADS_DOCBUILD
fi
# Set -j and -l for make (though -l is probably ignored by Sage)
export MAKEFLAGS="-j $SAGE_NUM_THREADS -l $((CPUTHREADS - 1)).8"
export MAKEFLAGS_DOCBUILD="-j $SAGE_NUM_THREADS_DOCBUILD -l $((CPUTHREADS_DOCBUILD - 1)).8"
