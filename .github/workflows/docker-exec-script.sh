#!/bin/sh -x
if [ $# -lt 3 ]; then
    echo >&2 "usage: docker-exec-script.sh CONTAINER WORKDIR [VAR=VALUE...] SCRIPT"
    exit 1
fi
CONTAINER=$1
WORKDIR=$2
shift 2
(echo "cd \"$WORKDIR\"";
 while [ $# -gt 1 ]; do
     echo "export \"$1\""
     shift
 done;
 cat "$1") | docker exec -i $CONTAINER bash -ex
