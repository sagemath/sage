#!/bin/sh

OFFSET=1024
# This script reads from stdin and prints to stdout as long as a the output
# does not exceed a certain number of bytes. When reading an EOF it prints the
# last $OFFSET lines if they have not been printed normally already.
# This script expects one argument, the number of bytes.

# Heavily inspired by a simlar strategy in make, https://stackoverflow.com/a/44849696/812379.

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

stdbuf -i0 -o0 -e0 awk -v limit=$1 -v firstMissingNR=-1 -v offset=$OFFSET -v bytes=0 \
'{
  if (bytes < limit) {
    # this probably gets multi-byte characters wrong, but that probably does
    # not matter for our purposes. (We add 1 for a UNIX newline.)
    bytes += length($0) + 1;
    print;
  } else {
    if (firstMissingNR == -1){
      print "[…output truncated…]";
      firstMissingNR = NR;
    }
    a[NR] = $0;
    delete a[NR-offset];
    printf "." > "/dev/stderr"
    if (/^sage-logger/){
      print > "/dev/stderr"
    }
  }
}
END {
  if (firstMissingNR != -1) {
    print "" > "/dev/stderr";
    for(i = NR-offset+1 > firstMissingNR ? NR-offset-1 : firstMissingNR; i<=NR ; i++){ print a[i]; }
  }
}
'

