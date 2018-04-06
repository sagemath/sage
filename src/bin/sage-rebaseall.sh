#!/bin/dash

# Authors:
# * Gary Zablackis <gzabl@yahoo.com>
# * Dmitrii Pasechnik <dimpase@gmail.com>
# * Jean-Pierre Flori <jean-pierre.flori@ssi.gouv.fr>
#
# Rebase all dlls in the SAGE_LOCAL directory (and its subdirectories)
# as well as the ones already stored in the system database,
# and update the database.
# This system-wide database is located in '/etc/rebase.db.i386' and
# includes the Cygwin system dlls.
#
# Invoke this script from a dash shell after going to the SAGE_LOCAL directory.
# Ensure that no other Cygwin processes are currently running.
# Note that you need write access to the system-wide rebase database
# (which usually means admin rights).

DIR=$(dirname "$(readlink -f "$0")")
exec "$DIR"/sage-rebase.sh --all $@
