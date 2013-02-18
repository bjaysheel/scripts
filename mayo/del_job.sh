#!/bin/sh

set -x

qstat -s $1 | awk '{ if ($0 ~ /^[0-9]/) { print $0; } }' | cut -f 1 -d ' ' | xargs qdel
