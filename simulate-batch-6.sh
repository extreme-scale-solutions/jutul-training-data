#!/bin/bash

export GKSwstype=nul

DATADIR=$1
if test "x$DATADIR" = x; then
    echo "Usage: $0 <datadir>"
    exit 1
fi

for THING in $DATADIR/*; do if test -f $THING/after40177days.png; then true; else echo $THING; fi; done | xargs -n64 -P6 julia --project scripts/simulate-batch.jl
