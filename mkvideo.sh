#!/bin/sh

if test x$1 = x; then
    echo "Usage: $0 <datadir> [<datadir> ...]"
    exit 1
fi

for DATADIR in "$@"; do
    if ! test -d $DATADIR; then
        echo "Datadir $DATADIR does not exist, or is not a dir?"
        exit 1
    fi
    if test x`find $DATADIR -type f -name 'after*.png' | wc -l` = x0; then
        echo "No png files in this dir, skipping."
        continue
    fi
    SEED=`basename $DATADIR`
    if test -f $DATADIR/$SEED.mp4; then
        echo "mp4 file already exists, skipping."
        continue
    fi
    ffmpeg -y -framerate 5 -pattern_type glob -i "$DATADIR/after*.png" -c:v libx264 -pix_fmt yuv420p $DATADIR/$SEED.mp4
done
