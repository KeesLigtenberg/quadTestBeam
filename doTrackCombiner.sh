#!/bin/bash

. ~/init_ilcsoftv02-00-01.sh

if [[ -n $FOLDER ]]; then
	cd $FOLDER
fi

echo "running job" > doCombine.out

#ALIGN=1
if [[ -z $ALIGN ]]; then
root -b >> doCombine.out <<EOF
.X ../TrackCombiner/combineTracks.h+("`ls run???_events.root`", "`ls 201810*_M26_TELESCOPE_combined.root`")
EOF
else
echo "doing alignment" > doCombine.out
root -b >> doCombine.out <<EOF
.L ../TrackCombiner/combineTracks.h+
combineTracksAndUpdateAlignment("`ls run???_events.root`", "`ls 201810*_M26_TELESCOPE_combined.root`")
EOF
fi

echo "finished job" >> doCombine.out
