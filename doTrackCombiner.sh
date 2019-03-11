#!/bin/bash

. ~/init_ilcsoftv02-00-01.sh

if [[ -n $FOLDER ]]; then
	cd $FOLDER
fi

echo "runing job" > doCombine.out

root -b >> doCombine.out <<EOF
.X ../TrackCombiner/combineTracks.h+("`ls run???_events.root`", "`ls 201810*_M26_TELESCOPE_combined.root`")
EOF


echo "finished job" >> doCombine.out
