#!/bin/bash

root -l <<EOF
.L /project/lepcol/users/cligtenb/quadTestBeam/TrackCombiner/combineTracks.h++
EOF

QUEUE="long"
qsub -q $QUEUE -W group_list=atlas -v FOLDER="/project/lepcol/users/cligtenb/quadTestBeam/run668" doTrackCombiner.sh
qsub -q $QUEUE -W group_list=atlas -v FOLDER="/project/lepcol/users/cligtenb/quadTestBeam/run672" doTrackCombiner.sh
qsub -q $QUEUE -W group_list=atlas -v FOLDER="/project/lepcol/users/cligtenb/quadTestBeam/run676" doTrackCombiner.sh

