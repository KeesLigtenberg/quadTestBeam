#!/bin/bash

. ~/init_ilcsoftv02-00-01.sh

if [[ -n ${1} ]]; then
	RUNNUMBER=$1
fi

cd /project/lepcol/users/cligtenb/quadTestBeam/eventBuilder/quadWithSource/

root -b <<EOF
.X processSource.h+("run${RUNNUMBER}_frames.root","run${RUNNUMBER}_clusters_droph.root")
EOF

