#!/bin/bash

root -b 	<<EOF
.L /project/lepcol/users/cligtenb/quadTestBeam/eventBuilder/quadWithSource/processSource.h++
EOF

QUEUE="short"

#for i in {779..787}; do
#	qsub -q $QUEUE -W group_list=atlas -v RUNNUMBER=$i processSource.sh
#done
for i in {770..777}; do
	qsub -q $QUEUE -W group_list=atlas -v RUNNUMBER=$i processSource.sh
done
#for i in {818..827}; do
#	qsub -q $QUEUE -W group_list=atlas -v RUNNUMBER=$i processSource.sh
#done

