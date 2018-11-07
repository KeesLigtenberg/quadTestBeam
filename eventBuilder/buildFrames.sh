#!/bin/bash

INPUTFILE=$1
OUTPUTFILE=${INPUTFILE%.root}_frames.root

root -l <<EOF
.L makeHistograms.h+
makeHistograms("${INPUTFILE}", "${OUTPUTFILE}");
EOF
