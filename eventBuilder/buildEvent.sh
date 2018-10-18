#!/bin/bash

INPUTFILE=$1
OUTPUTFILE=${INPUTFILE%.root}_events.root

root -l <<EOF
.L buildEvent.h+
convertToTree("${INPUTFILE}", "${OUTPUTFILE}");
EOF
