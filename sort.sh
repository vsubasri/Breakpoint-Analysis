#!/bin/bash
set -o nounset -o pipefail -o errexit

INPUTDIR=$1/*		
OUTPUTDIR=$2/

if [ ! -d $OUTPUTDIR ]; then
	mkdir "$OUTPUTDIR"
fi

for file in $INPUTDIR
do
	echo "$OUTPUTDIR$(basename $file)"
	sort -k1,1 -k2,2n "$file" > "$OUTPUTDIR$(basename $file)"
done



