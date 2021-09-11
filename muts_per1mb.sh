#!/bin/bash

cp $file $file.tmp
sed '$ d' $file.tmp > $file
rm -f $file.tmp
sed -e 's/^/chr/' $file > $file
bedtools coverage -a 1mb.window -b $file -counts > intersected
