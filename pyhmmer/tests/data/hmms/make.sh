#!/bin/sh

for x in txt/* ; do
    prefix=$(basename $x .hmm)
    hmmpress $x && mv $x.h3* -t db
    hmmconvert -b $x > bin/$prefix.h3m
    hmmconvert -2 $x > txt2/$prefix.hmm2
done
