#!/bin/csh

/home/astro/phulbz/scripts/splitter.py $1 $2

lcurve

lroche $1 $2 nfile=1 output=d.out scale=yes \\
lroche a.mod $2 nfile=1 output=a.out scale=no \\
lroche b.mod $2 nfile=1 output=b.out scale=no \\
lroche c.mod $2 nfile=1 output=c.out scale=no \\

/home/astro/phulbz/scripts/joiner.py $1 $2
