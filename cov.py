#!/usr/bin/env python
from __future__ import print_function
import sys

if len(sys.argv) < 3:
    print("usage: {} input.bed output.txt".format(sys.argv[0]))
    exit(1)

total = 0

with open(sys.argv[1]) as f,\
        open(sys.argv[2], 'w') as o:
    for line in f:
        array = line.split()
        chr = array[0]
        start = int(array[1])
        end = int(array[2])
        depth = array[4]
        total = total + (end - start)

        while (start < end):
            o.write(chr + "\t" + str(start) + "\t" + depth)
            o.write("\n")
            start = start + 1

print(total)
