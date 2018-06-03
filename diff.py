#coding=utf-8

"""
The script tells whether two binary files in DOUBLE format are identical or not

Binary will probably be slightly different then due to floating point non-associativity.
ASCII files out to 5-10 decimals should still be identical though.

"""

import sys
import numpy as np

if len(sys.argv) < 3:
    print('Usage: python', sys.argv[0], 'file1 file2')
    exit()

# Input file name
fname1 = sys.argv[1]
fname2 = sys.argv[2]

# Opens input file
data1 = np.fromfile(fname1)
data2 = np.fromfile(fname2)

delta = 1.0e-10
N = len(data1)

if N != len(data2):
    print(fname1, 'and', fname2, 'are not identical')
    print('They are not of the same size:', str(len(data1)), 'VS', str(len(data2)))
    exit()

for i in range(N):
    if abs(data1[i] - data2[i]) > delta:
        print(fname1, 'and', fname2, 'are not identical at', str(i))
        print(str(data1[i]), '', str(data2[i]))
        exit()

#print(fname1, " and ", fname2, "are identical")