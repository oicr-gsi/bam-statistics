#!/usr/bin/python

## my upstream target intervals were of size 1000 each.


import sys

# use stdin if stdin is full
if not sys.stdin.isatty():
	indexFile = sys.stdin

#otherwise, read from input
else:
	try:
		input_file = sys.argv[1]
	except IndexError:
		message = 'need filename as first argument if stdin is not full'
		raise IndexError(message)
	else:
		indexFile = open(input_file, 'rU')


count = 999;
start = 0;

for line in indexFile:

        if start==0:
                print count
                start = int(line)
                count=count-1
        else:
                val = int(line)
                if start == val - 1:
                        print count
                        count=count-1
                        start = int(line)
                else:
                        print 999
                        count = 998
                        start = int(line)
