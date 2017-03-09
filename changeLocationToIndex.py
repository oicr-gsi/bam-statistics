#!/usr/bin/python

import sys

with open(sys.argv[1]) as indexFile:

        count = 119;
        start = 0;

        while True:
                line1 = indexFile.readline()

                if start==0:
                        print count
                        start = int(line1)
                        count=count-1

                else:
                        val = int(line1)
                        if start == val - 1:
                                print count
                                count=count-1
                                start = int(line1)


                        else:
                                print 119
                                count = 118
                                start = int(line1)
