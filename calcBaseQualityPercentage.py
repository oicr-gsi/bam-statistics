#! /usr/bin/env python
###############################################################################
import sys
import math



inFile = open(sys.argv[1],'r')



for line in inFile:
        data = line.strip().split('\t')
        bp1 = data[0]
        bp2 = data[2]
        bp3 = data[4]
        bp4 = data[6]
        bp5 = data[8]
        bp6 = data[10]

        baseOne = data[1].upper()
        bases2 = data[3].upper()
        bases3 = data[5].upper()
        bases4 = data[7].upper()
        bases5 = data[9].upper()
        bases6 = data[11].upper()

        
        types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}

        i = 0
        qscore = 0
        while i < len(baseOne):
                base = baseOne[i]
                ahead = i+1

                if base == '.' or base == ',':
                        ## print 'oookay'
                        if ahead < len(baseOne) and baseOne[ahead] != '^':
                                ## print 'Yaaay'
                                qscore += 1
                                i += 1

                        elif ahead < len(baseOne) and baseOne[ahead] == '^':
                                ## print 'just once'
                                i += 3
                        else:
                                qscore += 1
                
                                i += 1

                
                else:
                        ## print 'Noooo'
                        i += 1

        # Second data:
        j = 0
        while j < len(bases2):
                base = bases2[j]
                ahead = j+1

                if base == '.' or base == ',':
                        ## print 'oookay'
                        if ahead < len(bases2) and bases2[ahead] != '^':
                                ## print 'Yaaay'
                                qscore += 1
                                j += 1

                        elif ahead < len(bases2) and bases2[ahead] == '^':
                                ## print 'just once'
                                j += 3
                        else:
                                qscore += 1
                
                                j += 1

                
                else:
                        j += 1

        # Third data:
        k = 0
        while k < len(bases3):
                base = bases3[k]
                ahead = k+1

                if base == '.' or base == ',':
                        ## print 'oookay'
                        if ahead < len(bases3) and bases3[ahead] != '^':
                                ## print 'Yaaay'
                                qscore += 1
                                k += 1

                        elif ahead < len(bases3) and bases3[ahead] == '^':
                                ## print 'just once'
                                k += 3
                        else:
                                qscore += 1
                
                                k += 1

                
                else:
                        k += 1

        # Fourth data:
        l = 0
        while l < len(bases4):
                base = bases4[l]
                ahead = l+1

                if base == '.' or base == ',':
                        ## print 'ooolay'
                        if ahead < len(bases4) and bases4[ahead] != '^':
                                ## print 'Yaaay'
                                qscore += 1
                                l += 1

                        elif ahead < len(bases4) and bases4[ahead] == '^':
                                ## print 'just once'
                                l += 3
                        else:
                                qscore += 1
                
                                l += 1

                
                else:
                        l += 1

        # Fifth data:
        m = 0
        while m < len(bases5):
                base = bases5[m]
                ahead = m+1

                if base == '.' or base == ',':
                        ## print 'ooolay'
                        if ahead < len(bases5) and bases5[ahead] != '^':
                                ## print 'Yaaay'
                                qscore += 1
                                m += 1

                        elif ahead < len(bases5) and bases5[ahead] == '^':
                                ## print 'just once'
                                m += 3
                        else:
                                qscore += 1
                
                                m += 1

                
                else:
                        m += 1

        # Sixth data:
        n = 0
        while n < len(bases6):
                base = bases6[n]
                ahead = n+1

                if base == '.' or base == ',':
                        ## print 'ooolay'
                        if ahead < len(bases6) and bases6[ahead] != '^':
                                ## print 'Yaaay'
                                qscore += 1
                                n += 1

                        elif ahead < len(bases6) and bases6[ahead] == '^':
                                ## print 'just once'
                                n += 3
                        else:
                                qscore += 1
                
                                n += 1

                
                else:
                        n += 1



        total = float(bp1) + float(bp2) + float(bp3) + float(bp4) + float(bp5) + float(bp6)

        if total == 0:
                print 0
        else:
                val = qscore/total * 100
                print val       		