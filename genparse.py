#!/usr/bin/python3

"""
    Generic CSV / data parsing code.

    TODO: add a safe_execute utility?
          https://stackoverflow.com/questions/36671077/one-line-exception-handling/36671208

    TODO: still having a problem with multi-line strings -- need to strip them in BS

"""

__author__  = "Sheila M Reynolds"
__version__ = "0.0.1"
__status__  = "Prototype"

import argparse
import ast
import csv
import json
import os
import pandas as pd
import string
import sys
import time
import warnings

from bs4 import BeautifulSoup as bs
from collections import Counter

##----------------------------------------------------------------------------------------------------

## BeautifulSoup issues a stern warning if the input text looks like a URL, but I don't care...
warnings.filterwarnings ( "ignore", category=UserWarning, module='bs4' )

##----------------------------------------------------------------------------------------------------

def estimateFileSize(fileName):

    print ("in estimateFileSize ... <%s>" % fileName)
    print (" ")

    chunkSizeBytes = 333333333
    chunkSizeBytes = 111111111

    fileSize = os.path.getsize(fileName)
    if ( fileSize < chunkSizeBytes ):
        chunkSizeBytes = fileSize
    elif ( fileSize < 5*chunkSizeBytes ):
        chunkSizeBytes = fileSize//5

    fh = open ( fileName, "r+" )
    lineBuf = fh.readlines(chunkSizeBytes)
    fh.close()

    nBuf = len(lineBuf)
    rowCountEst = int ( nBuf * float(fileSize) / float(chunkSizeBytes) )

    print ( fileSize, chunkSizeBytes, nBuf, rowCountEst )

    chunkSizeRows = float(chunkSizeBytes) / (float(fileSize)/float(rowCountEst)) 
    if ( chunkSizeRows < 0.90*rowCountEst ):
        if (   chunkSizeRows > 10000000 ): rndN = 1000000
        elif ( chunkSizeRows >  1000000 ): rndN =  100000
        elif ( chunkSizeRows >   100000 ): rndN =   10000
        else: rndN = 1000
        chunkSizeRows = rndN * int ( 1 + chunkSizeRows/rndN )
        numChunks = rowCountEst//chunkSizeRows + 1
    else:
        chunkSizeRows = -1
        numChunks = 1

    print ("    file size = %d MB" % (fileSize//1000000))
    print ("    row count estimate = %d" % (rowCountEst))
    if ( chunkSizeRows > 0 ):
        print ("    recommended chunk size = %d" % (chunkSizeRows))
        print ("    estimated number of chunks = %d" % (numChunks))
    else:
        print ("    file can be read in one chunk")
        chunkSizeRows = nBuf + 1000

    return ( fileSize, rowCountEst, chunkSizeRows, numChunks )

##----------------------------------------------------------------------------------------------------

def isTimeStamp(s):

    if not isinstance(s,str): return (False)
    if len(s) == 0: return (False)

    ## 2016-02-01T16:15:04.934Z
    ## 012345678901234567890123
    try:
        if ( s[-1] != 'Z' ): return (False)
    except:
        pass

    try:
        if ( s[4] != '-' ): return (False)
    except:
        pass

    try:
        if ( s[7] != '-' ): return (False)
    except:
        pass

    try:
        if ( s[10] != 'T' ): return (False)
    except:
        pass

    return (True)

##----------------------------------------------------------------------------------------------------

def isHexBlob(s):
    if s is None: return (False)
    if not isinstance(s,str): return (False)
    if len(s) == 0: return (False)
    if len(s) < 128: return (False)
    for c in s:
        if c not in string.hexdigits: return (False)
    return (True)

##----------------------------------------------------------------------------------------------------

def findFirstTwoX(s,x):
    i1 = s.find(x)
    if ( i1 > 0 and i1 < 10 ):
        i2 = s.find(x,i1+1)
        if ( i2 > 0 and i2 < 10 ):
            ## print ( s, i1, i2 )
            return ( i1, i2 )

    return ( 0, 0 )

##----------------------------------------------------------------------------------------------------
## TODO: I could potentially switch over to python-dateutil to deal with this ...

def handleDateStr(s):
    
    ## print ( " (a) in handleDateStr ... %s " % str(s)[:32] )

    ## types that we do not want to mess with:
    if s is None: return (s)
    if not isinstance(s,str): return(s)
    if len(s) == 0: return (s)
    if isHexBlob(s): return (s)
    if isTimeStamp(s): return (s)

    if ( s.find('://') > 0 ): return (s)

    ## print ( " (b) in handleDateStr ... %s " % s[:32] )

    ## we are looking at strings that start with something like:
    ##      2019-08-25
    ## or   2019/08/25
    ## or   8/25/2019
    ## or   8/5/2019
    ## or   8/5/19
        
    if ( len(s) < 6 ): return (s)

    i1, i2 = findFirstTwoX ( s, '-' )
    if ( i2 == 0 ):
        i1, i2 = findFirstTwoX ( s, '/' )
    if ( i2 == 0 ): return (s)

    ## print ( i1, i2, s[:30] )

    if ( (i2-i1) != 3 ): return (s)

    p1 = s[:i1]
    p2 = s[i1+1:i2]
    p3 = s[i2+1:]
    p4 = ''

    ## first let's test that the first two substrings 
    ## convert to integers w/o any problem ...
    try:
        n1 = int(p1)
        n2 = int(p2)
    except:
        return (s)

    ## the third substring might fail if it's too long
    ## or has something else coming after the date ...
    try:
        n3 = int(p3)
    except:
        if len(p3) < 2: return (s)
        i3 = p3.find(' ')
        if ( i3 < 0 ): i3 = p3.find('T')
        p4 = p3[i3:]
        p3 = p3[:i3]
        ## print ( " split into <%s> and <%s> " % (p3, p4) )

    ## print ( p1, p2, p3, p4 )
    ## print ( n1, n2, n3 )

    ## 012345678901234567890
    ##     -  -
    if ( i1==4 and i2==7 ):
        s = p1 + '-' + p2 + '-' + p3 + p4
        ## print ( "     returning <%s> " % s )
        return (s)

    if ( i1==2 and i2 ==5 ):
        ## we have either 08/25/19 or 08/25/2019
        ## and we want to convert to 2019-08-25
        if ( len(p3) == 2 ): p3 = '20' + p3
        s = p3 + '-' + p1 + '-' + p2 + p4
        ## print ( "     returning <%s> " % s )
        return (s)

    if ( i1==1 and i2 ==4 ):
        ## we have either 8/25/19 or 8/25/2019
        ## and we want to convert to 2019-08-25
        if ( len(p3) == 2 ): p3 = '20' + p3
        s = p3 + '-' + '0' + p1 + '-' + p2 + p4
        ## print ( "     returning <%s> " % s )
        return (s)

    print ( " failed to handle date string ??? !!! ", s )
    sys.exit(-1)

##----------------------------------------------------------------------------------------------------

def quickStrip(s):

    ## print ( " in quickStrip ... <%s> " % s )

    t = ''
    for u in s.splitlines():
        ## print ( "        <%s> " % u )
        if ( len(u) > 0 ):
            if ( len(t) > 0 ):
                t += ' ' + u
            else:
                t += u

    tabSplit = t.split('\t')
    t = ''
    for u in tabSplit:
        ## print ( "        <%s> " % u )
        if ( len(u) > 0 ):
            if ( len(t) > 0 ):
                t += ' ' + u
            else:
                t += u

    ## print ( "     returning ... <%s> " % t )
    ## if ( len(t) > 30 ): sys.exit(-1)

    return(t)

##----------------------------------------------------------------------------------------------------

def handleMultiLineStr(s):

    ## types that we do not want to mess with:
    if s is None: return (s)
    if not isinstance(s,str): return(s)
    if isHexBlob(s): return (s)

    ## I had written a bunch of hacky code and then realized that I could just use BeautifulSoup!!!
    if ( s.find('<')>= 0 or s.find('\\')>=0 ):
        ## print (" trying to make soup ... ", len(s))
        soup = bs(s)
        return ( soup.get_text(" ",strip=True) )
    else:
        return ( quickStrip(s) )

##----------------------------------------------------------------------------------------------------

def getMostCommonStrings ( ctrData, n ):

    ## cTmp will be a list of tuples
    cTmp = ctrData.most_common(n+1)

    cOut = []
    for t in cTmp:
        if ( t[0] != '' ):
            cOut += [ t ]

    ## print ( type(cTmp) )
    ## print ( cTmp )
    ## print ( cOut[:n] )

    return ( cOut[:n] )

##----------------------------------------------------------------------------------------------------

class DataTable:

    # Class Attributes
    pass

    # Initializer / Instance Attributes
    def __init__(self, dataFileName):

        try:
            self.dataFileName = dataFileName
            ( fileSize, rowCountEst, chunkSizeRows, numChunks ) = estimateFileSize(dataFileName)

            self.fileSize = fileSize
            self.rowCountEst = rowCountEst
            self.chunkSizeRows = chunkSizeRows
            self.numChunks = numChunks
            self.ithChunk = 0

            self.logFileName = dataFileName[:-4] + '.log'
            print (" opening log file <%s> " % self.logFileName)
            self.logFh = open ( self.logFileName, 'w' )

            hdrString = "colName\tithChunk\tcolType\tlen(colData)\tnumNull\tfracNull\tnumNotNull\tnumUniqVals\tcommonStr\tcommonLen\tcommonN\tminLen\tmaxLen\tavgLen"
            self.logFh.write("%s\n" % hdrString)

            ## print (" Trying to read input CSV from <%s> " % dataFileName)
            ## self.df = pd.read_csv(dataFileName, low_memory=False, skipinitialspace=True, keep_default_na=True, na_values=['[]','{}','na',':'])
            
            print (" getting csv iterator ... ")
            self.csvIter = pd.read_csv(dataFileName, dtype=str, low_memory=False, iterator=True, 
                                       skipinitialspace=True, keep_default_na=True, 
                                       na_values=['[]','{}','na',':','EditMe'])
            print (" --> got iterator! ")

        except:
            print (" Ooops -- handling an exception (?) in DataTable __init__ function ")
            print (" Assuming FATAL ERROR ")
            print ( sys.exc_info() )
            sys.exit(-1)


    @staticmethod
    def cleanOneCol(x):
        ## print ( " in cleanOneCol ... ", x.dtype )
    
        if ( x.dtype=='object' or x.dtype=='str' ):

            try:
                x = x.str.strip()
    
                x = x.str.replace('\btrue\b', 'True',  case=False,regex=True)
                x = x.str.replace('\bfalse\b','False', case=False,regex=True)
                x = x.str.replace('\byes\b',  'Yes',   case=False,regex=True)
                x = x.str.replace('\bno\b',   'No',    case=False,regex=True)
    
                x = x.str.replace('\x0d', '', case=True, regex=False)
            except:
                print ("         WARNING: in cleanOneCol, failed in the x.str.* block ")
    
            x = x.apply ( handleDateStr )
            x = x.apply ( handleMultiLineStr )
    
        return ( x )

    def getNextChunk(self):
        try:
            print (" in getNextChunk ... ", self.chunkSizeRows)
            self.chunkDf = self.csvIter.get_chunk(self.chunkSizeRows)
            print (" --> back from get_chunk()")
            print ( self.chunkDf.shape )
            self.ithChunk += 1
            return ( True )
        except StopIteration:
            print (" --> StopIteration exception: end of iteration")
            return ( False )
        except:
            print (" other error ???")
            return ( False )


    def examineChunk(self):
        ## this method will do the 'cleaning' of the data ...

        ## first, we apply the cleanOneCol to each column in the dataframe
        trim1 = self.chunkDf.apply(self.cleanOneCol)
        self.chunkDf = trim1

        ## and now we can git a bit further ...
        colNames = self.chunkDf.columns.array
        colTypes = self.chunkDf.dtypes
        numCol = len(colNames)
        numRows = len(self.chunkDf)

        for iCol in range(max(numCol,10)):

            ## get the name and the type for this column ...
            colName = colNames[iCol]
            colType = colTypes[iCol]

            print ("\n\n ------------------------------------------------------------------------------------ \n\n")

            ## and then get the data contained in this column
            colData = self.chunkDf[colName]

            ## now let's have a look at this data ...
            ## how many null values?  how many non-null values?
            numNull = colData.isnull().sum()
            fracNull = float(numNull) / float(numRows)
            numNotNull = numRows - numNull

            ## use the collections.Counter to find the unique values and their counts
            ctrData = Counter ( colData.dropna() )
            uniqDict = dict ( ctrData )
            numUniqVals = len(uniqDict)
            commonStrings = getMostCommonStrings ( ctrData, 5 ) 

            print ("====> ", colName, colType, len(colData), numNull, numNotNull, numUniqVals)

            outString = colName + "\t" + str(self.ithChunk) + "\t" + str(colType) + "\t" + str(len(colData)) + "\t" 
            outString += str(numNull) + "\t" + str(fracNull) + "\t" + str(numNotNull) + "\t" + str(numUniqVals) + "\t"

            ## get the single most common value
            ## and write it and and its length, and its count to the outString

            if ( numUniqVals > 0 ):
                c1 = commonStrings[0]
                try:
                    outString += str(c1[0])[:64] + "\t" + str(len(str(c1[0]))) + "\t" + str(c1[1]) + "\t"
                except:
                    outString += str(c1[0])[:64] + "\t" + "" + "\t" + str(c1[1]) + "\t"
            else:
                outString += "\t\t\t"

            ## let's also figure out how long these strings tend to be ...
            minLen = 99999999
            maxLen = 0
            sumLen = 0
            numU = 0
            for u in list(uniqDict.keys()): 
                uLen = len(str(u))
                sumLen += uLen
                numU += 1
                if ( maxLen < uLen ): 
                    maxLen = uLen
                    maxStr = str(u)
                if ( minLen > uLen ):
                    minLen = uLen
                    minStr = str(u)

            if ( numU > 0 ):
                avgLen = sumLen//numU
            else:
                avgLen = -1
                minLen = -1
                maxLen = -1

            outString += str(minLen) + "\t" + str(maxLen) + "\t" + str(avgLen) 
            self.logFh.write("%s\n" % outString)

            if ( maxLen > 0 ):
                print ( "LONGEST string representation is %d characters" % maxLen )
            if ( maxLen > 1024 ): 
                print ( "    THIS IS VERY LONG and may need to be truncated ... " )
                print ( "     --> first 128 chars: ", maxStr[:128] )
                print ( "     --> last  128 chars: ", maxStr[-128:] )
                print ( "  " )
            if ( avgLen > 8192 ):
                print ( "     AVERAGE length is %d ... field should probably be ommitted ??? " % avgLen )
    
            if ( numUniqVals == 0 ):
                print ("ALWAYS empty or null!")
            elif ( numUniqVals == numNotNull ):
                print ("ALL of the non-null values are UNIQUE!")
                if ( avgLen < 1024 ):
                    print ("    eg: ")
                    print ("        ", commonStrings[0][0])
                    try:
                        print ("        ", commonStrings[1][0])
                    except:
                        pass
                    try:
                        print ("        ", commonStrings[2][0])
                    except:
                        pass

            else:

                if ( numNull > (0.90*numRows) ):
                    print ("More than 90% of the values are NULL!")
                elif ( numNull > (0.50*numRows) ):
                    print ("More than half of the values are NULL!")
                elif ( numNull < (0.10*numRows) ):
                    print ("Fewer than 10% of the values are NULL!")

            if ( numUniqVals == 0 ):
                pass

            elif ( numUniqVals <= 6 ):
                print ("unique values: ", ' -;- '.join(str(x)[:80] for x in uniqDict))
            else:

                ## let's get a subset of the uniqDict -- only those values that are never repeated:
                nrDict = uniqDict
                print ( type(nrDict) )
                nrKeys = list(nrDict.keys())
                nrKeys.sort()
                for key in nrKeys:
                    if ( nrDict[key] > 1 ): del nrDict[key]
                print ( "number of unique non-null values: ", len(uniqDict) )
                print ( "number of non-null AND non-repeating values (occur only once each): ", len(nrDict) )

                ## pull out the single most common value ...
                ( aVal, aCount ) = ctrData.most_common(1)[0]
                print ( aCount )
                print ( type(aVal), len(str(aVal)) )
                ## print ( aVal )

                if isinstance(aVal,float): continue
                if isinstance(aVal,bool):  continue
                if isHexBlob(aVal): 
                    print ( "     --> this looks like a Hex Blob ... " )
                    continue

                ## and let's see if 'evaluating' this string produces anything
                ## different ...
                try:
                    ## print ( " trying out literal_eval ... " )
                    eVal1 = ast.literal_eval(aVal)
                    if ( eVal1 != aVal ):
                        if isinstance(eVal1,list):
                            print ( "     --> result of evaluation is a LIST of length %d " % len(eVal1) )
                            if ( len(eVal1) == 1 ):
                                eVal2 = ast.literal_eval(eVal1[0])

                                if ( eVal2 != eVal1[0] ):
                                    if isinstance(eVal2,list):
                                        print ( "     --> next result of evaluation is a LIST of length %d " % len(eVal2) )
                                    elif isinstance(eVal2,dict):
                                        print ( "     --> result of evaluation is a DICT of length %d " % len(eVal2) )

                                print ( "     2nd literal_eval : \n", type(eVal2), "\n", eVal2 )
                                if isinstance(eVal2,dict):
                                    print ( "     --> result of evaluation is a DICT of length %d " % len(eVal2) )
                                    eKeys = list(eVal2.keys())
                                    nKeys = len(eKeys)
                                    print ( "         dict has %d keys " % nKeys, eKeys )
                                    eKeys.sort()
                                    aKey = eKeys[nKeys//2]
                                    print ( "         for example ... %s: %s " % ( aKey, eVal2[aKey] ) )

                        elif isinstance(eVal1,dict):
                            print ( "     --> result of evaluation is a DICT of length %d " % len(eVal1) )
                            eKeys = list(eVal1.keys())
                            nKeys = len(eKeys)
                            print ( "     dict has %d keys " % nKeys, eKeys )
                            eKeys.sort()
                            aKey = eKeys[nKeys//2]
                            print ( "     for example ... %s: %s " % ( aKey, eVal1[aKey] ) )
                   
                except:
                    ## print ( "    failed to run literal_eval on this string <%s> " % str(aVal) )
                    ## print ( "    trying to load as JSON instead? " )
                    if ( aVal.find('{')>=0 or aVal.find('[')>=0 ):
                        try:
                            jVal1 = json.loads(str(aVal))
                            bVal = json.dumps(jVal1, indent=4) 
                            if ( bVal != aVal ):
                                print ( " YAY!!! JSON!!! " )
                                print ( bVal )
                                print ( " " )
                        except:
                            ## print ( "         also failed to parse as JSON ... " )
                            pass
                            
                ## now the THREE most common values ...
                if ( aCount > 1 ):
                    c3 = ctrData.most_common(3)
                    print ( "3 most common values: \n", c3[0], "\n", c3[1], "\n", c3[2], "\n" )


##----------------------------------------------------------------------------------------------------

def main (args):

    myData = DataTable(args.inputFile)

    done = False
    numChunks = 0
    while not done:
        r = myData.getNextChunk()
        if ( r ): 
            numChunks += 1
            print (" ***************** ")
            print (" **  CHUNK #%3d ** " % numChunks)
            print (" ***************** ")
            myData.examineChunk()
        else:
            done = True

    print (" ")
    print (" TOTAL number of chunks handled: %d " % numChunks)

    sys.exit(-1)


##----------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    t0 = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument ('-f', '--inputFileName', action='store', help='input CSV file', 
                         required=True, dest='inputFile', type=str)
    args = parser.parse_args()

    main (args)

##----------------------------------------------------------------------------------------------------
