#!/usr/local/bin/python3

"""
    Generic CSV / data parsing code.

    TODO: isolate the 'clean string' code into a simple function, so that I can call it
          after extracting deeper strings in lists/dicts using the eval thingy

    TODO: add a safe_execute utility?
          https://stackoverflow.com/questions/36671077/one-line-exception-handling/36671208

"""

__author__  = "Sheila M Reynolds"
__version__ = "0.0.1"
__status__  = "Prototype"

import argparse
import ast
import csv
import json
import pandas as pd
import string
import sys
import time

from collections import Counter

##----------------------------------------------------------------------------------------------------

def isHexBlob(s):
    if s is None: return ( False)
    if len(s) == 0: return ( False )
    for c in s:
        if c not in string.hexdigits: return ( False )
    return ( True )

##----------------------------------------------------------------------------------------------------

def findFirstTwoX(s,x):
    i1 = s.find(x)
    if ( i1 > 0 and i1 < 10 ):
        i2 = s.find(x,i1+1)
        if ( i2 > 0 and i2 < 10 ):
            ## print ( s, i1, i2 )
            return ( i1, i2 )

    return ( 0, 0 )


def handleDateStr(s):
    
    ## types that we do not want to mess with:
    if s is None: return (s)
    if isinstance(s,float): return(s)
    if isHexBlob(s): return (s)

    if ( s.find('://') > 0 ): return (s)

    ## print ( " in handleDateStr ... %s " % s )

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

    try:
        n1 = int(p1)
        n2 = int(p2)
        n3 = int(p3)
    except:
        return (s)

    ## print ( p1, p2, p3 )
    ## print ( n1, n2, n3 )

    if ( i1==4 and i2==7 ): return (s)

    if ( i1==2 and i2 ==5 ):
        ## we have either 08/25/19 or 08/25/2019
        ## and we want to convert to 2019-08-25
        if ( len(p3) == 2 ): p3 = '20' + p3
        s = p3 + '-' + p1 + '-' + p2
        ## print ( "new date string: %s" % s )
        return (s)

    if ( i1==1 and i2 ==4 ):
        ## we have either 8/25/19 or 8/25/2019
        ## and we want to convert to 2019-08-25
        if ( len(p3) == 2 ): p3 = '20' + p3
        s = p3 + '-' + '0' + p1 + '-' + p2
        ## print ( "new date string: %s" % s )
        return (s)

    print ( " failed to handle date string ??? !!! ", s )
    sys.exit(-1)

##----------------------------------------------------------------------------------------------------

def handleMultiLineStr(s):

    ## types that we do not want to mess with:
    if s is None: return (s)
    if isinstance(s,float): return(s)
    if isHexBlob(s): return (s)

    ## ok, now lets try to handle this...
    ## print ( " in handleMultiLineStr ... \n", s )
    sKeep = s

    ## first we're dumping non-ASCII characters ...
    ## if ( len(s) > 30 ): print ( s )
    t = ''
    for c in s:
        ## print ( s )
        ## print ( c )
        ## print ( ord(c) )
        if ( ord(c) >= 32 and ord(c) <= 126 ): 
            t += c
        else:
            t += ' '

    s = t
    ## print ( " after step 1 :", s )

    ## now we want to get rid of HTML breaks ...
    splitr = [ '<br/>', '<br>', '  ' ]
    ## print ( type(s) )
    ## print ( s )
    for c in splitr:
        ## print ( "     testing for <%s> " % c )
        done = 0
        while not done:
            ii = s.find(c)
            if ( ii >= 0 ):
                ## print ( " first part: <%s>   second part: <%s> " % ( s[:ii], s[ii+len(c):] ) )
                t = s[:ii] + ' ' + s[ii+len(c):]
                s = t
            else:
                done = 1

    ## if ( len(s) > 30 ): print ( " returning: ", s )
    ## if ( s == sKeep ): print ( " NO change " )
    ## if ( s != sKeep ): print ( " SOME change " )

    return ( s )
            
##----------------------------------------------------------------------------------------------------

class DataTable:

    # Class Attributes
    pass

    # Initializer / Instance Attributes
    def __init__(self, dataFileName):
        try:
            print (" Trying to read input CSV from <%s> " % dataFileName)
            self.df = pd.read_csv(dataFileName, skipinitialspace=True, keep_default_na=True, na_values=['[]','{}','na',':'])
            print (" Success!  ")
            print (self.df.shape)
            print (self.df.dtypes.value_counts())
            print (self.df.info())
            print ("\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n")
        except:
            print (" Ooops -- handling an excpetion (?) in DataTable __init__ function ")

    @staticmethod
    def cleanOneCol(x):
        ## print ( " in cleanOneCol ... ", x.dtype )
    
        if ( x.dtype=='object' or x.dtype=='str' ):
            x = x.str.strip()

            x = x.str.replace('\btrue\b', 'True',  case=False,regex=True)
            x = x.str.replace('\bfalse\b','False', case=False,regex=True)
            x = x.str.replace('\byes\b',  'Yes',   case=False,regex=True)
            x = x.str.replace('\bno\b',   'No',    case=False,regex=True)

            x = x.str.replace('\x0d', '', case=True, regex=False)
    
            x = x.apply ( handleDateStr )
            x = x.apply ( handleMultiLineStr )
    
        return ( x )

    def cleanDataFrame(self):
        ## this method will do the 'cleaning' of the data ...

        ## first, we apply the cleanOneCol to each column in the dataframe
        trim1 = self.df.apply(self.cleanOneCol)
        self.df = trim1

        ## and now we can git a bit further ...
        colNames = self.df.columns.array
        colTypes = self.df.dtypes
        numCol = len(colNames)
        numRows = len(self.df)

        for iCol in range(max(numCol,10)):

            ## get the name and the type for this column ...
            colName = colNames[iCol]
            colType = colTypes[iCol]

            print ("\n\n ------------------------------------------------------------------------------------ \n\n")

            ## and then get the data contained in this column
            colData = self.df[colName]

            ## now let's have a look at this data ...
            ## how many null values?  how many non-null values?
            numNull = colData.isnull().sum()
            numNotNull = numRows - numNull

            ## use the collections.Counter to find the unique values and their counts
            ctrData = Counter ( colData.dropna() )
            uniqDict = dict ( ctrData )
            numUniqVals = len(uniqDict)

            print ("====> ", colName, colType, len(colData), numNull, numNotNull, numUniqVals)

            ## let's also figure out how long these strings tend to be ...
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

            if ( numU > 0 ):
                avgLen = sumLen//numU
            else:
                avgLen = -1

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
                c3 = ctrData.most_common(3)
                if ( avgLen < 1024 ):
                    print ("    eg: ")
                    print ("        ", c3[0][0])
                    try:
                        print ("        ", c3[1][0])
                    except:
                        pass
                    try:
                        print ("        ", c3[2][0])
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
                    try:
                        jVal1 = json.loads(str(aVal))
                        print ( " YAY!!! JSON!!! " )
                        print ( json.dumps(jVal1, indent=4) )
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

    myData.cleanDataFrame()

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
