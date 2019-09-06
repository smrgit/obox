#!/usr/bin/python3

import  argparse
import  csv
import  sys
import  time

##------------------------------------------------------------------------------

def cleanToken ( a ):

  double_chars = [ '||', ' |', '| ' ]
  for d in double_chars:
    if ( a.find(d) >= 0 ): 
      b = a.replace ( d, '|' )
      a = b
    if ( a.startswith(d) ): a = a[2:]
    if ( a.endswith(d) ): a = a[:-2]

  double_chars = [ '  ' ]
  for d in double_chars:
    if ( a.find(d) >= 0 ): 
      b = a.replace ( d, ' ' )
      a = b
    if ( a.startswith(d) ): a = a[2:]
    if ( a.endswith(d) ): a = a[:-2]


  spec_chars = [ '\\', '\r', '\n', '\t' ]
  for s in spec_chars:
    if ( a.find(s) >= 0 ): 
      b = a.replace ( s, '' )    
      a = b
    if ( a.startswith(s) ): a = a[1:]
    if ( a.endswith(s) ): a = a[:-1]

  if ( a.startswith(' ') ): a = a[1:]
  if ( a.endswith(' ') ): a = a[:-1]

  if ( a == '{}' ): a = ''
  if ( a == '[]' ): a = ''
  if ( a.lower() == 'na' ): a = ''
  if ( a.lower() == 'n/a' ): a = ''
  if ( a.lower() == 'unspecified' ): a = ''
  if ( a == ':' ):  a = ''

  return ( a )

##------------------------------------------------------------------------------

def getBasicInfo ( foLog ):

    numTokens = 14

    firstField = ''
    fieldNames = []
    numRows = 0

    foLog.seek(0)
    lineNo = 0

    for aLine in foLog:

        if ( aLine.startswith("colName") ): continue

        lineNo += 1
        if ( ord(aLine[-1]) == 10 ): aLine = aLine[:-1]
        tokenList = aLine.split('\t')
        if ( len(tokenList) != 14 ):
            print ( " unexpected number of tokens (%d) at line #%d ??? " % (len(tokenList), lineNo) )
            print ( tokenList )
        else:
            aName = tokenList[0]
            if ( firstField == '' ): firstField = aName
            if ( aName not in fieldNames ): fieldNames += [ aName ]
            numChunks = int ( tokenList[1] )
            if ( aName == firstField ):
                numRows += int ( tokenList[3] )

    return ( fieldNames, numChunks, numRows )

##------------------------------------------------------------------------------

def getInfo ( foLog, fieldNames, numChunks ):

    foLog.seek(0)
    lineNo = 0

    d = {}
    for aName in fieldNames:
        d[aName] = [0] * numChunks

    for aLine in foLog:

        if ( aLine.startswith("colName") ): continue

        lineNo += 1
        if ( ord(aLine[-1]) == 10 ): aLine = aLine[:-1]
        tokenList = aLine.split('\t')
        if ( len(tokenList) != 14 ):
            print ( " unexpected number of tokens (%d) at line #%d ??? " % (len(tokenList), lineNo) )
            print ( tokenList )
        else:
            aName = tokenList[0]
            ithChunk = int ( tokenList[1] )
            if ( ithChunk <= numChunks ): 
                d[aName][ithChunk-1] = tokenList

    return ( d )

##------------------------------------------------------------------------------

def examineInfo ( d ):

    fieldNames = list(d.keys())
    fieldNames.sort()

    ##  colName  ithChunk   colType   len(colData)  numNull  fracNull   numNotNull   numUniqVals   commonStr                             commonLen  commonN  minLen  maxLen  avgLen
    ## ['_id',   '24',      'object', '240000',     '0',     '0.0',     '240000',    '240000',     'ObjectId(585aeff950b20e44120d558e)', '34',      '1',     '34',   '34',   '34\n']
    ##  0        1          2         3             4        5          6            7             8                                     9          10       11      12      13

    for aField in fieldNames:

        print ( " " )
        print ( "%s : " % aField )

        numChunks = len(d[aField])

        ## initialize various fields and counters ...
        typeList = []
        totRows = 0
        totNull = 0
        totUniq = 0
        cmnList = []
        cmnCount = 0
        minLen = -1
        maxLen = -1
        avgLen = 0
        avgN = 0

        for ii in range(numChunks):

            ## first check for a consistent type
            try:
                typeToken = d[aField][ii][2] 
                if ( typeToken not in typeList ): typeList += [ typeToken ]
            except:
                pass

            ## count up all of the data examined
            try:
                n = int ( d[aField][ii][3] )
                totRows += n
            except:
                pass

            ## count up all of the null values found
            try:
                n = int ( d[aField][ii][4] )
                totNull += n
            except:
                pass


            ## count up the number of unique values found
            try:
                n = int ( d[aField][ii][7] )
                totUniq += n
            except:
                pass

            ## get some common values (as strings) as well as counts and lengths
            try:
                aString = d[aField][ii][8] 
                if ( len(aString) > 0 ):
                    if ( aString not in cmnList ): cmnList += [ aString ]
                n = int ( d[aField][ii][10] )
                if ( len(aString) == 0 ):
                    if ( n > 0 ):
                        print ( " why is this happening ??? <%s> %d " % (aString, n) )
                        print ( aField, ii, n, d[aField][ii] )
                        sys.exit(-1)
                cmnCount += n
                n = int ( d[aField][ii][11] )
                if ( minLen == -1 ): minLen = n
                if ( minLen > n ): minLen = n
                n = int ( d[aField][ii][12] )
                if ( maxLen == -1 ): maxLen = n
                if ( maxLen > n ): maxLen = n
                n = int ( d[aField][ii][13] )
                if ( n > 0 ):
                    avgLen += n
                    avgN += 1
            except:
                pass


        ## done looping ...

        if ( avgN > 0 ): avgLen = avgLen / avgN

        if ( totNull == totRows ):
            print ( "         ALWAYS null !!! " )
        else:
            print ( "         types found : ", typeList )
            print ( "         totRows=%d   totNull=%d  fracNull=%f " % ( totRows, totNull, float(totNull)/float(totRows) ) )
            print ( "         number of unique values: ", totUniq )
            print ( "         common strings: ", cmnList )
            print ( "         total occurrences: ", cmnCount )
            if ( avgN > 0 ):
                print ( "         range of string lengths: (%d,%d) with avg=%d" % ( minLen, maxLen, avgLen ) )


##------------------------------------------------------------------------------

def main ( args ):

    print ( args )
  
    try:
        foLog = open ( args.logFile, 'r' )
    except:
        print ( " Failed to open log file <{}> ... EXITING ".format ( args.logFile ) )
        sys.exit(-1)
  
    ( fieldNames, numChunks, numRows ) = getBasicInfo ( foLog )
    print ( len(fieldNames), numChunks, numRows )
    print ( fieldNames )

    d = getInfo ( foLog, fieldNames, numChunks )
    print ( len(d) )

    examineInfo ( d )

  

##------------------------------------------------------------------------------

if __name__ == '__main__':

  t0 = time.time()

  parser = argparse.ArgumentParser()

  ## the first two arguments are required -- to fully specify the BQ table of interest
  parser.add_argument ( '-f',  '--inputLogFileName',  action='store', help='input log file (output of genparse)', required=True, dest='logFile', type=str )
  args = parser.parse_args()

  main ( args )

  t1 = time.time()

  ## print ( ' --> time taken in seconds: {dt}'.format(dt=(t1-t0)) )

##------------------------------------------------------------------------------
