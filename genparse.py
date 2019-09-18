#!/usr/bin/python3

"""
    Generic CSV / data parsing code.

    TODO: add a safe_execute utility?
          https://stackoverflow.com/questions/36671077/one-line-exception-handling/36671208

"""

__author__  = "Sheila M Reynolds"
__version__ = "0.0.2"
__status__  = "Prototype"

import argparse
import ast
import csv
import dateutil
import distutils.util as util 
import json
import logging
import os
import pandas as pd
import re
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

    logging.debug("in estimateFileSize ... <%s>" % fileName)

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
## BigQuery data types = ['string','bytes','integer','float','boolean','record','timestamp']

def getFieldType ( fConfig ):

    ## print (" ... in getFieldType ... ", fConfig)

    if ( 'bqType' in fConfig ):
        return ( fConfig['bqType'] )

    for aKey in fConfig:
        if ( aKey.lower().find("type") > 0 ): 
            if ( fConfig[aKey] == 'boolean' ): return ( 'boolean' )
            if ( fConfig[aKey] == 'datetime' ): return ( 'timestamp' )
            if ( fConfig[aKey] == 'genomic_locus' ): return ( 'string' )
            if ( fConfig[aKey] == 'chromosome' ): return ( 'string' )
            print ( "UHOH in getFieldType: what should I do with this ??? ", aKey, fConfig[aKey] )

    return ( 'string' )

##----------------------------------------------------------------------------------------------------
## optional descriptive text for this field that can go into the BigQuery JSON schema file

def getFieldDescription ( fConfig ):
    if ( 'description' in fConfig ):
        return ( fConfig['description'] )
    else:
        return ('')

##----------------------------------------------------------------------------------------------------
## TODO: get rid of this ???

def isHexBlob(s):

    logging.debug("in isHexBlob ... <%s>" % s[:128])

    if s is None: return (False)
    if not isinstance(s,str): return (False)
    if len(s) == 0: return (False)
    if len(s) < 128: return (False)
    for c in s:
        if c not in string.hexdigits: return (False)
    return (True)

##----------------------------------------------------------------------------------------------------

def stripBlanks(s):

    if ( str(s) == 'nan' ): return ('')
    if ( not isinstance(s,str) ): return (s)

    return ( s.strip() )

##----------------------------------------------------------------------------------------------------
## strip the specified prefix (if present) from the input string ...

def stripPrefixFromStr(s,p):

    ## print ("             in stripPrefixFromStr ... ")
    if ( str(s) == 'nan' ): return ('')

    if ( not isinstance(s,str) ):
        print (" UHOH ??? in stripPrefixFromStr but don't have a string ??? ", s, p)

    if ( s.startswith(p) ):
        np = len(p)
        t = s[len(p):]
        return (t)

    return (s)

##----------------------------------------------------------------------------------------------------
## strip the specified suffix (if present) from the input string ...

def stripSuffixFromStr(s,p):

    ## print ("             in stripSuffixFromStr ... ")
    if ( str(s) == 'nan' ): return ('')

    if ( s.endswith(p) ):
        np = len(p)
        t = s[:-len(p)]
        return (t)

    return (s)

##----------------------------------------------------------------------------------------------------
## handle the 'standard' types ...

def handleStandardType(s,p):

    if ( s == '' ): return (s)
    if ( str(s) == 'nan' ): return ('')

    standardTypeList = [ "boolean", "datetime", "date", "integer", "float", "string" ]

    ## print ("             in handleStandardType ... ", s, p)
    if ( p not in standardTypeList ):
        logging.warning('invalid standard type ' + p)
        print ( " UHOH -- invalid standard type " )
        return (s)

    if ( p=="boolean" ):
        ## try:
        if ( 1 ):
            ## print (" input: <%s> " % s.strip().lower() )
            t = str ( bool ( util.strtobool(s.strip().lower()) ) )
            ## print (" got back ", t )
            return (t)
        ## except:
        ##     print (" FAILED TO interpret/cast to boolean ??? <%s> " % s.strip().lower() )
        ##     sys.exit(-1)

    elif ( p=="datetime" ):
        ##   input: 2018-12-18T15:41:29.554Z
        ##   Python ISO format:  2018-12-18T15:41:29.554000+00:00
        try:
            t = dateutil.parser.parse(s.strip())
            ## print (" input: <%s> " % s.strip() )
            ## print (" ISO format: ", t.isoformat() )
            u = str(t.isoformat() )[:23] + 'Z'
            ## print (" --> <%s> " % u )
            return (u)
        except:
            print (" UHOH FAILED to interpret as datetime ??? ", s )

    elif ( p=="date" ):
        try:
            t = dateutil.parser.parse(s.strip())
            ## print (" input: <%s> " % s.strip() )
            ## print (" ISO format: ", t.isoformat() )
            u = str(t.isoformat() )[:10] 
            ## print (" --> <%s> " % u )
            return (u)
        except:
            if ( s == "00/00/0000" ): return ('')
            print (" UHOH FAILED to interpret as date ??? ", s )

    elif ( p=="integer" ):
        try:
            t = int(s.strip())
            u = str(t)
            return (u)
        except:
            print (" UHOH FAILED to interpret as integer ??? ", s )

    elif ( p=="float" ):
        try:
            t = float(s.strip())
            u = str(t)
            return (u)
        except:
            print (" UHOH FAILED to interpret as float ??? ", s )

    elif ( p=="string" ):
        ## nothing really needs to be done, this 'type' is just being
        ## included for completeness...
        return (s)

    else:
        print (" UHOH  TODO: implement handling for additional pyTypes ... ", p )
        sys.exit(-1)

##----------------------------------------------------------------------------------------------------
## TODO: add custom 'types' for the g.* c.* and p.* notation
##       also the ability to map single-letter amino acid abbreviations to 3-letter abbrev's
def handleCustomType(s,p):

    chrList = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
               '11','12','13','14','15','16','17','18','19', '20',
               '21','22', 'X', 'Y', 'M' ]

    if ( str(s) == 'nan' ): return ('')
    if ( s == '' ): return (s)

    customTypeList = [ "genomic_locus", "chromosome", ]

    ## print ("             in handleCustomType ... ", s, p)

    if ( p not in customTypeList ):
        logging.warning('invalid custom type ' + p)
        print (" UHOH ... invalid custom type ", p)
        return (s)

    if ( p=="genomic_locus" ):
        ## expecting strings like chrX:1234567 or chr1:98765432-98765437

        t = s
        if ( str(t) == 'nan' ): return ('')

        if ( not t.startswith('chr') ): t = 'chr' + t
        u = t.split(':')
        while ( '' in u ): u.remove('')
        ## print(u)

        v = u[0][3:]
        if ( v not in chrList ):
            print ( " INVALID chromosome ??? ", u )
            print (" --> re-setting entire string to blank ", s, t, u)
            return ('')

        ## is this something simple like chrX:1234567 ?
        if ( len(u) == 2 ):
            try:
                ipos = int(u[1])
                ## good to go!
                return (t)

            ## if it's not, then is it something like chr17:112233-123456 ?
            except:
                if ( u[1].find('-') > 0 ):
                    v = u[1].split('-')
                    if ( len(v) == 2 ):
                        try:
                            ipos = int(v[0])
                            jpos = int(v[1])
                            ## good to go!
                            return (t)
                        except:
                            print (" (c) UHOH more complicated ??? ", s, t, u )
                print (" (a) UHOH more complicated ??? ", s, t, u )
                print ("     --> returning (%s) " % t )
                ## oh well, just return what we have ...
                return (t)

        else:

            ## if thre is no ':' then what exactly do we have ???
            ## it's also possible that we have just something like 'chr17:' ...
            if ( len(u) == 1 ):
                t = u[0][3:]
                if t in chrList:
                    ## if its just the chromosome, it's not really a
                    ## valid 'locus' but we'll return it anyway ...
                    return ( u[0] )

                ## otherwise, let's dump it
                print (" UHOH BAD genomic locus ??? !!! <%s> " % t)
                print (" --> re-setting to blank")
                return ('')

            else:

                ## one common error is that the string looks like
                ## this:  chr3:chr3:987654
                ## in which case we want to remove the leading 'chr3:'
                if ( u[1].startswith('chr') ):
                    if ( u[0] == u[1] ):
                        u.remove(u[0])
                        if ( len(u) == 2 ):
                            try:
                                ipos = int(u[1])
                                return ( u[0] + ':' + u[1] )
                            except:
                                print (" (c) UHOH new error type ??? ", t)
                    else:
                        print (" UHOH BAD genomic locus ??? !!! <%s> " % t)
                print (" (b) UHOH more complicated ??? ", t )
                print ("     --> returning (%s) " % t )
                ## oh well, just return what we have ...
                return (t)

        print ( " UHOH ... CRASHING in handleCustomType ... ", s, t )
        sys.exit(-1)

    if ( p=="chromosome" ):
        if ( not s.startswith('chr') ): s = 'chr' + s
        t = s[3:]
        if t not in chrList:
            print ( " UHOH INVALID chromosome ??? ", s, t )
            print (" --> re-setting to blank")
            return ('')
        else:
            ## good to go!
            return (s)

    print (" UHOH TODO: write more code in handleCustomType !!! ", s)
    return (s)

##----------------------------------------------------------------------------------------------------
def reMatchTest(s,r):

    if ( s == '' ): return (s)
    if ( str(s) == 'nan' ): return ('')
    ## print ("             in reMatchTest ... ", s, r)

    t = re.match(r,s)
    if ( t ):
        ## print ("                 --> TRUE ")
        return (s)
    else:
        ## print ("                 --> FALSE ")
        return ("RE-MATCH-FAIL|"+s)

##----------------------------------------------------------------------------------------------------

def enforceAllowed(s,aList,mDict):

    if ( str(s) == 'nan' ): return ('')
    ## print (" in enforceAllowed ... ", aList, mDict)

    for m in mDict:
        b = m.lower()
        t = s.lower()
        if ( t == b ): 
            ## print ("     --> applying mapping", s, m, mDict[m])
            s = mDict[m]

    for a in aList:
        t = s.lower()
        b = a.lower()
        if ( t == b ): 
            ## print ("     --> found allowed match ", s, a)
            return (s)

    print (" UHOH failed to match to allowed strings !!! ", s, aList )
    sys.exit(-1)

##----------------------------------------------------------------------------------------------------

def string2list(s):

    if ( str(s) == 'nan' ): return ('')
    if ( s == '' ): return ('') 

    ## print (" in string2list ... ")
    ## print (" >>> %s <<< " % s )

    sList = []
    for u in s.splitlines():
        u = u.strip()
        if ( len(u) > 0 ): sList += [ u ]

    if ( len(sList) < 1 ): return ('')

    ## print ("     --> returning : ", len(sList), str(sList) )
    return ( str(sList) )

##----------------------------------------------------------------------------------------------------

def quickStrip(s):

    if ( str(s) == 'nan' ): return ('')

    logging.debug("in quickStrip ... <%s>" % str(s)[:88])
    ## print (" in quickStrip ... <%s> " % s )

    t = ''
    for u in s.splitlines():
        ## print ("        <%s> " % u )
        if ( len(u) > 0 ):
            if ( len(t) > 0 ):
                t += ' ' + u
            else:
                t += u

    tabSplit = t.split('\t')
    t = ''
    for u in tabSplit:
        ## print ("        <%s> " % u )
        if ( len(u) > 0 ):
            if ( len(t) > 0 ):
                t += ' ' + u
            else:
                t += u

    ## print ("     returning ... <%s> " % t )
    ## if ( len(t) > 30 ): sys.exit(-1)

    return (t)

##----------------------------------------------------------------------------------------------------

def getMostCommonStrings ( ctrData, n ):

    logging.debug("in getMostCommonStrings ... ")
    ## print (" in getMostCommonStrings ... ", n)

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
    def __init__(self, dataFileName, configFileName):

        ## initialize Logging setup
        ## set filemode='w' to overwrite with new logging file each time
        ##              'a' for appending
        self.logFileName = dataFileName[:-4] + '.logging'
        logging.basicConfig(format='%(asctime)s  %(levelname)s:%(message)s', 
                            datefmt='%Y/%m/%d %I:%M:%S %p',
                            level=logging.DEBUG, 
                            filename=self.logFileName, filemode='w')

        logging.debug("in DataTable.__init__ method")
        self.loadConfig(configFileName)
        self.initDataLoad(dataFileName)
        self.initDLogFile()


    # load configuration information from JSON file
    def loadConfig(self, configFileName):
        logging.debug("in DataTable.loadConfig method " + configFileName )
        try:
            with open(configFileName) as config_file:
                logging.debug(" opened config file " + configFileName)
                self.config = json.load(config_file)
            print(self.config)
        except:
            logging.critical('Failed to read config info from file ' + configFileName)
            self.config = {}
            sys.exit(-1)

    # get set up to read input CSV data file using an iterator
    def initDataLoad(self, dataFileName):
        logging.debug("in DataTable.initDataLoad method " + dataFileName )
        try:
            self.dataFileName = dataFileName
            ( fileSize, rowCountEst, chunkSizeRows, numChunks ) = estimateFileSize(dataFileName)

            self.fileSize = fileSize
            self.rowCountEst = rowCountEst
            self.chunkSizeRows = chunkSizeRows
            self.numChunks = numChunks
            self.ithChunk = 0

            ## print (" Trying to read input CSV from <%s> " % dataFileName)
            ## self.df = pd.read_csv(dataFileName, low_memory=False, skipinitialspace=True, keep_default_na=True, na_values=['[]','{}','na',':'])
            
            print (" getting csv iterator ... ")
            self.csvIter = pd.read_csv(dataFileName, dtype=str, low_memory=False, iterator=True, 
                                       skipinitialspace=True, keep_default_na=True, 
                                       na_values=['[]','{}','na',':','::','EditMe', 'Please edit me'])
            logging.info(' --> got CSV file iterator')

        except:
            logging.critical('Failed to open data file ' + dataFileName)
            print ( sys.exc_info() )
            sys.exit(-1)

    # initialize the data-log file (which is different from the 'logging' file)
    def initDLogFile(self):
            logging.debug("in DataTable.initDLogFile method")
            self.dlogFileName = self.dataFileName[:-4] + '.dlog'
            print (" opening log file <%s> " % self.dlogFileName)
            self.dlogFh = open ( self.dlogFileName, 'w' )

            logHdrString = "colName\tithChunk\tcolType\tlen(colData)\tnumNull\tfracNull\tnumNotNull\tnumUniqVals\tcommonStr\tcommonLen\tcommonN\tminLen\tmaxLen\tavgLen"
            self.dlogFh.write("%s\n" % logHdrString)

    def initOutput(self, outputFile, schemaFile):
        if ( outputFile ): self.outputFile = outputFile
        if ( schemaFile ): self.schemaFile = schemaFile
        self.fhOut = open ( self.outputFile, 'w' )
        self.emptyOutput = True
        ## self.outputSep = '\t'
        self.outputSep = ','

    def handleOneCol(self, x):
        logging.debug("in DataTable.handleOneCol method ... name is " + x.name + " and data-type is " + str(x.dtype))
        print("in DataTable.handleOneCol method ... name is " + x.name + " and data-type is " + str(x.dtype))

        # get configuration parameters for this column
        try:
            xConfig = self.config[x.name]
        except:
            xConfig = {}

        print (" in handleOneCol ... ", json.dumps(xConfig,indent=4) )

        # eliminate blank strings no matter what ...
        x = x.apply ( stripBlanks )

        # stripPrefix
        if ( "stripPrefix" in xConfig ):
            p = xConfig["stripPrefix"]
            ## print ("     --> prefix string : <%s> " % p )
            ## print ("     --> calling x.apply ... " )
            x = x.apply ( stripPrefixFromStr, args=(p,) )

        # stripSuffix
        if ( "stripSuffix" in xConfig ):
            p = xConfig["stripSuffix"]
            ## print ("     --> suffix string : <%s> " % p )
            ## print ("     --> calling x.apply ... " )
            x = x.apply ( stripSuffixFromStr, args=(p,) )

        # reMatch
        if ( "reMatch" in xConfig ):
            r = xConfig["reMatch"]
            ## print ("     --> re match string : <%s> " % r )
            x = x.apply ( reMatchTest, args=(r,) )

        # "built-in" handling for 'standard' Python data types ...
        if ( "pyType" in xConfig ):
            p = xConfig["pyType"]
            ## print ("     --> python type : <%s> " % p )
            x = x.apply ( handleStandardType, args=(p,) )

        # "custom" handling for 'unusual' data types ...
        if ( "customType" in xConfig ):
            p = xConfig["customType"]
            ## print ("     --> custom type : <%s> " % p )
            x = x.apply ( handleCustomType, args=(p,) )

        # other special handling ...
        if ( "stripString" in xConfig ):
            if ( xConfig['stripString'] == "True" ):
                x = x.apply ( quickStrip )

        if ( "string2list" in xConfig ):
            if ( xConfig['string2list'] == "True" ):
                x = x.apply ( string2list )

        if ( "allowed" in xConfig ):
            aList = xConfig['allowed']
            if ( "mappings" in xConfig ):
                mDict = xConfig['mappings']
            else:
                mDict = {}
            if ( len(aList) > 0 ):
                x = x.apply ( enforceAllowed, args=(aList,mDict,) )

        try:
            print ("         --> returning ", x[0], x[-1] )
        except:
            pass

        return (x)
    

    def getNextChunk(self):
        logging.debug("in DataTable.getNextChunk method")
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
            print (" UHOH other error ???")
            return ( False )


    def processChunk(self):
        ## this method processes one chunk of the input dataframe ...

        ## before we do anything else, we can drop any of the fields
        ## that have been flagged as to-be-dropped in the input config
        colNames = self.chunkDf.columns.array
        dropList = []
        for aName in colNames:
            if ( aName in self.config ):
                if ( "dropField" in self.config[aName] ):
                    if ( self.config[aName]["dropField"].lower() == "true" ):
                        dropList += [ aName ]
        if ( len(dropList) > 0 ):
            print (" DROPPING one or more columns ... ", dropList )
            print ( self.chunkDf.shape )
            self.chunkDf = self.chunkDf.drop (columns=dropList)
            print ( self.chunkDf.shape )
            colNames = self.chunkDf.columns.array
            print (" ")
       
        ## first, we apply the handleOneCol to each column in the dataframe
        trim1 = self.chunkDf.apply(self.handleOneCol)
        self.chunkDf = trim1

        ## and now we can git a bit further ...
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
            if ( '' in uniqDict ): del uniqDict['']
            print ("         uniqDict: ", str(uniqDict)[:188])
            numUniqVals = len(uniqDict)
            print ("         numUniqVals: ", numUniqVals)
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
            self.dlogFh.write("%s\n" % outString)

            if ( maxLen > 0 ):
                print ("LONGEST string representation is %d characters" % maxLen )
            if ( maxLen > 1024 ): 
                print ("    THIS IS VERY LONG and may need to be truncated ... " )
                print ("     --> first 128 chars: ", maxStr[:128] )
                print ("     --> last  128 chars: ", maxStr[-128:] )
                print ("  " )
            if ( avgLen > 8192 ):
                print ("     AVERAGE length is %d ... field should probably be ommitted ??? " % avgLen )
    
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
                ## print ( type(nrDict) )
                nrKeys = list(nrDict.keys())
                nrKeys.sort()
                for key in nrKeys:
                    if ( nrDict[key] > 1 ): del nrDict[key]
                print ("number of unique non-null values: ", len(uniqDict) )
                print ("number of non-null AND non-repeating values (occur only once each): ", len(nrDict) )

                ## pull out the single most common value ...
                ( aVal, aCount ) = ctrData.most_common(1)[0]
                print ( aCount )
                print ( type(aVal), len(str(aVal)) )
                ## print ( aVal )

                if isinstance(aVal,float): continue
                if isinstance(aVal,bool):  continue
                if isHexBlob(aVal): 
                    print ("     --> this looks like a Hex Blob ... " )
                    continue

                ## and let's see if 'evaluating' this string produces anything
                ## different ...
                try:
                    ## print (" trying out literal_eval ... " )
                    eVal1 = ast.literal_eval(aVal)
                    if ( eVal1 != aVal ):
                        if isinstance(eVal1,list):
                            print ("     --> result of evaluation is a LIST of length %d " % len(eVal1) )
                            if ( len(eVal1) == 1 ):
                                eVal2 = ast.literal_eval(eVal1[0])

                                if ( eVal2 != eVal1[0] ):
                                    if isinstance(eVal2,list):
                                        print ("     --> next result of evaluation is a LIST of length %d " % len(eVal2) )
                                    elif isinstance(eVal2,dict):
                                        print ("     --> result of evaluation is a DICT of length %d " % len(eVal2) )

                                print ("     2nd literal_eval : \n", type(eVal2), "\n", eVal2 )
                                if isinstance(eVal2,dict):
                                    print ("     --> result of evaluation is a DICT of length %d " % len(eVal2) )
                                    eKeys = list(eVal2.keys())
                                    nKeys = len(eKeys)
                                    print ("         dict has %d keys " % nKeys, eKeys )
                                    eKeys.sort()
                                    aKey = eKeys[nKeys//2]
                                    print ("         for example ... %s: %s " % ( aKey, eVal2[aKey] ) )

                        elif isinstance(eVal1,dict):
                            print ("     --> result of evaluation is a DICT of length %d " % len(eVal1) )
                            eKeys = list(eVal1.keys())
                            nKeys = len(eKeys)
                            print ("     dict has %d keys " % nKeys, eKeys )
                            eKeys.sort()
                            aKey = eKeys[nKeys//2]
                            print ("     for example ... %s: %s " % ( aKey, eVal1[aKey] ) )
                   
                except:
                    ## print ("    failed to run literal_eval on this string <%s> " % str(aVal) )
                    ## print ("    trying to load as JSON instead? " )
                    if ( aVal.find('{')>=0 or aVal.find('[')>=0 ):
                        try:
                            jVal1 = json.loads(str(aVal))
                            bVal = json.dumps(jVal1, indent=4) 
                            if ( bVal != aVal ):
                                print (" YAY!!! JSON!!! " )
                                print ( bVal )
                                print (" " )
                        except:
                            ## print ("         also failed to parse as JSON ... " )
                            pass
                            
                ## now the THREE most common values ...
                if ( aCount > 1 ):
                    c3 = ctrData.most_common(3)
                    print ("3 most common values: \n", c3[0], "\n", c3[1], "\n", c3[2], "\n" )


    def writeOutChunk(self):

        if ( self.emptyOutput ):
            self.chunkDf.to_csv ( self.fhOut, index=False, sep=self.outputSep )
            self.emptyOutput = False
        else:
            self.chunkDf.to_csv ( self.fhOut, index=False, sep=self.outputSep, header=False )

    def writeJsonSchema(self):

        ## BigQuery data types = ['string','bytes','integer','float','boolean','record','timestamp']

        ## one row of the JSON schema file should look like this:
        ## {"name": "_id", "type": "string", "mode": "nullable", "description": "<add description here>"},

        print ( "in writeJsonSchema ... ", self.schemaFile )

        allCols = list ( self.chunkDf.columns )

        with open(self.schemaFile, 'w') as schema_file:
            schema_file.write('[\n')
            numF = len(allCols)
            for iF in range(numF):
                fName = allCols[iF]
                print ( fName, self.config[fName] )
                fType = getFieldType ( self.config[fName] )
                fDesc = getFieldDescription ( self.config[fName] )
                oneLine = '    {"name": "' + fName + '", "type": "' + fType \
                        + '", "mode": "nullable", "description": "' + fDesc \
                        + '"}'
                if ( iF < (numF-1) ): oneLine += ','
                schema_file.write(oneLine+'\n')
            schema_file.write(']')
                



##----------------------------------------------------------------------------------------------------

def main (args):

    myData = DataTable(args.inputFile,args.configFile)

    if ( args.outputFile ): myData.initOutput ( args.outputFile, args.schemaFile )

    done = False
    numChunks = 0
    while not done:
        r = myData.getNextChunk()
        if ( r ): 
            numChunks += 1
            print (" ***************** ")
            print (" **  CHUNK #%3d ** " % numChunks)
            print (" ***************** ")
            myData.processChunk()
            if ( args.outputFile ): myData.writeOutChunk()
        else:
            done = True

    print (" ")
    print (" TOTAL number of chunks handled: %d " % numChunks)

    if ( args.schemaFile ): myData.writeJsonSchema()



##----------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    t0 = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument ('-f', '--inputFileName', action='store', help='input CSV data file', 
                         required=True, dest='inputFile', type=str)
    parser.add_argument ('-c', '--configFileName', action='store', help='input JSON config file', 
                         required=True, dest='configFile', type=str)
    parser.add_argument ('-o', '--outputFileName', action='store', help='output CSV data file', 
                         required=False, dest='outputFile', type=str)
    parser.add_argument ('-s', '--schemaFileName', action='store', help='output JSON schema file', 
                         required=False, dest='schemaFile', type=str)
    args = parser.parse_args()

    main (args)

##----------------------------------------------------------------------------------------------------
