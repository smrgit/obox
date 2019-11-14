#### #!/usr/bin/python3

"""
    Generic CSV / data parsing code.
"""

__author__ = "Sheila M Reynolds"
__version__ = "0.0.2"
__status__ = "Prototype"

import argparse
import ast
import csv
import datetime
import distutils.util as util
import json
import logging
import os
import re
import string
import sys
import time
import warnings
from collections import Counter

import dateutil
import pandas as pd

from Bio.SeqUtils import seq3
from biothings_client import get_client
from bs4 import BeautifulSoup as bs

# ----------------------------------------------------------------------------------------------------

# BeautifulSoup issues a stern warning if the input text looks like a URL, but I don't care...
warnings.filterwarnings("ignore", category=UserWarning, module='bs4')

gene_client = get_client('gene')
geneSymbolPattern = re.compile("[A-Z0-9]+")

# ----------------------------------------------------------------------------------------------------


def estimateFileSize(fileName):

    logging.debug("in estimateFileSize ... <%s>" % fileName)

    chunkSizeBytes = 333333333
    chunkSizeBytes = 111111111

    fileSize = os.path.getsize(fileName)
    if (fileSize < chunkSizeBytes):
        chunkSizeBytes = fileSize
    elif (fileSize < 5*chunkSizeBytes):
        chunkSizeBytes = fileSize//5

    fh = open(fileName, "r+")
    lineBuf = fh.readlines(chunkSizeBytes)
    fh.close()

    nBuf = len(lineBuf)
    rowCountEst = int(nBuf * float(fileSize) / float(chunkSizeBytes))

    print(fileSize, chunkSizeBytes, nBuf, rowCountEst)

    chunkSizeRows = float(chunkSizeBytes) / \
        (float(fileSize)/float(rowCountEst))
    if (chunkSizeRows < 0.90*rowCountEst):
        if (chunkSizeRows > 10000000):
            rndN = 1000000
        elif (chunkSizeRows > 1000000):
            rndN = 100000
        elif (chunkSizeRows > 100000):
            rndN = 10000
        else:
            rndN = 1000
        chunkSizeRows = rndN * int(1 + chunkSizeRows/rndN)
        numChunks = rowCountEst//chunkSizeRows + 1
    else:
        chunkSizeRows = -1
        numChunks = 1

    print("    file size = %d MB" % (fileSize//1000000))
    print("    row count estimate = %d" % (rowCountEst))
    if (chunkSizeRows > 0):
        print("    recommended chunk size = %d" % (chunkSizeRows))
        print("    estimated number of chunks = %d" % (numChunks))
    else:
        print("    file can be read in one chunk")
        chunkSizeRows = nBuf + 1000

    return (fileSize, rowCountEst, chunkSizeRows, numChunks)

# ----------------------------------------------------------------------------------------------------
# BigQuery data types = ['string','bytes','integer','float','boolean','record','timestamp']


def getFieldType(fConfig):

    ## print (" ... in getFieldType ... ", fConfig)

    if ('bqType' in fConfig):
        return (fConfig['bqType'])

    for aKey in fConfig:
        if (aKey.lower().find("type") > 0):
            if (fConfig[aKey] == 'boolean'):
                return ('boolean')
            if (fConfig[aKey] == 'chromosome'):
                return ('string')
            if (fConfig[aKey] == 'date'):
                return ('string')
            if (fConfig[aKey] == 'datetime'):
                return ('timestamp')
            if (fConfig[aKey] == 'float'):
                return ('float')
            if (fConfig[aKey] == 'genomic_locus'):
                return ('string')
            if (fConfig[aKey] == 'integer'):
                return ('integer')
            if (fConfig[aKey] == 'string'):
                return ('string')
            print("UHOH in getFieldType: what should I do with this ??? ",
                  aKey, fConfig[aKey])

    return ('string')

# ----------------------------------------------------------------------------------------------------
# optional descriptive text for this field that can go into the BigQuery JSON schema file


def getFieldDescription(fConfig):
    if ('description' in fConfig):
        return (fConfig['description'])
    else:
        return ('')

# ----------------------------------------------------------------------------------------------------
# TODO: get rid of this ???


def isHexBlob(s):

    logging.debug("in isHexBlob ... <%s>" % s[:128])

    if s is None:
        return (False)
    if not isinstance(s, str):
        return (False)
    if len(s) == 0:
        return (False)
    if len(s) < 128:
        return (False)
    for c in s:
        if c not in string.hexdigits:
            return (False)
    return (True)

# ----------------------------------------------------------------------------------------------------


def stripBlanks(s):

    if (str(s) == 'nan'): return ('')
    if (s == ''): return ('')
    if (not isinstance(s, str)): return (s)

    return (s.strip().strip(u'\u200b'))

# ----------------------------------------------------------------------------------------------------
# strip the specified prefix (if present) from the input string ...


def stripPrefixFromStr(s, p):

    ## print ("             in stripPrefixFromStr ... ")
    if (str(s) == 'nan'):
        return ('')

    if (not isinstance(s, str)):
        print(" UHOH ??? in stripPrefixFromStr but don't have a string ??? ", s, p)

    if (s.startswith(p)):
        np = len(p)
        t = s[len(p):]
        return (t)

    return (s)

# ----------------------------------------------------------------------------------------------------
# strip the specified suffix (if present) from the input string ...


def stripSuffixFromStr(s, p):

    ## print ("             in stripSuffixFromStr ... ")
    if (str(s) == 'nan'):
        return ('')

    if (s.endswith(p)):
        np = len(p)
        t = s[:-len(p)]
        return (t)

    return (s)

# ----------------------------------------------------------------------------------------------------

def hackDate(s):

    if ( s == "00/00/00" ): return ( s )
    if ( s == "00/00/0000" ): return ( s )

    today = datetime.date.today()
    nowYear = int ( today.timetuple()[0] )

    i1 = s.find('/')
    if ( i1 > 0 ):
        i2 = s.find('/',i1+1)
        ## print ( s, i1, i2 )
        if ( i2 > i1 ):
            iYear = int ( s[i2+1:] )
            if ( iYear == 0 ):
                print ( " UHOH year 0 ??? ", s )
                return ( s )
            elif ( iYear < 100 ):
                iYear += 2000
                if ( iYear > nowYear ): iYear -= 100
                ns = s[:i2+1] + str(iYear)
                ## print ( " --> returning <%s> " % ns )
                return ( ns )
            elif ( iYear >= 1850 ):
                return ( s )

    print ( " UHOH ... did something go wrong in hackDate ??? <%s> " % s )
    print ( i1, i2 )
    sys.exit(-1)

# ----------------------------------------------------------------------------------------------------
# handle the 'standard' types ...


def handleStandardType(s, p):

    if (s == ''):
        return (s)
    if (str(s) == 'nan'):
        return ('')

    standardTypeList = ["boolean", "datetime",
                        "date", "integer", "float", "string"]

    ## print ("             in handleStandardType ... ", s, p)
    if (p not in standardTypeList):
        logging.warning('invalid standard type ' + p)
        print(" UHOH -- invalid standard type ")
        return (s)

    if (p == "boolean"):
        try:
            ## print (" input: <%s> " % s.strip().lower() )
            t = str(bool(util.strtobool(s.strip().lower())))
            ## print (" got back ", t )
            return (t)
        except:
            print (" UHOH -- FAILED TO interpret/cast to boolean ??? <%s> " % s.strip().lower() )
            sys.exit(-1)

    elif (p == "datetime"):
        # input: 2018-12-18T15:41:29.554Z
        # Python ISO format:  2018-12-18T15:41:29.554000+00:00
        # 2018-12-18T15:41:29+00:00
        try:
            t = dateutil.parser.parse(s.strip())
            ## print (" input: <%s> " % s.strip() )
            ## print (" ISO format: ", t.isoformat() )
            u = str(t.isoformat())[:23]

            # the ISO format apparently is a little bit different if
            # the fractional part of the seconds is zero
            if (u.endswith('+00:')):
                u = u[:-4]

            ## print (" --> returning DATETIME <%s> " % u )
            return (u)
        except:
            print(" UHOH FAILED to interpret as datetime ??? ", s)

    elif (p == "date"):
        ## KNOWN BUG: if the input date is like 06/18/15
        ## there isn't any real way to know whether this should be 2015 or 1915 ...
        if ( s.find('/') > 0 ): s = hackDate(s)
        try:
            t = dateutil.parser.parse(s.strip())
            ## print (" input: <%s> " % s.strip() )
            ## print (" ISO format: ", t.isoformat() )
            # keep just the first 10 characters: YYYY-MM-DD
            u = str(t.isoformat())[:10]
            print (" --> returning DATE <%s> (derived from input string %s) " % ( u, s ) )
            return (u)
        except:
            if (s == "00/00/0000"): return ('')
            print (" UHOH FAILED to interpret as date ??? ", s )

    elif (p == "integer"):
        try:
            t = int(s.strip())
            u = str(t)
            return (u)
        except:
            print(" UHOH FAILED to interpret as integer ??? ", s)

    elif (p == "float"):
        try:
            t = float(s.strip())
            u = str(t)
            return (u)
        except:
            print(" UHOH FAILED to interpret as float ??? ", s)

    elif (p == "string"):
        # nothing really needs to be done, this 'type' is just being
        # included for completeness...
        return (s)

    else:
        print(" UHOH  TODO: implement handling for additional pyTypes ... ", p)
        sys.exit(-1)

# ----------------------------------------------------------------------------------------------------

def handleCustomType(s, p):

    chrList = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
               '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
               '21', '22', 'X', 'Y', 'M']

    if (str(s) == 'nan'):
        return ('')
    if (s == ''):
        return (s)

    customTypeList = ["genomic_locus", "chromosome", ]

    ## print ("             in handleCustomType ... ", s, p)

    if (p not in customTypeList):
        logging.warning('invalid custom type ' + p)
        print(" UHOH ... invalid custom type ", p)
        return (s)

    if (p == "genomic_locus"):
        # expecting strings like chrX:1234567 or chr1:98765432-98765437

        t = s
        if (str(t) == 'nan'):
            return ('')

        if (not t.startswith('chr')):
            t = 'chr' + t
        u = t.split(':')
        while ('' in u):
            u.remove('')
        # print(u)

        v = u[0][3:]
        if (v not in chrList):
            print(" INVALID chromosome ??? ", u)
            print(" --> re-setting entire string to blank ", s, t, u)
            return ('')

        # is this something simple like chrX:1234567 ?
        if (len(u) == 2):
            try:
                ipos = int(u[1])
                # good to go!
                return (t)

            # if it's not, then is it something like chr17:112233-123456 ?
            except:
                if (u[1].find('-') > 0):
                    v = u[1].split('-')
                    if (len(v) == 2):
                        try:
                            ipos = int(v[0])
                            jpos = int(v[1])
                            # good to go!
                            return (t)
                        except:
                            print(" (c) UHOH more complicated ??? ", s, t, u)
                print(" (a) UHOH more complicated ??? ", s, t, u)
                print("     --> returning (%s) " % t)
                # oh well, just return what we have ...
                return (t)

        else:

            # if thre is no ':' then what exactly do we have ???
            # it's also possible that we have just something like 'chr17:' ...
            if (len(u) == 1):
                t = u[0][3:]
                if t in chrList:
                    # if its just the chromosome, it's not really a
                    # valid 'locus' but we'll return it anyway ...
                    return (u[0])

                # otherwise, let's dump it
                print(" UHOH BAD genomic locus ??? !!! <%s> " % t)
                print(" --> re-setting to blank")
                return ('')

            else:

                # one common error is that the string looks like
                # this:  chr3:chr3:987654
                # in which case we want to remove the leading 'chr3:'
                if (u[1].startswith('chr')):
                    if (u[0] == u[1]):
                        u.remove(u[0])
                        if (len(u) == 2):
                            try:
                                ipos = int(u[1])
                                return (u[0] + ':' + u[1])
                            except:
                                print(" (c) UHOH new error type ??? ", t)
                    else:
                        print(" UHOH BAD genomic locus ??? !!! <%s> " % t)
                print(" (b) UHOH more complicated ??? ", t)
                print("     --> returning (%s) " % t)
                # oh well, just return what we have ...
                return (t)

        print(" UHOH ... CRASHING in handleCustomType ... ", s, t)
        sys.exit(-1)

    if (p == "chromosome"):
        if (not s.startswith('chr')):
            s = 'chr' + s
        t = s[3:]
        if t not in chrList:
            print(" UHOH INVALID chromosome ??? ", s, t)
            print(" --> re-setting to blank")
            return ('')
        else:
            # good to go!
            return (s)

    print(" UHOH TODO: write more code in handleCustomType !!! ", s)
    return (s)

# ----------------------------------------------------------------------------------------------------


def reMatchTest(s, r):

    if (s == ''):
        return (s)
    if (str(s) == 'nan'):
        return ('')
    ## print ("             in reMatchTest ... ", s, r)

    t = re.match(r, s)
    if (t):
        ## print ("                 --> TRUE ")
        return (s)
    else:
        ## print ("                 --> FALSE ")
        return ("RE-MATCH-FAIL|"+s)

# ----------------------------------------------------------------------------------------------------


def enforceAllowed(s, aList, mDict):

    if (str(s) == 'nan'):
        return ('')
    if (s == ''):
        return ('')
    ## print (" in enforceAllowed ... ", aList, mDict)

    for m in mDict:
        b = m.lower()
        t = s.lower()
        if (t == b):
            ## print ("     --> applying mapping", s, m, mDict[m])
            s = mDict[m]

    for a in aList:
        t = s.lower()
        b = a.lower()
        if (t == b):
            ## print ("     --> found allowed match ", s, a)
            return (s)

    print(" UHOH failed to match to allowed strings !!! ", s, aList)
    sys.exit(-1)

##----------------------------------------------------------------------------------------------------
## this function takes a string that may have line-delimiters and returns it as a list
## of strings split by line ...

def string2list(s):

    if (str(s) == 'nan'):
        return ('')
    if (s == ''):
        return ('')

    print (" in string2list ... ")
    print (" >>> %s <<< " % s )

    ## total hack to fix a really messed up input string
    if ( s == "MET c.3029C>T: p.T1010IMET c.3029C>T: p.T1010I" ):
        s = "MET c.3029C>T: p.T1010I"

    sList = []
    for u in s.splitlines():
        u = stripBlanks(u)
        if (len(u) > 0):
            sList += [u]

    if (len(sList) < 1):
        return ('')

    print ("     --> returning from string2list : ", len(sList), str(sList) )
    return ( str(sList) )

##----------------------------------------------------------------------------------------------------

def uniqueStringsFromX ( x ):

    if ( x is None ): return ( x )
    try:
        if ( len(x) == 0 ): return ( x )
    except:
        return ( x )

    ## print ( " XXXX ... ", x )

    if ( isinstance(x,list) ):
        nx = []
        for a in x:
            if len(a) == 0: continue
            ## print ( type(a), a )
            if ( isinstance(a,tuple) or isinstance(a,list) ):
                for b in a:
                    if len(b) == 0: continue
                    if b not in nx: nx += [ b ]
            else:
                if a not in nx: nx += [ a ]
        ## print ( " --> returning : ", nx )
        return ( nx )
    else:
        print ( " XXXX ... need to handle something that is not a list ??? " )
        print ( type(x) )
        print ( x )
        sys.exit(-1)

    return ( x )

##----------------------------------------------------------------------------------------------------

def handleKnownTerms ( s ):

    knownTerms = { 'amplification':-1, 'splice site':-1, 'promoter':-1, 'fusion':-1, \
                   'no clinvar entry for this variant':-1, \
                   'no clinvar entry':-1, 'clinvar reports as':-1, \
                   'no database entries for this variant':-1, \
                   'not expected in this sample':-1, \
                   'expected in this sample':-1, \
                   'low level variant with no cosmic or dbsnp entries':-1, \
                   'see previous analysis':-1, \
                   '(clinically relevant)':-1, \
                   'clinically relevant':-1, \
                   'germline':-1, 'exon':-1, \
                   'likely benign':-1, \
                   'likely_pathogenic':-1, 'likely pathogenic':-1, \
                   'pathogenic':-1, \
                   'unknown_significance':-1, 'unknown significance':-1, \
                   'tumor mutational burden':-1, 'TMB':-1, \
                   'microsatellite instability':-1, 'microsatellite':-1, 'MSI':-1 }

    lowerTerms = [x.lower() for x in list ( knownTerms.keys() )]

    ## first, we just need to "find" all the known terms ... 
    ## note that we are assuming that each known terms only appears ONCE
    ## 12nov2019: we need to find and remove the string(s), so that any of their
    ##            substrings don't also wind up getting 'found' ...

    foundTerms = []
    nn = 0
    t = s
    for k in knownTerms:
        knownTerms[k] = t.lower().find(k.lower())
        if ( knownTerms[k] >= 0 ):
            foundTerms += [k]
            tt = t[:knownTerms[k]] + t[knownTerms[k]+len(k):]
            t = tt
            nn += 1
    print ( " --> found %d terms ... " % nn, foundTerms )

    nn = 0
    ix = []
    for k in foundTerms:
        knownTerms[k] = s.lower().find(k.lower())
        if ( knownTerms[k] >= 0 ):
            print ( " term %s FOUND at %d ... " % ( k, knownTerms[k] ), s )
            nn += 1
            ix += [ knownTerms[k], knownTerms[k]+len(k) ]

    ## if we did NOT find any of these terms, we're done!
    if ( nn == 0 ): return ( [s] )

    ## otherwise, we need to do some more customized string-handling ...
    print ( " %d terms found in one alteration string " % nn )

    ## now let's figure out where we have to split this string ...
    ix = list(set(ix))
    ix.sort()
    print ( " split points : ", len(ix), ix )

    logString = "substrings after splitting : "
    ss = [''] * (len(ix)+1)
    print ( ss )
    for js in range(len(ix)+1):
        print ( js, ix, s, ss )
        if ( js == 0 ):
            ss[js] = stripBlanks ( s[:ix[js]] )
        elif ( js == len(ix) ):
            ss[js] = stripBlanks ( s[ix[js-1]:] )
        else:
            ss[js] = stripBlanks ( s[ix[js-1]:ix[js]] )
        logString += " <%s> " % ss[js]

    print ( logString )
    
    ## make sure that we don't have any blank substrings ...
    ## also strip off special characters at the start or end of substrings
    ## and split strings at ':' and ' ' unless they match a 'known term'
    ## ( yeah, I know this is hideous )

    us = []
    for u in ss:
        if ( len(u) == 0 ): continue
        if ( u[0] in [':','-'] ):
            u = stripBlanks ( u[1:] )
        if ( len(u) == 0 ): continue
        if ( u[-1] in [':','-'] ):
            u = stripBlanks ( u[:-1] )
        if ( len(u) == 0 ): continue
        v = u.split(':')
        for w in v:
            if ( w.lower() not in knownTerms and w[0]!='(' and w[-1]!=')' ):
                z = w.split(' ')
                for y in z:
                    if ( len(y) == 0 ): continue
                    if ( y not in us ): us += [ y ]
            else:
                if ( len(w) == 0 ): continue
                if ( w not in us ): us += [ w ]
    print ( " --> unique substrings : ", us )

    ## next up, there are 3 known terms that require special handling:
    ##     exon
    ##     microsatellite
    ##     tumor mutational burden

    for js in range(len(us)):
        u = us[js]
        if ( u.lower() == 'exon' ):
            print ( js, u )
            ns = us[js] + '=' + us[js+1]
            print ( ns )
            us.remove(us[js])
            us.remove(us[js])
            print ( us )
            us = us[:js] + [ns] + us[js:]
            print ( " --> after handling exon : ", us )
            break

    for js in range(len(us)):
        u = us[js]
        if ( u.lower() == 'microsatellite' ):
            if ( us[js+1].lower() == 'instability' ): us.remove(us[js+1])
            ns = 'MSI=' + us[js+1]
            us.remove(us[js])
            us.remove(us[js])
            us = us[:js] + [ns] + us[js:]
            print ( " --> after handling MSI : ", us )
            break

    for js in range(len(us)):
        u = us[js]
        if ( u.lower() == 'msi' ):
            if ( us[js+1].lower() == 'instability' ): us.remove(us[js+1])
            ns = 'MSI=' + us[js+1]
            us.remove(us[js])
            us.remove(us[js])
            us = us[:js] + [ns] + us[js:]
            print ( " --> after handling MSI : ", us )
            break

    for js in range(len(us)):
        u = us[js]
        if ( u.lower() == 'tumor mutational burden' ):
            ns = 'TMB=' + us[js+1]
            us.remove(us[js])
            us.remove(us[js])
            us = us[:js] + [ns] + us[js:]
            print ( " --> after handling TMB : ", us )
            break

    for js in range(len(us)):
        u = us[js]
        if ( u.lower() == 'tmb' ):
            ns = 'TMB=' + us[js+1]
            us.remove(us[js])
            us.remove(us[js])
            us = us[:js] + [ns] + us[js:]
            print ( " --> after handling TMB : ", us )
            break

    ## we also need to eliminate the blank from strings like 'splice site':

    for js in range(len(us)):
        u = us[js]
        if ( u.lower() == 'splice site' ):
            ns = 'splice_site'
            us = us[:js] + [ns] + us[js+1:]
            print ( " --> after handling splice site : ", us )
            break

    for js in range(len(us)):
        u = us[js]
        if ( u.lower() == 'likely pathogenic' ):
            ns = 'likely_pathogenic'
            us = us[:js] + [ns] + us[js+1:]
            print ( " --> after handling likely pathogenic : ", us )
            break

    for js in range(len(us)):
        u = us[js]
        if ( u.lower() == 'unknown significance' ):
            ns = 'unknown_significance'
            us = us[:js] + [ns] + us[js+1:]
            print ( " --> after handling unknown significance : ", us )
            break

    return ( us )

##----------------------------------------------------------------------------------------------------

def noCommonWords ( s ):

    cw = ['the', 'but', 'not', 'was', 'is', 'be', 'may', 'in', 'it', 'be', 'at', \
          'of', 'from', 'with', 'only', 'like' ]

    print ( " checking for common words ... <%s> " % s )

    t = s.lower()
    for w in cw:
        iw = t.find(w)
        if ( iw == 0 ): 
            try:
                if ( t[len(w)]==' ' ): 
                    print ( "     FOUND CW ", iw, w )
                    return ( False )
            except:
                pass
        if ( iw > 0 ):
            try:
                if ( t[iw-1]==' ' ):
                    if ( t[iw+len(w)]==' ' ):
                        print ( "     FOUND CW ", iw, w )
                        return ( False )
            except:
                pass

    print ( "     NO CW FOUND " )
    return ( True )

##----------------------------------------------------------------------------------------------------

def lookForHGVSabbrevs ( s ):

    if ( s == '' ): return ( False )

    abbrevs = ['c.', 'g.', 'r.', 'p.', 'm.', 'n.']

    t = s
    ## but not inside parenthetical statements ...
    while ( t.find ('(') >= 0 ):
        i1 = t.find('(')
        i2 = t.find(')', i1)
        if ( i2 > i1 ):
            t = t[:i1] + t[i2+1:]
        else:
            t = t[:i1]

    if ( t == '' ): return ( False )

    ## also remove a final '.'
    if ( t[-1] == '.' ): t = t[:-1]

    for a in abbrevs:
        skipFlag = False
        ia = t.find(a)
        if ( ia >= 0 ):
            ## is there another lower-case letter preceding this?
            ja = ia - 1
            if ( ja >= 0 ):
                ka = ord(t[ja])
                print ( ja, t, t[ja], ka )
                if ( ka >= 97 and ka <= 122 ):
                    skipFlag = True

            if ( not skipFlag ): return ( True )

    return ( False )

##----------------------------------------------------------------------------------------------------

def removeCertainBlanks ( u ):

    if ( u == '' ): return ( '' )

    d = [' =', '= ', ' :', ': ']
    t = u
    for e in d:
        print ( " e = <%s> " % e )
        while ( t.find(e) >= 0 ):
            ii = t.find(e)
            print ( "     ii = %d " % ii )
            try:
                if ( t[ii] == ' ' ): t = t[:ii] + t[ii+1:]
            except:
                continue
            try:
                if ( t[ii+1] == ' ' ): t = t[:ii+1] + t[ii+2:]
            except:
                continue

    if ( t == '' ): return ( '' )
    if ( t[0] in ['=',':'] ):
        t = t[1:]

    if ( t == '' ): return ( '' )
    if ( t[-1] in ['=',':'] ):
        t = t[:-1]

    t = stripBlanks ( t )
    return ( t )

##----------------------------------------------------------------------------------------------------
## example:
##     we start with:                 ASXL1 c.2083C>T: p.Q695X
##     after "first pass" we have:  ['ASXL1 c.2083C>T: p.Q695X']
##     after "second pass" we have: ['ASXL1', 'c.2083C>T', 'p.Q695X']

def handleFreeAlterationText ( u ):

    print ( " ... in handleFreeAlterationText ... ", type(u), len(u), u )

    u = removeCertainBlanks ( u )

    ## there are a few key words and phrases that we need to look for ...
    ## but maybe only when there is something like a coding- or protein-variant ???
    if ( lookForHGVSabbrevs(u) and noCommonWords(u) ):
        v = handleKnownTerms ( u )
    else:
        print ( " skipping the handleKnownTerms part ... " )
        v = [ u ]
    print ( " ==> first pass produces : ", v )

    ## is there anything else we need to do ??? probably is ...
    vv = []
    for w in v:
        if ( lookForHGVSabbrevs(w) and noCommonWords(w) ):
            z = uniqueStringsFromX ( re.split('\s|:', w) )
        else:
            z = [ w ]
        vv += z

    print ( " ==> second pass produces : ", vv )

    ## evidently sometimes we get here and the 'c.' and/or the 'p,'
    ## have been separated from the rest of the HGVS description ...
    if ( 'c.' in vv ):
        if ( 1 ):
            ic = vv.index('c.')
            print ( "     found 'c.' at %d " % ic, vv )
            pref = vv[:ic]
            try:
                suff = vv[ic+2:]
            except:
                suff = []
            print ( "     ", pref, suff )
            try:
                if ( vv[ic+1][0] in [ '*', 'A', 'C', 'G', 'T', '-', '1', '2', '3', '4', '5', '6', '7', '8', '9' ] ):
                    ns = vv[ic] + vv[ic+1]
                else:
                    ns = vv[ic+1]
            except:
                ns = ''
            print ( "     ns = <%s> " % ns )
            vv = pref + [ns] + suff
            print ( " FIXED : ", vv )
        ## except:
        ##     print ( " FAILED to fix a problem ??? ", vv )

    if ( 'p.' in vv ):
        if ( 1 ):
            ic = vv.index('p.')
            print ( "     found 'p.' at %d " % ic, vv )
            if ( ic == len(vv)-1 ):
                print ( "          but it's at the end, so ditching ... " )
                vv = vv[:-1]
            else:
                pref = vv[:ic]
                try:
                    suff = vv[ic+2:]
                except:
                    suff = []
                print ( "     ", pref, suff )
                try:
                    if ( ord(vv[ic+1][0])>=65 and ord(vv[ic+1][0])<=90 ):
                        ns = vv[ic] + vv[ic+1]
                    else:
                        ns = vv[ic+1]
                except:
                    ns = ''
                print ( "     ns = <%s> " % ns )
                vv = pref + [ns] + suff
                print ( " FIXED : ", vv )
        ## except:
        ##     print ( " FAILED to fix a problem ??? ", vv )

    return ( vv )

##----------------------------------------------------------------------------------------------------

def findLongestInteger(s):

    if ( s[0] == '*' ):
        print ( " in findLongestInteger ... ignoring leading * ... " )
        s = s[1:]

    ii = len(s)
    done = False
    while not done:
        try:
            ipos = int(s[:ii])
            done = True
            return ( ipos )
        except:
            ii -= 1
            ## print ( " in findLongestInteger ... nothing worked !!! ", s )
            if ( ii == 0 ): return ( None )

##----------------------------------------------------------------------------------------------------
## one problematic example seems to be c.GG34_35TT

def handleRefPosVar ( s ):

    ## not really sure if * should be considered a valid nucleotide ... ???
    validN = ['A', 'C', 'G', 'T', 'N', '*']

    try:
        j1 = 0
        while ( s[j1] in validN ): j1 += 1
        j2 = len(s) - 1
        while ( s[j2] in validN ): j2 -= 1

        print ( s, j1, j2 )
        ref = s[:j1]
        pstr = s[j1:j2+1]
        var = s[j2+1:]
        print ( " in handlRefPosVar ... ", ref, pstr, var )

    except:
        return ( s )

    try:
        ipos = int ( pstr )
        print ( "     --> integer position : ", ipos )
    except:
        print ( " UHOH ... first attempt at getting integer position failed ??? ", pstr )
        try:
            ipos = int ( pstr[:-1] )
            print ( "     --> second attempt: ", ipos )
        except:
            print ( " UHOH ... still not working in handlRefPosVar ??? ", pstr )
            if ( pstr.find('_') > 0 ):
                print ( "     ... AAAAAH ... need to deal with underscore too ... " )
            sys.exit(-1)

    return ( 'c.' + str(ipos) + ref + '>' + var )

##----------------------------------------------------------------------------------------------------

def validateCodingAlt(w):

    ## 12nov2019: this function was bombing on c.*2336_*2339del  ... FIXED

    ## not really sure if * should be considered a valid nucleotide ... ???
    validN = ['A', 'C', 'G', 'T', 'N', '*']
    validD = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '*']

    origW = w

    if ( w.startswith('c.') ):
        ## validate the rest of this string ...
        print ( " validating coding alteration ... <%s> " % w )

        ## very hacky, but seeing this kind of string:
        ##         c.[146_153del; c.172_175del]:
        ## getting split into:
        ##         c.[146_153del;
        ##     and c.172_175del]:
        ## or similar :(
        yucky = ['[', ']', ';', ':' ]
        for y in yucky:
            iy = w.find(y)
            if ( iy == 0 ):
                w = w[1:]
            elif ( iy > 0 ):
                w = w[:iy] + w[iy+1:]
        if ( w != origW ): 
            print ( "     un-yucked ... <%s> to <%s> " % ( origW, w ) )
            ## sys.exit(-1)

        ## 't' substring starts after the 'c.' prefix
        t = w[2:]

        ## we are expecting something like c.694C>T but sometimes get instead c.C694T
        ## or even c.GG1044>TT 
        ## or c.*2321A>G
        if ( t[0] not in validD ):
            print ( " NBNB!!! non-standard coding alteration string ??? <%s> " % w )
            nw = handleRefPosVar ( t )
            print ( " --> got back %s from handleRefPosVar " % nw )
            return ( nw )

        hasUnderscore = False
        if ( w.find('_') >= 0 ):
            ## one problematic example seems to be c.GG34_35TT
            hasUnderscore = True
            print ( " NBNB!!! coding alteration string has UNDERSCORE !!! %s " % w )

            iu = w.find('_')
            wp1 = w[2:iu]
            wp2 = w[iu+1:]
            print ( "     two parts : %s and %s " % ( wp1, wp2 ) )
            ipos1 = findLongestInteger(wp1)
            ipos2 = findLongestInteger(wp2)
            if ( ipos1 is not None and ipos2 is not None ):
                return ( w )
            else:
                print ( " UHOH ??? still some kind of problem ??? %s " % w )
            sys.exit(-1)

        ipos = findLongestInteger(t)
        print ( "     --> got back longest integer %d " % ipos )
        if ( ipos is not None ): 
            if ( ipos < 0 ): print ( " note ... in validateCodingAlt ... negative position found ... ", ipos, w )
            if ( ipos == 0 ): 
                print ( " WHAT ? there is no position 0 !!! ??? ", ipos, w )
                sys.exit(-1)
            return (w)

        if ( 1 ):
            j1 = 0
            while ( t[j1] in validN ): j1 += 1
            j1 -= 1

            j2 = len(t)-1
            while ( t[j2] in validN ): j2 -= 1
            j2 += 1

            if ( t[j1] in validN and t[j2] in validN ):
                ipos = findLongestInteger(t[j1+1:])
                if ( ipos <= 0 ):
                    print ( " INVALID coding alteration string ??? !!! ", w )
                    sys.exit(-1)
                else:
                    nw = 'c.' + str(ipos) + t[:j1+1] + '>' + t[j2:]
                    print ( " FIXED: %s --> %s " % ( w, nw ) )
                    return(nw)
    else:
        print ( " INVALID coding alteration string ??? ", w )
        sys.exit(-1)


##----------------------------------------------------------------------------------------------------
##    alanine       - Ala - A 
##    cysteine      - Cys - C 
##    aspartic acid - Asp - D 
##    glutamic acid - Glu - E 
##    phenylalanine - Phe - F 
##    glycine       - Gly - G 
##    histidine     - His - H 
##    isoleucine    - Ile - I 
##    lysine        - Lys - K 
##    leucine       - Leu - L 
##    methionine    - Met - M 
##    asparagine    - Asn - N 
##    proline       - Pro - P 
##    glutamine     - Gln - Q 
##    arginine      - Arg - R 
##    serine        - Ser - S 
##    threonine     - Thr - T 
##    valine        - Val - V 
##    tryptophan    - Trp - W 
##    tyrosine      - Tyr - Y 
##   (termination)  - Ter - *

standard_aa_names = [ "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", 
                      "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", 
                      "Trp", "Xaa", "Tyr", "Glx", "Ter" ] 

aa1 = "ACDEFGHIKLMNPQRSTVWXYZ*" 
aa3 = standard_aa_names 
aa3_upper = [x.upper() for x in standard_aa_names]

special_terms = [ "del", "dup", "fsTer", "fs", "ins", "delins" ]

def threeLetterAAs(s):

    ## only attempt to edit this string if it starts with p.
    if ( not s.startswith('p.') ): return(s)

    os = s
    if ( s.startswith('p.(') and s.endswith(')') ): 
        inParens = True
        s = s[3:-1]
    else:
        inParens = False
        s = s[2:]

    print ( " 3AA : in threeLetterAAs ... <%s> " % os )

    ii = 0
    done = False
    while not done:

        ## first lets make sure we're not starting with one 
        ## of the 'special' terms ...
        for x in special_terms:
            if ( s[ii:].startswith(x) ): 
                ii += len(x)
                continue
        if ( ii >= len(s) ): 
            done = True
            break

        ## next we grab the next 3 letters and see if they are in aa3_upper
        t = s[ii:min(ii+3,len(s))]
        ## print ( ii, done, t, s )
        if ( t.upper() in aa3_upper ): 
            ## print ( "     --> %s is a valid 3-letter abbreviation " % t )
            ii += 3

        ## if they are not, they try for a single-character match ...
        else:
            if ( s[ii] in aa1 ):
                ## print ( "     --> %s is a valid 1-letter abbreviation " % s[ii] )
                ns = standard_aa_names[aa1.find(s[ii])]
                ## print ( "     --> subbing in %s for %s " % ( ns, s[ii] ) )
                s = s[:ii] + ns + s[ii+1:]
                ## print ( "     --> got : ", s )
                ii += 3
            ## sometimes it seems like 'x' is used instead of 'X'
            ## so we'll allow it ...
            elif ( s[ii] == 'x' ):
                ns = standard_aa_names[aa1.find('X')]
                s = s[:ii] + ns + s[ii+1:]
                ii += 3
            else:
                ii += 1
        if ( ii >= len(s) ): done = True

    ## final check because I saw something weird ...
    okFlag = True
    if ( s[-3:] not in standard_aa_names ):
        if ( s[-3:] not in special_terms ):
            if ( s[-2:] not in special_terms ):
                if ( s != 'Met1?' and s != '?' ):
                    okFlag = False
                    iTer = s.find('fsTer')
                    if ( iTer >= 0 ):
                        try:
                            print ( iTer, s[iTer+5:], s )
                            jTer = int ( s[iTer+5:] )
                            okFlag = True
                        except:
                            okFlag = False

    ## one last hack-y test ...
    if ( not okFlag ):
        if ( s[-4:-1] in standard_aa_names ): s = s[:-1]
        okFlag = True

    if ( not okFlag ): print ( " UHOH THIS LOOKS WRONG ???   <%s> " % s )

    if ( inParens ):
        s = 'p.(' + s + ')'
    else:
        s = 'p.' + s

    print ( "               returning ... <%s> " % s )
    return(s)

##----------------------------------------------------------------------------------------------------

def checkHGVS(s):

    if ( str(s) == 'nan' ): return ('')
    if ( s == '' ): return ('')

    print ( "\n\n" )
    print ( "===> in checkHGVS function !!! ", type(s), s )

    ## first we need to revert from the string representation of a list
    ## to the underlying list ...
    if ( len(s) > 16384 ): 
        print ( "SVGH!   ACK ACK ACK ... string is really really LONG ... ", len(s) )
        print ( s[:80], "\n", s[-80:], "\n\n" )
        print ( " --> returning as-is ... " )
        return (s)

    newS = ''

    t = ast.literal_eval(s)
    print ( "SVGH!         literal_eval produces: ", type(t), len(t), t )
    if ( t == ['None'] ): return ( '' )

    for u in t:
        print ( "SVGH!             looping within t ... ", type(u), u )
        ## print ( "SVGH! <<<%s>>> " % u )
        v = handleFreeAlterationText ( u )
        print ( "SVGH!                 back from splitting u : ", type(v), len(v), v )

        ## now we are going to assemble a dict of information ...
        print ( "SVGH!             next we want to build up uDict ... " )
        uDict = {}

        for w in v:
            print ( "SVGH!                 looping within v ... ", type(w), w )
            if ( len(w) < 1 ): continue
            ## print ( "SVGH!         --> ", w, len(w) )
            gotGene = False
            if ( geneSymbolPattern.fullmatch(w) ):
                ## print ( "SVGH!         this looks like a gene symbol !!! " )
                try:
                    ## print ( "SVGH! calling gene_client.query !!! " )
                    r = gene_client.query ( w, species='human', size=1 )
                    if ( len(r['hits']) == 1 ):
                        h = r['hits'][0]
                        print ( "SVGH!     got : ", h )
                        ### for k in ['symbol', 'entrezgene', 'taxid', 'name' ]:
                        for k in ['symbol', 'entrezgene' ]:
                            if k in h: uDict[k] = h[k]
                        print ( "SVGH!     uDict so far : ", uDict )
                        gotGene = True
                    else:
                        print ( "SVGH! no hit ??? or multiple hits !!! " )
                        print ( r )
                except:
                    print ( "SVGH! error when calling gene_client.query ??? UHOH API request FAILED ??? !!! " )

            ## if we were successful in identifying a gene symbol, then 
            ## there is no need to dig further ...
            if ( gotGene ): continue

            ## from https://varnomen.hgvs.org/bg-material/refseq/
            if ( w.startswith('g.') ):
                uDict['genomic'] = w
            elif ( w.startswith('m.') ):
                uDict['mito'] = w
            elif ( w.startswith('c.') ):
                uDict['coding'] = validateCodingAlt(w)
            elif ( w.startswith('n.') ):
                uDict['non-coding'] = w
            elif ( w.startswith('r.') ):
                uDict['rna'] = w
            elif ( w.startswith('p.') ):
                uDict['protein'] = threeLetterAAs(w)
            else:
                ## there may be other bits of information that we still get here with ...
                if ( len(w) > 1 ):
                    if ( 'other' not in uDict ): uDict['other'] = []
                    uDict['other'] += [ w ]
                    print ( "SVGH! dumping into uDict other catch-all ... ", len(w), w )
                else:
                    print ( "SVGH: ditching this entirely ", len(w), w )
                ## print ( len(re.findall(r'\d+',w)), re.findall(r'\d+',w) )
                ## print ( len(re.findall(r'\D+',w)), re.findall(r'\D+',w) )
                ## print ( "SVGH! " )

        print ( "SVGH! built up uDict : ", uDict )
        newS += str(uDict)

    ## final string to be returned ...
    print ( " FINAL HGVS-parsed & re-built string : ", len(newS), newS )

    return ( newS )

##----------------------------------------------------------------------------------------------------

def quickStrip(s):

    if (str(s) == 'nan'):
        return ('')

    logging.debug("in quickStrip ... <%s>" % str(s)[:88])
    ## print (" in quickStrip ... <%s> " % s )

    t = ''
    for u in s.splitlines():
        ## print ("        <%s> " % u )
        if (len(u) > 0):
            if (len(t) > 0):
                t += ' ' + u
            else:
                t += u

    tabSplit = t.split('\t')
    t = ''
    for u in tabSplit:
        ## print ("        <%s> " % u )
        if (len(u) > 0):
            if (len(t) > 0):
                t += ' ' + u
            else:
                t += u

    ## print ("     returning ... <%s> " % t )
    ## if ( len(t) > 30 ): sys.exit(-1)

    return (t)

# ----------------------------------------------------------------------------------------------------


def getMostCommonStrings(ctrData, n):

    logging.debug("in getMostCommonStrings ... ")
    ## print (" in getMostCommonStrings ... ", n)

    # cTmp will be a list of tuples
    cTmp = ctrData.most_common(n+1)

    cOut = []
    for t in cTmp:
        if (t[0] != ''):
            cOut += [t]

    ## print ( type(cTmp) )
    ## print ( cTmp )
    ## print ( cOut[:n] )

    return (cOut[:n])

# ----------------------------------------------------------------------------------------------------


class DataTable:

    # Class Attributes
    pass

    # Initializer / Instance Attributes
    def __init__(self, dataFileName, configFileName):

        # initialize Logging setup
        # set filemode='w' to overwrite with new logging file each time
        # 'a' for appending
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
        logging.debug("in DataTable.loadConfig method " + configFileName)
        try:
            with open(configFileName) as config_file:
                logging.debug(" opened config file " + configFileName)
                self.config = json.load(config_file)
            print(self.config)
        except:
            logging.critical(
                'Failed to read config info from file ' + configFileName)
            self.config = {}
            sys.exit(-1)

    # get set up to read input CSV data file using an iterator
    def initDataLoad(self, dataFileName):
        logging.debug("in DataTable.initDataLoad method " + dataFileName)
        try:
            self.dataFileName = dataFileName
            (fileSize, rowCountEst, chunkSizeRows,
             numChunks) = estimateFileSize(dataFileName)

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
                                       na_values=['[]','{}','na','NA',':','::','[""]',"['']",
                                                  'EditMe','Please edit me',
                                                  'None.','None found','Not determined','Unknown.'])
            logging.info(' --> got CSV file iterator')

        except:
            logging.critical('Failed to open data file ' + dataFileName)
            print(sys.exc_info())
            sys.exit(-1)

    # initialize the data-log file (which is different from the 'logging' file)
    def initDLogFile(self):
        logging.debug("in DataTable.initDLogFile method")
        self.dlogFileName = self.dataFileName[:-4] + '.dlog'
        print(" opening log file <%s> " % self.dlogFileName)
        self.dlogFh = open(self.dlogFileName, 'w')

        logHdrString = "colName\tithChunk\tcolType\tlen(colData)\tnumNull\tfracNull\tnumNotNull\tnumUniqVals\tcommonStr\tcommonLen\tcommonN\tminLen\tmaxLen\tavgLen"
        self.dlogFh.write("%s\n" % logHdrString)

    def initOutput(self, outputFile, schemaFile):
        if (outputFile):
            self.outputFile = outputFile
        if (schemaFile):
            self.schemaFile = schemaFile
        self.fhOut = open(self.outputFile, 'w')
        self.emptyOutput = True
        ## self.outputSep = '\t'
        self.outputSep = ','

    def handleOneCol(self, x):
        logging.debug("in DataTable.handleOneCol method ... name is " +
                      x.name + " and data-type is " + str(x.dtype))
        print("in DataTable.handleOneCol method ... name is " +
              x.name + " and data-type is " + str(x.dtype))

        # get configuration parameters for this column
        try:
            xConfig = self.config[x.name]
        except:
            xConfig = {}

        print(" in handleOneCol ... ", json.dumps(xConfig, indent=4))

        # eliminate blank strings no matter what ...
        x = x.apply(stripBlanks)

        # stripPrefix
        if ("stripPrefix" in xConfig):
            p = xConfig["stripPrefix"]
            ## print ("     --> prefix string : <%s> " % p )
            ## print ("     --> calling x.apply ... " )
            x = x.apply(stripPrefixFromStr, args=(p,))

        # stripSuffix
        if ("stripSuffix" in xConfig):
            p = xConfig["stripSuffix"]
            ## print ("     --> suffix string : <%s> " % p )
            ## print ("     --> calling x.apply ... " )
            x = x.apply(stripSuffixFromStr, args=(p,))

        # reMatch
        if ("reMatch" in xConfig):
            r = xConfig["reMatch"]
            ## print ("     --> re match string : <%s> " % r )
            x = x.apply(reMatchTest, args=(r,))

        # "built-in" handling for 'standard' Python data types ...
        if ("pyType" in xConfig):
            p = xConfig["pyType"]
            ## print ("     --> python type : <%s> " % p )
            x = x.apply(handleStandardType, args=(p,))

        # "custom" handling for 'unusual' data types ...
        if ("customType" in xConfig):
            p = xConfig["customType"]
            ## print ("     --> custom type : <%s> " % p )
            x = x.apply(handleCustomType, args=(p,))

        # other special handling ...
        if ("stripString" in xConfig):
            if (xConfig['stripString'] == "True"):
                x = x.apply(quickStrip)

        if ( "string2list" in xConfig ):
            if ( xConfig['string2list'] == "True" ):
                x = x.apply ( string2list )
            if ( "HGVS" in xConfig ):
                if ( xConfig['HGVS'] == "True" ):
                    x = x.apply ( checkHGVS )

        if ("allowed" in xConfig):
            aList = xConfig['allowed']
            if ("mappings" in xConfig):
                mDict = xConfig['mappings']
            else:
                mDict = {}
            if (len(aList) > 0):
                x = x.apply(enforceAllowed, args=(aList, mDict,))

        try:
            print("         --> returning ", x[0], x[-1])
        except:
            pass

        return (x)

    def getNextChunk(self):
        logging.debug("in DataTable.getNextChunk method")
        try:
            print(" in getNextChunk ... ", self.chunkSizeRows)
            self.chunkDf = self.csvIter.get_chunk(self.chunkSizeRows)
            print(" --> back from get_chunk()")
            print(self.chunkDf.shape)
            self.ithChunk += 1
            return (True)
        except StopIteration:
            print(" --> StopIteration exception: end of iteration")
            return (False)
        except:
            print(" UHOH other error ???")
            return (False)

    def processChunk(self):
        # this method processes one chunk of the input dataframe ...

        # before we do anything else, we can drop any of the fields
        # that have been flagged as to-be-dropped in the input config
        colNames = self.chunkDf.columns.array
        dropList = []
        for aName in colNames:
            if (aName in self.config):
                if ("dropField" in self.config[aName]):
                    if (self.config[aName]["dropField"].lower() == "true"):
                        dropList += [aName]
        if (len(dropList) > 0):
            print(" DROPPING one or more columns ... ", dropList)
            print(self.chunkDf.shape)
            self.chunkDf = self.chunkDf.drop(columns=dropList)
            print(self.chunkDf.shape)
            colNames = self.chunkDf.columns.array
            print(" ")

        # first, we apply the handleOneCol to each column in the dataframe
        trim1 = self.chunkDf.apply(self.handleOneCol)
        self.chunkDf = trim1

        # and now we can git a bit further ...
        colTypes = self.chunkDf.dtypes
        numCol = len(colNames)
        numRows = len(self.chunkDf)

        for iCol in range(max(numCol, 10)):

            # get the name and the type for this column ...
            colName = colNames[iCol]
            colType = colTypes[iCol]

            print(
                "\n\n ------------------------------------------------------------------------------------ \n\n")

            # and then get the data contained in this column
            colData = self.chunkDf[colName]

            # now let's have a look at this data ...
            # how many null values?  how many non-null values?
            numNull = colData.isnull().sum()
            fracNull = float(numNull) / float(numRows)
            numNotNull = numRows - numNull

            # use the collections.Counter to find the unique values and their counts
            ctrData = Counter(colData.dropna())
            uniqDict = dict(ctrData)
            if ('' in uniqDict):
                del uniqDict['']
            print("         uniqDict: ", str(uniqDict)[:188])
            numUniqVals = len(uniqDict)
            print("         numUniqVals: ", numUniqVals)
            commonStrings = getMostCommonStrings(ctrData, 5)

            print("====> ", colName, colType, len(colData),
                  numNull, numNotNull, numUniqVals)

            outString = colName + "\t" + \
                str(self.ithChunk) + "\t" + str(colType) + \
                "\t" + str(len(colData)) + "\t"
            outString += str(numNull) + "\t" + str(fracNull) + \
                "\t" + str(numNotNull) + "\t" + str(numUniqVals) + "\t"

            # get the single most common value
            # and write it and and its length, and its count to the outString

            if (numUniqVals > 0):
                c1 = commonStrings[0]
                try:
                    outString += str(c1[0])[:64] + "\t" + \
                        str(len(str(c1[0]))) + "\t" + str(c1[1]) + "\t"
                except:
                    outString += str(c1[0])[:64] + "\t" + \
                        "" + "\t" + str(c1[1]) + "\t"
            else:
                outString += "\t\t\t"

            # let's also figure out how long these strings tend to be ...
            minLen = 99999999
            maxLen = 0
            sumLen = 0
            numU = 0
            for u in list(uniqDict.keys()):
                uLen = len(str(u))
                sumLen += uLen
                numU += 1
                if (maxLen < uLen):
                    maxLen = uLen
                    maxStr = str(u)
                if (minLen > uLen):
                    minLen = uLen
                    minStr = str(u)

            if (numU > 0):
                avgLen = sumLen//numU
            else:
                avgLen = -1
                minLen = -1
                maxLen = -1

            outString += str(minLen) + "\t" + str(maxLen) + "\t" + str(avgLen)
            self.dlogFh.write("%s\n" % outString)

            if (maxLen > 0):
                print("LONGEST string representation is %d characters" % maxLen)
            if (maxLen > 1024):
                print("    THIS IS VERY LONG and may need to be truncated ... ")
                print("     --> first 128 chars: ", maxStr[:128])
                print("     --> last  128 chars: ", maxStr[-128:])
                print("  ")
            if (avgLen > 8192):
                print(
                    "     AVERAGE length is %d ... field should probably be ommitted ??? " % avgLen)

            if (numUniqVals == 0):
                print("ALWAYS empty or null!")
            elif (numUniqVals == numNotNull):
                print("ALL of the non-null values are UNIQUE!")
                if (avgLen < 1024):
                    print("    eg: ")
                    print("        ", commonStrings[0][0])
                    try:
                        print("        ", commonStrings[1][0])
                    except:
                        pass
                    try:
                        print("        ", commonStrings[2][0])
                    except:
                        pass

            else:

                if (numNull > (0.90*numRows)):
                    print("More than 90% of the values are NULL!")
                elif (numNull > (0.50*numRows)):
                    print("More than half of the values are NULL!")
                elif (numNull < (0.10*numRows)):
                    print("Fewer than 10% of the values are NULL!")

            if (numUniqVals == 0):
                pass

            elif (numUniqVals <= 6):
                print("unique values: ",
                      ' -;- '.join(str(x)[:80] for x in uniqDict))
            else:

                # let's get a subset of the uniqDict -- only those values that are never repeated:
                nrDict = uniqDict
                ## print ( type(nrDict) )
                nrKeys = list(nrDict.keys())
                nrKeys.sort()
                for key in nrKeys:
                    if (nrDict[key] > 1):
                        del nrDict[key]
                print("number of unique non-null values: ", len(uniqDict))
                print(
                    "number of non-null AND non-repeating values (occur only once each): ", len(nrDict))

                # pull out the single most common value ...
                (aVal, aCount) = ctrData.most_common(1)[0]
                print(aCount)
                print(type(aVal), len(str(aVal)))
                ## print ( aVal )

                if isinstance(aVal, float):
                    continue
                if isinstance(aVal, bool):
                    continue
                if isHexBlob(aVal):
                    print("     --> this looks like a Hex Blob ... ")
                    continue

                # and let's see if 'evaluating' this string produces anything
                # different ...
                try:
                    ## print (" trying out literal_eval ... " )
                    eVal1 = ast.literal_eval(aVal)
                    if (eVal1 != aVal):
                        if isinstance(eVal1, list):
                            print(
                                "     --> result of evaluation is a LIST of length %d " % len(eVal1))
                            if (len(eVal1) == 1):
                                eVal2 = ast.literal_eval(eVal1[0])

                                if (eVal2 != eVal1[0]):
                                    if isinstance(eVal2, list):
                                        print(
                                            "     --> next result of evaluation is a LIST of length %d " % len(eVal2))
                                    elif isinstance(eVal2, dict):
                                        print(
                                            "     --> result of evaluation is a DICT of length %d " % len(eVal2))

                                print("     2nd literal_eval : \n",
                                      type(eVal2), "\n", eVal2)
                                if isinstance(eVal2, dict):
                                    print(
                                        "     --> result of evaluation is a DICT of length %d " % len(eVal2))
                                    eKeys = list(eVal2.keys())
                                    nKeys = len(eKeys)
                                    print("         dict has %d keys " %
                                          nKeys, eKeys)
                                    eKeys.sort()
                                    aKey = eKeys[nKeys//2]
                                    print("         for example ... %s: %s " %
                                          (aKey, eVal2[aKey]))

                        elif isinstance(eVal1, dict):
                            print(
                                "     --> result of evaluation is a DICT of length %d " % len(eVal1))
                            eKeys = list(eVal1.keys())
                            nKeys = len(eKeys)
                            print("     dict has %d keys " % nKeys, eKeys)
                            eKeys.sort()
                            aKey = eKeys[nKeys//2]
                            print("     for example ... %s: %s " %
                                  (aKey, eVal1[aKey]))

                except:
                    ## print ("    failed to run literal_eval on this string <%s> " % str(aVal) )
                    ## print ("    trying to load as JSON instead? " )
                    if (aVal.find('{') >= 0 or aVal.find('[') >= 0):
                        try:
                            jVal1 = json.loads(str(aVal))
                            bVal = json.dumps(jVal1, indent=4)
                            if (bVal != aVal):
                                print(" YAY!!! JSON!!! ")
                                print(bVal)
                                print(" ")
                        except:
                            ## print ("         also failed to parse as JSON ... " )
                            pass

                # now the THREE most common values ...
                if (aCount > 1):
                    c3 = ctrData.most_common(3)
                    print("3 most common values: \n",
                          c3[0], "\n", c3[1], "\n", c3[2], "\n")

    def writeOutChunk(self):

        if (self.emptyOutput):
            self.chunkDf.to_csv(self.fhOut, index=False, sep=self.outputSep)
            self.emptyOutput = False
        else:
            self.chunkDf.to_csv(self.fhOut, index=False,
                                sep=self.outputSep, header=False)

    def writeJsonSchema(self):

        # BigQuery data types = ['string','bytes','integer','float','boolean','record','timestamp']

        # one row of the JSON schema file should look like this:
        ## {"name": "_id", "type": "string", "mode": "nullable", "description": "<add description here>"},

        print("in writeJsonSchema ... ", self.schemaFile)

        allCols = list(self.chunkDf.columns)

        with open(self.schemaFile, 'w') as schema_file:
            schema_file.write('[\n')
            numF = len(allCols)
            for iF in range(numF):
                fName = allCols[iF]
                print(fName, self.config[fName])
                fType = getFieldType(self.config[fName])
                fDesc = getFieldDescription(self.config[fName])
                oneLine = '    {"name": "' + fName + '", "type": "' + fType \
                    + '", "mode": "nullable", "description": "' + fDesc \
                    + '"}'
                if (iF < (numF-1)):
                    oneLine += ','
                schema_file.write(oneLine+'\n')
            schema_file.write(']')


# ----------------------------------------------------------------------------------------------------

def main(args):

    myData = DataTable(args.inputFile, args.configFile)

    if (args.outputFile):
        myData.initOutput(args.outputFile, args.schemaFile)

    done = False
    numChunks = 0
    while not done:
        r = myData.getNextChunk()
        if (r):
            numChunks += 1
            print(" ***************** ")
            print(" **  CHUNK #%3d ** " % numChunks)
            print(" ***************** ")
            myData.processChunk()
            if (args.outputFile):
                myData.writeOutChunk()
        else:
            done = True

    print(" ")
    print(" TOTAL number of chunks handled: %d " % numChunks)

    if (args.schemaFile):
        myData.writeJsonSchema()


# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    t0 = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--inputFileName', action='store', help='input CSV data file',
                        required=True, dest='inputFile', type=str)
    parser.add_argument('-c', '--configFileName', action='store', help='input JSON config file',
                        required=True, dest='configFile', type=str)
    parser.add_argument('-o', '--outputFileName', action='store', help='output CSV data file',
                        required=False, dest='outputFile', type=str)
    parser.add_argument('-s', '--schemaFileName', action='store', help='output JSON schema file',
                        required=False, dest='schemaFile', type=str)
    args = parser.parse_args()

    main(args)

# ----------------------------------------------------------------------------------------------------
