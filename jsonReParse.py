

"""
    JSON parsing ... at least initially for Watson-for-Genomics reports ... 
"""

__author__ = "Sheila M Reynolds"
__version__ = "0.0.1"
__status__ = "Prototype"

import argparse
import json
import sys
import time
import uuid

# ----------------------------------------------------------------------------------------------------

inputFilename = ''
outFh = None
maxDepth = 0

ithGen = 0
genList = []

naStrings = ['[]', '{}', 'na', 'NA', ':', '::', '[""]', "['']", \
             'EditMe', 'Please edit me', \
             'None', 'None.', 'None found', 'Not determined', \
             'Unknown', 'Unknown.']

skipLabels = [ 'LegalNotice' ]

# ----------------------------------------------------------------------------------------------------

def quickStrip(s):

    if (str(s) == 'nan'):
        return ('')

    t = ''
    for u in s.splitlines():
        if (len(u) > 0):
            if (len(t) > 0):
                t += ' ' + u.strip()
            else:
                t += u.strip()

    tabSplit = t.split('\t')
    t = ''
    for u in tabSplit:
        if (len(u) > 0):
            if (len(t) > 0):
                t += ' ' + u.strip()
            else:
                t += u.strip()

    t = t.strip()

    if ( t.find(',') >= 0 ): 

        if ( t.find('"') >= 0 ): t = t.replace('"',"'")

        if ( t.find('"') < 0 ):
            t = '"' + t + '"'
        else:
            print ( " FATAL ERROR with commas and quotes in a string !!! " )
            print ( " <%s> " % t )

    return (t)

# ----------------------------------------------------------------------------------------------------

def isContainerNode(key,value):

    cFlag = True
    print ( " DDDD : in isContainerNode ... ", key, str(type(value)) )

    ## if this is a dict, we consider that a 'container'
    if isinstance(value,dict): return ( cFlag )

    ## if this is a list, we consider that a 'container'
    if isinstance(value,list): return ( cFlag )

    ## this may not matter too much, but if there isn't going to 
    ## actually be any 'data' coming from this (key,value) pair
    ## then we can also return 'True'

    try:
        textString = quickStrip(str(value))
    except:
        print ( " DDDD :         what now ??? ", type(value), value )
        print ( " ERROR !!! BAILING !!! " )
        sys.exit(-1)

    if ( len(textString)==0 or textString in naStrings ): return ( cFlag ) 

    cFlag = False
    return ( cFlag )

# ----------------------------------------------------------------------------------------------------

def isContainer2Node(key,value):

    cFlag = True
    print ( " DDDD : in isContainer2Node ... ", key, str(type(value)) )

    ## if this is a dict, check if the next layer is a container ...
    if isinstance(value,dict): 
        for k,v in value.items():
            return ( isContainerNode(k,v) )

    ## if this is a list, check if the next layer is a container ...
    if isinstance(value,list):
        for v in value:
            return ( isContainerNode(None,v) )

    cFlag = False
    return ( cFlag )

# ----------------------------------------------------------------------------------------------------

def handleNodeData(key,value,gpID,pID):

    print ( " HND : in handleNodeData ... ", gpID, pID, key, str(type(value)) )

    ## if the incoming 'thing' (aka node) is a dict, then we go the next level
    ## since in JSON, 'dicts' are what define 'levels' 
    ## (or at least that's our working interpretation)
    if isinstance(value,dict):
        dLen = len(value)
        print ( " HND :     this value is a dict of length %d ... " % dLen )

        if ( dLen > 1 ):
            eID = uuid.uuid4()
            print ( " HND :     generating a new UUID ... ", eID )
            print ( " HND :    calling handleNextLevel ... ", pID, eID )
            handleNextLevel(key,value,pID,eID)
        else:
            eID = None
            print ( " HND :     NOT generating a new UUID at this time ... " )
            print ( " HND :    calling handleNextLevel ... ", gpID, pID )
            handleNextLevel(key,value,gpID,pID)

    ## if, on the other hand, the incoming 'thing' is a list, then we loop over
    ## the things in the list and recursively call this function ...
    elif isinstance(value,list):
        vLen = len(value)
        print ( " HND :     this value is a list of length %d ... " % vLen )
        try:
            print ( " HND :        --> list of ", str(type(value[0])) )
        except:
            pass

        ## loop over contents of list:
        for a in value:

            print ( " HND : calling handleNodeData ... ", gpID, pID )
            handleNodeData(key,a,gpID,pID)

    ## and if it is neither a dict nor a list, then it must be something fairly
    ## simple and we can just handle it right now
    else:

        try:
            textString = quickStrip(str(value))
        except:
            print ( " HND :         what now ??? ", type(value), value )
            print ( " ERROR !!! BAILING !!! " )
            sys.exit(-1)

        if ( len(textString) < 1 ): 
            print ( " HND :        zero length textString ??? !!! " )
        if ( textString in naStrings ): 
            print ( " HND :        textString is in naStrings ", textString )
        if ( len(textString)==0 or textString in naStrings ): return() 

        longLabel = ''
        if ( len(genList) > 0 ): longLabel = ' > '.join(genList) 
        longLabel = longLabel.strip()

        parentLabel = ''
        elementLabel = ''
        try:
            elementLabel = genList[-1]
            parentLabel = genList[-2]
        except:
            pass

        if ( parentLabel in skipLabels or elementLabel in skipLabels ): 
            print ( " HND :         label in skipLabels ", parentLabel, elementLabel )
            return()

        print ( " HND :     writing string value {%s = %s} " % ( longLabel, textString ) )

        outString = inputFilename + '\t' + str(ithGen) + '\t' \
                  + str(gpID) + '\t' + str(pID) + '\t' \
                  + longLabel + '\t' + parentLabel + '\t' + elementLabel + '\t' \
                  + textString

        outFh.write ( "%s\n" % outString )
        print ( " HND : %s " % outString )
        print ( " \n\n" )

    print ( " HND : leaving handleNodeData ... " )
    return()

# ----------------------------------------------------------------------------------------------------
# inspired by:
# https://stackoverflow.com/questions/28194703/recursive-xml-parsing-python-using-elementtree

def handleNextLevel(key,value,gpID,pID):

    global ithGen, genList

    print ( " HNL : in handleNextLevel ... ", ithGen, genList, key, str(type(value)) )

    ## at this point we do one of two main things: if it is NOT a dict, we handle
    ## this node's 'data' ... if it IS a dict, then we loop over its items() ...

    if ( not isinstance(value,dict) ):

        ## we should never get here w/o value being of type dict !!!
        print ( " DO WE EVER ACTUALLY GET HERE ??? !!! " )
        sys.exit(-1)

    else:

        ## YES this IS a dict ...
        bkpGen = ithGen
        bkpList = genList

        ## loop over all *next* (key,value) pairs beneath the incoming 'node':
        for nkey,nvalue in value.items():
            print ( " HNL :     --> key = <%s> : value is %s " % ( nkey, str(type(nvalue)) ) )
    
            ithGen = bkpGen + 1
            genList = bkpList + [ nkey ]

            print ( " HNL : calling handleNodeData ... ", gpID, pID )
            handleNodeData(nkey,nvalue,gpID,pID)

        ## back up ...
        ithGen = bkpGen
        genList = bkpList
        
    return()

# ----------------------------------------------------------------------------------------------------

def main(args):

    global inputFilename, outFh, maxDepth

    ## print ( args )

    if ( not args.inputFile.endswith(".json") ):
        print ( " is this a JSON file ??? <%s> " % args.inputFile )
        print ( " ERROR !!! BAILING !!! " )
        sys.exit(-1)

    inputFilename = args.inputFile

    allowedModes = [ 'w', 'w+', 'a', 'a+' ]
    if ( args.writeMode not in allowedModes ):
        print ( " ERROR: invalid writeMode <%s> " % args.writeMode )
        print ( " ERROR !!! BAILING !!! " )
        sys.exit(-1)

    try:
        outFh = open ( args.outputFile, args.writeMode )
    except:
        print ( " FAILED to open output file <%s> with <%s> " % ( args.outputFile, args.writeMode ) )
        print ( " ERROR !!! BAILING !!! " )
        sys.exit(-1)

    maxDepth = args.maxDepth

    if ( args.colLabels.lower().startswith('t') ):
        outString = 'inputJSONfile\tjsonDepth\tparentID\tbundleID\tfull_label\tparent_label\telement_label\telement_text'
        outFh.write ( "%s\n" % outString )

    try:
        inFh = open ( inputFilename, 'r' )
        jdata = json.load(inFh)
        #XML xmlp = ET.XMLParser(encoding="latin-1")
        #XML tree = ET.parse(inputFilename,parser=xmlp)
    except:
        print ( " FAILED to parse input JSON file ??? <%s> " % inputFilename )
        print ( " ERROR !!! BAILING !!! " )
        sys.exit(-1)
    
    ## handle any information that exists at the root level...
    ## the first two parameters in handleNodeData are a (key,value) 
    ## but at the 'root' there is no 'key' and the 'value' is the entire dict...
    root = jdata
    print ( " CCCC : calling handleNodeData from main ... ", None, None )
    handleNodeData(None,root,None,None)

    if ( 0 ):
        ## and then go to the next level (this is the recursive function)
        print ( " CCCC : calling handleNextLevel from main ... ", None, None )
        handleNextLevel(None,root,None,None)

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    t0 = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--inputFileName', action='store', help='input JSON data file',
                        required=True, dest='inputFile', type=str)
    parser.add_argument('-o', '--outputFileName', action='store', help='output TSV data file',
                        required=True, dest='outputFile', type=str)
    parser.add_argument('-w', '--writeMode', action='store', help='output TSV data file',
                        required=False, dest='writeMode', type=str, default='w+')
    parser.add_argument('-d', '--maxDepth', action='store', help='maximum depth to output TSV',
                        required=False, dest='maxDepth', type=int, default=99)
    parser.add_argument('-l', '--colLabels', action='store', help='write column labels to output TSV',
                        required=False, dest='colLabels', type=str, default="True")
    parser.add_argument('-s', '--schemaFileName', action='store', help='output JSON schema file',
                        required=False, dest='schemaFile', type=str)
    args = parser.parse_args()

    main(args)

    print ( "\n\n\n" )

# ----------------------------------------------------------------------------------------------------
