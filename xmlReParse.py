

"""
    XML parsing ... at least initially for N-of-One reports ... 
"""

__author__ = "Sheila M Reynolds"
__version__ = "0.0.1"
__status__ = "Prototype"

import argparse
import sys
import time
import uuid

import xml.etree.ElementTree as ET

# ----------------------------------------------------------------------------------------------------

inputFilename = ''
outFh = None
maxDepth = 0

ithGen = 0
genList = []

naStrings = ['[]','{}','na','NA',':','::','[""]',"['']",'EditMe','Please edit me','None','None.','None found','Not determined','Unknown.']

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

def handleNodeData(node,gpID,pID):

    if ( len(node.attrib) > 0 ):

        ## print ( " found attributes : ", ithGen, node.tag, len(node.attrib), node.attrib )
        ## need to handle the attributes:

        ## NB (!!!) attributes need to have their own 'bundle' identifier that is separate
        ## from the other 'contents' of this node...
        aID = uuid.uuid4()
        for a in node.attrib:

            qs = quickStrip(str(node.attrib[a]))
            if ( len(qs)==0 or qs in naStrings ): 
                skipFlag = True
            else:

                longLabel = ''
                if ( len(genList) > 0 ): longLabel = ' > '.join(genList) + ' > '
                longLabel += node.tag + ' > ' + a 

                outString = inputFilename + '\t' + str(ithGen) + '\t' \
                          + str(pID) + '\t' + str(aID) + '\t' \
                          + longLabel + '\t' \
                          + node.tag + '\t' \
                          + a + '\t' + qs 

                outFh.write ( "%s\n" % outString )

    qs = quickStrip(str(node.text))
    if ( len(qs)==0 or qs in naStrings ):
        skipFlag = True
    else:

        longLabel = ''
        if ( len(genList) > 0 ): longLabel = ' > '.join(genList) + ' > '
        longLabel += node.tag 

        outString = inputFilename + '\t' + str(ithGen) + '\t' \
                  + str(gpID) + '\t' + str(pID) + '\t' \
                  + longLabel + '\t' \
                  + genList[-1] + '\t' \
                  + node.tag + '\t' + qs

        outFh.write ( "%s\n" % outString )

# ----------------------------------------------------------------------------------------------------
# inspired by:
# https://stackoverflow.com/questions/28194703/recursive-xml-parsing-python-using-elementtree

def handleNextLevel(root,gpID,pID):

    global ithGen, genList

    ithGen += 1
    genList += [ root.tag ]

    ## loop over all children of the incoming node:
    for elem in root:

        if ( ithGen <= maxDepth ):
            handleNodeData(elem,gpID,pID)

        eID = uuid.uuid4()
        handleNextLevel(elem,pID,eID)

    ## back up...
    ithGen -= 1
    genList = genList[:-1]

# ----------------------------------------------------------------------------------------------------

def main(args):

    global inputFilename, outFh, maxDepth

    ## print ( args )

    if ( not args.inputFile.endswith(".xml") ):
        print ( " is this an XML file ??? <%s> " % args.inputFile )
        sys.exit(-1)

    inputFilename = args.inputFile

    allowedModes = [ 'w', 'w+', 'a', 'a+' ]
    if ( args.writeMode not in allowedModes ):
        print ( " ERROR: invalid writeMode <%s> " % args.writeMode )
        sys.exit(-1)

    try:
        outFh = open ( args.outputFile, args.writeMode )
    except:
        print ( " FAILED to open output file <%s> with <%s> " % ( args.outputFile, args.writeMode ) )
        sys.exit(-1)

    maxDepth = args.maxDepth

    if ( args.colLabels.lower().startswith('t') ):
        outString = 'inputXMLfile\txmlDepth\tparentID\tbundleID\tfull_label\tparent_label\telement_label\telement_text'
        outFh.write ( "%s\n" % outString )

    try:
        xmlp = ET.XMLParser(encoding="latin-1")
        tree = ET.parse(inputFilename,parser=xmlp)
    except:
        print ( " FAILED to parse input XML file ??? <%s> " % inputFilename )
        sys.exit(-1)
    
    ## get the root of the tree
    root = tree.getroot()
    rootID = uuid.uuid4()

    handleNodeData(root,rootID,rootID)

    ## start by actually handling the root itself and its attributes, and then start 
    ## digging down ...
    handleNextLevel(root,rootID,rootID)

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    t0 = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--inputFileName', action='store', help='input XML data file',
                        required=True, dest='inputFile', type=str)
    parser.add_argument('-o', '--outputFileName', action='store', help='output TSV data file',
                        required=False, dest='outputFile', type=str)
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

# ----------------------------------------------------------------------------------------------------
