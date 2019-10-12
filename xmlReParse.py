

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

def isContainerNode(node):

    ## print ( " CCCC : in isContainerNode ... " )

    ## we'll start with an assumption that it IS just a container node ...
    cFlag = True

    if ( len(node.attrib) > 0 ):

        for a in node.attrib:

            textString = quickStrip(str(node.attrib[a]))
            if ( len(textString)==0 or textString in naStrings ): 
                skipFlag = True
            else:
                cFlag = False

    textString = quickStrip(str(node.text))
    if ( len(textString)==0 or textString in naStrings ):
        skipFlag = True
    else:
        cFlag = False

    ## print ( " CCCC : leaving isContainerNode ... ", cFlag )

    return ( cFlag )

# ----------------------------------------------------------------------------------------------------

def handleNodeData(node,pID,eID):

    ## print ( " DDDD : in handleNodeData ... ", pID, eID )

    ## keep track of whether or not we have written something out for this node ...
    outFlag = False

    if ( len(node.attrib) > 0 ):

        ## print ( " found attributes : ", ithGen, node.tag, len(node.attrib), node.attrib )
        ## need to handle the attributes:

        ## NB (!!!) attributes need to have their own 'bundle' identifier that is separate
        ## from the other 'contents' of this node...
        for a in node.attrib:

            textString = quickStrip(str(node.attrib[a]))
            if ( len(textString)==0 or textString in naStrings ): 
                skipFlag = True
                ## print ( " DDDD :     skipping attribute <%s> " % a )
            else:

                ## write out the 'attrib' contents of this particular node
                longLabel = ''
                if ( len(genList) > 0 ): longLabel = ' > '.join(genList) + ' > '
                longLabel += node.tag + ' > ' + a 
                longLabel = longLabel.strip()

                ## print ( " DDDD :     writing attribute {%s = %s} " % ( longLabel, textString ) )

                outString = inputFilename + '\t' + str(ithGen) + '\t' \
                          + str(pID) + '\t' + str(eID) + '\t' \
                          + longLabel + '\t' \
                          + node.tag + '\t' \
                          + a + '\t' + textString

                outFh.write ( "%s\n" % outString )
                outFlag = True

    textString = quickStrip(str(node.text))
    if ( len(textString)==0 or textString in naStrings ):
        ## if this node has no contents, we will still need to write 
        ## out a line with a 'null' where the text would have been
        ## so that will be able to track the bundle identifiers
        textString = ''

    ## print ( " DDDD :     textString = %s " % textString )

    ## write out the text contents of this particular node
    longLabel = ''
    if ( len(genList) > 0 ): longLabel = ' > '.join(genList) + ' > '
    longLabel += node.tag 
    longLabel = longLabel.strip()

    parentLabel = ''
    if ( len(genList) > 0 ): parentLabel = genList[-1]

    if ( (textString!='') or (not outFlag) ):
        ## print ( " DDDD :     writing text {%s = %s} " % ( longLabel, textString ) )
        outString = inputFilename + '\t' + str(ithGen) + '\t' \
                  + str(pID) + '\t' + str(eID) + '\t' \
                  + longLabel + '\t' \
                  + parentLabel + '\t' \
                  + node.tag + '\t' + textString
        outFh.write ( "%s\n" % outString )

    ## print ( " DDDD : leaving handleNodeData ... " )

# ----------------------------------------------------------------------------------------------------
# inspired by:
# https://stackoverflow.com/questions/28194703/recursive-xml-parsing-python-using-elementtree

def handleNextLevel(root,pID):

    global ithGen, genList

    ithGen += 1
    genList += [ root.tag.strip() ]

    ## print ( " RRRR : in handleNextLevel ... ", ithGen, genList )

    ## loop over all children of the incoming node:
    for elem in root:

        ## is this node just a "container" node with no information???
        if ( isContainerNode(elem) ):
            if ( ithGen <= maxDepth ):
                ## print ( " RRRR :    moving on to next level ... " )
                handleNextLevel(elem,pID)

        else:

            ## generate a new identifier for this node
            eID = uuid.uuid4()
            ## print ( " RRRR :    generating new identifier ... ", eID )

            handleNodeData(elem,pID,eID)

            ## if we haven't gone too deep yet...
            if ( ithGen <= maxDepth ):
                handleNextLevel(elem,eID)

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

    ## handle any information that exists at the root level ...
    handleNodeData(root,None,rootID)

    ## and then go to the next level (this is the recursive function)
    handleNextLevel(root,rootID)

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
