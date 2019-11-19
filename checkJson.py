

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

def writeJsonString2File ( jString, fName ):

    try:
        jdata = json.loads(jString)
        print ( " --> SUCCESS! parsed string as JSON " )

        ## ok, since we know this is good JSON, we want to write this out 
        print ( fName )
        outFh = open ( fName, 'w' )
        outFh.write ( "%s\n" % jString )
        outFh.close()

        print ( " --> output file written ", fName )
        print ( "\n\n" )

        return ( 1 )

    except:
        print ( " --> FAILED in writeJsonString2File ", len(jString), fName )
        return ( 0 )

# ----------------------------------------------------------------------------------------------------

def main(args):

    if ( not args.inputFile.endswith(".json") ):
        print ( " is this a JSON file ??? <%s> " % args.inputFile )
        print ( " ERROR !!! BAILING !!! " )
        sys.exit(-1)

    print ( "\n\n" )
    print ( " in checkJson ... <%s> " % args.inputFile )

    inputFilename = args.inputFile
    numOutFiles = 0
    numGarbled = 0

    ## first, we try to read and parse the entire file as a single JSON object
    try:
        inFh = open ( inputFilename, 'rb' )
        allData = inFh.read()
        try:
            allData = allData.decode('latin')
        except:
            print ( " FAILED TO DECODE ??? " )
        inFh.close()

        jdata = json.loads(allData)
        print ( " --> SUCCESS! input was parsed on first attempt with json.loads() " )

        ## ok, since we know this is good JSON, we want to write this out 
        outputFilename = inputFilename[:-5] + '.check.json'
        ## print ( outputFilename )
        outFh = open ( outputFilename, 'w' )
        outFh.write ( "%s\n" % allData )
        outFh.close()
        numOutFiles += 1
        print ( " Number of output files: ", numOutFiles, len(allData) )
        return 

    except:
        print ( " first try: FAILED to parse input JSON file ??? <%s> " % inputFilename )
        inFh.close()

    ## if that fails, then we assume that there are actually multiple concatenated
    ## JSON objects...

    ## now read the entire file into a string
    done = False
    inFh = open ( inputFilename, 'rb' )

    allData = inFh.read()
    try:
        allData = allData.decode('latin')
    except:
        print ( " FAILED TO DECODE ??? " )

    print ( " allData: ", len(allData), allData[:40], " ... ", allData[-40:] )

    openCu  = 0
    closeCu = 0
    openSq  = 0
    closeSq = 0

    tdString = '"testDetails":'

    iiStart = 0
    iiStop = 0

    for ii in range(len(allData)):

        if ( allData[ii:ii+len(tdString)] == tdString ): 
            if ( (ii-iiStart) > 3 ):
                print ( " UNEXPECTED FIND ... <%s> at ii=%d " % ( tdString, ii ) )
                print ( " delta = ", ii-iiStart )
                print ( ii, openCu, closeCu, openSq, closeSq, iiStart, iiStop )
                print ( allData[max(ii-64,0):min(ii+64,len(allData))] )
                print ( " " )
                if ( (openCu - closeCu) != 1 ):
                    print ( " --> curly brace imbalance ... hmmmmm ... " )

                    iiOld = iiStart

                    print ( " --> HACKING iiStart (and open/closeCu counts) ... " )
                    print ( " --> (re)setting iiStart from %d to %d " % ( iiStart, ii-1 ) )
                    iiStart = ii - 1
                    openCu  = 1
                    closeCu = 0

                    print ( " LOST (garbled) DATA: " )
                    print ( allData[iiOld:iiStart] )
                    print ( " " )
                    numGarbled += 1

                    ## print ( " A: ", allData[max(ii-16,0):min(ii+16,len(allData))] )
                    ## print ( " B: %s ... %s" % ( allData[iiStart:iiStart+40], allData[ii-40:ii] ) )
                    ## print ( "\n\n" )

        otherFlag = False
        if ( allData[ii] == '{' ): 
            openCu  += 1
            if ( allData[iiStart] != '{' ): 
                print ( " --> (re)setting iiStart from %d to %d " % ( iiStart, ii ) )
                print ( "        ... %s ... " % allData[ii-16:ii+17] )
                iiStart = ii
        elif ( allData[ii] == '}' ): closeCu += 1
        elif ( allData[ii] == '[' ): openSq  += 1
        elif ( allData[ii] == ']' ): closeSq += 1
        elif ( allData[ii] == '\n' or allData[ii] == '\r' ):

            print ( " --> newline character found at %d ... iiStart=%d " % ( ii, iiStart ) )

            iiStop = ii 
            tmpString = allData[iiStart:iiStop].strip()

            if ( len(tmpString) > 2 ):

                jPass = False

                print ( " len(tmpString): ", len(tmpString), iiStart, iiStop )
                print ( " --> newline character was found ... trying to process JSON ... " )
                ## print ( " C: newline character found !!! ", ii )
                ## print ( " D: ", allData[max(ii-64,0):min(ii+64,len(allData))] )
                ## print ( " E: openCu=%d  and closeCu=%d " % ( openCu, closeCu ) )
                ## print ( " F: iiStart=%d and iiStop=%d and ii=%d " % ( iiStart, iiStop, ii ) )
                ## print ( " G: NOW WHAT ??? " )
                try:
                    jdata = json.loads(tmpString)
                    print ( "         --> SUCCESS !!!! (%d-%d) " % ( iiStart, iiStop ) )
                    jPass = True
                except:
                    print ( " json.loads FAILED ... ", openCu, closeCu, openSq, closeSq )
                    print ( " len(tmpString) = ", len(tmpString) )
                    if ( len(tmpString) > 130 ):
                        print ( " %s ... %s " % ( tmpString[:64], tmpString[-64:] ) )
                    else:
                        print ( tmpString )

                    if ( tmpString[-1] == ']' ): 
                        try:
                            print ( "         --> appending curly brace ", openCu, closeCu )
                            tmpString = tmpString + '}'
                            jdata = json.loads(tmpString)
                            print ( "         --> SUCCESS !!!! (%d-%d) " % ( iiStart, iiStop ) )
                            jPass = True
                        except:
                            ## print ( " PERSISTENT FAILURE ??? " )
                            print ( " keep on keepin on ... ", openCu, closeCu, openSq, closeSq )
                            ## sys.exit(-1)
                    else:
                        print ( " NEED TO DO SOMETHING DIFFERENT HERE I THINK " )
                        print ( " --> maybe we just keep going and collect more JSON ??? " )
                        print ( " final character in tmpString : <%s> " % tmpString[-1] )
                        print ( len(tmpString) )
                        if ( len(tmpString) > 130 ):
                            print ( " %s ... %s " % ( tmpString[:64], tmpString[-64:] ) )
                        else:
                            print ( tmpString )
                        print ( " " )

    
                if ( jPass ):
                    fName = inputFilename[:-5] + '.check.' + str(iiStart) + '_' + str(iiStop) + '.json'
                    print ( " --> calling writeJsonString2File ... ", iiStart, iiStop )
                    r = writeJsonString2File ( tmpString, fName )
                    if ( r == 0 ):
                        print ( " THIS SHOULD HAVE WORKED ??? !!! ... calling this garbled " )
                        numGarbled += 1
                    numOutFiles += r
                    openCu = 0
                    closeCu = 0
                    print ( " --> (re)setting iiStart from %d to %d " % ( iiStart, ii+1 ) )
                    iiStart = ii + 1
                else:
                    print ( " jPass is FALSE ... just keep going (?) " )

                otherFlag = True

        else: 
            otherFlag = True

        if not otherFlag: 
            ## print ( ii, openCu, closeCu, openSq, closeSq )
            if ( openCu > 0 ):
                if ( openCu == closeCu ): 
                    print ( " matching curly braces !!! " )

                    ## test if this segment can be parsed as JSON...
                    iiStop = ii + 1
                    jString = allData[iiStart:iiStop]

                    jPass = False

                    try:
                        jdata = json.loads(jString)
                        print ( "         --> SUCCESS !!! (%d-%d) " % ( iiStart, iiStop ) )
                        jPass = True
                    except:
                        print ( "         --> FAILED to parse jString ??? ", iiStart, iiStop ) 


                    if ( jPass ):
                        fName = inputFilename[:-5] + '.check.' + str(iiStart) + '_' + str(iiStop) + '.json'
                        print ( " --> calling writeJsonString2File ... ", iiStart, iiStop )
                        r = writeJsonString2File ( jString, fName )
                        numOutFiles += r
                        if ( r == 0 ):
                            print ( " THIS SHOULD HAVE WORKED ??? !!! ... calling this garbled " )
                            numGarbled += 1
                        openCu  = 0
                        closeCu = 0
                        print ( " --> (re)setting iiStart from %d to %d " % ( iiStart, ii+1 ) )
                        iiStart = ii+1
                    else:
                        print ( " jPass is FALSE ... " )

            if ( openSq > 0 ):
                if ( openSq == closeSq ): 
                    ## print ( " matching square braces !!! " )
                    openSq  = 0
                    closeSq = 0

    if ( numOutFiles < 1 ):
        print ( " NO output files ??? !!! " )

    print ( " Number of output files = ", numOutFiles, " number garbled = ", numGarbled, iiStop, ii )
    print ( " EXITING " )

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    t0 = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--inputFileName', action='store', help='input JSON data file',
                        required=True, dest='inputFile', type=str)
    args = parser.parse_args()

    main(args)

    print ( "\n\n\n" )

# ----------------------------------------------------------------------------------------------------
