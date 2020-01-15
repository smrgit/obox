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

def parseHeaderLine ( aLine ):

    tokens = aLine.split()
    outFilename = tokens[1] + ".tsv"
    fhOut = open ( outFilename, 'w' )

    i1 = aLine.find("(")
    i2 = aLine.find(")")

    csv = aLine[i1+1:i2]
    tokens = csv.split(", ")
    print ( tokens )

    nTokens = len(tokens)

    fhOut.write ( '\t'.join(tokens) )
    fhOut.write ( '\n' )

    return ( nTokens, tokens, fhOut )

##------------------------------------------------------------------------------

def parseDataLine ( aLine, nExp ):

    tsv = aLine[:-1]
    tokens = tsv.split("\t")
    nTokens = len(tokens)

    if ( nTokens != nExp ):
        print ( " ERROR ... unexpected number of tokens ??? !!! " )
        print ( nExp, nTokens )
        print ( tokens )
        sys.exit(-1)

    newList = []
    for t in tokens:
        if ( len(t) == 2 ):
            if ( ord(t[0]) == 92 ):
                if ( ord(t[1]) == 78 ):
                    t = ''

        newList += [t]

    return ( newList )

##------------------------------------------------------------------------------

def parseSQLdumpFile ( fhDmp ):

    fhDmp.seek(0)
    lineNo = 0

    startTable = False

    for aLine in fhDmp:
        if ( aLine.startswith("COPY ") ):
            startTable = True
            print ( " here we are at the start of a data block ... " )
            ( nTokens, tokens, fhOut ) = parseHeaderLine ( aLine )
            numLines = 0
        elif ( ord(aLine[0]) == 92 ):
            print ( ord(aLine[1]) )
            print ( " here we are at the end of the data ... ", numLines )
            print ( aLine ) 
            fhOut.close()
            startTable = False
        else:
            if ( startTable ):
                tokens = parseDataLine ( aLine, nTokens )
                fhOut.write ( '\t'.join(tokens) )
                fhOut.write ( '\n' )
                numLines += 1

##------------------------------------------------------------------------------

def main ( args ):

    print ( args )
  
    try:
        fhDmp = open ( args.sqlDmpFile, 'r' )
    except:
        print ( " Failed to open log file <{}> ... EXITING ".format ( args.sqlDmpFile ) )
        sys.exit(-1)

    parseSQLdumpFile ( fhDmp )

    print ( " " )
    print ( " DONE!!! " )
    print ( " " )
  
##------------------------------------------------------------------------------

if __name__ == '__main__':

  t0 = time.time()

  parser = argparse.ArgumentParser()

  ## the first argument is required
  parser.add_argument ( '-f',  '--inputSQLdumpFileName',  action='store', help='input log file (output of genparse)', required=True, dest='sqlDmpFile', type=str )
  args = parser.parse_args()

  main ( args )

  t1 = time.time()

  ## print ( ' --> time taken in seconds: {dt}'.format(dt=(t1-t0)) )

##------------------------------------------------------------------------------
