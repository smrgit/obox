#### #!/usr/bin/python3

"""
    General purpose CSV parsing and data-munging code ...
"""

__author__ = "Sheila M Reynolds"
__version__ = "0.0.3"
__status__ = "Prototype"

import argparse
import ast
import csv
import datetime
import distutils.util as util
import json
import logging
import os
import pprint
import re
import requests
import string
import sys
import time
import warnings
from collections import Counter

import dateutil
import pandas as pd

# ----------------------------------------------------------------------------------------------------

geneSymbolPattern = re.compile("[A-Z0-9-]+")

# ----------------------------------------------------------------------------------------------------

rspDB = {}
rstDB = {}

def localRefSeqDatabase ( s  ):

    global rspDB, rstDB

    if ( len(rspDB) == 0 ):
        print ( " --> opening local Ensembl identifier database files ...  " )

        ## fh = open ( "/Users/sheila/Data/RefSeq/protein_id.map.tsv2" )
        fh = open ( "/home/sheila/RefSeq/protein_id.map.tsv2" )
        for aLine in fh:
            aLine = aLine.strip()
            tList = aLine.split()
            if ( len(tList) != 2 ): continue
            if ( tList[1] == "gene" ): continue
            gene_symbol = tList[1]
            protein_id = tList[0]
            rspDB[protein_id] = gene_symbol
        fh.close()

        ## fh = open ( "/Users/sheila/Data/RefSeq/transcript_id.map.tsv2" )
        fh = open ( "/home/sheila/RefSeq/transcript_id.map.tsv2" )
        for aLine in fh:
            aLine = aLine.strip()
            tList = aLine.split()
            if ( len(tList) != 2 ): continue
            if ( tList[1] == "gene" ): continue
            gene_symbol = tList[1]
            transcript_id = tList[0]
            rstDB[transcript_id] = gene_symbol
        fh.close()

    print ( " looking up identifier <%s> in current RefSeq databases " % ( s ) )
    if ( s[2] != "_" ):
        print ( "     INVALID RefSeq identifier ??? ", s )
        return ( None )

    if ( s.find('.') >= 0 ):
        s2 = s.split('.')[0]
    else:
        s2 = s

    if ( s2 in rspDB ): return ( rspDB[s] )
    if ( s2 in rstDB ): return ( rstDB[s] )

    print ( "     not found ... " )
    return ( None )

# ----------------------------------------------------------------------------------------------------

ensgDB = {}
enspDB = {}
enstDB = {}
ensNamesDB = {}

def localEnsemblDatabase ( s ):

    global ensgDB, enspDB, enstDB, ensNamesDB

    if ( len(ensgDB) == 0 ):
        print ( " --> opening local Ensembl identifier database files ...  " )

        ## fh = open ( "/Users/sheila/Data/Ensembl/GTF/gene_id.map.tsv" )
        fh = open ( "/home/sheila/Ensembl/GTF/gene_id.map.tsv" )
        for aLine in fh:
            aLine = aLine.strip()
            tList = aLine.split()
            if ( len(tList) != 2 ): continue
            if ( tList[1] == "gene_name" ): continue
            gene_symbol = tList[1]
            gene_id = tList[0]
            ensgDB[gene_id] = gene_symbol
            ensNamesDB[gene_symbol] = gene_id
        fh.close()

        ## fh = open ( "/Users/sheila/Data/Ensembl/GTF/protein_id.map.tsv" )
        fh = open ( "/home/sheila/Ensembl/GTF/protein_id.map.tsv" )
        for aLine in fh:
            aLine = aLine.strip()
            tList = aLine.split()
            if ( len(tList) != 2 ): continue
            if ( tList[1] == "gene_name" ): continue
            gene_symbol = tList[1]
            protein_id = tList[0]
            enspDB[protein_id] = gene_symbol
        fh.close()

        ## fh = open ( "/Users/sheila/Data/Ensembl/GTF/transcript_id.map.tsv" )
        fh = open ( "/home/sheila/Ensembl/GTF/transcript_id.map.tsv" )
        for aLine in fh:
            aLine = aLine.strip()
            tList = aLine.split()
            if ( len(tList) != 2 ): continue
            if ( tList[1] == "gene_name" ): continue
            gene_symbol = tList[1]
            transcript_id = tList[0]
            enstDB[transcript_id] = gene_symbol
        fh.close()

    print ( " looking up identifier <%s> in current Ensembl databases " % ( s ) )

    if ( s.startswith("ENS") ):
        if ( s[3] == "G" ):
            if ( s in ensgDB ): return ( ensgDB[s] )
        if ( s[3] == "P" ):
            if ( s in enspDB ): return ( enspDB[s] )
        if ( s[3] == "T" ):
            if ( s in enstDB ): return ( enstDB[s] )

    if ( s in ensNamesDB ): return ( ensNamesDB[s] )

    print ( "     not found in any Ensembl maps/lists ... ", s )
    return ( None )

# ----------------------------------------------------------------------------------------------------

hgncDB = {}

def localHGNCdatabase ( s ):

    global hgncDB

    if ( len(hgncDB) == 0 ):
        print ( " --> opening local HGNC database file ... " )
        ## fh = open( "/Users/sheila/Data/RefSeq/HGNC_genes.tsv", 'r' )
        fh = open( "/home/sheila/RefSeq/HGNC_genes.tsv", 'r' )
        for aLine in fh:
            aLine = aLine.strip()
            tList = aLine.split()
            if ( len(tList) != 3 ): continue
            if ( tList[0] == "gene_symbol" ): continue
            gene_symbol = tList[0]
            HGNC_ID = tList[1]
            NCBI_geneID = tList[2]
            hgncDB[gene_symbol] = tList
        fh.close()
        print ( " --> have internal HGNC database of length %d " % len(hgncDB) )

    print ( " looking up gene <%s> in current HGNC database (%d) " % ( s, len(hgncDB) ) )
    if ( s in hgncDB ):
        print ( "      FOUND !!! ", hgncDB[s] )
        return ( hgncDB[s] )
    else:
        print ( "      NOT found in local HGNC database ... " )
        return ( [] )

# ----------------------------------------------------------------------------------------------------

def lookupHGNC ( gene ):

    print ( " in lookupHGNC ... ", gene )
    return ( localHGNCdatabase ( gene ) )

# ----------------------------------------------------------------------------------------------------

vepInfoDB = {}

def localVEPdatabase ( s, d ):

    global vepInfoDB

    if ( len(vepInfoDB) == 0 ):
        print ( " --> opening local VEP database file ... " )
        ## fh = open( "/Users/sheila/Data/Ensembl/VEP/vepInfoDB.GRCh37.json", 'r' )
        fh = open( "/home/sheila/Ensembl/VEP/vepInfoDB.GRCh37.json", 'r' )
        for aLine in fh:
            print ( aLine )
            aLine = aLine.strip()
            print ( aLine )
            ## td = eval(aLine)
            td = json.loads(aLine)
            print ( td )
            if ( 'vep_input_hgvs' in td ):
                vep_id = td['vep_input_hgvs']
                ##  do we already have this in our local db/dict?
                if ( vep_id in vepInfoDB ):
                    print ( " already have this <%s> ... need to check dates ... " % vep_id )
                    print ( " existing ... date : ", vepInfoDB[vep_id]['date_created'] )
                    print ( " new ........ date : ", td['date_created'] )
                    existing_date = vepInfoDB[vep_id]['date_created']
                    new_date = td['date_created']
                    if ( existing_date == new_date ):
                        print ( "     --> two entries with same date ... ", vep_id )
                    elif ( existing_date > new_date ):
                        print ( "     --> seems like we should keep the existing one ... " )
                    else:
                        print ( "     --> seems like new one is newer ... " )
                        vepInfoDB[vep_id] = td
                else:
                    vepInfoDB[vep_id] = td
        fh.close()
        print ( " --> have internal VEP database of length %d " % len(vepInfoDB) )

    if ( len(d) == 0 ):
        print ( " looking up variant <%s> in current VEP database (%d) " % ( s, len(vepInfoDB) ) )
        if ( s in vepInfoDB ):
            print ( "      FOUND !!! ", vepInfoDB[s] )
            return ( vepInfoDB[s] )
        else:
            print ( "      NOT found in local VEP database ... " )
            return ( {} )

    if ( len(d) > 0 ):
        if ( s in vepInfoDB ):
            print ( " <%s> is already in local VEP database ??? " % s )
            y = vepInfoDB[s]

            ## remove 'date_created' key so that is not part of the comparison ...
            ## https://stackoverflow.com/questions/11277432/how-to-remove-a-key-from-a-python-dictionary
            y.pop('date_created', None)

            if ( 0 ):
                if ( y != d ):
                    print ( "     --> information is different ??? UHOH ??? !!! " )
                    print ( y )
                    print ( d )
                    sys.exit(-1)

            else:
                for yk in y:
                    if ( y[yk] is not None ):
                        if ( yk in d ):
                            if ( str(y[yk]) != str(d[yk]) ):
                                print ( "     --> different information ??? !!! ", yk, y[yk], d[yk] )
                                print ( "     BAILING BAILING BAILING !!! " )
                                print ( "     y : ", y )
                                print ( "     d : ", d )
                                sys.exit(-1)

        else:
            print ( " ADDING new information to database !!! ", s, d )
            vepInfoDB[s] = d
            print ( "     --> database length is now ", len(vepInfoDB) )

# ----------------------------------------------------------------------------------------------------

def tryVEPagain ( gene, alt, errString ):

    ## right now the only type of error that we're handling is a message of this form:
    ## "Reference allele extracted from $reference:$start-$end ($refseq_allele) does not match reference allele given by HGVS notation $hgvs ($ref_allele)"

    if ( errString.find("Reference allele extracted from") < 0 ): return({})
    if ( errString.find("does not match reference allele given by HGVS notation") < 0 ): return({})

    print ( " " )
    print ( " VEP_AGAIN: in tryVEPagain ... ", gene, alt )
    print ( " VEP_AGAIN:     errString = ", errString )
    print ( " " )

    try:

        op1 = errString.find('(',0)
        cp1 = errString.find(')',op1)
        op2 = errString.find('(',cp1+1)
        cp2 = errString.find(')',op2)
    
        ## print ( op1, cp1, op2, cp2 )
        ## print ( errString[op1+1:cp1] )
        ## print ( errString[op2+1:cp2] )

        s1 = errString[op1+1:cp1]
        s2 = errString[op2+1:cp2]

        is1 = alt.find(s1)
        is2 = alt.find(s2)

        ## print ( " --> ", s1, is1, s2, is2 )

        if ( is1 < 0 or is2 < 0 ): 
            print ( " VEP_AGAIN:    failed to find either parenthetical sequence in the error message  ??? " )
            print ( "               ", s1, s2 )
            return({})
            
        if ( is2 < is1 ):
            new_alt = alt[:is2] + s1 + alt[is2+len(s2):is1] + s2 + alt[is1+len(s1):]
            print ( " VEP_AGAIN:   --> (a) new_alt = ", new_alt )
            return ( callEnsemblVEP ( gene, new_alt ) )
        elif ( is2 > is1 ):
            new_alt = alt[:is1] + s2 + alt[is1+len(s1):is2] + s1 + alt[is2+len(s2):]
            print ( " VEP_AGAIN:   --> (b) new_alt = ", new_alt )
            return ( callEnsemblVEP ( gene, new_alt ) )
        else:
            print ( " VEP_AGAIN:  in tryVEPagain ...  this does not make  sense ??? !!! " )
            return({})

    except:
        print ( " VEP_AGAIN:  in tryVEPagain except clause ... " )

    print ( " VEP_AGAIN:  --> returning nothing from tryVEPagain ... " )
    return({})

# ----------------------------------------------------------------------------------------------------

def getFirstIntegerFromString ( s ):

    print ( "  need to get integer from this string : ", s )

    ks = 0
    while ( ord(s[ks]) < 48 or ord(s[ks]) > 57 ): ks += 1

    js = ks+1
    while ( ord(s[js]) >= 48 and ord(s[js]) <= 57 ): js += 1

    print ( ks, js, s )

    try:
        kVal = int ( s[ks:js] )
        return ( kVal )
    except:
        return ( -999999 )

# ----------------------------------------------------------------------------------------------------

def getGenomicCoords ( gene, alt ):

    print ( " in getGenomicCoords ... ", gene, alt )

    mapBases = {}
    mapBases['A'] = 'T'
    mapBases['T'] = 'A'
    mapBases['C'] = 'G'
    mapBases['G'] = 'C'

    if ( gene == "TERT" and alt[:3] == "c.-" ):
        try:
            m = findLongestInteger ( alt[3:] )
            p = 1295104 + m
            gString = '5:g.' + str(p)

            lenm = len(str(m))
            print ( " --> need to flip : ", alt[6:] )
            ii = 6
            while ii < len(alt):
                if ( alt[ii] in mapBases ):
                    gString += mapBases[alt[ii]]
                else:
                    gString += alt[ii]
                ii += 1

            print ( " --> gString = ", gString )

        except:
            print ( " FAILED TO MAP TERT PROMOTER VARIANT ??? !!! " )
            gString = ''

    else:
        print ( " in getGenomicCoords ... do not know this gene ", gene, alt )
        gString = ''

    return ( gString )

# ----------------------------------------------------------------------------------------------------
##  in callEnsemblVEP ...  RB1 c.107_109delinsCG
##  in callEnsemblVEP ...  RB1 p.Asp36Alafs

def callEnsemblVEP ( gene, alt ):

    keyList = ['gene_symbol', 'gene_symbol_source', 'hgnc_id', 'gene_id', \
               'transcript_id', 'codons', 'cds_start', 'cds_end', \
               'amino_acids', 'protein_start', 'protein_end', 'impact', \
               'consequence_terms', \
               'polyphen_score', 'polyphen_prediction', \
               'sift_score', 'sift_prediction']

    ordered_consequence_terms = \
        ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', \
         'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', \
         'transcript_amplification', 'inframe_insertion', 'inframe_deletion', \
         'missense_variant', 'protein_altering_variant', 'splice_region_variant', \
         'incomplete_terminal_codon_variant', 'start_retained_variant', \
         'stop_retained_variant', 'synonymous_variant', 'coding_sequence_variant', \
         'mature_miRNA_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', \
         'non_coding_transcript_exon_variant', 'intron_variant', 'NMD_transcript_variant', \
         'non_coding_transcript_variant', 'upstream_gene_variant', 'downstream_gene_variant', \
         'TFBS_ablation', 'TFBS_amplification', 'TF_binding_site_variant', \
         'regulatory_region_ablation', 'regulatory_region_amplification', 'feature_elongation', \
         'regulatory_region_variant', 'feature_truncation', 'intergenic_variant']
                                

    print ( " in callEnsemblVEP ... ", gene, alt )

    xd = localVEPdatabase ( gene + ":" + alt, {} )
    if ( len(xd) > 0 ): return ( xd )

    print ( " need to go ping Ensembl REST API ... " )

    if ( 0 ):
        ## here we need to trap things like TERT promoter variants, eg c.-124C>T 
        ## has to be translated to g.1295228G>A
        if ( alt.startswith("c.-") ):
            gString = getGenomicCoords ( gene, alt )
            ext = "/vep/human/hgvs/" + gString + "?"
        else:
            ext = "/vep/human/hgvs/" + gene + ":" + alt + "?"
    else:
        ext = "/vep/human/hgvs/" + gene + ":" + alt + "?"

    server = "https://grch37.rest.ensembl.org"
    try:
        print  ( " calling VEP REST API ...  ", server+ext )
        r = requests.get ( server+ext, headers={"Content-Type":"application/json"}, timeout=360 )
    except: 
        print ( " VEP request error ? ", server+ext )
        print ( " --> returning nothing (a) " )
        return ({})

    ## 20-apr-2020: want to be able to deal with errors that come back from VEP that look like this:
    ##      --> API request FAILED ...  
    ##      https://grch37.rest.ensembl.org/vep/human/hgvs/CUX1:c.2907+1G>A?
    ##          {'error': "Unable to parse HGVS notation 'CUX1:c.2907+1G>A':  : Reference "
    ##           'allele extracted from CUX1:101845452-101845452 (A) does not match '
    ##           'reference allele given by HGVS notation CUX1:c.2907+1G>A (G)'}

    ## errString = Unable to parse HGVS notation 'CUX1:c.2907+1G>A':  : Reference allele extracted from CUX1:101845452-101845452 (A) does not match reference allele given by HGVS notation CUX1:c.2907+1G>A (G)

    ## assuming we got *something* back fromthe API call...
    if not r.ok:
        print ( "     --> API request FAILED ... ", server+ext )
        ## pprint.pprint(r.json())
        try:
            print ( r.json() )
        except:
            print ( " can't even dump the returned result ??? oh well... " )

        ## lets try to trap things here like TERT promoter variants etc
        if ( alt.startswith("c.-") ):
            gString = getGenomicCoords ( gene, alt )
            ext = "/vep/human/hgvs/" + gString + "?"
        else:
            ext = "/vep/human/hgvs/" + gene + ":" + alt + "?"

        try:
            print  ( " calling VEP REST API **AGAIN** ...  ", server+ext )
            r = requests.get ( server+ext, headers={"Content-Type":"application/json"}, timeout=360 )
        except: 
            print ( " VEP request error ? ", server+ext )
            print ( " --> returning nothing (c) " )
            return ({})

        if not r.ok:
            print ( " --> returning nothing (d) " )
            return ( {} )

        ## taking this back OUT ... just a few hours after implementing it ;)
        if ( 0 ):
            try:
                rj = r.json()
                errString = rj['error']
                print ( errString )
                if ( errString.find('Reference allele extracted from') >= 0 ):
                    return ( tryVEPagain ( gene, alt, errString ) )
            except:
                print ( " failed to get error string out of returned JSON ??? " )

    if not r.ok:
        print ( " --> returning nothing (e) " )
        return ( {} )

    else:
        decoded = r.json()
        if ( len(decoded) != 1 ):
            print ( "     --> API ERROR ??? expecting list of length 1 ??? ", len(decoded) )
        else:
            d = decoded[0]
            ## pprint.pprint(d)

            ## TODO: need to extract the first integer that I can from the 'alt' string ...
            firstInteger = getFirstIntegerFromString ( alt )
            print ( " got back firstInteger = ", firstInteger )

            xd = {}
            xd['vep_input_hgvs'] = d['input']
            xd['assembly_name'] = d['assembly_name']
            xd['seq_region_name'] = d['seq_region_name']
            xd['seq_region_start'] = d['start']
            xd['seq_region_end'] = d['end']
            xd['strand'] = d['strand']
            xd['most_severe_consequence'] = d['most_severe_consequence']

            if ( 'transcript_consequences' not in d ):
                print ( "     --> NOTHING much returned by VEP ... ", d )
                print ( " WHY AM I CALLING localVEPdatabase HERE ??? !!! (a) " )
                localVEPdatabase ( d['input'], xd )
                print ( "     --> returning what I have so far (c) ... ", len(xd) )
                return ( xd )

            print ( " transcript_consequences : ", len(d['transcript_consequences']) )
            print ( d['transcript_consequences'] )

            transcriptID_list = []
            for t in d['transcript_consequences']:
                try:
                    if ( alt[0] == 'c' ):
                        if ( t['cds_start'] == firstInteger ):
                            print ( "  CDS_START MATCH ON FIRST INTEGER !!! " )
                            print ( t )
                            print ( " ---- " )
                            ## if ( t['gene_symbol_source'] != 'HGNC' ): continue
                            if ( t['biotype'] != 'protein_coding' ): continue
                            transcriptID_list += [ t['transcript_id'] ]
                    if ( alt[0] == 'p' ):
                        if ( t['protein_start'] == firstInteger ):
                            print ( "  PROTEIN_START MATCH ON FIRST INTEGER !!! " )
                            print ( t )
                            print ( " ---- " )
                            ## if ( t['gene_symbol_source'] != 'HGNC' ): continue
                            if ( t['biotype'] != 'protein_coding' ): continue
                            transcriptID_list += [ t['transcript_id'] ]
                except:
                    ## if we get here it might be because we have a promoter variant...
                    try:
                        if ( alt[0]  == 'c' ):
                            if ( t['distance'] == firstInteger ):
                                print ( " DISTANCE MATCH ON FIRST INTEGER !!! " )
                                print ( t )
                                print  ( " ---- " )
                                if ( t['biotype'] != 'protein_coding' ): continue
                                transcriptID_list += [ t['transcript_id'] ]
                    except:
                        if ( len(transcriptID_list) == 0 ):
                            print ( "  --> hit EXCEPT clause ... is this bad ??? ", alt, t )

            if ( len(transcriptID_list) > 0 ):
                print ( " IDs found : ", len(transcriptID_list), transcriptID_list )
            else:
                print ( " NO transcript IDs found ??? !!!! " )

            print  ( " at line 335 ... " )
            if ( len(transcriptID_list) == 0 ):

                if ( d['most_severe_consequence'] == '?' ):
                    print ( "     --> NOTHING much returned by VEP ... ", d )
                    print ( " WHY AM I CALLING localVEPdatabase HERE ??? !!! (b) " )
                    localVEPdatabase ( d['input'], xd )
                    print ( "     --> returning what I have so far (d) ... ", len(xd) )
                    return ( xd )
    
                ## it is possible that the 'most_severe_consequence' never occurs in a
                ## 'protein_coding' biotype ... so we might need to change how we're
                ## doing this ...
                found_terms = []
                for t in d['transcript_consequences']:
                    if ( t['gene_symbol_source'] != 'HGNC' ): continue
                    if ( t['biotype'] != 'protein_coding' ): continue
                    for e in t['consequence_terms']:
                        if e not in found_terms: found_terms += [ e ]
    
                if ( xd['most_severe_consequence'] not in found_terms ):
                    ## then we want to find the most_severe_consequence in a protein_coding
                    ## biotype if it exists ...
                    if ( len(found_terms) == 0 ):
                        print ( "     --> NOTHING much returned by VEP ... ", d )
                        print ( " WHY AM I CALLING localVEPdatabase HERE ??? !!! (c) " )
                        ## in this call, d['input'] is the HGVS input string, eg: 5:g.1295228G>A
                        ## and xd is what just came back from the VEP call ...
                        localVEPdatabase ( d['input'], xd )
                        print ( "     --> returning what I have so far (e) ... ", len(xd) )
                        return ( xd )
                    iMin = 999
                    for e in found_terms:
                        ii = ordered_consequence_terms.index(e)
                        if ( ii < iMin ): iMin = ii
                    print ( " --> resetting most_severe_consequence: %s to %s " % \
                            ( xd['most_severe_consequence'], ordered_consequence_terms[iMin] ) )
                    xd['most_severe_consequence'] = ordered_consequence_terms[iMin]
    
                ## ok, now we also try and choose the *worst* Polyphen and SIFT scores
                maxScore = -1.
                for t in d['transcript_consequences']:
                    if ( t['gene_symbol_source'] != 'HGNC' ): continue
                    if ( t['biotype'] != 'protein_coding' ): continue
                    if ( d['most_severe_consequence'] in t['consequence_terms'] ):
                        ## if both scores are available, then:
                        try:
                            avgScore = ( t['polyphen_score'] + 1.0 - t['sift_score'] ) / 2.
                        except:
                            ## otherwise, try just using SIFT ...
                            try:
                                avgScore = ( 1.0 - t['sift_score'] )
                            except:
                                try:
                                    avgScore = t['polyphen_score']
                                except:
                                    ## print ( " --> no score information available at all " )
                                    avgScore = -2. 
                        ## print ( " avgScore = ", avgScore )
                        if ( avgScore > maxScore ): maxScore = avgScore
            
                if ( maxScore >= 0. ):
                    scoreFlag = 1
                    print ( " maxScore = ", maxScore )
                    maxScore -= 0.0001
                else:
                    scoreFlag = 0
    
                for t in d['transcript_consequences']:
                    if ( t['gene_symbol_source'] != 'HGNC' ): continue
                    if ( t['biotype'] != 'protein_coding' ): continue
                    if ( d['most_severe_consequence'] in t['consequence_terms'] ):
                        if ( scoreFlag ):
                            ## if both scores are available, then:
                            try:
                                avgScore = ( t['polyphen_score'] + 1.0 - t['sift_score'] ) / 2.
                            except:
                                ## otherwise, try just using SIFT ...
                                try:
                                    avgScore = ( 1.0 - t['sift_score'] )
                                except:
                                    try:
                                        avgScore = t['polyphen_score']
                                    except:
                                        ## print ( " --> no score information available at all " )
                                        avgScore = -2. 
                        else:
                            avgScore = 1.
    
                        if ( avgScore > maxScore ):
                            for k in keyList: 
                                if k in t: 
                                    xd[k] = t[k]
                            ## set this high so that this if clause can only occur once!
                            maxScore = 99.

                if ( 'hgnc_id' in xd ):
                    try:
                        intID = int(xd['hgnc_id'])
                        xd['hgnc_id'] = "HGNC:" + str(intID)
                    except:
                        pass

                print ( " WHY AM I CALLING localVEPdatabase HERE ??? !!! (d) " )
                localVEPdatabase ( d['input'], xd )
                print ( "     --> returning what I have so far (f) ... ", len(xd) )
                return ( xd )

            print  ( " at line 428 ... " )
            if ( len(transcriptID_list) == 1 ):
                pickT = transcriptID_list[0]
            elif ( len(transcriptID_list) > 1 ):
                minInt = 9999999999
                for t in transcriptID_list:
                    intID = int ( t[4:] )
                    if ( intID < minInt ):
                        pickT = t
                        minInt = intID

            ## now go back and pick out the chosen transcript ...
            for t in d['transcript_consequences']:
                if ( t['transcript_id'] == pickT ):
                    print ( " WHAT DO I WANT OUT OF THIS ??? " )
                    print ( "     currenty keyList : ", keyList )
                    print ( "     available information : ", t )
                    print ( " ---- " )
                    print ( " ---- " )
                    for k in keyList: 
                        if k in t: 
                            xd[k] = t[k]

                    ## do we have a coding_change or protein_change string ?!
                    try:
                        if ( alt.startswith('c.') ):
                            if ( 'coding_change' not in xd ):
                                xd['coding_change'] = alt
                            else:
                                if ( xd['coding_change'] != alt ):
                                    print ( " HOW CAN THIS BE ??? !!!! ", xd['coding_change'], alt )
                                    print ( xd )
                                    print ( " " )
                        if ( alt.startswith('p.') ):
                            if ( 'protein_change' not in xd ):
                                xd['protein_change'] = alt
                            else:
                                if ( xd['protein_change'] != alt ):
                                    print ( " HOW CAN THIS BE ??? !!!! ", xd['protein_change'], alt )
                                    print ( xd )
                                    print ( " " )
                    except:
                        pass
           
                    if ( 'hgnc_id' in xd ):
                        try:
                            intID = int(xd['hgnc_id'])
                            xd['hgnc_id'] = "HGNC:" + str(intID)
                        except:
                            pass

                    print ( " WHY AM I CALLING localVEPdatabase HERE ??? !!! (e) " )
                    localVEPdatabase ( d['input'], xd )
                    print ( "     --> returning what I have so far (g) ... ", len(xd) )
                    return ( xd )

            print ( " " )
            print ( "  we should never get here should we ??? " )

    if ( 0 ):
        print ( " " )
        print ( " will return this dict : " )
        pprint.pprint(xd)

    print  ( " at line 553 ... " )
    print ( " WHY AM I CALLING localVEPdatabase HERE ??? !!! (f) " )
    localVEPdatabase ( d['input'], xd )
    print ( "     --> returning what I have so far (h) ... ", len(xd) )
    return ( xd )

# ----------------------------------------------------------------------------------------------------

def myToJson ( fhOut, chunkDf ):

    ## chunkDf.to_json ( fhOut, orient='records', lines=True )

    j1 = chunkDf.to_json(orient='records',lines=True)
    ## print ( " got j1 ... " )
    ## print ( "     type   : ", type(j1) )
    ## print ( "     length : ", len(j1) )

    rows = j1.splitlines(keepends=False)
    ## print ( len(rows) )

    for r in rows:
        ## print ( r )
        try:
            d = eval(r)
        except:
            try:
                d = json.loads(r)
            except:
                print ( " WHY DID THIS HAPPEN in myToJson ??? ", len(rows) )
                print ( r )
                sys.exit(-1)

        if ( not isinstance(d,dict) ):
            print ( " FATAL ERROR in myToJson ??? !!! " )
            print ( r )
            print ( d )
            sys.exit(-1)

        print ( len(d), d )
        sd = {}
        ## 11-apr-2020: added sorted() below
        for k in sorted(d):
            if ( isinstance(d[k],str) ): d[k] = d[k].strip()
            try:
                if ( len(d[k]) == 0 ):
                    ## print ( " skipping key <%s> " % k )
                    pass
                else:
                    sd[k] = d[k]
            except:
                try:
                    if ( d[k] is None ):
                        pass
                except:
                    print ( " WTH is going on here ??? in myToJson ... " )
                    print ( len(d), d )
                    print ( k )
                    print ( d[k] )
                    sys.exit(-1)

        ## if ( len(sd) < len(d) ): print ( len(sd), sd )

        sss = json.dumps(sd,sort_keys=True)
        ## sanity check that there should be either an _id field or an id field
        if ( sss.find('_id": ') < 0 and sss.find('"id": ') < 0 and sss.find('_print": ') < 0 ):
            print ( " " )
            print ( " WHAT IS GOING ON HERE expecting a row identifier ??? !!! " )
            print ( "   r : ", r )
            print ( "   d : ", d )
            print ( "  sd : ", sd )
            print ( " sss : ", sss )
            print ( " " )
            sys.exit(-1)
        
        fhOut.write ( json.dumps(sd,sort_keys=True)+"\n" )

    return ()

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

    ## 04-mar-2020: adding some additional stripping of leading/trailing
    ## punctuation ...
    try:
        cList = [',', '.', ';', ':']
        if ( s[0] in cList ): s = s[1:]
        if ( s[-1] in cList ): s = s[:-1]
    except:
        pass

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

            ## 09-mar-2020: returning a string now
            return ( str(u) )
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
# aList is the list of allowed values, for example: [ 'F', 'M', 'U' ]
# and mDict is the mapping dictionary, eg:
# {'f': 'F', 'm': 'M', 'male': 'M', 'female': 'F', 'unknown': 'U', 'unspecified': 'U'}
# which specifies how other possible inputs (eg 'f', 'm', 'male', 'femal', 'unknown', 'unspecified')
# should be mapped

def enforceAllowed(s, aList, mDict):

    ## start by just returning blank strings as appropriate
    if (str(s) == 'nan'): return ('')
    if (s == ''): return ('')

    t = s.lower()
    print (" in enforceAllowed ... ", aList, mDict, "    input: ", s, t)

    ## loop over the allowed list -- if we already match (case-insensitive)
    ## then we're good to go!
    for a in aList:
        b = a.lower()
        if (t == b):
            print ("     --> found allowed match ", s, a)
            return (a)

    ## loop over the input mapping-dictionary and get all of the keys
    for m in mDict:
        b = m.lower()
        if (t == b):
            print ("     --> applying mapping", s, m, mDict[m])
            return ( mDict[m] )

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

    if ( 0 ):
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

    print ( " in uniqueStringsFromX ... ", x )

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

def hasKnownTerms ( s ):

    knownTerms = { 'amplification':-1, 'deletion':-1, \
                   'rearrangement':-1, 'intronic':-1, \
                   'splice site':-1, 'promoter':-1, 'fusion':-1, \
                   'skipping':-1, 'translocation':-1, \
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
                   'tumor mutational burden':-1, 'TMB':-1, 'burden':-1, \
                   'microsatellite instability':-1, 'microsatellite':-1, 'MSI':-1, 'CNV':-1, 'ITD':-1 }

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

    return ( foundTerms )

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
                   'tumor mutational burden':-1, 'TMB':-1, 'burden':-1, \
                   'microsatellite instability':-1, 'microsatellite':-1, 'MSI':-1, 'CNV':-1, 'ITD':-1 }

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

    ## print ( " checking for common words ... <%s> " % s )

    t = s.lower()
    for w in cw:
        iw = t.find(w)
        if ( iw == 0 ): 
            try:
                if ( t[len(w)]==' ' ): 
                    ## print ( "     FOUND CW ", iw, w )
                    return ( False )
            except:
                pass
        if ( iw > 0 ):
            try:
                if ( t[iw-1]==' ' ):
                    if ( t[iw+len(w)]==' ' ):
                        ## print ( "     FOUND CW ", iw, w )
                        return ( False )
            except:
                pass

    ## print ( "     NO CW FOUND " )
    return ( True )

##----------------------------------------------------------------------------------------------------

def getAccessionAndVersion ( s ):

    accession = 'N/A'
    version = 'N/A'

    ## it's possible that there's a c. or p. string embedded ...
    ## if there is we need to strip that off ... 
    ## note that we are assuming that it is at the end / comes *after* the identifier
    ic = s.find('c.')
    if ( ic > 0 ): 
        if ( s[ic-1] == ':' ): s = s[:ic-1]
    ip = s.find('p.')
    if ( ip > 0 ): 
        if ( s[ip-1] == ':' ): s = s[:ip-1]

    ## even without a clear c. or p. string, there might still be a ':' delimiter
    ix = s.find(':')
    if ( ix > 0 ): s = s[:ix]

    ## 02-jun-2020 ran into a bug here when the input string looked like this: NM_002116.7.4
    iDot = s.find('.')

    if ( iDot > 0 ):
        try:
            testInt1 = int ( s[:iDot] )
            accession = s[:iDot]
        except:
            print ( " (a) failed to extract integer from <%s> ??? " % s[:iDot] )
            print ( " --> trying again ... " )
            try:
                accession = str ( findLongestInteger ( s[:iDot] ) )
            except:
                print ( " giving up ... " )
                accession = ''
                version = ''
                return ( accession, version )

        try:
            testInt2 = int ( s[iDot+1:] )
            version = s[iDot+1:]
        except:
            iu = s[iDot+1:].find('_')
            if ( iu > 0 ):
                try:
                    testInt2 = int ( s[iDot+1:iDot+1+iu] )
                    version = s[iDot+1:iDot+1+iu]
                except:
                    print ( " (b) failed to extract integer from <%s> ??? " % s[iDot+1:] )
                    print ( " --> trying again ... " )
                    try:
                        version = str ( findLongestInteger ( s[iDot+1:] ) )
                    except:
                        print ( " giving up ... " )
                        version = ''
            else:
                version = str ( findLongestInteger ( s[iDot+1:] ) )
                return ( accession, version )

        return ( accession, version )

    else:
        try:
            testInt1 = int ( s )
            accession = s
        except:
            print ( " failed to extract integer from <%s> ??? " % s )
            sys.exit(-1)

        version = ''
        return ( accession, version )

##----------------------------------------------------------------------------------------------------
## RefSeq accession prefixes 
## https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.ensembl_accession_numbers_and_mole/?report=objectonly

def lookForIdentifiers ( u ):

    ## if the second letter is C, G, T, W, or Z it is a genomic identifier
    ## if the second letter is M it is a protein-coding mRNA transcript
    ## if the second letter is R it is a non-protein-coding RNA transcript
    ## if the second letter is P it is a protein identifier
    refseq_prefixes = ['AC', 'NC', 'NG', 'NT', 'NW', 'NZ', 'NM', 'NR', 'XM', 'XR', 'AP', \
                       'NP', 'YP', 'XP', 'WP' ]
    ensembl_prefixes = ['ENSE', 'ENSFM', 'ENSG', 'ENSP', 'ENSR', 'ENST' ]

    refSeqIdList = []
    ensemblIdList = []

    print ( " in lookForIdentifiers ... <%s> " % u )

    uList = u.split()
    nuList = []
    for u in uList:
        if ( u.find(':') > 0 ):
            tu = u.split(':')
            for t in tu:
                if ( t not in nuList ): nuList += [ t ]
        else:
            if ( u not in nuList ): nuList += [ u ]

    uList = nuList
    print ( "     --> split like this : ", uList )

    for v in uList:
        if ( len(v) < 4 ): continue

        ## first: look for NM_ or NP_ style identifiers...
        if ( v[0]=='N' and v[2]=='_' ):
            prefix = v[:2]
            if ( prefix in refseq_prefixes ):
                print ( " looks like RefSeq identifier: ", v )
                ( refseq_accession, refseq_version ) = getAccessionAndVersion ( v[3:] )
                gene_from_refseqID = localRefSeqDatabase ( prefix + '_' + refseq_accession )
                if ( gene_from_refseqID is not None ):
                    print ( " --> got gene from refseq ID ... ", prefix+'_'+refseq_accession, gene_from_refseqID )
                full_id = prefix + '_' + refseq_accession
                if ( refseq_version != '' ): full_id += '.' + refseq_version
                if ( full_id not in refSeqIdList ): refSeqIdList += [ full_id ]
            else:
                print ( " sort of looks like RefSeq, but ??? ", v )

        ## next: look for ENS* style identifiers...
        if ( v[:3] == 'ENS' ):
            prefix = v[:4]
            if ( prefix in ensembl_prefixes ):
                print ( " looks like Ensembl identifier: ", v )
                ( ensembl_accession, ensembl_version ) = getAccessionAndVersion ( v[4:] )
                gene_from_ensemblID = localEnsemblDatabase ( prefix+ensembl_accession )
                if ( gene_from_ensemblID is not None ):
                    print ( " --> got gene from ensembl ID ... ", prefix+ensembl_accession, gene_from_ensemblID )
                full_id = prefix + ensembl_accession
                if ( ensembl_version != '' ): full_id += '.' + ensembl_version
                if ( full_id not in ensemblIdList ): ensemblIdList += [ full_id ]
            else:
                print ( " sort of looks like Ensembl, but ??? ", v )

        ## continue...

    if ( len(refSeqIdList)>0 or len(ensemblIdList)>0 ): 
        print ( " --> list of identifiers found: ", refSeqIdList, ensemblIdList )
        print ( "     from: <%s> " % u )
        try:
            print ( "         gene from RefSeq ID : ", gene_from_refseqID )
        except:
            pass
        try:
            print ( "         gene from Ensembl ID : ", gene_from_ensemblID )
        except:
            pass

    idDict = {}
    if ( 0 ):
        if ( len(refSeqIdList) > 0 ):
            idDict['refseq'] = refSeqIdList
        if ( len(ensemblIdList) > 0 ):
            idDict['ensembl'] = ensemblIdList
    else:

        ## handle RefSeq identifiers ...
        if ( len(refSeqIdList) > 0 ):
            for id in refSeqIdList:
                if ( gene_from_refseqID is not None ): idDict['refseq_gene'] = gene_from_refseqID
                if ( id[1] in ['C','G','T','W','Z'] ):
                    if ( 'refseq_genomic_id' in idDict ):
                        if ( idDict['refseq_genomic_id'] > id ):
                            idDict['refseq_genomic_id'] = id
                    else:
                        idDict['refseq_genomic_id'] = id
                elif ( id[1] in ['M','R'] ):
                    if ( 'refseq_transcript_id' in idDict ):
                        if ( idDict['refseq_transcript_id'] > id ):
                            idDict['refseq_transcript_id'] = id
                    else:
                        idDict['refseq_transcript_id'] = id
                elif ( id[1] in ['P'] ):
                    if ( 'refseq_protein_id' in idDict ):
                        if ( idDict['refseq_protein_id'] > id ):
                            idDict['refseq_protein_id'] = id
                    else:
                        idDict['refseq_protein_id'] = id
                else:
                    print ( " SHOULD NEVER GET HERE NO ??? ", refSeqIdList, id )

        ## handle Ensembl identifiers ...
        if ( len(ensemblIdList) > 0 ):
            if ( gene_from_ensemblID is not None ): idDict['ensembl_gene'] = gene_from_ensemblID
            for id in ensemblIdList:
                if ( id.startswith('ENSG') ):
                    if ( 'ensembl_gene_id' in idDict ):
                        if ( idDict['ensembl_gene_id'] > id ):
                            idDict['ensembl_gene_id'] = id
                    else:
                        idDict['ensembl_gene_id'] = id
                elif ( id.startswith('ENST') ):
                    if ( 'ensembl_transcript_id' in idDict ):
                        if ( idDict['ensembl_transcript_id'] > id ):
                            idDict['ensembl_transcript_id'] = id
                    else:
                        idDict['ensembl_transcript_id'] = id
                elif ( id.startswith('ENSP') ):
                    if ( 'ensembl_protein_id' in idDict ):
                        if ( idDict['ensembl_protein_id'] > id ):
                            idDict['ensembl_protein_id'] = id
                    else:
                        idDict['ensembl_protein_id'] = id
                else:
                    print ( " SHOULD NEVER GET HERE NO ??? ", ensemblIdList, id )

    return ( idDict )

##----------------------------------------------------------------------------------------------------
## returns a list of dicts ...

def lookForGeneSymbols ( u ):

    shortGeneSymbols = [ 'AR', 'C2', 'C3', 'C5', 'C6', 'C7', 'C9', 'CP', 'CS', \
                         'F2', 'F3', 'F5', 'F7', 'F8', 'F9', 'FH', 'GC', 'GK', \
                         'HP', 'HR', 'IK', 'KL', 'KY', 'MB', \
                         'PC', 'SI', 'TF', 'TH', 'XG', 'XK', 'XM', 'XS' ]

    geneFixes = {}
    geneFixes['AC004943.2'] = 'RP5-991G20.1'
    geneFixes['AC055811.2'] = 'RP11-45M22.4'
    geneFixes['AC007780.1'] = 'RP11-120M18.2'
    geneFixes['AC078883.1'] = 'AC093818.1'
    geneFixes['AL139393.2'] = 'RP3-428L16.2'
    geneFixes['AL359922.1'] = 'RP11-145E5.5'
    geneFixes['AP003108.2'] = 'RP11-286N22.8'
    geneFixes['C11ORF30']   = 'EMSY'
    geneFixes['DNMT']       = 'DNMT3A'
    geneFixes['DNMT3']      = 'DNMT3A'
    geneFixes['FAM175A']    = 'ABRAXAS1'
    geneFixes['FAM46C']     = 'TENT5C'
    geneFixes['FGFR']       = 'FGFR1'
    geneFixes['GPR124']     = 'ADGRA2'
    geneFixes['H3F3C']      = 'H3-5'
    geneFixes['HIST1H1C']   = 'H1-2'
    geneFixes['HIST1H3A']   = 'H3C1'
    geneFixes['HIST1H3B']   = 'H3C2'
    geneFixes['HIST1H3J']   = 'H3C12'
    geneFixes['HIST2H3D']   = 'H3C13'
    geneFixes['HIST1H3G']   = 'H3C8'
    geneFixes['HIST3H3']    = 'H3-4'
    geneFixes['H3F3A']      = 'H3-3A'
    geneFixes['H3F3B']      = 'H3-3B'
    geneFixes['MRE11A']     = 'MRE11'
    geneFixes['MYCL1']      = 'MYCL'
##  geneFixes['PARK2']      = 'PRKN'
    geneFixes['PDPK2P']     = 'PDPK2'
    geneFixes['PIK3C']      = 'PIK3CA'
    geneFixes['TOMM70A']    = 'TOMM70'
    geneFixes['TCEB1']      = 'ELOC'
    geneFixes['TP5']        = 'TP53'
    geneFixes['TXLNGY']     = 'TXLNG2P'
    geneFixes['WHSC1']      = 'NSD2'
    geneFixes['WHSC1L1']    = 'NSD3'
    geneFixes['WISP3']      = 'CCN6'

    geneList = []
    symList = []

    try:
        if ( u == '' ): return ( geneList )
    except:
        if ( len(u) == 0 ): return ( geneList )

    print ( " in lookForGeneSymbols ... ", len(u), type(u), u )

    if ( 1 ):

        ## we used to use:      re.split('\W+', u)
        ## but there are quite a few gene symbols that have hyphens in them, so we really
        ## need to first consider that as a possibility ...

        ## so we'll try a couple of approaches and then build the list that we are 
        ## going to test, in the order we want to test...

        vList0 = u.split()
        vList1 = re.findall('(?=\S*[-])([A-Z0-9-]+)', u)
        vList2 = re.split('\W+', u)

        if ( len(vList0) > 0 ): print ( "     vList0 : ", vList0 )
        if ( len(vList1) > 0 ): print ( "     vList1 : ", vList1 )
        if ( len(vList2) > 0 ): print ( "     vList2 : ", vList2 )

        vListP = []
        for v in vList0+vList1+vList2:
            if ( len(v) < 2 ): continue
            if ( len(v)==2 and (v not in shortGeneSymbols) ): continue
            if ( v.lower().startswith('egfrv') ): 
                vListP += [ 'EGFR' ]
            else:
                if ( geneSymbolPattern.fullmatch(v) ): 
                    ## in case this is actually just an integer, try to cast it...
                    try:
                        iv = int(v)
                    except:
                        if ( v not in vListP ): vListP += [ v ]
                else:
                    if ( v.find('.') > 0 ):
                        try:
                            fv = float(v)
                        except:
                            if ( v not in vListP ): vListP += [ v ]

        if ( len(vListP) > 0 ): print ( "     vListP : ", vListP )
        
        print ( " here at line 1363 ... " )
        for v in vListP:
            print ( " planning to test: <%s> " % ( v ) )
            if ( v == 'NM' ): sys.exit(-1)

            ## but if we have already found a valid gene symbol that has this
            ## as a prefix, then we skip...
            alreadyFound = False
            print ( "     --> checking if already found ... ", symList )
            for y in symList:
                if ( y.startswith(v) or y.endswith(v) ): 
                    print ( "     >>> GENE SYMBOL ALREADY FOUND ", y, v )
                    alreadyFound = True

            if ( alreadyFound ): 
                print ( " should be skipping past next block ... " )
                continue

            print ( " not skipping past ... ??? " )
            ## EGFR hack ...
            if ( v != 'EGFR' ):
                ie = v.find('EGFR')
                if ( ie >= 0 and len(v)>(ie+4) ):
                    if ( v != 'EGFR-AS1' ):
                        if ( v[ie+4] == 'v' or v[ie+4] == 'V' ): v = 'EGFR'
                        print ( " --> special EGFR hack to produce : <%s> " % v )

            ## fix some obsolete gene symbols ...
            if ( v in geneFixes ):
                print ( "     --> FIXING gene symbol from %s to %s " % ( v, geneFixes[v] ) )
                v = geneFixes[v]

            ## first check if this is an HGNC gene symbol ...
            try:
                print ( "         calling lookupHGNC with ", v )
                gInfo = lookupHGNC ( v )
                if ( len(gInfo) > 0 ):
                    if ( gInfo[0] not in symList ):
                        gd = {}
                        gd['symbol'] = gInfo[0]
                        gd['HGNC_ID'] = gInfo[1]
                        gd['NCBI_geneID'] = gInfo[2]
                        if ( gd not in geneList ): geneList += [ gd ] 
                        if ( gInfo[0] not in symList ): symList += [ gInfo[0] ]
            except:
                pass

            ## if not, maybe try for Ensembl ???
            if ( len(symList) == 0 ):
                print ( "         no luck with lookupHGNC ... trying Ensembl instead ... ", v )
                s = localEnsemblDatabase ( v )
                if ( s is not None ):
                    t = localEnsemblDatabase ( s )
                    print ( " HEREHERE !!! ", v, s, t )
                    if ( s.startswith("ENS") and (t is not None) ):
                        gd = {}
                        gd['symbol'] = t
                        if ( gd not in geneList ): geneList += [ gd ] 
                        if ( t not in symList ): symList += [ t ]
                    else:
                        gd = {}
                        gd['symbol'] = s
                        if ( gd not in geneList ): geneList += [ gd ] 
                        if ( s not in symList ): symList += [ s ]

        ## BUT...
        ## if we wind up with two gene symbols AND there is a hyphen between the two, maybe
        ## we should also check to see if that is also thought of as a gene? eg BORCS8-MEF2B
        if ( len(geneList) == 2 ):
            print ( " --> because we got TWO gene symbols, we are checking whether this hyphenated gene might itself be a gene ... " )
            print ( "     u : <%s> " % u, type(u) )
            gA = geneList[0]['symbol']
            gB = geneList[1]['symbol']
            iH = u.find('-')
            print ( "     ", gA, gB, iH )
            if ( iH > 0 ):
                iA = u.find(gA)
                iB = u.find(gB)
                if ( iA < 0 or iB < 0 ):
                    print ( " HOW DID THIS HAPPEN ??? why didn't we find the genes ??? ", gA, gB, u )
                else:
                    testG = ''
                    if ( iA < iB ):
                        if ( iA<iH and iH<iB ): testG = gA + '-' + gB
                    if ( iB < iA ):
                        if ( iB<iH and iH<iA ): testG = gB + '-' + gA
                    if ( testG != '' ):
                        try:
                            print ( "         calling lookupHGNC with ", testG )
                            gInfo = lookupHGNC ( v )
                            if ( len(gInfo) > 0 ):
                                if ( gInfo[0] not in symList ):
                                    gd = {}
                                    gd['symbol'] = gInfo[0]
                                    gd['HGNC_ID'] = gInfo[1]
                                    gd['NCBI_geneID'] = gInfo[2]
                                    if ( gd not in geneList ): geneList += [ gd ] 
                                    if ( gInfo[0] not in symList ): symList += [ gInfo[0] ]
                        except:
                            pass

        ## OR
        ## if we wind up with NO genes because, eg HLA-A got truncated to HLA ...
        ## HACKING THIS OUT (12:30pm Thu 05-mar-2020)
        if ( len(geneList) == 0 and False ):
            print ( " --> because we got NO gene symbols, we are trying again, allowing hyphens ... ", u )
            for v in u.split():
                if ( v.find('-') < 0 ): continue
                if ( v.startswith('c.') or v.startswith('p.') ): continue
                if ( len(v) < 2 ): continue
                if ( len(v)==2 and (v not in shortGeneSymbols) ): continue
                if ( geneSymbolPattern.fullmatch(v) ):
                    print ( " from split got this string with a hyphen: <%s> " % ( v ) )
    
                    ## fix some obsolete gene symbols ...
                    if ( v in geneFixes ):
                        print ( "     --> FIXING gene symbol from %s to %s " % ( v, geneFixes[v] ) )
                        v = geneFixes[v]
    
                    try:
                        print ( "         calling lookupHGNC with ", v )
                        gInfo = lookupHGNC ( v )
                        if ( len(gInfo) > 0 ):
                            if ( gInfo[0] not in symList ):
                                gd = {}
                                gd['symbol'] = gInfo[0]
                                gd['HGNC_ID'] = gInfo[1]
                                gd['NCBI_geneID'] = gInfo[2]
                                if ( gd not in geneList ): geneList += [ gd ] 
                                if ( gInfo[0] not in symList ): symList += [ gInfo[0] ]
                    except:
                        pass

    if ( len(geneList) == 0 ):
        print ( " leaving lookForGeneSymbols ... never found a valid gene symbol ... ", u )

    print ( " in lookForGeneSymbols ... returning: ", geneList )

    return ( geneList )

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
## looking at this in detail on 2020-04-24 ... seems like it is causing problems ...

def removeCertainBlanks ( u ):

    u = u.strip()
    if ( u == '' ): return ( '' )

    t = u
    print ( " hacking blanks ... ", t )

    ## actually we want to replace any of these with just a space
    for c in [':', '[', ']', '(', ')']:
        while ( t.find(c) >= 0 ):
            ii = t.find(c)
            t = t[:ii].strip() + ' ' + t[ii+1:].strip()
            t = t.strip()
        print ( "     --> ", t )

    ## however, if there is a ' =' I think we need to remove the space
    ## before the equal sign ...
    while ( t.find(' =') >= 0 ):
        ii = t.find(' =')
        t = t[:ii].strip() + t[ii+1:].strip()
        t = t.strip()
    print ( "     --> ", t )

    ## this doesn't belong here, but ...
    ## need to make sure there is a blank before a 'c.'
    ii = 0
    while ( t.find('c.',ii) > 0 ):
        ii = t.find('c.',ii)
        print ( "     --> ", ii, t )
        if ( t[ii-1] != ' ' ):
            jj = len(t[:ii].strip()) + 3
            t = t[:ii].strip() + ' ' + t[ii:].strip()
            t = t.strip()
            ii = t.find('p.',jj) + 2
        else:
            ii = t.find('c.',ii) + 2

    ## need to make sure there is a blank before a 'p.'
    ii = 0
    while ( t.find('p.',ii) > 0 ):
        ii = t.find('p.',ii)
        print ( "     --> ", ii, t )
        if ( t[ii-1] != ' ' ):
            jj = len(t[:ii].strip()) + 3
            t = t[:ii].strip() + ' ' + t[ii:].strip()
            t = t.strip()
            ii = t.find('p.',jj) + 2
        else:
            ii = t.find('p.',ii) + 2

    ## on other hand, we need to make sure there is no blank
    ## AFTER 'c.' or 'p.' ...
    while ( t.find('c. ') > 0 ):
        ii = t.find('c. ')
        t = t[:ii+2].strip() + t[ii+3:].strip()
    print ( "     --> ", t )
        
    while ( t.find('p. ') > 0 ):
        ii = t.find('p. ')
        t = t[:ii+2].strip() + t[ii+3:].strip()
    print ( "     --> ", t )
        
    if ( t != u ):
        print ( "     --> changes made by removeCertainBlanks: ", u, " --> ", t )

    return ( t )
    
##----------------------------------------------------------------------------------------------------

def removeGenesFromString ( s, foundGenes ):

    print ( " !!! in removeGenesFromString ", s, foundGenes )

    cList = [',', '.', ';', ':']

    for g in foundGenes:
        try:
            y = g['symbol']
            while ( s.find(y) >= 0 ):
                iy = s.find(y)
                s1 = s[:iy].strip()
                if ( len(s1) > 0  and  s1[0] in cList ): s1 = s1[1:]
                if ( len(s1) > 0  and  s1[-1] in cList ): s1 = s1[:-1]
                s2 = s[iy+len(y):].strip()
                if ( len(s2) > 0  and  s2[0] in cList ): s2 = s2[1:]
                if ( len(s2) > 0  and  s2[-1] in cList ): s2 = s2[:-1]
                s = s1 + ' ' + s2
                s = s.strip()
        except:
            pass

    print ( " !!! --> returning : <%s> " % s )
    return ( s )

##----------------------------------------------------------------------------------------------------
## we come in here with an input string s and a dict of 'found' identifiers, eg:
## {'refseq_protein_id': 'NP_004976.2', 'refseq_transcript_id': 'NM_004985.4'}

def removeIdentifiersFromString ( s, foundIDs ):

    print ( " !!! in removeIdentifiersFromString ", s, foundIDs )

    cList = [',', '.', ';', ':']

    ## outer loop is over the keys (eg refseq_protein_id or refseq_transcript_id)
    for aKey in foundIDs:
        try:
            y = foundIDs[aKey]
            while ( s.find(y) >= 0 ):
                iy = s.find(y)
                s1 = s[:iy].strip()
                if ( len(s1) > 0  and  s1[0] in cList ): s1 = s1[1:]
                if ( len(s1) > 0  and  s1[-1] in cList ): s1 = s1[:-1]
                s2 = s[iy+len(y):].strip()
                if ( len(s2) > 0  and  s2[0] in cList ): s2 = s2[1:]
                if ( len(s2) > 0  and  s2[-1] in cList ): s2 = s2[:-1]
                s = s1 + ' ' + s2
                s = s.strip()
        except:
            pass

    print ( " !!! --> returning : <%s> " % s )
    return ( s )

#### OLD CODE -- KEEPING TEMPORARILY ...
####  def removeIdentifiersFromString ( s, foundIDs ):
####  
####      print ( " !!! in removeIdentifiersFromString ", s, foundIDs )
####  
####      cList = [',', '.', ';', ':']
####  
####      ## outer loop is over the keys (eg refseq_protein_id or refseq_transcript_id)
####      for aKey in foundIDs:
####          try:
####              aList = foundIDs[aKey]
####              for y in aList:
####                  while ( s.find(y) >= 0 ):
####                      iy = s.find(y)
####                      s1 = s[:iy].strip()
####                      if ( len(s1) > 0  and  s1[0] in cList ): s1 = s1[1:]
####                      if ( len(s1) > 0  and  s1[-1] in cList ): s1 = s1[:-1]
####                      s2 = s[iy+len(y):].strip()
####                      if ( len(s2) > 0  and  s2[0] in cList ): s2 = s2[1:]
####                      if ( len(s2) > 0  and  s2[-1] in cList ): s2 = s2[:-1]
####                      s = s1 + ' ' + s2
####                      s = s.strip()
####          except:
####              pass
####  
####      print ( " !!! --> returning : <%s> " % s )
####      return ( s )

##----------------------------------------------------------------------------------------------------
## TODO: may also want to try and remove duplicates from string, eg:
## TERT TERT NM_198253.2:c.-146C>T promoter c.-146C>T promoter

def trimWhiteSpace ( s ):

    s = s.strip()
    if ( len(s) == 0 ): return ( '' )

    t = s[0]
    for ii in range(1,len(s)):
        if ( s[ii]==' ' and t[-1]==' ' ):
            continue
        else:
            t += s[ii]

    while ( t.find(':') > 0 ):
        jj = t.find(':')
        t = t[:jj] + ' ' + t[jj+1:]

    u = t.split()
    uu = []
    for v in u:
        if ( v.startswith("na_for_") ):
            print ( "     --> ditching the na_for_ string ... ", v )
            pass
        elif ( v.startswith("AF=1000") ):
            print ( "     --> ditching the AF=1000 string ... ", v )
            pass
        else:
            if ( v not in uu ): uu += [ v ]

    t = ''
    for v in uu:
        t += v + ' '

    t = t.strip()
    return ( t )

##----------------------------------------------------------------------------------------------------
## aGene can look like this: {'symbol': 'SH2B3', 'HGNC_ID': '29605', 'NCBI_geneID': '10019'}

def fixListOrder ( z, aGene ):

    print ( " in fixListOrder ... " )
    print ( z )
    print ( aGene )
    if ( 'symbol' in aGene ):
        gene_symbol = aGene['symbol']
    else:
        gene_symbol = "NO*KNOWN*GENE"

    gString = ''
    cString = ''
    pString = ''
    nString = ''
    rString = ''

    otherz = []

    for ii in range(len(z)):
        if ( z[ii] == gene_symbol ):   
            gString = z[ii]
        elif ( z[ii].startswith("c.") ): 
            cString = z[ii]
        elif ( z[ii].startswith("p.") ): 
            pString = z[ii]
        elif ( z[ii].startswith("n.") ): 
            nString = z[ii]
        elif ( z[ii].startswith("r.") ): 
            rString = z[ii]
        else:
            print ( " in fixListOrder ... what is this ??? ", z[ii], ii, z )
            otherz += [ z[ii] ]

    print ( " what I found : ", gString, cString, pString, nString, rString )

    newz = []
    if ( gString != '' ): newz += [ gString ]
    if ( cString != '' ): newz += [ cString ]
    if ( pString != '' ): newz += [ pString ]
    if ( nString != '' ): newz += [ nString ]
    if ( rString != '' ): newz += [ rString ]

    newz += otherz

    print ( " --> returning from fixListOrder : ", newz )

    return ( newz )

##----------------------------------------------------------------------------------------------------
## example:
##     we start with:                 ASXL1 c.2083C>T: p.Q695X
##     after "first pass" we have:  ['ASXL1 c.2083C>T: p.Q695X']
##     after "second pass" we have: ['ASXL1', 'c.2083C>T', 'p.Q695X']

def handleFreeAlterationText ( u, vepFlag ):

    ## special hack
    if ( u.lower() == "not available" ): return ( {} )

    ## special hack for: 'TERT promoter -124C>T'
    if ( u.find(' -124C>T') >= 0 ):
        i1 = u.find(' -124C>T')
        nu = u[:i1] + ' c.-124C>T'
        if ( len(u) > (i1+7) ): nu += u[i1+8:]
        print ( " --> special hack changing {%s} to {%s} " % ( u, nu ) )
        u = nu

    ## special hack for the "EML4" typo as "ELM4" ... 
    if ( u.find('ELM4') >= 0 ):
        i1 = u.find('ELM4')
        nu = u[:i1] + 'EML4'
        if ( len(u) > (i1+4) ): nu += u[i1+4:]
        print ( " --> special hack changing {%s} to {%s} " % ( u, nu ) )
        u = nu

    print ( " ... in handleFreeAlterationText ... ", type(u), len(u), u, vepFlag )

    ## initialize the "return dict" ...
    rd = {}

    ## this first step just removes blanks before and after '=' and ':'
    u = removeCertainBlanks ( u )

    ## note that this is called after literal_eval so there should not be any '\n'
    ## characters in the input string ...
    if ( u.find('\n') >= 0 ): 
        print ( " in handleFreeAlterationText ... NOT expecting EOL character ??? !!! " )
        print ( " NOT EXITING ... YET ... (1) " )
        ## sys.exit(-1)

    ## check for AF= ...
    (sRest, afD) = extractAF(u)

    ## reset string u and dict rd ...
    u = sRest
    rd = afD

    ## now check if there are any 'common' words in this string ...
    noCW = noCommonWords(u)
    if ( not noCW ):
        print ( " in handleFreeAlterationText ... HAS common words ... " )
    else:
        print ( " in handleFreeAlterationText ... NO common words found ... " )

    ## now check for various 'known terms' ...
    knownTerms = hasKnownTerms(u)
    numKnownTerms = len(knownTerms)
    if ( numKnownTerms > 0 ):
        print ( " in handleFreeAlterationText ...  HAS %d known terms ... " % numKnownTerms )
        print ( knownTerms )
    else:
        print ( " in handleFreeAlterationText ...  has NO known terms ... " )

    ## now lets see if we find some gene symbols ...
    foundGenes = lookForGeneSymbols(u)
    if ( len(foundGenes) > 0 ): print ( " found gene symbols : ", foundGenes )

    ## and lets see if we find any identifiers, eg 'NM_001754.4' or 'NP_001745.2' or 'ENS...'
    foundIds = lookForIdentifiers(u)

    if ( len(foundIds) > 0 ): 

        print ( " found identifiers : ", foundIds )
        rd.update ( foundIds )

        ## and then we also want to remove them from the string...
        u = removeIdentifiersFromString ( u, foundIds )

        ## and if we don't have a 'foundGene', we might be able to pull them out of 
        ## what came back from lookForIdentifiers ...
        if ( len(foundGenes) == 0 ):
            gd = {}
            if ( 'ensembl_gene' in foundIds ):
                gd['symbol'] = foundIds['ensembl_gene']
            elif ( 'refseq_gene' in foundIds ):
                gd['symbol'] = foundIds['refseq_gene'] 
            try:
                if ( gd['symbol'] is not None ):
                    foundGenes = [ gd ]
            except:
                pass
            if ( len(foundGenes) > 0 ):
                print ( " managed to rescue the lack of a gene symbol !!! ", foundGenes )

    ## check to see if it 'looks like' HGVS ...
    looksLikeHGVS = lookForHGVSabbrevs(u)

    ## example 'bad' string: ARc.1385_1420del:p.Gly462_Gly473del 
    if ( looksLikeHGVS and len(foundGenes)==0 ):
        ic = u.find('c.')
        if ( ic > 0 ):
            if ( u[ic-1] != ' ' ):
              new_u = u[:ic] + ' ' + u[ic:]
              print ( " change string from <%s> to <%s> ??? " % ( u, new_u ) )

              new_foundGenes = lookForGeneSymbols(new_u)
              if ( len(new_foundGenes) >= 0 ):
                print ( " --> YAY !!! now able to find a gene !!! ", new_foundGenes )
                u = new_u
                foundGenes = new_foundGenes


    ## ----
    ## if we found any valid gene symbols, what we do may depend on how many we found ...
    if ( len(foundGenes) == 0 ):
        if ( looksLikeHGVS ):
            print ( " NO GENE SYMBOLS FOUND ??? !!! ", u )
        pass

    elif ( len(foundGenes) == 1 ):
        print ( "     --> found ONE gene symbol! " )

        ## then check if it looks like there are HGVS variation description terms in this string ...
        if ( looksLikeHGVS ):
            print ( " in handleFreeAlterationText ... looks like HGVS ... " )
            z = uniqueStringsFromX ( re.split('\s|:', u) )
            print ( z )
            ## ideally at this point we have something like: ['TET2', 'c.2578C>T', 'p.Q860*']
            z = fixListOrder ( z, foundGenes[0] )

            ## >>>>
            ## special hack ...

            if ( len(z) == 3 ):
                print ( " >>>> special hack #1 ... " )

                ## first check if we have a c.* string and p.* string, then all good
                if ( z[1].startswith('c.') and len(z[1])>2 ):
                    if ( z[2].startswith('p.') and len(z[2])>2 ):
                        ## good
                        pass
                    else:
                        if ( z[2] not in ['promoter','intronic',"5'UTR"] ):
                            print ( " --> adding p. prefix to ", z[2] )
                            z[2] = 'p.' + z[2]

                ## if we just have c but not c. ...
                elif ( z[1].startswith('c') and len(z[1])>2 ):
                    print ( " --> inserting a . into: ", z[1] )
                    z[1] = 'c.' + z[1][1:]
                    if ( z[2].startswith('p.') and len(z[2])>2 ):
                        ## good
                        pass

                else:
                    ## if we get here, then z[1] doesn't start with either c or c.
                    if ( z[2].startswith('p.') and len(z[2])>2 ):
                        print ( " --> inserting a c. in front of ", z[1] )
                        z[1] = 'c.' + z[1]

                    elif ( z[2].startswith('p') and len(z[2])>2 ):
                        print ( " --> inserting a . into: ", z[2] )
                        z[2] = 'p.' + z[2][1:]
        
                        print ( " --> inserting a c. in front of ", z[1] )
                        z[1] = 'c.' + z[1]

                        pass
                print ( "         --> result of special hack: ", z )
            ## <<<<

            ## >>>>
            ## 26-apr-2020 I look at this below and I have NO clue what this is 
            ## even trying to do ...
            ## --> OH ... ok, it's trying to merge strings if somehow the c. or the p.
            ## get separated from the information that comes next ...
            ## special hack ...
            if ( len(z) >= 3 ):
                print ( " >>>> special hack #2 ... " )
                newz = []
                for ii in range(len(z)):
                    if ( z[ii] == 'c.' ):
                        try:
                            newz += [ z[ii]+z[ii+1] ]
                            z[ii+1] = ''
                        except:
                            pass
                    elif ( z[ii] == 'p.' ):
                        try:
                            newz += [ z[ii]+z[ii+1] ]
                            z[ii+1] = ''
                        except:
                            pass
                    else:
                        newz += [ z[ii] ]
                print ( "         --> result of special hack: ", newz )
                z = newz
            ## <<<<

            ## try and verify coding and protein alterations ...
            cDot = ''
            pDot = ''
            for ii in range(len(z)):
                s = z[ii]
                if ( s.startswith('c.') ): 
                    t = validateCodingAlt ( s )
                    if ( z[ii] != t ):
                        print ( "     --> changed from {%s} to {%s} " % ( s, t ) )
                        z[ii] = t
                    if ( cDot == '' ): 
                        if ( t != '' ):
                            cDot = t 
                    else:
                        if ( t != '' ):
                            if ( cDot != t ):
                                print ( " --> UHOH ??? DIFFERENT coding alterations in same string ??? !!! ", cDot, t )
                                print ( "     --> arbitrarily keeping only the first one: ", cDot )

                if ( s.startswith('p.') ): 
                    t = threeLetterAAs ( s )
                    if ( z[ii] != t ):
                        print ( "     --> changed from {%s} to {%s} " % ( s, t ) )
                        z[ii] = t
                    if ( pDot == '' ): 
                        if ( t != '' ):
                            pDot = t 
                    else:
                        if ( t != '' ):
                            if ( pDot != t ):
                                print ( " --> UHOH ??? DIFFERENT coding alterations in same string ??? !!! ", pDot, t )
                                print ( "     --> arbitrarily keeping only the first one: ", pDot )

            ## ok, at this point we know we have one gene-symbol in foundGenes[0]['symbol']
            ## and we also have some HGVS abbreviations ... if we have 'c.' or 'p.' we'll 
            ## try to use either of those to call VEP ...
            vepInfo = {}
            if ( vepFlag == "withVEP" ):
                if ( cDot != '' ):
                    vepInfo = callEnsemblVEP ( foundGenes[0]['symbol'], cDot )
                    if ( len(vepInfo) == 0 ):
                        if ( 'ensembl_gene' in foundIds ):
                            if ( foundIds['ensembl_gene'] != foundGenes[0]['symbol'] ):
                                 print ( " --> trying again with ensembl_gene ... " )
                                 vepInfo = callEnsemblVEP ( foundIds['ensembl_gene'], cDot )
                    if ( len(vepInfo) == 0 ):
                        if ( 'refseq_gene' in foundIds ):
                            if ( foundIds['refseq_gene'] != foundGenes[0]['symbol'] ):
                                 print ( " --> trying again with refseq_gene ... " )
                                 vepInfo = callEnsemblVEP ( foundIds['refseq_gene'], cDot )
    
                if ( len(vepInfo) == 0 ):
                    if ( pDot != '' ):
                        vepInfo = callEnsemblVEP ( foundGenes[0]['symbol'], pDot )
                        if ( len(vepInfo) == 0 ):
                            if ( 'ensembl_gene' in foundIds ):
                                if ( foundIds['ensembl_gene'] != foundGenes[0]['symbol'] ):
                                    print ( " --> trying again with ensembl_gene ... " )
                                    vepInfo = callEnsemblVEP ( foundIds['ensembl_gene'], pDot )
                        if ( len(vepInfo) == 0 ):
                            if ( 'refseq_gene' in foundIds ):
                                if ( foundIds['refseq_gene'] != foundGenes[0]['symbol'] ):
                                    print ( " --> trying again with refseq_gene ... " )
                                    vepInfo = callEnsemblVEP ( foundIds['refseq_gene'], pDot )

            if ( len(vepInfo) > 0 ):
                print ( " got VEP info : ", vepInfo )

                ## we don't need to store the details of the VEP information, just the hgvs key
                rd['vep_input_hgvs'] = vepInfo['vep_input_hgvs']
                rd.update ( {'biomarker_type':'small_variant', 'gene_symbol':foundGenes[0]['symbol']} )
                if ( cDot != '' ): rd['coding_change'] = cDot
                if ( pDot != '' ): rd['protein_change'] = pDot

                ## if there are any "known" terms, add those in now:
                if ( numKnownTerms > 0 ):
                    for kt in knownTerms:
                        rd[kt] = 'True'

                print ( " " )
                print ( " --> returning: ", rd )
                return ( rd )

            else:
                print ( " did NOT get any VEP info ??? " )

                rd.update ( {'biomarker_type':'small_variant', 'gene_symbol':foundGenes[0]['symbol']} )
                if ( cDot != '' ): rd['coding_change'] = cDot
                if ( pDot != '' ): rd['protein_change'] = pDot

                ## if there are any "known" terms, add those in now:
                if ( numKnownTerms > 0 ):
                    for kt in knownTerms:
                        rd[kt] = 'True'

                print ( " " ) 
                print ( " --> returning: ", rd )
                return ( rd )

        else:
            print ( " in handleFreeAlterationText ... does NOT look like HGVS ... " )

            ## handle special case of an amplification
            if ( u.lower().find('amplification') >= 0  or  u.lower().find('amplificaton') >= 0 ):
                rd.update ( {'biomarker_type':'CNV', \
                             'gene_symbol':foundGenes[0]['symbol'], 'variant_name':'amplification'} )
                return ( rd )

            ## handle special case of a translocation
            if ( u.lower().find('translocation') >= 0 ):
                print ( " translocation string : <%s> " % u )
                rd.update ( {'biomarker_type':'complex_variant', 'gene_symbol':foundGenes[0]['symbol'] } )
                uTokens = u.split()
                print ( uTokens )
                vname = ''
                for t in uTokens:
                    if ( t.lower() == 'translocation' ): continue
                    if ( t == rd['gene_symbol'] ): continue
                    vname += t + ' '
                if ( len(vname) > 1 ):
                    print ( " vname : <%s> " % vname )
                    rd['variant_name'] = vname.strip()
                    print ( " --> returning: ", rd )
                    return ( rd )
                else:
                    print ( " WHAT ELSE CAN I DO HERE ??? " )
                    print ( " NOT EXITING ... YET ... (2) " )
                    ## sys.exit(-1)

            ## handle special case of a rearrangement
            if ( u.lower().find('rearrangement') >= 0 ):
                print ( " rearrangement string : <%s> " % u )
                rd.update ( {'biomarker_type':'complex_variant', 'gene_symbol':foundGenes[0]['symbol'] } )
                uTokens = u.split()
                print ( uTokens )
                vname = ''
                for t in uTokens:
                    if ( t.lower() == 'rearrangement' ): continue
                    if ( t == rd['gene_symbol'] ): continue
                    vname += t + ' '
                if ( len(vname) > 1 ):
                    print ( " vname : <%s> " % vname )
                    rd['variant_name'] = vname.strip()
                    print ( " --> returning: ", rd )
                    return ( rd )
                else:
                    print ( " WHAT ELSE CAN I DO HERE ??? " )
                    print ( " NOT EXITING ... YET ... (2b) " )
                    ## sys.exit(-1)

            ## handle special case of a splice site
            if ( u.lower().find('splice site') >= 0 ):
                print ( " splice site string : <%s> " % u )
                rd.update ( {'biomarker_type':'complex_variant', 'gene_symbol':foundGenes[0]['symbol'] } )
                uTokens = u.split()
                print ( uTokens )
                vname = ''
                for t in uTokens:
                    if ( t.lower() == 'splice' ): continue
                    if ( t.lower() == 'site' ): continue
                    if ( t == rd['gene_symbol'] ): continue
                    vname += t + ' '
                if ( len(vname) > 1 ):
                    print ( " vname : <%s> " % vname )
                    rd['variant_name'] = vname.strip()
                    print ( " --> returning: ", rd )
                    return ( rd )
                else:
                    print ( " WHAT ELSE CAN I DO HERE ??? " )
                    print ( " NOT EXITING ... YET ... (3) " )
                    ## sys.exit(-1)

            ## handle special case of EGFR vIII variant
            if ( foundGenes[0]['symbol'] == 'EGFR' ):
                if ( u.lower().find('viii') >= 0 ):
                    rd.update ( {'biomarker_type':'complex_variant', 'gene_symbol':'EGFR', 'variant_name':'vIII'} )
                    return ( rd )

            ## handle special case of single-gene fusion ... which we will call a 'complex_variant' 
            ## rather than a 'fusion' ...
            if ( numKnownTerms > 0 ):
                print ( numKnownTerms, knownTerms )
                if ( 'fusion' in knownTerms ):
                    uTokens = u.split()
                    rd.update ( {'biomarker_type':'complex_variant', 'gene_symbol':foundGenes[0]['symbol']} )
                    w = ''
                    for v in uTokens:
                        if ( v != foundGenes[0]['symbol'] ): w += v + ' '
                    rd.update ( {'variant_name': w.strip()} )
                    return ( rd )
                    
            ## handle exon-skipping similar to above fusion ...
            if ( numKnownTerms > 0 ):
                print ( numKnownTerms, knownTerms )
                if ( 'skipping' in knownTerms ):
                    uTokens = u.split()
                    rd.update ( {'biomarker_type':'complex_variant', 'gene_symbol':foundGenes[0]['symbol']} )
                    w = ''
                    for v in uTokens:
                        if ( v != foundGenes[0]['symbol'] ): w += v + ' '
                    rd.update ( {'variant_name': w.strip()} )
                    return ( rd )
                    
            ## handle exon deletion similar to above skipping ...
            if ( numKnownTerms > 0 ):
                print ( numKnownTerms, knownTerms )
                if ( 'deletion' in knownTerms ):
                    uTokens = u.split()
                    rd.update ( {'biomarker_type':'complex_variant', 'gene_symbol':foundGenes[0]['symbol']} )
                    w = ''
                    for v in uTokens:
                        if ( v != foundGenes[0]['symbol'] ): w += v + ' '
                    rd.update ( {'variant_name': w.strip()} )
                    return ( rd )

            ## handle a "CNV" ... 
            if ( numKnownTerms > 0 ):
                if ( 'CNV' in knownTerms ):
                    u = removeGenesFromString ( u, foundGenes )
                    u = trimWhiteSpace(u)
                    rd.update ( {'biomarker_type':'CNV', 'gene_symbol':foundGenes[0]['symbol']} )
                    if ( len(u) > 1 ): rd['free_text'] = u
                    return ( rd )

            ## handle a "rearrangement" ... 
            if ( numKnownTerms > 0 ):
                if ( 'rearrangement' in knownTerms ):
                    u = removeGenesFromString ( u, foundGenes )
                    u = trimWhiteSpace(u)
                    rd.update ( {'biomarker_type':'complex_variant', 'gene_symbol':foundGenes[0]['symbol']} )
                    if ( len(u) > 1 ): rd['free_text'] = u
                    return ( rd )

            ## handle a "ITD" (internal tandem duplication) ... 
            if ( numKnownTerms > 0 ):
                if ( 'ITD' in knownTerms ):
                    u = removeGenesFromString ( u, foundGenes )
                    u = trimWhiteSpace(u)
                    rd.update ( {'biomarker_type':'complex_variant', 'gene_symbol':foundGenes[0]['symbol']} )
                    if ( len(u) > 1 ): rd['free_text'] = u
                    return ( rd )

            ## what should we do with any remaining knownTerms ???
            if ( numKnownTerms > 0 ):
                print ( " got this far ... what about knownTerms ? ", numKnownTerms, knownTerms )
                    
            ## if we get this far, we'll just remove any "found" gene symbols from
            ## our text string and continue ...
            rd['gene_symbol'] = foundGenes[0]['symbol']
            u = removeGenesFromString ( u, foundGenes )
            u = trimWhiteSpace(u)
            if ( len(u) > 1 ): rd['free_text'] = u
            print ( " RETURNING from handleFreeAlterationText ... ", rd )
            return ( rd )

    elif ( len(foundGenes) == 2 ):
        gA = foundGenes[0]['symbol']
        gB = foundGenes[1]['symbol']
        print ( "     --> found TWO gene symbols! ", gA, gB )

        ## if we do not have any knownTerms (in particular not the word "fusion"), maybe we can *assume*
        ## that this is a gene fusion if there is a '-' between these two gene symbols...
        if ( 'fusion' not in knownTerms ):
            iA = u.find(gA)
            iB = u.find(gB)
            iH = u.find('-')
            if ( iH > 0 ):
                if ( iA < iB ):
                    if ( iH<iB and iH>iA ): knownTerms += ['fusion']
                    if ( iH<iA and iH>iB ): knownTerms += ['fusion']
            numKnownTerms = len(knownTerms)

        if ( numKnownTerms > 0 ):
            print ( numKnownTerms, knownTerms )
            if ( 'fusion' in knownTerms ):

                u0 = ''
                u1 = ''

                uTokens = u.split('-')
                if ( len(uTokens) == 2 ):
                    u0 = uTokens[0].split()[-1]
                    u1 = uTokens[1].split()[0]
                    if ( u0[0] == '(' ): u0 = u0[1:]
                    if ( u1[-1] == ')' ): u1 = u1[:-1]
               
                if ( ( u0==gA or u0==gB ) and ( u1==gA or u1==gB ) ):
                    ## yay!
                    pass
                else:
                    iA = u.find(gA)
                    iB = u.find(gB)
                    if ( iA < iB ):
                        u2 = u[iA:iB+len(gB)]
                        u0 = gA
                        u1 = gB
                    else:
                        u2 = u[iB:iA+len(gA)]
                        u0 = gB
                        u1 = gA

                    print ( " trying to figure out fusion string ... ", u, u2, u0, u1 )
                
                ## removed 'biomarker_type':'3'
                if ( u0==gA and u1==gB ):
                    return ( {'biomarker_type':'fusion', 'a_gene_symbol':gA, 'b_gene_symbol':gB} )
                elif ( u0==gB and u1==gA ):
                    return ( {'biomarker_type':'fusion', 'a_gene_symbol':gB, 'b_gene_symbol':gA} )
                elif ( u0==gA ):
                    return ( {'biomarker_type':'fusion', 'a_gene_symbol':gA, 'b_gene_symbol':gB} )
                elif ( u0==gB ):
                    return ( {'biomarker_type':'fusion', 'a_gene_symbol':gB, 'b_gene_symbol':gA} )
                elif ( u1==gA ):
                    return ( {'biomarker_type':'fusion', 'a_gene_symbol':gB, 'b_gene_symbol':gA} )
                elif ( u1==gB ):
                    return ( {'biomarker_type':'fusion', 'a_gene_symbol':gA, 'b_gene_symbol':gB} )

                else:
                    print ( " problem figuring out fusion ??? ", u0, u1, uTokens )
                    print ( " NOT EXITING ... YET ... (6) " )
                    ## sys.exit(-1)

        print ( " (b) WHAT SHOULD I DO NOW ??? " )
        print ( " NOT EXITING ... YET ... (7) " )
        ## sys.exit(-1)
        ## TO DO NEXT ...

    elif ( len(foundGenes) > 2 ):
        ## one example is: pathogenic:SLC34A2-ROS1;GOPC fusion ... how to interpert something like this???
        print ( "     --> found more than TWO gene symbols ??? ", foundGenes )
        print ( " (c) WHAT SHOULD I DO NOW ??? " )
        print ( " NOT EXITING ... YET ... (8) " )
        ## sys.exit(-1)

    ## ----
    ## if we found any 'known terms' ...
    if ( numKnownTerms > 0 ):
        print ( numKnownTerms, knownTerms )

        for kt in knownTerms:

            ## first we handle TMB 
            ## removed 'biomarker_type':'2'
            if ( kt.lower()=='tmb' or kt.lower().find('burden') >= 0 ):
                print ( " have TMB string ... what next ? ", u )
                if ( u.lower().find('high') >= 0 ):
                    return ( {'biomarker_type':'TMB', 'biomarker_string_value':'high'} )
                elif ( u.lower().find('low') >= 0 ):
                    print ( " TMB low " )
                    return ( {'biomarker_type':'TMB', 'biomarker_string_value':'low'} )
                elif ( u.lower().find('intermed') >= 0 ):
                    print ( " TMB intermediate " )
                    return ( {'biomarker_type':'TMB', 'biomarker_string_value':'intermediate'} )
                elif ( u.lower().find('interm') >= 0 ):
                    return ( {'biomarker_type':'TMB', 'biomarker_string_value':'intermediate'} )
                else:
                    print ( " unclear TMB ??? " )
                    print ( " NOT EXITING ... YET ... (9) " )
                    ## sys.exit(-1)

            ## and then MSI
            ## removed 'biomarker_type':'4'
            elif ( kt.lower()=='msi' or kt.lower().find('micro') >= 0 ):
                print ( " have MSI string ... what next ? ", u )
                if ( u.lower().find('high') >= 0 ):
                    return ( {'biomarker_type':'MSI', 'biomarker_string_value':'high'} )
                elif ( u.lower().find('stable') >= 0 ):
                    return ( {'biomarker_type':'MSI', 'biomarker_string_value':'stable'} )
                elif ( u.lower().find('low') >= 0 ):
                    return ( {'biomarker_type':'MSI', 'biomarker_string_value':'low'} )
                else:
                    ## maybe there's a floating point value we can extract?
                    fVal = -999.
                    ux = u.split()
                    for x in ux:
                        try:
                            fVal = float(x)
                        except:
                            pass
                    if ( fVal >= 0. ): return ( {'biomarker_type':'MSI', 'biomarker_numeric_value':fVal} )

                    print ( " other MSI value ??? " )
                    print ( " NOT EXITING ... YET ... (10) " )
                    ## sys.exit(-1)

            elif ( kt.lower()=='fusion' ):
                print ( " have unknown FUSION string ... what next ? ", u )
                print ( " foundGenes : ", foundGenes )
                ## u = removeGenesFromString ( u, foundGenes )
                u = trimWhiteSpace(u)
                return ( {'biomarker_type':'fusion', 'free_text':u} )

            else:
                print ( " have other string .... what is it ??? ", u )
 
        ## no special handling ...
        print ( " --> have NO special handling for this term: <%s> " % kt )

    ## ----
    ## at this point we just return whatever we might have and anything that is left
    ## as 'free_text' ...
    u = trimWhiteSpace(u)
    if ( len(u) > 1 ): rd['free_text'] = u
    print ( " RETURNING from handleFreeAlterationText ... ", rd )
    return ( rd )


##----------------------------------------------------------------------------------------------------

def findLongestInteger(s):

    if ( s[0] == '*' ):
        print ( " in findLongestInteger ... ignoring leading * ... " )
        s = s[1:]

    ## we assume that the input string looks something like this: 1234abcXYZ
    ## if it looked like this: 12abc34, it would just return 12

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
                print ( "     oh, well, for now just returning what we started with ", s )
                return ( s )
            sys.exit(-1)

    return ( 'c.' + str(ipos) + ref + '>' + var )

##----------------------------------------------------------------------------------------------------

def validateCodingAlt(w):

    ## 12nov2019: this function was bombing on c.*2336_*2339del  ... FIXED

    ## not really sure if * should be considered a valid nucleotide ... ???
    validN = ['A', 'C', 'G', 'T', 'N', '*']
    validD = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '*']

    print ( " in validateCodingAlt ... ", w )

    origW = w

    if ( w == 'c' ): return ( '' )
    if ( w == 'c.' ): return ( '' )

    ## sometimes it starts with 'c.c.' ...
    if ( w.startswith('c.c.') ):
        w = w[2:]

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

    if ( s == 'p' ): return ( '' )
    if ( s == 'p.' ): return ( '' )

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
## leaving most of this code in place, but actually not returning the AF value because
## DavidB says that this information should be obtained elsewhere anyway ...

def extractAF(s):

    cList = [',', '.', ';', ':']

    if ( s.find('AF=') >= 0 ):

        print ( " extracting allele frequency value from ", s )

        ## try:
        if ( 1 ):

            ## first we extract the floating point allele frequency
            js = s.find('AF=')
            t1 = s[js+3:]
            t2 = t1.split()
            t3 = t2[0].strip()
            try:
                afValue = float ( t3 )
                print ( "     >>> extracted floating point value %f from <%s> " % ( afValue, t3 ) )
            except:
                print ( "     >>> failed to extract floating point value from <%s> " % ( t3 ) )
                afValue = 1000.

            if ( 0 ):
                ## now we look at the actual value
                if ( afValue > 999. ):
                    ## we will discard this value
                    print ( "         --> discarding this, assuming it is just a flag of some sort " )
                    pass
                elif ( afValue > 100. or afValue < 0.01 ):
                    ## this is unexpected
                    print ( "         --> UNEXPECTED allele frequency ??? ", afValue )
                elif ( afValue > 1. ):
                    ## AF between 1. and 100 -- assume it is a percentage  
                    print ( "         --> allele frequency assumed to be a PERCENTAGE value ", afValue )
                elif ( afValue < 1. ):
                    ## AF value less than 1 ... assumed to be a fraction
                    afValue *= 100.
                    print ( "         --> multiplied AF by 100 to get: ", afValue )
                elif ( afValue == 1. ):
                    ## AF value is exactly 1 ... what to do?
                    print ( "          --> AF is exactly 1 ??? what to do ??? ", afValue )
                    sys.exit(-1)
                    
                ## create output afD
                if ( afValue >=1. and afValue < 100. ):
                    afD = {'allele_frequency': afValue}
                else:
                    afD = {}

            ## now we want to strip out that part of the string and returning what is left ...
            s0 = s[:js].strip()
            try:
                if ( s0[0] in cList ): s0 = s0[1:]
            except:
                pass
            try:
                if ( s0[-1] in cList ): s0 = s0[:-1]
            except:
                pass

            if ( len(t2) > 1 ):
                for jt in range(1,len(t2)):
                    s0 += ' ' + t2[jt].strip()

            ## what is this doing ???
            try:
                if ( s0[-1] in cList ): s0 = s0[:-1]
                s0.strip()
            except:
                print ( " oops, in except clause!!! ??? ", s0, cList )
                pass

            print ( "     >>> planning to return remaining string : <%s> " % s0 )
            if ( s0.find('AF=') >= 0 ): sys.exit(-1)

            sRest = s0

        ## except:
        ##    print ( "     >>> failed to extract floating point value from <%s> " % s )

        ## return ( sRest, afD )
        return ( s, {} )

    else:
        return ( s, {} )

##----------------------------------------------------------------------------------------------------

def checkHGVS(s, vepFlag, p):

    if ( str(s) == 'nan' ): return ('')
    if ( s == '' ): return ('')

    print ( "\n\n" )
    print ( "===> in checkHGVS function !!! ", type(s), s, type(p), p )

    ## first we need to revert from the string representation of a list
    ## to the underlying list ...

    ## 25-mar-2020 removing this check ... or at least changing to just a warning ...
    if ( len(s) > 16384 ): 
        print ( "SVGH!   ACK ACK ACK ... string is really really LONG ... ", len(s) )
        print ( s[:80], "\n", s[-80:], "\n\n" )
        if ( 1 ):
            print ( " --> used to bail at this point, but now just issuing a WARNING ... " )
        else:
            print ( " --> returning as-is ... " )
            return (s)

    ## converting what was just a string to a list of strings (09-mar-2020)
    sList = []

    t = ast.literal_eval(s)
    print ( "SVGH!         literal_eval produces: ", type(t), len(t), t )
    if ( t == ['None'] ): return ( '' )

    for u in t:

        newS = ''

        print ( "SVGH!             looping within t ... ", type(u), u )
        ## print ( "SVGH! <<<%s>>> " % u )

        ## ONLY PLACE WHERE handleFreeAlterationText is called ...
        uDict = handleFreeAlterationText ( u, vepFlag )

        if ( 'free_text' in uDict ):
            print ( " NOTA BENE !!! wound up with some FREE_TEXT here ... ", uDict['free_text'] )
            ## sys.exit(-1)

        print ( " >>>> back from handleFreeAlterationText ... ", uDict )
        if ( p != '' ):
            print ( "     --> may need to add this additional key:value pair : ", p )
            pSplit = p.split(':')
            if ( len(pSplit) == 2 ):
                pKey = pSplit[0]
                pValue = pSplit[1]
            else:
                pKey = pSplit[0]
                pValue = "True"

            if ( pKey not in uDict ):
                uDict.update ( {pKey:pValue} )
                print ( "         --> now have this : ", uDict )
            else:
                print ( "         --> label already present !!! ", uDict[pKey] )
                if ( uDict[pKey] != pValue ):
                    print ( " WHAT NOW ???  config file specified <%s>:<%s> " % ( pKey, pValue) )

        newS = json.dumps(uDict,sort_keys=True)
        sList += [ newS ]

    ## final list of strings to be returned ...
    print ( " " )
    print ( " FINAL HGVS-parsed & re-built string : ", len(sList) )
    for ii in range(len(sList)):
        print ( "     ", ii, sList[ii] )
    print ( " " )

    ## NOW convert this list of JSON strings to a single STRING ...
    bigString = repr(sList)
    if ( len(bigString) > 48 ):
        print ( "     --> FINAL output (REPR) string ", len(bigString), bigString[:20], " ... ", bigString[-20:] )
    else:
        print ( "     --> FINAL output (REPR) string ", len(bigString), bigString )

    return ( sList )

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
            (fileSize, rowCountEst, chunkSizeRows, numChunks) = estimateFileSize(dataFileName)

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
                                       sep=',',
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
                    if ( "addTerm" in xConfig ):
                        p = xConfig["addTerm"]
                    else:
                        p = ''
                    x = x.apply (checkHGVS, args=("withVEP",p,))
            if ( "HGVS_noVEP" in xConfig ):
                if ( xConfig['HGVS_noVEP'] == "True" ):
                    if ( "addTerm" in xConfig ):
                        p = xConfig["addTerm"]
                    else:
                        p = ''
                    x = x.apply (checkHGVS, args=("noVEP",p,))

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
            print("     --> calling csvIter.get_chunk ... ")
            self.chunkDf = self.csvIter.get_chunk(self.chunkSizeRows)
            print(" --> back from get_chunk()")
            print(self.chunkDf.shape)
            self.ithChunk += 1
            return (True)
        except StopIteration:
            print(" --> StopIteration exception: end of iteration")
            return (False)
        except:
            print(" UHOH other type of error in getNextChunk ???")
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

        ### THERE USED TO BE A BUNCH OF CODE BELOW FOR FURTHER "DIGGING" INTO THE DATAFRAME
        ### TO FURTHER CHARACTERIZE IT, BUT I AM GOING TO GET RID OF ALL OF THAT ...
        ### --> deleted approx 218 lines

    def writeOutChunk(self):

        if ( 0 ):
            ## write output to CSV file ...
            if (self.emptyOutput):
                self.chunkDf.to_csv(self.fhOut, index=False, sep=self.outputSep)
                self.emptyOutput = False
            else:
                self.chunkDf.to_csv(self.fhOut, index=False,
                                    sep=self.outputSep, header=False)

        if ( 1 ):
            myToJson ( self.fhOut, self.chunkDf )

            if ( 0 ):
                ## write output to JSON file ...
                if (self.emptyOutput):
                    print ( " " )
                    print ( " WRITING OUTPUT AS JSON RECORDS " )
                    self.chunkDf.to_json ( self.fhOut, orient='records', lines=True )
                    self.emptyOutput = False
                else:
                    self.chunkDf.to_json ( self.fhOut, orient='records', lines=True )

        

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
    parser.add_argument('-f', '--inputFileName', action='store', help='input CSV data file (NB: this must be a CSV)',
                        required=True, dest='inputFile', type=str)
    parser.add_argument('-c', '--configFileName', action='store', help='input JSON config file',
                        required=True, dest='configFile', type=str)
    parser.add_argument('-o', '--outputFileName', action='store', help='output CSV/JSON data file',
                        required=False, dest='outputFile', type=str)
    parser.add_argument('-s', '--schemaFileName', action='store', help='output JSON schema file',
                        required=False, dest='schemaFile', type=str)
    args = parser.parse_args()

    main(args)

# ----------------------------------------------------------------------------------------------------
