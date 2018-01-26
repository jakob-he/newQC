#-*- coding: utf - 8 - *-
#!/usr/bin/env python3

"""
Created on Monday 25 11:23:04 2017

@author: Tori,Jakob

"""

import datetime as dt
from pathlib import Path
import json
import os
import shutil
import re  # Used to retrieve directoy information
import urllib3  # Used to change HGVS with mutalyzer
import re
import csv
import pandas as pd
import sys
import getopt
import hgvs.parser  # hgvs parsing and validation library
import hgvs.dataproviders.uta
import hgvs.variantmapper
import hgvs.validator
import hgvs.exceptions


# ===============================
# ===== helper functions ========
# ===============================


def purge(dir, pattern):
    """Deletes files in a directory which include a certain pattern"""
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


def copy(file, dir):
    """
    Copies file to the curated folder and changes the type for all entries in geneList to int

    Args:
        file (string) = path to file
        dir (string) = path to directory

    Returns:
        exists (boolean) = False if file is already in the directory
    """

    exists = my_file.is_file(dir + file)
    if exists:
        with open(os.path.join(dir, file), 'r') as json_data:
            data = json.load(json_data)
                for entry in data["geneList"]:
                    if "gene_id" in entry:
                        entry["gene_id"] = int(entry["gene_id"])
                    if "gene_omim_id" in entry:
                        entry["gene_omim_id"] = int(entry["gene_omim_id"])
        with open(dir + file, 'w') as outputfile:
            json.dump(data, outputfile)


def check_hgvs(parser, validator, hgvs, submitterteam, submitter, overview, file):
    """Checks the HGVS code. Prints error and the HGVS code if HGVS is not recognized and saves
    the error information in the overview directory

    Args:
        parser (hgvs.parser.Parser object): Used to parse the HGVS code
        validator (hgvs.validator.Validator object): Used to validate the parsed HGVS code
        hgvs (string): Input that is going to be tested
        submitterteam (string): Team of submitters in the overview dictionary
        submitter (string): Name of submitter
        overview (dictionary): Dictionary containing the information of previously parsed JSONs
        file (string): directory of json-file

    """
    try:
        hgvstest = validator.validate(parser.parse_hgvs_variant(hgvs))
    except Exception:
        if file not in overview[submitterteam]['incorrect JSONs']:
            overview[submitterteam]['incorrect JSONs'][file] = {}
        if 'falscher HGVS-Code' not in overview[submitterteam]['incorrect JSONs'][file]:
            overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'] = {}
            overview[submitterteam]['incorrect JSONs'][file]['submitter'] = submitter
            overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'][hgvs] = "Parsing error"
        else:
            overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'][hgvs] = "Parsing error"

    if hgvstest == False:
        if file not in overview[submitterteam]['incorrect JSONs']:
            overview[submitterteam]['incorrect JSONs'][file] = {}
        if 'falscher HGVS-Code' not in overview[submitterteam]['incorrect JSONs'][file]:
            overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'] = {}
            overview[submitterteam]['incorrect JSONs'][file]['submitter'] = submitter
            overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'][hgvs] = "internally incorrect HGVS"
        else:
            overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'][hgvs] = "internally incorrect HGVS"


def append_incorrect(overview, file, submitterteam, submitter, str):
    """Appends the input json-file to the list of incorrect files with the corresponding error"""

    if file not in overview[submitterteam]['incorrect JSONs']:
        overview[submitterteam]['incorrect JSONs'][file] = {}
        overview[submitterteam]['incorrect JSONs'][file]['submitter'] = submitter
    overview[submitterteam]['incorrect JSONs'][file][str] = True


def newentry(dictionary, submitterteam):
    """Adds a new entry for a submitterteam to a dictionary"""

    dictionary = {}
    dictionary[submitterteam] = {}
    dictionary[submitterteam]['team members'] = {}
    dictionary[submitterteam]['number of cases'] = 1
    dictionary[submitterteam]['VCFs'] = 0
    dictionary[submitterteam]['correct JSONs'] = {}
    dictionary[submitterteam]['correct JSONs']['number of correct jsons'] = 0  # []
    dictionary[submitterteam]['correct JSONs']['list of correct jsons'] = []
    dictionary[submitterteam]['incorrect JSONs'] = {}
    dictionary[submitterteam]['not monogenic'] = []
    return(dictionary)

def savetomultivcf(overview,jsoncurratedfolder):
    jsonlist = []
    jsonlist2 = []
    novcf = []
    withvcf = []

    for file in os.listdir(jsoncurratedfolder+"/"):
        if file endswith(".json"):
            with open(jsoncurratedfolder + '/' + file) as json_data:
                jsondata = json.load(json_data)
                        for submitter in overview.keys():
                            # and jsondata['vcf']=='noVCF':
                            if file in overview[submitter]['correct JSONs']['list of correct jsons']:
                                jsonlist.append(file)
                                if jsondata['vcf'] == 'noVCF':
                                    novcf.append(file)
                                else:
                                    withvcf.append(file)

    # multiVCF erstellen
    multivcf = pd.DataFrame(columns=[
                            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NM'])
    vcfcounter = 0
    x = 0


    for file in jsonlist:
        with open(jsoncurratedfolder + '/' + file) as json_data:
            d = json.load(json_data)
            caseID = d['case_id']
            hgvslist = []
            for mutation in d['genomicData']:
                if type(mutation)==dict:
                    if 'HGVS-code' in mutation['Mutations'].keys():
                        hgvs = mutation['Mutations']['HGVS-code']
                        hgvslist.append(hgvs)
                    elif 'Mutation 1' in mutation['Mutations'].keys():
                        for mutationnr, mutationdict in mutation['Mutations'].items():
                            if 'HGVS-code' in mutationdict.keys():
                                hgvs = mutationdict['HGVS-code']
                                hgvslist.append(hgvs)
                    else:
                        print 'no Mutation: ', file
                        continue
                    genotype = mutation['Test Information']['Genotype']
                    if genotype == 'Hemizygous':
                        genotype = '1'
                    elif genotype == 'Homozygous':
                        genotype = '1/1'
                    elif genotype == 'Heterozygous' or genotype == 'Compound Heterozygous':
                        genotype = '0/1'
                    else:
                        genotype = './1'
                    for hgvscode in hgvslist:
                        try:
                            hgvsparser = hgvs.parser.Parser()
                            parsedhgvs = hgvsparser.parse_hgvs_variant(hgvscode)
                            chrom =
                            chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(
                                str(hgvscode), genome, get_transcript=get_transcript)
                            if hgvscode in list(multivcf['NM']):
                                index = list(multivcf['NM']).index(hgvscode)
                                multivcf.set_value(index, caseID, genotype)
                                if caseID not in jsonlist2:
                                    jsonlist2.append(caseID)
                            else:
                                chromo = chrom.split('chr')[1]
                                multivcf.set_value(x, '#CHROM', str(chromo))
                                try:
                                    multivcf.set_value(x, 'sort', int(chromo))
                                except ValueError, e:
                                    multivcf.set_value(x, 'sort', 30)
                                multivcf.set_value(x, 'NM', hgvscode)
                                multivcf.set_value(x, 'POS', offset)
                                multivcf.set_value(x, 'ID', '.')
                                multivcf.set_value(x, 'REF', str(ref))
                                multivcf.set_value(x, 'ALT', str(alt))
                                multivcf.set_value(x, 'QUAL', '.')
                                multivcf.set_value(x, 'FILTER', '.')
                                multivcf.set_value(
                                    x, 'INFO', 'HGVS="' + hgvscode + '"')
                                multivcf.set_value(x, 'FORMAT', 'GT')
                                multivcf.set_value(x, caseID, genotype)
                                x = x + 1
                                if caseID not in jsonlist2:
                                    jsonlist2.append(caseID)
                        except ValueError, e:  # 'falsche' HGVS-Codes überspringen und anzeigen
                            print 'Error:', file, hgvs, e
                            continue

    # data_vcf sortieren
    print 'Sort DataFrame ...'
    multivcf = multivcf.sort_values(by=['POS'])
    multivcf = multivcf.reset_index(drop=True)

    # leere Felder füllen
    multivcf = multivcf.fillna(value='0/0')

    multivcf.to_csv(mVCF + ".tmp", sep='\t', index=False,
                    header=True, quoting=csv.QUOTE_NONE)

    with open(mVCF, 'w') as outfile:
        outfile.write(
            '##fileformat=VCFv4.1\n##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        with open(mVCF + ".tmp", 'r') as infile:
            for line in infile:
                outfile.write(line)

    os.remove(mVCF + ".tmp")

    # Liste der Cases im MultiVCF
    with open(configfile, 'w') as jsoninmultivcf:
        jsoninmultivcf.write('SINGLE_SAMPLES:\n')
        for item in novcf:
            item = item.split('.')[0]
            jsoninmultivcf.write(" - %s\n" % item)
        jsoninmultivcf.write('VCF_SAMPLES:\n')
        for item in withvcf:
            item = item.split('.')[0]
            jsoninmultivcf.write(" - %s\n" % item)


# ===============================
# ===== Inspection function =====
# ===============================

def checkjson(file, overview, checked = False):
    """Checks the json-file for missing entries and appends the
    information to a dictionary. Either the dictionary is already present
    or newly created.

    Args:
        file (string): Directory of the JSON-File
        overview (dictionary): Dictionary containing the information of previously parsed JSONs

    Returns:
        overview (dictionary): Changed dictionary
    """
    if file.endswith('.json'):
        with open(file) as json_data:
            d = json.load(json_data)

            #delete old entry for previously checked files
            if checked == True:
                try:
                    del d["submitter"]["user_team"]["incorrect JSONs"][file]
                except KeyError:
                    del d["submitter"]["user_name"]["incorrect JSONs"][file]

            # Catch Json without a submitter
            try:
                submitterteam = d['submitter']['user_team']
            except KeyError:
                return("No Submitter!")

            if submitterteam is None:
                submitterteam = d['submitter']['user_name']

            # Get submittername
            submitter = d['submitter']['user_name']

            # Define overview dictionary
            if checked == False:
                if submitterteam not in overview:
                    newentry(overview, submitterteam)
                else:
                    overview[submitterteam]['number of cases'] += 1

            # Check if team member is in the submitterteam
            if submitter not in overview[submitterteam]['team members']:
                overview[submitterteam]['team members'][submitter] = d['submitter']['user_email']

            # Check if there is a VCF
            if checked == False:
                vcfexists = d['vcf'] != 'noVCF'
                if vcfexists:
                    overview[submitterteam]['VCFs'] += 1

            # Get highest gestalt score
            if isinstance(d['geneList'], list):
                gscore_list = []
                for gene in d['geneList']:
                    gscore = gene['gestalt_score']
                    gscore_list.append(float(gscore))
            if not gscore_list:
                max_gscore = max(gscore_list)
            else:
                max_gscore = 0
                append_incorrect(overview, file, submitterteam,
                                 submitter, 'No gene list!')

            # Check for annotated features
            if not d['features']:
                append_incorrect(overview, file, submitterteam,
                                 submitter, 'No Features')
            # No image uploaded
            if max_gscore == 0:
                append_incorrect(overview, file, submitterteam,
                                 submitter, 'No Picture')

            # Check for the molecular diagnosis
            if d['selected_syndromes'] == 'Non selected':
                append_incorrect(overview, file, submitterteam,
                                 submitter, 'No Diagnosis listed')

            # HGVS-code check

            # Initialize HGVS parser and validator
            parser = hgvs.parser.Parser()
            validator = hgvs.validator.IntrinsicValidator(parser)

            if not d['genomicData']:
                append_incorrect(overview, file, submitterteam,
                                 submitter, 'No Mutations listed')
            for mutation in d['genomicData']:
                if not mutation['Mutations']:
                    append_incorrect(overview, file, submitterteam,
                                     submitter, 'One mutation not listed')
                    continue
                if 'HGVS-code' in mutation['Mutations']:
                    hgvs = mutation['Mutations']['HGVS-code']
                    check_hgvs(parser, validator, hgvs,
                               submitterteam, submitter, overview, file)

                # Check for compound heterozygous
                elif 'Mutation 1' in mutation['Mutations']:
                    hgvslist = []
                    for mutationnr, mutationdict in mutation['Mutations'].items():
                        if 'HGVS-code' in mutationdict:
                            hgvs = mutationdict['HGVS-code']
                            hgvslist.append(hgvs)
                            check_hgvs(parser, validator, hgvs, submitterteam,
                                       submitter, overview, file)

                # Check for monogenic
                elif mutation['Test Information']['Mutation Type'] == 'Monogenic':
                    if 'protein level' in mutation['Mutations']:
                        append_incorrect(overview, file, submitterteam,
                                         submitter, 'Mutation on Proteinlevel')
                    else:
                        if not vcfexists:
                            append_incorrect(overview, file, submitterteam,
                                             submitter, 'No HGVS-code listed')
                        else:
                            append_incorrect(overview, file, submitterteam,
                                             submitter, 'VCF listed but no HGVS-Code')
                else:
                    append_incorrect(
                        overview, file, submitterteam, submitter, 'No HGVS-Code listed (2)')
                if mutation['Test Information']['Mutation Type'] != 'Monogenic':
                    overview[submitterteam]['not monogenic'].append(file)

            # Check if the file is correct and append it to the list of correct jsons
            if file not in overview[submitterteam]['incorrect JSONs']:
                overview[submitterteam]['correct JSONs']['number of correct jsons'] += 1
                overview[submitterteam]['correct JSONs']['list of correct jsons'].append(
                    file)

        json.dump(overview, maindict)
        close(maindict)

    return overview


# ===============================
# == errorcorrection functions===
# ===============================

def errorcorrbydict(overview, errordict):
    """
    Checks every incorrect json file in the overview dictionary and tries to correct the HGVS-Codes with the errordictionary

    Args:
        overview (dictionary): dictionary containing the information of previously parsed JSONs
        errordict (string): path to errordictionary

    Returns:
        overview (dictionary): dictionary with corrected entries
        dicterrors (int): Value representing the number HGVS-Codes corrected with the errordict
    """

    dicterrors = 0
    with open(errordict) as json_data:
        transcript_errors = json.load(json_data)

        for submitter in overview:
            for jsonname in overview[submitter]['incorrect JSONs']:
                if unicode('falscher HGVS-Code') in overview[submitter]['incorrect JSONs'][jsonname]:
                    with open(jsoncurratedfolder + '/' + jsonname, 'r') as jsonfile:
                        jsondata = json.load(jsonfile)
                        for hgvs, error in overview[submitter]['incorrect JSONs'][jsonname]['falscher HGVS-Code'].items():
                            # hgvs=str(hgvs.replace(" ", "")) #Unicode-Fehlermeldung entgehen
                            if hgvs in transcript_errors:
                                for mutation in jsondata['genomicData']:
                                    if 'HGVS-code' in mutation['Mutations']:
                                        if mutation['Mutations']['HGVS-code'] == hgvs:
                                            mutation['Mutations']['HGVS-code'] = transcript_errors[hgvs]
                                            dicterrors = dicterrors + 1
                                    elif 'Mutation 1' in mutation['Mutations']:
                                        for multimut, description in mutation['Mutations'].items():
                                            if mutation['Mutations'][multimut]['HGVS-code'] == hgvs:
                                                mutation['Mutations'][multimut]['HGVS-code'] = transcript_errors[hgvs]
                                                dicterrors = dicterrors + 1
                        with open(jsoncurratedfolder + '/' + jsonname, 'w') as dicttojson:
                            json.dump(jsondata, dicttojson)
                        overview = checkjson(
                            jsoncurratedfolder + '/' + jsonname, overview, checked = True)
    return(overview, dicterrors)

    def errorcorrbymut(overview, errordict):
        """
        Tries to correct the HGVS-Code for every incorrect json file in the overview dictionary with the mutalyzer

        Args:
            overview (dictionary): dictionary containing the information of previously parsed JSONs
            errordict (string): path to errordictionary

        Returns:
            overview (dictionary): dictionary with corrected entries
            muterrors (int): Value representing the number corrected HGVS-Codes
        """
        muterrors = 0
        with open(errordict, 'r') as dicttojson:
            transcript_errors = json.load(dicttojson)

        new_hgvs = []
        corrected_transcripts = {}
        for submitter in overview:
            for jsonfile in overview[submitter]['incorrect JSONs']:
                if unicode('falscher HGVS-Code') in overview[submitter]['incorrect JSONs'][jsonfile]:
                    for key in overview[submitter]['incorrect JSONs'][jsonfile]['falscher HGVS-Code']:
                        if 'required' in overview[submitter]['incorrect JSONs'][jsonfile]['falscher HGVS-Code'][key]:
                            url = "https://mutalyzer.nl/position-converter?assembly_name_or_alias=GRCh37&description=" + key
                            try:
                                page = urllib2.urlopen(url)
                                data = page.read()
                            except Exceptio:
                                print 'Could not connect: ', jsonfile, key
                                data = 'empty'
                            if 'Found transcripts' in data:
                                data = data.split(
                                    'variant region</h4>')[1].split('<br></pre>')[0]
                                transcript = key.split('.')[0]
                                transcriptversion = key.split(
                                    ':')[0].split('.')[1]
                                positions = [m.start()
                                             for m in re.finditer(transcript, data)]
                                alt_transcripts = []
                                for position in positions:
                                    trans = data[position:].split(':')[0]
                                    trans_vers = trans.split(transcript)[1].split(':')[
                                        0].split('.')[1]
                                    if trans_vers > transcriptversion:
                                        transcriptversion = trans_vers
                                transcriptversion = transcript + '.' + \
                                    transcriptversion + ':' + key.split(':')[1]
                                with open(jsoncurratedfolder + '/' + jsonfile) as json_data:
                                    d = json.load(json_data)
                                    for mutation in d['genomicData']:
                                        if 'HGVS-code' in mutation['Mutations']:
                                            if mutation['Mutations']['HGVS-code'] == key:
                                                mutation['Mutations']['HGVS-code'] = transcriptversion
                                                # append_errordict()
                                            elif 'Mutation 1' in mutation['Mutations']:
                                                for multimut, description in mutation['Mutations'].items():
                                                    if mutation['Mutations'][multimut]['HGVS-code'] == key:
                                                        mutation['Mutations']['HGVS-code'] = transcriptversion
                                                        # append_errordict
                                transcript_errors[key] = transcriptversion
                                muterrors = muterrors + 1
                                new_hgvs.append(transcriptversion)
                                print jsonfile, 'vorher: ', key, 'nachher: ', transcriptversion
                                corrected_transcripts[key] = transcriptversion
                            elif 'We found these versions' in data:
                                # print jsonfile, data
                                newtranscript = data.split('We found these versions: ')[
                                    1].split('<p></p>')[0].split('</p>')[0]
                                newtranscript = newtranscript + \
                                    ':' + key.split(':')[1]
                                with open(jsoncurratedfolder + '/' + jsonfile) as json_data:
                                    d = json.load(json_data)
                                    for mutation in d['genomicData']:
                                        if 'HGVS-code' in mutation['Mutations']:
                                            if mutation['Mutations']['HGVS-code'] == key:
                                                mutation['Mutations']['HGVS-code'] = newtranscript
                                                # append_errordict()
                                            elif 'Mutation 1' in mutation['Mutations']:
                                                for multimut, description in mutation['Mutations'].items():
                                                    if mutation['Mutations'][multimut]['HGVS-code'] == key:
                                                        mutation['Mutations']['HGVS-code'] = newtranscript
                                                        # append_errordict()
                                transcript_errors[key] = newtranscript
                                muterrors = muterrors + 1
                                print jsonfile,  'vorher: ', key, 'nachher: ', newtranscript
                                new_hgvs.append(newtranscript)
                                corrected_transcripts[key] = newtranscript
                            elif 'could not be found in our database (or is not a transcript).' in data:
                                print 'no transcript found: ', submitter, jsonfile, key
                            else:
                                print 'please check: ', jsonfile, key

            with open(errordict, 'w') as dicttojson:
                json.dump(transcript_errors, dicttojson)
    return(overview, muterrors)

    def errorcorrmanual(overview, errordict):
        """
        Iterates through every incorrect json file of the overview dictionary and trie to correct the HGVS-Codes with user input

        Args:
            overview (dictionary): dictionary containing the information of previously parsed JSONs
            errordict (string): path to errordictionary

        Returns:
            overview (dictionary): dictionary with corrected entries
            manerros (int): Value representing the number manually corrected HGVS-Codes
        """
        manerrors = 0
        not_solved = {}

        for submitter in overview:
            for jsonname in overview[submitter]['incorrect JSONs']:
                if unicode('falscher HGVS-Code') in overview[submitter]['incorrect JSONs'][jsonname]:
                    with open(jsoncurratedfolder + '/' + jsonname, 'r') as jsonfile:
                        jsondata = json.load(jsonfile)
                        for hgvs, error in overview[submitter]['incorrect JSONs'][jsonname]['falscher HGVS-Code'].items():
                            print hgvs, error
                            for mutation in jsondata['genomicData']:
                                if 'HGVS-code' in mutation['Mutations']:
                                    if mutation['Mutations']['HGVS-code'] == hgvs:
                                        print '\n\n', jsonname, '\n', hgvs, error, '\n', mutation.items()
                                        right_hgvs = raw_input('right HGVS:')
                                        if right_hgvs == 'n':
                                            not_solved[jsonname] = hgvs
                                            continue
                                        else:
                                            mutation['Mutations']['HGVS-code'] = right_hgvs
                                            transcript_errors[hgvs] = right_hgvs
                                            manerrors = manerrors + 1
                                        with open(errordict, 'w') as dicttojson:
                                            json.dump(
                                                transcript_errors, dicttojson)
                                elif 'Mutation 1' in mutation['Mutations']:
                                    for multimut, description in mutation['Mutations'].items():
                                        if mutation['Mutations'][multimut]['HGVS-code'] == hgvs:
                                            print '\n\n', jsonname, '\n', hgvs, error, mutation.items()
                                            right_hgvs = raw_input(
                                                'right HGVS:')
                                            if right_hgvs == 'n':
                                                not_solved[jsonname] = hgvs
                                                continue
                                            else:
                                                mutation['Mutations'][multimut]['HGVS-code'] = right_hgvs
                                                transcript_errors[hgvs] = right_hgvs
                                                manerrors = manerrors + 1
                                            with open(errordict, 'w') as dicttojson:
                                                json.dump(
                                                    transcript_errors, dicttojson)
                            with open(jsoncurratedfolder + '/' + jsonname, 'w') as dicttojson:
                                json.dump(jsondata, dicttojson)
                            overview = checkjson(
                                jsoncurratedfolder + '/' + jsonname, overview, checked = true)
    return(overview, manerrors)


# ===============================
# ===== main script =============
# ===============================

def main():
    """
    main function with the following steps:
    1. Read all inputs
    2. Copy the input json file to the currated folder if not yet present
    3. Depending on user input execute error correction for all incorrect json files in overview
    4. If wanted by the user drop the information in overview to a vcf file
    5. Rearrange overview
    6. Save log file
    """
    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, "h::", ["help", "jsonfile =", "overview =", "mappedjsons=",
                                                 "log=", "debugfolder=", "jsoncurrated=", "errordict=", "config=", "vcf="])
    except getopt.GetoptError as e:
        print(e)
        print('JsonQC.py --jsonsoriginal --log ')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('jsonToTable.py -i <input-folder> -o <output-file.tsv>')
            sys.exit(1)
        elif opt in ("--jsonfile"):
            jsonfile = arg
        elif opt in ("--overview"):
            overviewpath = overview
        elif opt in ("--mappedjsons"):
            jsonsoriginal = arg
        elif opt in ("--log"):
            logfile = arg
        elif opt in ("--debugfolder"):
            debugfolder = arg
        elif opt in ("--jsoncurrated"):
            jsoncurratedfolder = arg
        elif opt in ("--errordict"):
            errordict = arg
        elif opt in ("--vcf"):
            mVCF = arg
        elif opt in ("--config"):
            configfile = arg

    print('Downloadfolder jsons:', jsonsoriginal)
    print('Debug folder:', debugfolder)
    print('JSON currated folder:', jsoncurratedfolder)
    print('Errordict:', errordict)
    print('ConfigFile:', configfile)
    print('Multi-VCF:', mVCF)

    # Create debugfolder if it does not exist
    if not os.path.exists(debugfolder):
        os.makedirs(debugfolder)

    # Initialize the time variables 
    now = dt.datetime.now()
    date = now.strftime('%Y-%m-%d')
    time = now.strftime('%H:%M:%S')

    # copy json to currated
    copy(jsonfile, jsoncurrated)

    # check for the overview file
    try:
        maindict = open(overviewpath, "a")
        overview = json.load(maindict)
    except IOError:
        overview = {}

    # append the information of the new json file to overview
    overview = checkjson(jsonfile, overview)

    # collect information for the logfile
    jsons = 0
    vcfs = 0
    submitter_count = 0
    cor_jsons = 0

    for submitter in overview:
        submitter_count = submitter_count + 1
        jsons = jsons + int(overview[submitter]['number of cases'])
        vcfs = vcfs + int(overview[submitter]['VCFs'])
        cor_jsons = cor_jsons + \
            int(overview[submitter]['correct JSONs']
                ['number of correct jsons'])

    # get user input for errorcorrection
    errorcorrbydict = input('Correct JSONs from errordict? (y OR n)')
    if(errorcorrbydict == "y"):
        overview, dicterrors = errorcorrbydict(overview, errordict)
    errorcorrbymut = input('Correct JSONs with mutalyzer? (y OR n)')
    if(errorcorrbymut == "y"):
        overview, muterrors = errorcorrbydḿut(overview, errordict)
    errorcorrmanual = input('Correct HGVS-Codes manually? (y OR n)')
    if(errorcorrmanual == "y"):
        overview, manerrors = errorcorrmanual(overview, errordict)


    # korrigierte JSONs in MultiVCF / dump corrected JSONs to multi-VCF
    savetovcf = input('Dump to multiVCF? (y OR n) ')
    if(savetovcf == "y"):
        savetomultivcf()

    # Ordnen der falschen JSONs nach Submitter / arrange incorrect JSONs per submitter
    wantto = raw_input('Rearrange result file? (y OR n) ')

    if wantto != 'n':
        submitter_cases = {}
        with open(debugfolder + '/result_' + date + '.json') as json_data:
            result = json.load(json_data)
            for submitterteam in result.keys():
                submitter_cases[submitterteam] = {}
                for submitter in result[submitterteam]['team members'].keys():
                    submitter_cases[submitterteam][submitter] = {}
                for jsonfile in result[submitterteam]['incorrect JSONs'].keys():
                    submitter = result[submitterteam]['incorrect JSONs'][jsonfile]['submitter']
                    indeldups = ['ins', 'del', 'dup']
                    submitter_cases[submitterteam][submitter][jsonfile] = {}
                    submitter_cases[submitterteam][submitter][jsonfile].update(
                        result[submitterteam]['incorrect JSONs'][jsonfile])
                    del submitter_cases[submitterteam][submitter][jsonfile]['submitter']
                    for error in result[submitterteam]['incorrect JSONs'][jsonfile].keys():
                        if 'falscher HGVS-Code' in error:
                            for key in result[submitterteam]['incorrect JSONs'][jsonfile]['falscher HGVS-Code'].keys():
                                if any(indeldup in key for indeldup in indeldups):
                                    del submitter_cases[submitterteam][submitter][jsonfile]['falscher HGVS-Code']
                    if len(submitter_cases[submitterteam][submitter][jsonfile]) == 0:
                        del submitter_cases[submitterteam][submitter][jsonfile]

        with open(debugfolder + '/resultpersubmitter' + date + '.json', 'w') as dicttojson:
            json.dump(submitter_cases, dicttojson)

    ## Dokumentation in Tabellenform

    qc_jsons = 0
    qc_vcfs = 0

    with open(debugfolder + '/result_' + date + '.json') as jsondata:
        result = json.load(jsondata)
        for submitter in result.keys():
            qc_jsons = qc_jsons + \
                int(result[submitter]['correct JSONs']['number of correct jsons'])
            for jsonfile in result[submitter]['correct JSONs']['list of correct jsons']:
                with open(jsoncurratedfolder + '/' + jsonfile) as json_data:
                    vcfin = json.load(json_data)
                    if vcfin['vcf'] != 'noVCF':
                        qc_vcfs = qc_vcfs + 1


    progress = pd.read_excel(open(logfile, 'rb'), sheetname=0)
    newest = [date, time, submitter_count, jsons, vcfs, cor_jsons,
              dicterrors, muterrors, manerror, qc_jsons, qc_vcfs]
    progress.loc[len(progress)] = newest
    progress.to_excel(logfile, index=False)



if __name__ == "__main__":
    main()
