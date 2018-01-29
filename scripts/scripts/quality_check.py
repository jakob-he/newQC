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
import hgvs
import hgvs.parser  # hgvs parsing and validation library
import hgvs.dataproviders.uta
import hgvs.assemblymapper
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


def parsehgvs(hgvs, parser):
    parsedhgvs = parser.parse_hgvs_variant(hgvs)
    chrom = str(int(parsedhgvs.ac.split(".")[0][-2:]))
    offset = parsedhgvs.posedit.pos.start.offset
    ref = parsedhgvs.posedit.ref
    alt = parsedhgvs.posedit.alt
    return(chrom, offset, ref, alt)


def copy(file, dir):
    """
    Copies file to the curated folder and changes the type for all entries in geneList to int

    Args:
        file (string) = path to file
        dir (string) = path to directory

    Returns:file
        exists (boolean) = False if file is already in the directory
    """
    jsonname= os.path.basename(os.path.normpath(file))
    present = Path(dir + "/" + jsonname).is_file()
    print("copying...",jsonname)
    if not present:
        with open(file, 'r') as json_data:
            data = json.load(json_data)
            for entry in data["geneList"]:
                if "gene_id" in entry:
                    entry["gene_id"] = int(entry["gene_id"])
                if "gene_omim_id" in entry:
                    entry["gene_omim_id"] = int(entry["gene_omim_id"])
        with open(dir + "/" + jsonname, 'w+') as outputfile:
            json.dump(data, outputfile)
    return(present)


def check_hgvs(parser, validator, hdp, hgvscode, jsonfile):
    """Checks the HGVS code. Prints error and the HGVS code if HGVS is not recognized and saves
    the error information in the overview directory

    Args:
        parser (hgvs.parser.Parser object): Used to parse the HGVS code
        validator (hgvs.validator.Validator object): Used to validate the parsed HGVS code
        hgvs (string): Input that is going to be tested

    Returns:
        correct (boolean): True if HGVS-Code is correct
    """
    try:
        variant = parser.parse_hgvs_variant(hgvscode)
    except hgvs.exceptions.HGVSError as e:
        if "char 1: expected a letter or digit" in e:
            return(errorcorrbymut(hgvscode,jsonfile))

    vm = hgvs.assemblymapper.AssemblyMapper(
        hdp, assembly_name='GRCh37', alt_aln_method='splign')
    var_g = vm.c_to_g(variant)
    try:
        correct = validator.validate(hp.parse_hgvs_variant(str(variant)))
    except Exception:
        return(False)

    return(True)


def savetomultivcf(jsonfile):
    jsonlist = []
    jsonlist2 = []
    novcf = []
    withvcf = []

    # create multivcf
    multivcf = pd.DataFrame(columns=[
                            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NM'])
    vcfcounter = 0
    x = 0

    with open(jsonfile) as json_data:
        d = json.load(json_data)
        if d['vcf'] == 'noVCF':
            novcf.append(file)
        else:
            withvcf.append(file)
        caseID = d['case_id']
        hgvslist = []
        for mutation in d['genomicData']:
            if type(mutation) == dict:
                if 'HGVS-code' in mutation['Mutations']:
                    hgvs = mutation['Mutations']['HGVS-code']
                    hgvslist.append(hgvs)
                elif 'Mutation 1' in mutation['Mutations']:
                    for mutationnr, mutationdict in mutation['Mutations'].items():
                        if 'HGVS-code' in mutationdict:
                            hgvs = mutationdict['HGVS-code']
                            hgvslist.append(hgvs)
                genotype = mutation['Test Information']['Genotype']
                if genotype == 'Hemizygous':
                    genotype = '1'
                elif genotype == 'Homozygous':
                    genotype = '1/1'
                elif genotype == 'Heterozygous' or genotype == 'Compound Heterozygous':
                    genotype = '0/1'
                else:
                    genotype = './1'
                hp = hgvs.parser.Parser()
                for hgvscode in hgvslist:
                    try:
                        chrom, offset, ref, alt = parsehgvs(hgvscode, hp)
                        if hgvscode in list(multivcf['NM']):
                            index = list(multivcf['NM']).index(hgvscode)
                            multivcf.set_value(index, caseID, genotype)
                            if caseID not in jsonlist2:
                                jsonlist2.append(caseID)
                        else:
                            multivcf.set_value(x, '#CHROM', chrom)
                            try:
                                multivcf.set_value(x, 'sort', int(chromo))
                            except ValueError:
                                multivcf.set_value(x, 'sort', 30)
                                multivcf.set_value(x, 'NM', hgvscode)
                                multivcf.set_value(x, 'POS', offset)
                                multivcf.set_value(x, 'ID', '.')
                                multivcf.set_value(x, 'REF', ref)
                                multivcf.set_value(x, 'ALT', alt)
                                multivcf.set_value(x, 'QUAL', '.')
                                multivcf.set_value(x, 'FILTER', '.')
                                multivcf.set_value(
                                    x, 'INFO', 'HGVS="' + hgvscode + '"')
                                multivcf.set_value(x, 'FORMAT', 'GT')
                                multivcf.set_value(x, caseID, genotype)
                                x=x + 1
                                if caseID not in jsonlist2:
                                    jsonlist2.append(caseID)
                    except ValueError:
                        print('Error:', file, hgvs)
                        continue

    # data_vcf sortieren
    print('Sort DataFrame ...')
    multivcf=multivcf.sort_values(by=['POS'])
    multivcf=multivcf.reset_index(drop=True)

    # leere Felder füllen
    multivcf=multivcf.fillna(value='0/0')

    with open(jsonfile + "vcf", 'w') as outfile:
        outfile.write(
            '##fileformat=VCFv4.1\n##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

    multivcf.to_csv(jsonfile + "vcf", sep='\t', index=False,
                    header=True, quoting=csv.QUOTE_NONE)



# ===============================
# ===== Inspection function =====
# ===============================

def checkjson(file, hdp):
    """Checks the json-file for missing values and the HGVS-Code.

    Args:
        file (string): directory of the JSON-File

    Returns:
        correct (boolean): True if json-file is correct
        error (string): error description if json file is incorrect
    """
    if file.endswith('.json'):
        with open(file) as json_data:
            d=json.load(json_data)

            # Catch Json without a submitter
            try:
                submitterteam=d['submitter']['team']
            except KeyError:
                return(False, "no Submitter!")

            if submitterteam is None:
                submitterteam=d['submitter']['name']

            # Get submittername
            submitter=d['submitter']['name']

            # Check if there is a VCF
            vcfexists=d['vcf'] != 'noVCF'

            # Get highest gestalt score
            if isinstance(d['geneList'], list):
                gscore_list=[]
                for gene in d['geneList']:
                    if "gestalt_score" in gene:
                        gscore=gene['gestalt_score']
                        gscore_list.append(float(gscore))
            if gscore_list:
                max_gscore=max(gscore_list)
            else:
                return(False, "no geneList")

            # Check for annotated features
            if not d['features']:
                return(False, "no features")
            # No image uploaded
            if max_gscore == 0:
                return(False, "no picture")

            # Check for the molecular diagnosis
            if not d["ranks"]:
                return(False, "no selected syndrome")

            # HGVS-code check

            # Initialize HGVS parser and validator
            parser=hgvs.parser.Parser()
            validator=hgvs.validator.Validator(hdp=hdp)

            if not d['genomicData']:
                return(False, "no genomicData")
            for mutation in d['genomicData']:
                if not mutation['Mutations']:
                    print("one Mutation not listed")
                    continue
                if 'HGVS-code' in mutation['Mutations']:
                    hgvscode=mutation['Mutations']['HGVS-code']
                    hgvscheck,errorstring=check_hgvs(parser, validator, hdp, hgvscode, jsonfile)
                    if not(hgvscheck):
                        return(hgvscheck,errorstring)

                # Check for compound heterozygous
                elif 'Mutation 1' in mutation['Mutations']:
                    hgvslist=[]
                    for mutationnr, mutationdict in mutation['Mutations'].items():
                        if 'HGVS-code' in mutationdict:
                            hgvscode=mutationdict['HGVS-code']
                            hgvslist.append(hgvscode)
                            if not (check_hgvs(parser, validator, hdp, hgvscode, jsonfile)):
                                return(False, "wrong HGVS-code")

                # Check for monogenic
                elif mutation['Test Information']['Mutation Type'] == 'Monogenic':
                    if 'protein level' in mutation['Mutations']:
                        return(False, "Mutation on Protein level")
                    else:
                        if not vcfexists:
                            return(False, "HGVS-Code not listed")
                        else:
                            return(False, "VCF present but no HGVS-Code")
                else:
                    return(False, "HGVS-Code not listed 2")
                if mutation['Test Information']['Mutation Type'] != 'Monogenic':
                    return(False, "not monogenic")

    return(True, "no error")


# ===============================
# == errorcorrection functions===
# ===============================

def errorcorrbydict(jsonfile, errordict, hdp):
    """
    Checks if json file can be corrected with errordictioanry

    Args:
        jsonfile (string): path to json file
        errordict (string): path to errordictionary

    Returns:
        correct (boolean): True if json-file has been corrected
    """

    with open(errordict) as json_data:
        transcript_errors=json.load(json_data)

        with open(jsonfile) as jsonfile:
            jsondata=json.load(jsonfile)
            for mutation in jsondata['genomicData']:
                if 'HGVS-code' in mutation['Mutations']:
                    if mutation['Mutations']['HGVS-code'] in transcript_errors:
                        mutation['Mutations']['HGVS-code']=transcript_errors[mutation['Mutations']['HGVS-code']]
                elif 'Mutation 1' in mutation['Mutations']:
                    for multimut, description in mutation['Mutations'].items():
                        if mutation['Mutations']['HGVS-code'] in transcript_errors:
                            mutation['Mutations']['HGVS-code']=transcript_errors[mutation['Mutations']['HGVS-code']]
        correct=checkjson(jsonfile, hdp)
    return(correct)


    def errorcorrbymut(hgvscode,jsonfile):
        """
        Tries to correct hgvs code without transcript with the mutalyzer

        Args:
            hgvscode (string): hgvscode with missing transcript

        Returns:
            correct(boolean): True if json file was corrected
            hgvscode(string): containing either the old hgvscode (if not corrected) or the corrected version
        """

        with open(jsonfile) as json_file:
            jsondata = json.load(json_file)
            url = "https://mutalyzer.nl/position-converter?assembly_name_or_alias=GRCh37&description=" + hgvscode
            try:
                page = urllib2.urlopen(url)
                data = page.read()
            except Exception:
                return(False)
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
                for mutation in d['genomicData']:
                    if 'HGVS-code' in mutation['Mutations']:
                        if mutation['Mutations']['HGVS-code'] == hgvscode:
                            mutation['Mutations']['HGVS-code'] = transcriptversion
                    elif 'Mutation 1' in mutation['Mutations'].keys():
                        for multimut, description in mutation['Mutations'].items():
                            if mutation['Mutations'][multimut]['HGVS-code'] == hgvscode:
                                mutation['Mutations']['HGVS-code'] = transcriptversion
                json.dump(jsondata,json_file)#TODO conitnue working on newhgvs + mutalyzer
                checkjson(jsonfile)
            elif 'We found these versions' in data:
                # print jsonfile, data
                newtranscript = data.split('We found these versions: ')[
                    1].split('<p></p>')[0].split('</p>')[0]
                newtranscript = newtranscript + \
                    ':' + key.split(':')[1]
                with open(jsoncurratedfolder + '/' + jsonfile) as json_data:
                    d = json.load(json_data)
                    for mutation in d['genomicData']:
                        if 'HGVS-code' in mutation['Mutations'].keys():
                            if mutation['Mutations']['HGVS-code'] == key:
                                mutation['Mutations']['HGVS-code'] = newtranscript
                                # append_errordict()
                            elif 'Mutation 1' in mutation['Mutations'].keys():
                                for multimut, description in mutation['Mutations'].items():
                                    if mutation['Mutations'][multimut]['HGVS-code'] == key:
                                        mutation['Mutations']['HGVS-code'] = newtranscript
                                        # append_errordict()
                transcript_errors[key] = newtranscript
                muterrors = muterrors + 1
                print("{} vorher: {} nachher: {}".format(jsonfile, key, newstranscript))
                new_hgvs.append(newtranscript)
                corrected_transcripts[key] = newtranscript
            elif 'could not be found in our database (or is not a transcript).' in data:
                print('no transcript found: {} {} {}'.format(submitter, jsonfile, key))
            else:
                print('please check: {}'.format(jsonfile, key))

            with open(errordict, 'w') as dicttojson:
                json.dump(transcript_errors, dicttojson)

            step = 'mutalyzer'
            overview = checkjsons(step, jsoncurratedfolder)

    def errorcorrmanual(jsonfile, errordict, hdp):
        """
        Tries to correct errors manually

        Args:
            jsonfile (string): path to jsonfile
            errordict (string): path to errordictionary

        Returns:
            correct (boolean): True if the json file has been corrected
            error (string): error description
        """
        with open(jsonfile, 'r') as jsonfile:
            jsondata=json.load(jsonfile)
            for mutation in jsondata['genomicData']:
                if 'HGVS-code' in mutation['Mutations']:
                    print(mutation['Mutations']['HGVS-code'])
                    right_hgvs=input('right HGVS (u for unknown):')
                    if right_hgvs == 'u':
                        continue
                    else:
                        mutation['Mutations']['HGVS-code']=right_hgvs
                        transcript_errors[hgvs]=right_hgvs
                        with open(errordict, 'w') as dicttojson:
                            json.dump(transcript_errors, dicttojson)
                elif 'Mutation 1' in mutation['Mutations']:
                    for multimut, description in mutation['Mutations'].items():
                        print(mutation['Mutations'][multimut]['HGVS-code'])
                        right_hgvs=input('right HGVS (u for unknown):')
                        if right_hgvs == 'u':
                            continue
                        else:
                            mutation['Mutations'][multimut]['HGVS-code']=right_hgvs
                            transcript_errors[hgvs]=right_hgvs
                        with open(errordict, 'w') as dicttojson:
                            json.dump(transcript_errors, dicttojson)
                with open(jsonfile, 'w') as dicttojson:
                    json.dump(jsondata, dicttojson)
                correct, error=checkjson(jsonfile, hdp)
    return(correct, error)

class JsonFile:
    '''Compare a json object againt the specified schema.
    '''
    old_schema = {
            'submitter': {
                'team':''
                ,'name':''
            }
            ,'vcf':''
            ,'geneList':[{'gestalt_score':''}]
            ,'features':''
            ,'ranks':''
            ,'genomicData':[
                {   'Mutations':{'HGVS-code':'','Mutation 1':''}
                    ,'Test Information':{'Mutation Type':''}
                }
            ]
    }

    schema = {'algo_deploy_version': '',
            'case_id': '',
            'detected_syndromes': [{'combined_score': '',
                'feature_score': '',
                'gestalt_score': '',
                'has_mask': '',
                'omim_id': '',
                'syndrome_name': ''}],
            'documents': [],
            'features': [],
            'genomic_entries': [],
            'selected_syndromes': [{'has_mask': '', 'omim_id': '', 'syndrome_name': ''}],
            'submitter': {'user_email': '', 'user_name': '', 'user_team': ''}}


    def __init__(self, path=None, inputdict=None):
        self.process = {
            'genomic_entries': [self.load_entries]
            }
        self.raw = None
        self.errors = None
        if path is not None:
            self.read_file(path)
        elif inputdict is not None:
            self.read_dict(inputdict)

    def load_entries(self, filename, path):
        return json.load(open(os.path.join(path,filename),'r'))

    def process(self):
        self.check()
        self.raw

    def read_file(self, inputpath):
        read_json = json.load(open(inputpath,'r'))
        self.raw = read_json

    def read_dict(self, inputdict):
        self.raw = inputdict

    def generate(self):
        return self.generate_schema(self.raw)

    def generate_schema(self, data):
        '''Recursively get schema from available raw data'''
        if isinstance(data, dict):
            return {k:self.generate_schema(v) for k,v in data.items()}
        elif isinstance(data, list):
            res = [self.generate_schema(v) for v in data]
            out = []
            for r in res:
                if isinstance(r, dict):
                    if set(r.keys()) not in [set(o) for o in out]:
                        out.append(r)
                elif isinstance(r, str):
                    if r not in out:
                        out.append(r)
                else:
                    out.append(r)
            return out
        else:
            return ''

    def filter_schema(self, schema):
        if isinstance(schema,dict):
            out = {}
            for k,v in schema.items():
                if isinstance(v, tuple):
                    if not v[0]:
                        out[k] = v
                else:
                    o = self.filter_schema(v)
                    if not(o == {} or o == []):
                        out[k] = o
            return out
        elif isinstance(schema,list):
            out = []
            for v in schema:
                if isinstance(v,tuple):
                    if not v[0]:
                        out.append(v)
                else:
                    o = self.filter_schema(v)
                    if not(o =={} or o ==[]):
                        out.append(o)
            return out
        else:
            print(schema)
            raise RuntimeError()


    def check(self, false_only=True):
        schema =  self.check_schema(self.schema, self.raw)
        self.error = schema
        if false_only:
            return self.filter_schema(schema)
        else:
            return schema

    def check_schema(self, schema, rawjs):
        '''Recursively check the specified schema against the provided schema
        '''
        if isinstance(schema, dict):
            if not isinstance(rawjs,dict):
                return (False, "Not a dict")
            else:
                out = {}
                for k,v in schema.items():
                    if k not in rawjs:
                        out[k] = (False, "No key")
                    else:
                        out[k] = self.check_schema(v,rawjs[k])
                return out
        elif isinstance(schema, list):
            if not isinstance(rawjs,list):
                return (False, "Not a list")
            else:
                if len(schema) == 0:
                    return (True, "")
                out = []
                for v in rawjs:
                    res = []
                    for s in schema:
                        res.append(self.check_schema(s,v))
                    if len(schema) == 1:
                        res = res[0]
                    out.append(res)
                return out
        elif hasattr(schema, '__call__'):
            return schema(rawjs)
        else:
            if rawjs is None:
                return (False, "No value")
            elif schema == '':
                return (True, "")
            elif schema == rawjs:
                return (True, "matches expected string")
            else:
                return (False,"No value")


# ===============================
# ===== main script =============
# ===============================

def main():
    """
    main function with the following steps:
    1. Read all inputs
    2. Copy the input json file to the currated folder if not yet present
    3. If json file is incorrect try error correction (by errordict,mutalyzer and manually)
    4. Create vcf for json file
    """

    inputdir='aws_dir/cases'
    outputdir='out'
    json_files = os.listdir(inputdir)
    testcase = json_files[0]
    key_sets = []
    for f in json_files:
        fp = os.path.join(inputdir,f)
        qc = JsonFile(path=fp)
        print(qc.check())

    argv=sys.argv[1:]

    try:
        opts, args=getopt.getopt(
            argv, "h::", ["help", "jsoncurratedfolder=", "pathtojsonfile=", "errordict="])
    except getopt.GetoptError as e:
        print(e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('jsonToTable.py -i <input-folder> -o <output-file.tsv>')
            sys.exit(1)
        elif opt in ("--pathtojsonfile"):
            jsonfile=arg
        elif opt in ("--errordict"):
            errordict=arg
        elif opt in ("--jsoncurratedfolder"):
            jsoncurratedfolder=arg

    print("path to json-file",jsonfile)
    print('Errordict:', errordict)
    print("jsoncurratedfolder",jsoncurratedfolder)

    # copy json to currated; if present the return is Truej
    present=copy(jsonfile, jsoncurratedfolder)
    if present:
        print("File already currated")
        return(0)


    # Initialize server for HGVS-code validation
    hdp=hgvs.dataproviders.uta.connect()

    # check if json file is correct
    correct, error=checkjson(jsonfile, hdp)
    print(correct,error)
    # error correction
    if(error == "HGVS Syntax error"):
        correct,error = errorcorrbymut(jsonfile,errordict)
    if error == "wrong HGVS-code":
        correct, error=errorcorrbydict(jsonfile, errordict, hdp)
    if(error == "wrong HGVS-code"):
        correct, error=errorcorrmanual(jsonfile, errordict, hdp)

    # save correct json to vcf
    if correct:
        savetomultivcf(jsonfile)
    else:
        print("error in json file: " + error)


if __name__ == "__main__":
    main()
