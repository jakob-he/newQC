from JsonFile import JsonFile
import pandas as pd
import hgvs
import argparse
import hgvs
import hgvs.parser  # hgvs parsing and validation library
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.validator
import hgvs.exceptions
import pdb
import csv

def parsehgvs(hgvscode, parser,hdp):
    hgvscode="NM_033419.3:c.861G>T"
    print(hgvscode)
    variant = parser.parse_hgvs_variant(hgvscode)
    vm = hgvs.assemblymapper.AssemblyMapper(
        hdp, assembly_name='GRCh37', alt_aln_method='splign')
    var_g = vm.c_to_g(variant)
    chrom = str(int(var_g.ac.split(".")[0][-2:]))
    offset = variant.posedit.pos.start.offset
    ref = variant.posedit.edit.ref
    alt = variant.posedit.edit.alt
    print(variant,var_g)
    print(chrom, offset, ref, alt)
    return(chrom, offset, ref, alt)


hdp = hgvs.dataproviders.uta.connect()
jsondata = JsonFile(
    '/home/jakobh/Dokumente/PEDIA/PEDIA-workflow-new_download/1_qualityCheck/json/mapped/16057.json').raw

vcffolder="./"
# create multivcf

multivcf = pd.DataFrame(columns=[
                        '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NM'])

caseID = jsondata['case_id']
hgvslist = []
rowcounter = 0
hp = hgvs.parser.Parser()
for mutation in jsondata['genomicData']:
    if 'HGVS-code' in mutation['Mutations']:
        hgvscode = mutation['Mutations']['HGVS-code']
        hgvslist.append(hgvscode)
    elif 'Mutation 1' in mutation['Mutations']:
        for mutationnr, mutationdict in mutation['Mutations'].items():
            if 'HGVS-code' in mutationdict:
                hgvscode = mutationdict['HGVS-code']
                hgvslist.append(hgvscode)
        genotype = mutation['Test Information']['Genotype']
        if genotype == 'Hemizygous':
            genotype = '1'
        elif genotype == 'Homozygous':
            genotype = '1/1'
        elif genotype == 'Heterozygous' or genotype == 'Compound Heterozygous':
            genotype = '0/1'
        else:
            genotype = './1'
    for hgvscode in set(hgvslist):
        try:
            chrom, offset, ref, alt = parsehgvs(hgvscode, hp,hdp)
            multivcf.set_value(rowcounter, '#CHROM', chrom)
            multivcf.set_value(rowcounter, 'NM', hgvscode)
            multivcf.set_value(rowcounter, 'POS', offset)
            multivcf.set_value(rowcounter, 'ID', '.')
            multivcf.set_value(rowcounter, 'REF', ref)
            multivcf.set_value(rowcounter, 'ALT', alt)
            multivcf.set_value(rowcounter, 'QUAL', '.')
            multivcf.set_value(rowcounter, 'FILTER', '.')
            multivcf.set_value(
                rowcounter, 'INFO', 'HGVS="' + hgvscode + '"')
            multivcf.set_value(rowcounter, 'FORMAT', 'GT')
            multivcf.set_value(rowcounter, caseID, genotype)
            rowcounter = rowcounter + 1
        except ValueError:
            print ('HGVS-Code not parsed:', caseID, hgvscode)
            continue

# data_vcf sortieren
multivcf = multivcf.sort_values(by=['#CHROM', "POS"])
multivcf = multivcf.reset_index(drop=True)
multivcf = multivcf.fillna(value='0/0')
print(multivcf)

with open(vcffolder + "/" + caseID + ".vcf", 'w+') as outfile:
    outfile.write(
        '##fileformat=VCFv4.1\n##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

multivcf.to_csv(vcffolder + "/" + caseID + ".vcf", sep='\t', index=False,
                header=True, quoting=csv.QUOTE_NONE)
