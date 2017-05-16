#!/usr/bin/env python3

# Written by Michael Olvera 2017
# parse_genbank.py takes in a SNP file from NCBI and a genbank file from NCBI and
# joins the two, allowing the new ganbank file to be loaded into SnapGene with SNP
# locations annotated. 
#
# ./parse_genbank.py GENE_SEQUENCE.gp SNP_LIST.xml -o output_file
#
# XML snp files can be aquired from https://www.ncbi.nlm.nih.gov/variation/view/.
# Genbank files can be downloaded from https://www.ncbi.nlm.nih.gov/gene/, after 
# selecting the 'GenBank' sequence option.

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature, ExactPosition
import xml.etree.ElementTree as ET
import re, os, argparse, sys


def getChrStartStop(genbank_object):
	""" 
	Gets the start and stop on the chromosome from a genbank file, assuming its in the
	"accessions" node.

	Most genbank files have a header that starts as below:

	LOCUS       NC_000010              31464 bp    DNA     linear   CON 06-JUN-2016
	DEFINITION  Homo sapiens chromosome 10, GRCh38.p7 Primary Assembly.
	ACCESSION   NC_000010 REGION: 119648701..119680164 GPC_000001302

	Where '119648701..119680164' is the chromosomal position of the genbank. This function should
	reliably extract that infomation.

	TODO: Add an option if the information is not there. 
	"""
	accession = genbank_object.annotations["accessions"]
	non_decimal = re.compile(r'[^\d.]+')
	# Within the Genbank File, grab the chromosomal locations from  
	g_range=[non_decimal.sub('',i) for i in accession if re.search('[0-9]+\.\.[0-9]+', i)] # genomic location
	start,end = g_range[0].split('..')
	return int(start), int(end)

def parseXMLforSNPs(xml_file):
	"""
	Parses the xml file containing SNPs and outputs a list containing tuples of the name of the SNP,
	the genomic location, and the MAF. 

	TODO Include a cutoff option for minimum MAF value. 
	"""
	tree = ET.parse(xml_file).getroot()
	snp_list = []
	for snp in tree.findall('./Variants/Variant'):
		snp_list.append( (snp.find('VariantID').text, \
			int(snp.find('Location').attrib['start']), float(snp.find('MinorAlleleFrequencies/MAF').text)))
	return snp_list

def addSNPFeatures(genbank, snp_list):
	"""
	Function to add SNP lists to the genbank file. The only qualifier sofar is the 'name', which
	merges the SNP name with the MAF for visibility on SnapGene.

	There is an unfortunate bug in SnapGene where '1bp long' features are automatically converted to
	2 base pairs. I will contact Snapgene to try to get the issue resolved, although know that
	the frount of the feature is its location.
	"""
	for snp in snp_list:
		l = ExactPosition(snp[1])
		location = FeatureLocation(l,l)
		name = snp[0] + '_' + '%.3f' %snp[2]
		snp_feature = SeqFeature(location,type='snp', id=name, qualifiers={'label':name})
		genbank.features.append(snp_feature)
	return genbank

def trimSNPs(snp_list, start, end):
	"""
	Trims the length of the SNP locations to reflect the location in the abridged genbank file.
	By default, SMP locations are based on chromosomal coordinates, while genbank is relative. 
	This will convert all SNPs to their proper locations while also removing SNPs outside the 
	bounds of the genbank file. 
	"""
	length = end - start
	name, loc, freq = [sep for sep in zip(*snp_list)]
	loc = tuple(snp - start for snp in loc)
	oor_indexes = oorSNPs(loc,length)

	filtered_snps = [(name[i],loc[i],freq[i]) for i, _ in enumerate(name) if i not in oor_indexes]
	return filtered_snps


def oorSNPs(snp_len_list, genbank_len):
	"""
	Works with trimSNPs() to remove SNP out of bounds from the genbank file. Outputs the number
	of SNPs excluded, and information about how to extend the genbank file if you wish to include
	the SNPs in the future. 
	"""
	oor_snps_minus = [(i, snp) for i, snp in enumerate(snp_len_list) if snp < 0]
	oor_spns_plus = [(i, snp) for i, snp in enumerate(snp_len_list) if snp > genbank_len]
	oor_snps = oor_spns_plus + oor_snps_minus
	if not oor_snps:
		return []

	index, loc = [i for i in zip(*oor_snps)]
	errorM = "\nWarning: %i SNPs out of range " %len(index) 
	if not oor_snps_minus:
		errorM += "(0, "
	else: errorM += "(%i, " %min(loc)
	if not oor_spns_plus:
		errorM += "+0)."
	else: errorM += "+%i).\n" %(max(loc)-genbank_len)

	print(errorM)
	return index
		
def checkAssembly(genbank, xml):
	"""
	Make sure the genome assembly files are the same. For the genbnak files, its
	usually burried in the comment section:

	COMMENT     REFSEQ INFORMATION: The reference sequence is identical to
            CM000672.2.
            On Feb 3, 2014 this sequence version replaced gi:224589801.
            Assembly Name: GRCh38.p7 Primary Assembly
            The DNA sequence is composed of genomic sequence, primarily
            finished clones that were sequenced as part of the Human Genome
            Project. PCR products and WGS shotgun sequence have been added
            where necessary to fill gaps or correct errors. All such additions
            are manually curated by GRC staff. For more information see:
            http://genomereference.org.

    Wherease for the XML file is in the VariationData/MetaData/AssemblyName node:

  <VariationData>
  	<Metadata>
	    <TaxID>9606</TaxID>
	    <OrganismScientificName>Homo sapiens</OrganismScientificName>
	    <AssemblyName>GRCh38.p7</AssemblyName>
	"""
	tree = ET.parse(xml).getroot()
	version = tree.find('./Metadata/AssemblyName').text

	return version in genbank.annotations['comment'] # True if assemblies are the same

def handleArgs():
	"""
	Handles taking in command line arguments. 
	"""
	parser = argparse.ArgumentParser(description='Take in a raw genbank file and adds SNP data. Output can be used in Snapgene.')
	parser.add_argument('genbank', type=str,
                    help='A genbank file (.gb).')
	parser.add_argument('snps_input',  type=str,
                    help='A list of SNPs in XML format.')
	parser.add_argument('-o','--output', type=str,
                    help='File name to output, without the .gb.')
	return parser.parse_args()

def main():

	args = handleArgs()
	genome = SeqIO.read(args.genbank,'genbank') # You MUST tell SeqIO what format is being read
	# Run a check to make sure assemblies are matching.
	if not checkAssembly(genome, args.snps_input):
		raise StandardError('Genome assembly files are not the same between SNP list and genbank file.')
	# Get start and end points. 
	start, end = getChrStartStop(genome)
	snps = parseXMLforSNPs(args.snps_input)
	# Check for SNP outside the range of the genbank file.
	snps = trimSNPs(snps, start, end)
	# Add new SNP features 
	genome = addSNPFeatures(genome, snps)
	with open(args.output + '.gb', 'w') as out:
		SeqIO.write(genome, out,'genbank')

if __name__ == '__main__':
	main()




