#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature, ExactPosition
import primer3
import argparse


class SingleNucleotidePolymorphism(object):
	"""
	A class to hold SNP information.
	"""
	def __init__(self, name, loc, maf):
		self.name = name
		self.loc = loc
		self.maf = maf
	def __str__(self):
		return "%s at %i with a MAF  of %f" %(self.name, self.loc, self.maf) 

class TLAfragment(object):
	"""
	A class to hold fragment infomration, SNP placement, primer location, and 
	primer scores for TLA.
	"""
	def __init__(self, seq, loc):
		self.seq = seq
		self.gc = self.getGC(self.seq)
		self.start = loc[0]
		self.end = loc[1]
		self.snps = []
		self.primers = {}
		self.score = []
	def __len__(self):
		return len(self.seq)
	def __str__(self):
		string = 'index: %s\t' %str((self.start,self.end)) + 'length: %i\t' %len(self.seq) + 'gc: %i\t' %self.gc + 'snps: %i\n' %len(self.snps)
		if len(self.seq) > 20:
			string += self.seq[:10] + ".." + self.seq[-10:]
		else: string += self.seq 
		return string
	def getGC(self, seq):
		return float((seq.count("G") + seq.count("C")))/len(seq) * 100.0 
	def addSNPs(self, snp):
		if not isinstance(snp, SingleNucleotidePolymorphism):
			raise('Error: Object not of class SingleNucleotidePolymorphism.')
		self.snps.append(snp)
	def addPrimers(self, primers):
		self.primers = primers


def basicDigest(seq, re_site,cut_index):
	"""
	Digests a DNA (as a string), returning TLAFragment objects that are digested
	along with coordination information relative to the initial DNA string.  
	"""
	fragments = []
	last_i = 0
	while True:
		i = seq.find(re_site)
		if i == -1: # Tail end of sequence
			fragments.append(TLAfragment(seq, (last_i, last_i+len(seq)))) 
			return fragments
		end = i+cut_index
		fragments.append(TLAfragment(seq[:end],(last_i, last_i+end)))
		seq = seq[i+cut_index:]
		last_i += end

def extractSNPfromGB(genbank_file):
	"""
	Extract SNP data from a genbank object, per the 'parse_genbank.py' specification. 
	The assumption being made is that snps are labeled as 'snp' in the genbank file, 
	and can be accessed as such. Also assumes labels ar in a "name_maf" format, per
	the 'parse_genbank.py' file.

	TODO: Add another funciton to handle SNP data straign from the XML 
	to return SNP objects.
	"""
	snp_list = [(feature.qualifiers['label'][0].split('_'),int(feature.location.start)) for feature in genbank_file.features if feature.type == 'snp']
	return [SingleNucleotidePolymorphism(snp[0][0],snp[1],float(snp[0][1])) for snp in snp_list]


def asignSNPS(list_of_fragmets, list_of_snps):
	"""
	Assigns snps from a list to their coorisponding genomic location on TLA fragments
	"""
	for snp in list_of_snps:
		for fragment in list_of_fragmets:
			if fragment.start <= snp.loc < fragment.end:
				fragment.addSNPs(snp)
				break
	return list_of_fragmets

def subsetTLA(list_of_fragments, min_length=180, max_length=1000, min_snps=0,
 min_gc=30, max_gc=70, flank_size=50):
	"""
	With SNP information added, subset a list of TLA fragments by cutoffs desirable for
	TLA. By default, a minimum length of 180, at least 1 SNP within the first or last 50
	nucleotides, and a %GC between 30 and 70. 

	This is the function that could use the most work, at it would be preferable to return 
	rankings or scores for each fragment, rather than arbitraty cutoffs. 

	TODO: Check for excessive k-mers or repetitive sequences ('TAATAATAA').
	"""
	filtered_fragments = []
	for frag in list_of_fragments:
		if min_length <= len(frag) <= max_length and len(frag.snps) >= min_snps \
		and min_gc <= frag.gc <= max_gc:
			if min_snps == 0:
				filtered_fragments.append(frag)
			else:
				for snp in frag.snps:
					if frag.start <= snp.loc <= (frag.start + flank_size) \
					or (frag.end - flank_size) <= snp.loc <= frag.end:
						filtered_fragments.append(frag)
						break
	return filtered_fragments

def addPrimers(list_of_fragments, flank_size=50):
	# (Idealy) define primers 
	# Not working right now, safe for later
	for frag in list_of_fragments:
		clipped = frag.seq[flank_size:-flank_size]
		l = len(clipped)
		frag.addPrimers(primer3.bindings.designPrimers(
	{
		'SEQUENCE_ID': 'MH1000',
		'SEQUENCE_TEMPLATE': clipped
	},
	{
		'PRIMER_OPT_SIZE': 20,
		'PRIMER_PICK_INTERNAL_OLIGO': 1,
		'PRIMER_INTERNAL_MAX_SELF_END': 8,
		'PRIMER_MIN_SIZE': 18,
		'PRIMER_MAX_SIZE': 22,
		'PRIMER_OPT_TM': 60.0,
		'PRIMER_MIN_TM': 57.0,
		'PRIMER_MAX_TM': 63.0,
		'PRIMER_MIN_GC': 20.0,
		'PRIMER_MAX_GC': 80.0,
		'PRIMER_MAX_POLY_X': 4,
		'PRIMER_PRODUCT_SIZE_RANGE': [[l*.2,l*.8],[l*.3,l*.7]],
	}))
	return list_of_fragments

def addTLAFeatures(genbank, fragment_list):
	"""
	Function to add SNP lists to the genbank file. The only qualifier sofar is the 'name', which
	merges the SNP name with the MAF for visibility on SnapGene.

	There is an unfortunate bug in SnapGene where '1bp long' features are automatically converted to
	2 base pairs. I will contact Snapgene to try to get the issue resolved, although know that
	the frount of the feature is its location.
	"""
	count = 0
	for frag in fragment_list:
		location = FeatureLocation(ExactPosition(frag.start),ExactPosition(frag.end))
		tla_feature = SeqFeature(location,type='tla', id='tla', qualifiers={'label':'TLA_Region_%i'%count})
		genbank.features.append(tla_feature)
		count += 1
	return genbank

def handleArgs():
	"""
	Handles taking in command line arguments. 
	"""
	parser = argparse.ArgumentParser(description='Identifies optimal primers for TLA.')
	parser.add_argument('genbank', type=str,
                    help='A genbank file (.gb).')
	parser.add_argument('output', type=str,
                    help='File name to output, without the .gb.')
	return parser.parse_args()

def main():
	args = handleArgs()
	genome = SeqIO.read(args.genbank,'genbank')
	dna = str(genome.seq)

	#	   0000000000111111111122222
	#	   0123456789012345678901234
	#dna = 'AAAAACATGTTTTCATGCCCATGAA'

	tla_fragments = basicDigest(dna,'CATG',4)

	snp_list = extractSNPfromGB(genome)
	tla_fragments = asignSNPS(tla_fragments, snp_list)

	tla_fragments = subsetTLA(tla_fragments)
	
	# Not working as expected. 
	#tla_fragments = addPrimers(tla_fragments)
	genome = addTLAFeatures(genome, tla_fragments)
	with open(args.output + '.gb', 'w') as out:
		SeqIO.write(genome, out,'genbank')

if __name__ == '__main__':
	main()

