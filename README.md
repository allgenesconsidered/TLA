# TLA
Finding TLA-optimized primers

Written by Michael Olvera 2017
parse_genbank.py takes in a SNP file from NCBI and a genbank file from NCBI and
joins the two, allowing the new ganbank file to be laded into SnapGene with SNP
locations annotated. 

./parse_genbank.py GENE_SEQUENCE.gp SNP_LIST.xml -o output_file

XML snp files can be aquired from https://www.ncbi.nlm.nih.gov/variation/view/.
Genbank files can be downloaded from https://www.ncbi.nlm.nih.gov/gene/, after 
selecting the 'GenBank' sequence option.

find_tla_primers.py will (eventually) identify optimum TLA primers from an input 
genomic locus containing SNP info. 
