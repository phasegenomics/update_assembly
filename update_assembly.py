# -*- coding: utf-8 -*-
'''
Created on June 25, 2018

@author: Aakash Sur

This script ingests one juicebox assembly file, and outputs two juicebox assembly 
files, one for each haplotype. At this point, the Falcon-Phase step to determine
phasing of contigs has been run, the final polishing step is finished, and all
that remains to be completed is generating the fasta files for each haplotype.
In order to do so, we must generate two assembly files so that the
juicebox_assembly_to_scaffolds.py program can use them to generate fasta files.

Copyright 2018, Phase Genomics Inc. All rights reserved.

The contents of this file are proprietary and private and are not intended for
distribution or use by any person or entity except Phase Genomics. You may not
use, modify, or distribute it in any fashion. You may not copy this file. You
may not describe the contents of this file to any other party.
'''

import argparse
import logging
import os

def setup():
	'''Return command line arguments

	    Returns:
	        parser.parse_args(): object with command line parsed arguments
	'''
	parser = argparse.ArgumentParser(
		description = 'Update Juicebox assembly file to use one haplotype'
	)
	parser.add_argument(
		'-a', 
		'--assembly',
		metavar = 'File', 
		required = True, 
		help = 'Assembly file name'
	)
	parser.add_argument(
		'-p', 
		'--phase',
		metavar = 'File', 
		required = True, 
		help = 'FALCON-Phase output file name'
	)

	parser.add_argument(
		'--sequence-length',
		type = int,  
		default = 9,
		help = 'Length of sequence names. Default is 9')

	parser.add_argument(
		'-l',
		'--logging',
		default = 'DEBUG',
		choices = ['WARNING', 'INFO', 'DEBUG'], 
		help = 'Logging level')
	return parser.parse_args()


def parse_phase(filename):
	'''Parse the phase information and produce two haplotypes.

    Args:
        filename(str): filename of phase information file

    Returns:
        list: two lists, one for each haplotype, each containing the set of
        contigs that belong in that phase
	'''
	logging.debug('Parsing phase information.')
	maternal = []
	paternal = []

	with open(filename) as infile:
		for line in infile:
			split = line.strip().split()
			maternal.append(split[1])
			paternal.append(split[2])

	return maternal, paternal

def write_assembly(filename, output, length, haplotype):
	'''Produce a new Juicebox assembly file with a single haplotype.

		The contig name always appears in the beginning, and with Falcon these 
		are of predictable length: e.g.
		">000021F_1 1 6822245"
		">000005F_1:::fragment_2:::debris 43 50000"
		This means we can extract the the contig name from the 1:10 slice, 
		i.e. "000021F_1" for the first example. We then check if the name 
		appears in our haplotype list, if not, then it must be the other haplotig, 
		which according to falcon-phase naming convention means swapping the 
		ending number, i.e. "000021F_1" -> "000021F_0". After determining the 
		correct name, we recreate essentially the original assembly file by just 
		changing the names and preserving the remaining text. 

    Args:
        filename(str): name of input Juicebox assembly file to be read
        output(str): name of output Juicebox assembly file to be written
        haplotype(list): list of contigs in that particular phase

    Returns:
        None

    Effect:
    	Writes a file based on the input Juicebox assembly file and make the 
    	relevant swaps so that all the contigs in the assembly file belong to 
    	one phase. 
	'''
	logging.debug('Writing new Juicebox assembly file.')
	infile = open(filename)
	outfile = open(output, 'w')

	for line in infile:
		# See description for text parsing logic. 
		if '>' in line:
			name = line[1:1+length]
			suffix = line[1+length:]

			if name in haplotype:
				outfile.write(line)
			else:
				if name[-1] == '0':
					name = '>{0}1'.format(name[:-1])
				else:
					name = '>{0}0'.format(name[:-1])
				
				outline = name + suffix
				outfile.write(outline)
		else:
			outfile.write(line)

	infile.close()
	outfile.close()	

def main():
	arguments = setup()
	logging.basicConfig(level = arguments.logging)
	maternal, paternal = parse_phase(arguments.phase)
	split = os.path.split(os.path.abspath(arguments.assembly))

	# Write twice, once for each haplotype.
	write_assembly(
		arguments.assembly, 
		os.path.join(split[0],'0-{0}'.format(split[1])), 
		arguments.sequence_length, 
		maternal
	)
	write_assembly(
		arguments.assembly, 
		os.path.join(split[0],'1-{0}'.format(split[1])),
		arguments.sequence_length, 
		paternal
	)

if __name__ == '__main__':
	main()