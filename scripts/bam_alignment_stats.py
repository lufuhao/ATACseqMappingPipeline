#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################# Requirements #########################
import getopt
import sys
import os
import re
import numpy as np
from pandas import Series, DataFrame
import pysam
debug_mode = 0
verbose_mode = 0
opt_bam_list=[]
opt_output='MyOut.stats'
opt_threads=1
opt_nuchr=0;

#@click.command()
#debug_mode = 0
#@click.version_option(version='v20201208')
#@click.option('--bam', '-b',  multiple=True, type=string, hide_input=False, required=True, help="BAM input")
#@click.option('--threads', '-t', default=1, nargs=1, type=int, hide_input=False, required=False, help="Number of threads")
#@click.option('--out', '-o', default="MyOut.stats", nargs=1, type=string, hide_input=False, required=False, help="Statistics Output")

def ArgsParser(argv):
	functionname = 'ArgsParse'
	partdate = '202001208'
	global opt_bam_list
	global opt_output
	global opt_threads
	
	try:
		opts, args = getopt.getopt(argv, "hi:o:t:", ["help", "bam=", "output=", "threads="]) 
	except getopt.GetoptError:
		print('Error: ***.py -i bam,bam2,bam3,... -o <output>')
		print('  or: ***.py --bam=<bam,bam2,bam3,...> --output=<output>')
		sys.exit(2)

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print('***.py -i bam,bam2,bam3,... -o <output>')
			print('or: ***.py --bam=<bam,bam2,bam3,...> --output=<output>')
			sys.exit(1)
		elif opt in ("-i", "--bam"):
			opt_bams = arg
			bamarray=opt_bams.split(',')
			for indbam in bamarray:
				opt_bam_list.append(indbam)
		elif opt in ("-o", "--output"):
			opt_output = arg
		elif opt in ("-t", "--threads"):
			opt_threads = arg
	print('-----------------------------------------------------------------------')
	print(opts)
	print(args)
	print('Input BAMs: ', opt_bams)
	print (opt_bam_list)
	print('Ouput: ', opt_output)
	print('threads: ', opt_threads)
	print('-----------------------------------------------------------------------')


def getBamRef (bamfiles):
	global opt_threads
	ChrsDict={}
	if debug_mode==1 or verbose_mode==1:
		print (bamfiles)
	for bamfile in bamfiles:
		if debug_mode==1 or verbose_mode==1:
			print ('-@', opt_threads, '-H', bamfile)
		bamline = pysam.view('-@', str(opt_threads), '-H', bamfile).split("\n")
		for headline in bamline:
			if debug_mode==1 or verbose_mode==1:
				print("Test: " + headline) 
			if re.match(r'^@SQ', headline):
				if debug_mode==1 or verbose_mode==1:
					print ("Test: " + headline)
				pattern=re.compile(r'SN:(\S+)\s+LN:(\d+)')
				result=pattern.findall(headline)
				if debug_mode==1 or verbose_mode==1:
					print (result)
					print (result[0])
					print (result[0][0] + "\t" + str(result[0][1]))
				var_chr=result[0][0]
				var_len=result[0][1]
				if var_chr in ChrsDict:
					if var_len!=ChrsDict[var_chr]:
						sys.stderr.write ("Error: same sequences but different length: "+ var_chr)
				else:
					ChrsDict[var_chr]=var_len
	return ChrsDict





if __name__=='__main__':
	ArgsParser(sys.argv[1:])
	
	Chrs=list(getBamRef(opt_bam_list).keys())
	
	if debug_mode==1 or verbose_mode==1:
		print (Chrs)
		sys.exit(0)
	reads_count = DataFrame(np.ones((len(Chrs),len(opt_bam_list))),index=Chrs,columns=opt_bam_list)

	for indbam in opt_bam_list:
		if indbam[-4:].lower() == '.bam':
			bamf = pysam.Samfile(indbam, "rb" )
			## Is there an indexed .bai for the BAM? Check.
			try:
				for entry in bamf.fetch():
					codes = map(lambda x: x[0],entry.cigar)
					break
			except Exception:
				### Make BAM Indexv lciv9df8scivx 
				print ('.'),
				bam_dir = str(indbam)
				#On Windows, this indexing step will fail if the __init__ pysam file line 51 is not set to - catch_stdout = False
				pysam.index(indbam)
				bamf = pysam.Samfile(indbam, "rb" ) 
		reads_count[indbam] = [int(pysam.view("-@", str(opt_threads), "-c", indbam, Chr).strip()) for Chr in Chrs]
	if debug_mode==1 or verbose_mode==1:
		print ("reads_count: ")
		print (reads_count)
	reads_count = reads_count.T
	if debug_mode==1 or verbose_mode==1:
		print ("reads_count: ")
		print (reads_count)
		print (reads_count.iloc[:,:len(opt_bam_list)])
		print (reads_count.iloc[:,len(opt_bam_list):])
	
	reads_count['nuc_ratio'] = np.sum(reads_count.iloc[:,:len(Chrs)],axis=1) / np.sum(reads_count.iloc[:,:],axis=1)
	reads_count['org_ratio'] = np.sum(reads_count.iloc[:,len(Chrs):],axis=1) / np.sum(reads_count.iloc[:,:],axis=1)
	print(reads_count)
	reads_count.to_csv(opt_output, sep='\t', encoding='utf-8')
	sys.exit(0)
