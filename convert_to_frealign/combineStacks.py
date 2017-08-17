#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will combine all the little stacks from lmbfgs into a single imagic stack
#
#import modules, this could probably stand to be pruned
#

# import system / house keeping modules
import os
import subprocess
import sys
import optparse
from global_def import *
from sparx import *
#
#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -i <stack_list> ")
	parser.add_option("-i",dest="input",type="string",metavar="FILE",
		help="list of particle stacks, one per line (easy to generate with ls *.mrcs > stacklist.txt)")	

	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit()
	params={}
	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params
#=========================
def checkExistence(f_name):
	if not os.path.exists(f_name):
		print "Error: file %s does not exist.  Exiting..."%f_name
		sys.exit()
#=========================
if __name__ == "__main__":
	params=setupParserOptions()

	infile=params['input']
	checkExistence(infile)

	f = open(infile,'r')
	stacks = f.readlines()

	a = EMData()
	combined_stack_index = 0

	for stack in stacks:
		stack=stack.strip()
		checkExistence(stack)
		imn = EMUtil.get_image_count(stack)
		for i in xrange(imn):
			a.read_image(stack,i)
			a.write_image('combined_stack.hed',combined_stack_index)
			combined_stack_index += 1
		

