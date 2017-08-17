#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will generate half stacks from the CTF parameter files output by a previous run of HalfStacksfromMulti.py
# This is useful for reproducing previous work, or promulgating half stack assignments to different versions of the data,
# i.e. data with different binning.
# Note that this will only work for one model at a time, but stacks for multiple models can be generated independently in series
# From the respective CTF parameter files

# import system / house keeping modules
import os
import sys
import optparse
import subprocess
import linecache

# import eman / sparx specific modules for dealing with images
from EMAN2 import *
from sparx import *
from global_def import *

# to write mrcs
from pyami import mrc

#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog --whitestack <phase-flipped stack> --blackstack <un-CTF corrected stack> --ctf_a <ctf from half a> --ctf_b <ctf from half b> --apix <apix>")
	parser.add_option("--whitestack",dest="whitestack",type="string",metavar="FILE",
		help="phaseflipped particle stack for EMAN2 gold-standard FSC refinement")
	parser.add_option("--blackstack",dest="blackstack",type="string",metavar="FILE",
		help="black (un-CTF corrected) particle stack for FREALIGN gold-standard FSC refinement")
	parser.add_option("--ctf_a",dest="ctf_a",type="string",metavar="FILE",
		help="CTF file for half stack a output by HalfStacksfromMulti.py")
	parser.add_option("--ctf_b",dest="ctf_b",type="string",metavar="FILE",
		help="CTF file for half stack b output by HalfStacksfromMulti.py")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",
                help="angstroms per pixel")
	

	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 5:
		parser.print_help()
		sys.exit()
	params={}
	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params
#=========================
def checkConflicts(params):
	if not params['whitestack']:
		print "\nWarning: no EMAN2 stack specified\n"
	elif not os.path.exists(params['whitestack']):
		print "\nError: EMAN2 stack file '%s' does not exist\n" % params['whitestack']
		sys.exit()
	if not params['blackstack']:
		print "\nWarning: no FREALIGN stack specified\n"
	elif not os.path.exists(params['blackstack']):
		print "\nError: FREALIGN stack file '%s' does not exist\n" % params['blackstack']
	if not params['ctf_a']:
		print "\nError: ctf parameter file a not specified"
		sys.exit()
	if not os.path.isfile(params['ctf_a']):
		print "\nError: ctf parameter file '%s' does not exist\n" % params['ctf_a']
		sys.exit()
	if not params['ctf_b']:
		print "\nError: ctf parameter file b not specified"
		sys.exit()
	elif not os.path.isfile(params['ctf_b']):
		print "\nError: ctf parameter file '%s' does not exist\n" % params['ctf_b']
		sys.exit()
#=========================
def makeStacksfromCTF(ctfname,whiteoutname,blackoutname,params):
	f=open(ctfname,'r')
	
	#initialize an EM image object for reading particles
	whiteparticle=EMData()
	blackparticle=EMData()

	#get relevant variables from params
	apix=float(params['apix'])
	whitestack_name=params['whitestack']
	blackstack_name=params['blackstack']
		
	i=0
	for line in f:
		l =  line.split()
		# CTF parameters start with 1, stack starts with zero, so there is an offset of 1
		particle=int(l[0])-1
		whiteparticle.read_image(whitestack_name,particle)
		blackparticle.read_image(blackstack_name,particle)
		
		#this code is taken from Richard Hall's uphead.py script
		#it sets EMAN2 attributes / parameters in header for a fresh alignment
		whiteparticle.set_attr_dict({'active':1})
		t2 = Transform({"type":"spider","phi":0,"theta":0,"psi":0})
		whiteparticle.set_attr("xform.projection", t2)
		whiteparticle.set_attr("apix_x",apix )

		whiteparticle.write_image(whiteoutname,i)
		blackparticle.write_image(blackoutname,i,EMUtil.ImageType.IMAGE_IMAGIC)
		i+=1
	f.close()

#=========================
def mrcHeader(box,nump):
	#copied verbatim from EMAN2toFREALIGN.py
        header = mrc.newHeader()
        mrc.updateHeaderDefaults(header)
        header['nx']=box
        header['ny']=box
        header['nz']=nump
        header['mode']=2
        header['mx']=box
        header['my']=box
        header['mz']=nump
        header['xlen']=box
        header['ylen']=box
        header['zlen']=nump
        header['amin']=0.0
        header['amax']=0.0
        header['amean']=0.0
        header['rms']=0.0
        header['xorigin']=0.0
        header['yorigin']=0.0
        header['zorigin']=0.0
        hbytes = mrc.makeHeaderData(header)
        return hbytes	
#=========================
def recombineHalvesFre(stack_a,stack_b,size):
	combined_stack = "FRE_recombined_stack"

	a_part = EMData()
	b_part = EMData()

	a_part.read_image(stack_a,0)
	nx = a_part.get_xsize()
	
	final_stack_size=2*size
	header = mrcHeader(nx,final_stack_size)
	
	position = 0

	for i in xrange(size):
		a_part.read_image(stack_a,i)
		b_part.read_image(stack_b,i)

		#mirror images
		a_part.process_inplace("xform.mirror",{"axis":"y"})
		b_part.process_inplace("xform.mirror",{"axis":"y"})

		#write in sequence
		a_part.write_image(combined_stack+".img",position,EMUtil.ImageType.IMAGE_IMAGIC)
		position+=1
		b_part.write_image(combined_stack+".img",position,EMUtil.ImageType.IMAGE_IMAGIC)
		position+=1
		
	tmphf = "tmpheaderfile.bin"
	h = open(tmphf,'wb')
	h.write(header)
	h.close()

	mrc_command = "cat %s %s.img > %s.mrc"%(tmphf,combined_stack,combined_stack)
	subprocess.Popen(mrc_command,shell=True).wait()

	os.remove(tmphf)
	os.remove(combined_stack+".img")
	os.remove(combined_stack+".hed")

#=========================
if __name__ == "__main__":
	params=setupParserOptions()
	checkConflicts(params)
	ctf_a=params['ctf_a']
	ctf_b=params['ctf_b']
	ctfs=(ctf_a,ctf_b)
	halves=('a','b')
	
	for half in xrange(0,2):
		whiteoutname="stack_half_%s.hdf"%halves[half]
		blackoutname="FRE_stack_half_%s.img"%halves[half]
		makeStacksfromCTF(ctfs[half],whiteoutname,blackoutname,params)
	
	f=open(ctf_a,'r')
	f_lines=f.readlines()
	size=len(f_lines)
	f.close()
	f_lines=None
	stack_a="FRE_stack_half_a.img"
	stack_b="FRE_stack_half_b.img"
	recombineHalvesFre(stack_a,stack_b,size)

