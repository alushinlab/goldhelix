#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will split a stack output from appion into half stacks stacks
# Based on the output parameters from a multi-model EMAN2 / sparx refinement
# for gold-standard FSC refinement of each model
#
# It will also split per-particle FSC parameters into the same groups for FREALIGN
#

#
#import modules, this could probably stand to be pruned
#

# import system / house keeping modules
import os
from optparse import OptionParser
import sys
import optparse
import glob
import subprocess
import linecache
import struct
import shutil

# import eman / sparx specific modules for dealing with images
from EMAN2 import *
from sparx import *
from global_def import *

# import random for shuffling
import random
#
# to write mrcs
from pyami import mrc

#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog --whitestack <phase-flipped stack> --blackstack <un-CTF corrected stack> --param <multimodel_params> --ctf <perpart_ctf> --apix <apix>")
	parser.add_option("--whitestack",dest="whitestack",type="string",metavar="FILE",
		help="phaseflipped particle stack for EMAN2 gold-standard FSC refinement")
	parser.add_option("--blackstack",dest="blackstack",type="string",metavar="FILE",
		help="black (un-CTF corrected) particle stack for FREALIGN gold-standard FSC refinement")
	parser.add_option("--param",dest="param",type="string",metavar="FILE",
		help="EMAN2 output parameter file from sorting EMAN2 run")
	parser.add_option("--ctf",dest="ctf",type="string",metavar="FILE",
		help="per-particle CTF information file from APPION")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",
                help="angstroms per pixel")
	

	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 4:
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
	if not params['param']:
		print "\nError: no EMAN2 parameter file specified"
		sys.exit()
	if not os.path.isfile(params['param']):
		print "\nError: EMAN2 parameter file '%s' does not exist\n" % params['param']
		sys.exit()
	if not params['ctf']:
		print "\nError: no CTF parameter file specified"
		sys.exit()
	elif not os.path.isfile(params['ctf']):
		print "\nError: Appion CTF parameter file '%s' does not exist\n" % params['ctf']
		sys.exit()

#=========================
def getEMANPath():	
	### get the imagicroot directory	
	emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip() 
	
	if emanpath:
		emanpath = emanpath.replace("EMAN2DIR=","")
	if os.path.exists(emanpath):
		return emanpath
	print "EMAN2 was not found, make sure it is in your path"
	sys.exit()

#=========================
def getNumModels(params):
	# find number of models included in reconstruction
	# taken from Gabe Lander's EMAN2FREALIGN script
	f=open(params['param'])
	mods = []
	for line in f:
		l = line.split()
		model=float(l[-1])
		if model == 888 :
			continue
		if model == 999 :
			model = 0
			mods.append(model)
		if model not in mods:
			mods.append(model)
	f.close()
	return len(mods)

#==========================
def splitFromMulti(params):
	# split stacks and perpart ctf into temp stacks before generating random halfs
	
	#get relevant variables from params
	parm=params['param']
	numMods = params['num']
			
	#initialize an empty list of lists that will keep track of which model each segment belongs to
	#and a list of lists that contains particles
	modAssignments=[]
	
	for y in range (numMods):
		modAssignments.append([])
		
	#open EMAN2 parameter file
	f=open(parm,'r')

	#intialize particle counter
	part_count = 0

	#iterate over lines in the eman2 output
	for line in f:
		l = line.split()
		
		#get model number
		model=float(l[-1])
		
		#check for single model or particle exclusion
		if model == 999 :
			model = 0
		if model == 888 :
			part_count += 1
			continue
		#add particle to the proper list for its model
		modAssignments[int(model)].append(part_count)
		part_count += 1	
					
	#return the lists of particles and ctf lines split by model
	
	return modAssignments
	
#=========================
def randomizeHalves(model_Assignments,model_number):
	#This function randomly assigns the particles and parameters for a model to equally sized half stacks
	#Since we need to match the parameters with the particles this is slightly cumbersome		
	#generate a list of integers with the length of the number of particles which match the model
	picklist=list(range(len(model_Assignments[model_number])))
		
	#randomly sort this list to determine assignments
	random.shuffle(picklist)
		
	#split this shuffled picklist into alternating halves
	first_half = picklist[::2]
	second_half = picklist[1::2]
	
	if len(first_half) > len(second_half):
		popped = first_half.pop()
		print "Odd number of particles; particle %s excluded to keep half stacks equal in size."%(str(popped))
	if len(second_half) > len(first_half):
		popped = second_half.pop()
		print "Odd number of particles; particle %s excluded to keep half stacks equal in size."%(str(popped))

	#initialize empty lists to accept picks from ctf parameters and stacks
	assignment_halves=[[],[]]
			
	#split with two loops, don't know a more elegant way to do this
	for first_h in range(len(first_half)):
		assignment_halves[0].append(model_Assignments[model][first_half[first_h]])
	for second_h in range(len(second_half)):
		assignment_halves[1].append(model_Assignments[model][second_half[second_h]])
	
	return assignment_halves

#========================
def writeHalfStacksandCTF(assigned_Halves,model_number,params):
	#This function writes out EMAN2 hdf files for the half stacks 
	
	#initialize an EM image object for reading particles
	whiteparticle=EMData()
	blackparticle=EMData()

	#get relevant variables from params
	apix=float(params['apix'])
	whitestack_name=params['whitestack']
	blackstack_name=params['blackstack']
	ctf_name=params['ctf']
	
	#two halves
	halves=('a','b')
	
	#for each half
	for half in range(len(halves)):
		#intialize file names for write
		whiteout_stack="stack_model_%02i_half_%s.hdf"%(model_number,halves[half])	
		blackout_stack="FRE_stack_model_%02i_half_%s.img"%(model_number,halves[half])

		out_ctf="ctf_model_%02i_half_%s"%(model_number,halves[half])
		
		ctf_output=open(out_ctf,'w')
	
		#for each particle
		for part in range(len(assigned_Halves[half])):
			#find the right particle in the list
			particletowrite=assigned_Halves[half][part]
			
			#the ctf line to write is this plus one
			ctflinetowrite=particletowrite+1

			whiteparticle.read_image(whitestack_name,particletowrite)
			blackparticle.read_image(blackstack_name,particletowrite)

			ctf_line=linecache.getline(ctf_name,ctflinetowrite)
			
			#this code is taken from Richard Hall's uphead.py script
			#it sets EMAN2 attributes / parameters in header for a fresh alignment
			whiteparticle.set_attr_dict({'active':1})
			t2 = Transform({"type":"spider","phi":0,"theta":0,"psi":0})
		        whiteparticle.set_attr("xform.projection", t2)
			whiteparticle.set_attr("apix_x",apix )
			
			whiteparticle.write_image(whiteout_stack,part)
			
			#write black particle out as well
			blackparticle.write_image(blackout_stack,part,EMUtil.ImageType.IMAGE_IMAGIC)

			#write CTF line
			ctf_output.write(ctf_line)
		ctf_output.close()
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
def recombineHalvesFre(assigned,mod):
	stack_a_name = "FRE_stack_model_%02i_half_a.img"%mod
	stack_b_name = "FRE_stack_model_%02i_half_b.img"%mod
	combined_stack = "FRE_recombined_stack_model_%02i"%mod

	a_part = EMData()
	b_part = EMData()

	finalstacksize=2*(len(assigned[0]))

	a_part.read_image(stack_a_name,0)
	nx = a_part.get_xsize()
	
	header = mrcHeader(nx,finalstacksize)
	
	position = 0

	for i in xrange(len(assigned[0])):
		a_part.read_image(stack_a_name,i)
		b_part.read_image(stack_b_name,i)

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

	getEMANPath()

	checkConflicts(params)
	params['num']=getNumModels(params)
	numMods = params['num']
	modelAssignments=splitFromMulti(params)
	
	for model in range (0,numMods):
		halves_assigned = randomizeHalves(modelAssignments,model)
		writeHalfStacksandCTF(halves_assigned,model,params)
		recombineHalvesFre(halves_assigned,model)	

	
