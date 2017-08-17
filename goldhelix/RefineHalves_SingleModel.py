#!/usr/bin/env python

# Author: Greg Alushin, galushin@rockefeller.edu

# This script will perform gold-standard FSC (i.e. independent half dataset) refinement 
# of a single model from a multimodel refinement
# Note that it is important to return to a highly-filtered initial model!!!!
#

#
#import modules, this could probably stand to be pruned
#

from __future__ import division

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
import time

# import eman / sparx specific modules for dealing with images
from EMAN2 import *
from sparx import *
from global_def import *

#========================
def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog --half_a <phase_flipped_half_stack_1> --half_b <phase_flipped_half_stack2> --init <initial_model> --apix <apix> --ang <angular steps> --dtrans <translational steps> --trans <translational range> --hpars <twist rise> --nproc <num processors>")

	parser.add_option("--half_a",dest="stack_a",type="string",metavar="FILE",
		help="half stack 1 (output from HalfStacksfromMulti.py)")
	parser.add_option("--half_b",dest="stack_b",type="string",metavar="FILE",
		help="half stack 2 (output from HalfStacksfromMulti.py)")
	parser.add_option("--init",dest="init_mod",type="string",metavar="FILE",
		help="initial model (highly filtered)")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",default=2.18,
                help="angstroms per pixel")
	parser.add_option("--ang",dest="angular_steps",type="string",metavar="STRING",default='-1 -1 30 30 10 10',
		help="angular step series, separated by spaces and surrounded by single quotes")
	parser.add_option("--delta",dest="delta",type="string",metavar="STRING",default='4 4 3 3 2 2',
		help="angular step of reference projections, separated by spaces and surrounded by single quotes")
	parser.add_option("--dtrans",dest="trans_steps",type="string",metavar="STRING",default='4 4 2 2 1 1',
		help="translational step series, separated by spaces and surrounded by single quotes")
	parser.add_option("--trans",dest="trans_ranges",type="string",metavar="STRING",default='16 16 8 8 4 4',
		help="translational range series, separated by spaces and surrounded by single quotes")
	parser.add_option("--hpars",dest="helical_params",type="string",metavar="STRING",default='-166.5 27.7',
		help="twist and rise, separated spaces and surrounded by single quotes. Left handed twist is negative")
	parser.add_option("--nprocs",dest="n_procs",type="int",metavar="INT",default=32,
		help="number of processors.  Must be an even integer multiple of 32 to run correctly on biowulf g72 nodes.")
	parser.add_option("--crit_pass",dest="crit_pass",type="float",metavar="FLOAT",default=0.5,
		help="criterion for choosing butterworth lowpass filter pass frequency from FSC, default 0.5")
	parser.add_option("--crit_stop",dest="crit_stop",type="float",metavar="FLOAT",default=0.3,
		help="criterion for choosing butterworth lowpass filter stop frequency from FSC, default 0.3")
	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit()
	params={}
	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params
#=========================
def checkFileExistence(params):
	if not os.path.isfile(params['stack_a']):
		print "\nError: stack file '%s' does not exist\n" % params['stack_a']
		sys.exit()
	if not os.path.isfile(params['stack_b']):
		print "\nError: stack file '%s' does not exist\n" % params['stack_b']
		sys.exit()
	if not os.path.isfile(params['init_mod']):
		print "\nError: 3D model '%s' does not exist\n" % params['init_mod']
		sys.exit()
#=========================
def checkNprocSanity(params):
	num_procs=float(params['n_procs'])
	num_nodes = num_procs / 32
	integer_test = num_nodes.is_integer()

	if integer_test is False:
		print "nproc must be divisible by 32 to run correctly on biowulf as this script was intended.\n"
		sys.exit()
#=========================
def loadEMANPaths():	
	eman2load="module load EMAN2/2.1"
	subprocess.Popen(eman2load,shell=True).wait() 
	
	emanload="source /data/Alushinlab/environment/EMAN1-1.9/eman.bashrc"	
	subprocess.Popen(emanload,shell=True).wait()
#=========================
def genClusterScript(runname,stackname,refname,nproc,outdir,xrad,t_s,delta,an,twist,rise,scriptout):
	#this module generates a run.sh style script which can be submitted to the biowulf cluster via PBS
	shebang="""#!/bin/bash"""

	pbscommand1="""#SBATCH --job-name='"""
	pbscommand2="""'"""
	#pbscommand3="""#PBS -k oe"""
	pbscommand4="""cd $SLURM_SUBMIT_DIR"""
	pbscommand5="""source /data/Alushinlab/environment/.bashrc"""
	#pbscommand6="""source /usr/local/EMAN1-1.9/eman.bashrc"""
	pbscommand7="""module load EMAN2/2.1"""

	executecommand1="""mpirun -np """
	executecommand2=""" /data/Alushinlab/lab_scripts/EMAN2/EMAN2_actin2/refine.py """
	executecommand3=""" """
	executecommand4=""" """
	executecommand5=""" --ou=120 --rs=1 --xr='"""
	executecommand6="""' --ts='"""
	executecommand7="""' --delta='"""
	executecommand8="""' --an='"""
	executecommand9="""' --snr=0.08 --maxit=1 --ref_a=S --cutoff=10 --pix_cutoff=10000 --MPI --olmask=80 --lmask=280 --protos='1' --hpars='"""
	executecommand10=""" """
	executecommand11="""' --hsearch='0.0 50.0' --oplane=10 --recon_pad=2 --full_output > out.log"""

	f=open(scriptout, 'w')
	f.write("%s\n%s%s%s\n\n%s\n\n%s\n\n%s\n\n"%(shebang,pbscommand1,runname,pbscommand2,pbscommand4,pbscommand5,pbscommand7))
	f.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n"%(executecommand1,str(nproc),executecommand2,stackname,executecommand3,refname,executecommand4,outdir,executecommand5,str(xrad),executecommand6,str(t_s),executecommand7,str(delta),executecommand8,str(an),executecommand9,str(twist),executecommand10,str(rise),executecommand11))
	f.close()

#===========================
def SubmitandWait(script_tosubmit,nproc,delta_step,ref_dir):
	num_nodes = int(nproc)
	
	cmd = "sbatch --ntasks=%s --exclusive --ntasks-per-core=1 --time=10-00:00:00  %s"%(str(num_nodes),script_tosubmit)
	print(cmd)
	subprocess.Popen(cmd,shell=True).wait()

	#get path to parameters output, which is last file written by refinement
	working_dir = os.getcwd()
	param_path="%s/%s/paramout_00%s_00"%(working_dir,ref_dir,str(delta_step))	
	
	#print param_path

	while not os.path.exists (param_path):
		time.sleep(10)
#============================
def GenHParsFile(twist,rise):
	hpars_file = "helical_params.spi"
	f = open(hpars_file,'w')
	f.write("; spi/spi\n")
	f.write("    1 2 %11.5f %11.5f\n"%(twist,rise))
	f.close()
#============================
def FindFSCcut(fsc_curve,cutoff):
	for i in range(len(fsc_curve[1])):
		if fsc_curve[1][i]<=cutoff:
			band=fsc_curve[0][i]
			break
	try:
		return (band)
	except NameError:
		print "Warning: FSC curve did not fall below %s cutoff."%(str(cutoff))
		return (False)
#============================
def PostRecon(vol_a,vol_b,mask,hpars,apix,round,pass_crit,stop_crit):
	#sum vols for calculating helical parameters
	cmd_sum = "proc3d %s %s sum=tempsum.mrc"%(vol_a,vol_b)
	subprocess.Popen(cmd_sum,shell=True).wait()
	
	#
	# calculate new helical parameters
	#

	#rotate temp volume and convert to spider
	cmd_rot = "proc3d tempsum.mrc tempsumrot.spi rot=90,0,0 spidersingle"
	subprocess.Popen(cmd_rot,shell=True).wait()

	#run hsearch_lorentz
	cmd_hsearch = "hsearch_lorentz tempsumrot.spi %s %.4f 0.0 50.0 0.03 0.03"%(hpars,apix)
	subprocess.Popen(cmd_hsearch,shell=True).wait()
	
	#
	# impose helical parameters on each half volume
	#
	
	#get most recent rise and twist: code taken from Gabe Lander's helical module
	f = open(hpars)
	lines = f.readlines()
	f.close()
	last_pars=lines[-1].strip().split()
	
	#get rise and twist
	rise=float(last_pars[-1])
	twist=float(last_pars[-2])

	#stupid sparx needs negative of twist, 2/3 of vol for helicizing
	helical_region = 0.67
	helical_region = float(helical_region)
	sparx_twist=-1*twist

	#open and helicize volumes
	
	a=EMData()
	a.read_image(vol_a)
	a=a.helicise(apix,rise,sparx_twist,helical_region)

	b=EMData()
	b.read_image(vol_b)
	b=b.helicise(apix,rise,sparx_twist,helical_region)

	a.write_image("%s_h.hdf"%vol_a[0:-4])
	b.write_image("%s_h.hdf"%vol_b[0:-4])	
	
	#calculate FSC
	fsc_filename=("gold_fsc_round_%i"%(round))
	
	#filter half volumes to use as reference for next round
	maskfile=EMData()
	maskfile.read_image(mask)
	w=1
	fsc_gold=fsc_mask(a,b,maskfile,w,fsc_filename)
	
	filt_pass=FindFSCcut(fsc_gold,pass_crit)
	filt_stop=FindFSCcut(fsc_gold,stop_crit)

	#also process sum volume for looking at

        fsc_point143=FindFSCcut(fsc_gold,0.143)
        filt143=float(apix/fsc_point143)

        print "FSC 0.143 between half volumes is %4.2f"%filt143
        mvcommand="mv tempsum.mrc sumvol_unfilt_round_%i.mrc"%round
        subprocess.Popen(mvcommand,shell=True).wait()

	#go with butterworth filter
	try:	
		a_filtered=filt_btwl(a,filt_pass,filt_stop,True)
		b_filtered=filt_btwl(b,filt_pass,filt_stop,True)
		a_filtered.write_image("%s_h_filt.hdf"%vol_a[0:-4])
		b_filtered.write_image("%s_h_filt.hdf"%vol_b[0:-4])
		print "Reference half volumes filtered at pass band: %4.2f stop band: %4.2f"%((apix/filt_pass),(apix/filt_stop))
		return 
	except:
		#still write unfiltered references if finding cutoffs fails
		print "Warning!!!! no filter applied, continuing with misnamed unfiltered references"
		a_filtered=a
		b_filtered=b
		a_filtered.write_image("%s_h_filt.hdf"%vol_a[0:-4])
       		b_filtered.write_image("%s_h_filt.hdf"%vol_b[0:-4])
		return	
#
#=========================
#
if __name__ == "__main__":

	#initial setup and sanity check
	params=setupParserOptions()
	#loadEMANPaths()
	checkFileExistence(params)
	checkNprocSanity(params)

	#unpack some parameters to make things simpler
	init_mod=str(params['init_mod'])
	apix=float(params['apix'])
	nprocs=int(params['n_procs'])
	init_model=str(params['init_mod'])
	crit_pass=float(params['crit_pass'])
	crit_stop=float(params['crit_stop'])

	angular_steps=str(params['angular_steps']).strip().split()
	trans_steps=str(params['trans_steps']).strip().split()
	trans_ranges=str(params['trans_ranges']).strip().split()
	deltas=str(params['delta']).strip().split()
	helical_params=str(params['helical_params']).strip().split()
	
	twist=float(helical_params[0])
	rise=float(helical_params[1])

	#generate helical parameters file
	GenHParsFile(twist,rise)
        hpars="helical_params.spi"
	
	#main loop
	
	for round in range(len(angular_steps)):
		print "Beginning round %i"%round
		#for simplicity's sake will currently run jobs in sequence rather than parallel

		halves=('a','b')
		
		#for the hell of it get current hpars to see if half vols give back the same helical parameters
		#once again use Gabe Lander's code
		f=open(hpars)
		lines=f.readlines()
		f.close()
		current_hpars=lines[-1].strip().split()
		
		for half in range(len(halves)):
			if round == 0:
				current_ref = init_mod
			else:
				current_ref = "model_%s_%02d_h_filt.hdf"%(halves[half],(round-1))	
			
			current_runname='run_%s_%02d'%(halves[half],round)
			current_scriptname='run_%s_%02d.sh'%(halves[half],round)
			current_outdir='out_%s_%02d'%(halves[half],round)
			current_stack=str(params['stack_%s'%(halves[half])])
                        
			#TESTING123TESTING123
			#for now make fake directory and files
			#fake_dir_cmd = "mkdir %s"%(current_outdir)
			#print fake_dir_cmd
			#subprocess.Popen(fake_dir_cmd,shell=True).wait()
			#relevant_output = "volNoSym_%03d_00.hdf"%(int(deltas[round]))
			#fake_vol_cmd = "touch %s/%s"%(current_outdir,relevant_output)
			#subprocess.Popen(fake_vol_cmd,shell=True).wait()
			###############################################

			#hpars 2 should be twist, hpars 3 should be rise
			
			genClusterScript(current_runname,current_stack,current_ref,nprocs,current_outdir,trans_ranges[round],trans_steps[round],deltas[round],angular_steps[round],float(current_hpars[2]),float(current_hpars[3]),current_scriptname)
			print "Submitting refinement for half %s"%halves[half]

			SubmitandWait(current_scriptname,nprocs,deltas[round],current_outdir)
			
			#link unsymmetrized output volume back to main directory for post-processing
			relevant_output = "volNoSym_%03d_00.hdf"%int(deltas[round])
			link_cmd = "ln -s %s/%s model_%s_%02d.hdf"%(current_outdir,relevant_output,halves[half],round)
			subprocess.Popen(link_cmd,shell=True).wait()
			
			#link mask back to main directory, need not do every time but won't hurt
			link_mask_cmd = "ln -s %s/mask3D_cyl.mrc ."%(current_outdir)
			subprocess.Popen(link_mask_cmd,shell=True).wait()
	
		#post processing
		volume_a="model_a_%02d.hdf"%round
		volume_b="model_b_%02d.hdf"%round
		mask="mask3D_cyl.mrc"

		PostRecon(volume_a,volume_b,mask,hpars,apix,round,crit_pass,crit_stop)
