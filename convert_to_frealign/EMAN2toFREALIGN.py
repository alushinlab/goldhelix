#!/usr/bin/env python

import os,sys
from optparse import OptionParser
import glob
import subprocess
import linecache
import struct
import shutil
from pyami import mrc

def setupParserOptions():
	parser = OptionParser()
	parser.set_usage("%prog -f <stack> -p <parameter> -c <ctf>")
	parser.add_option("-f",dest="stack",type="string",metavar="FILE",
		help="raw particle stack (black particles) - if not specified, only parameter files will be created, no new stack")
	parser.add_option("-p",dest="param",type="string",metavar="FILE",
		help="EMAN2 output parameter file")
	parser.add_option("-c",dest="ctf",type="string",metavar="FILE",
		help="per-particle CTF information file from APPION (optional)")
	parser.add_option("--mag",dest="mag",type="float", metavar="FLOAT", default=10000,
		help="actual magnification of images (default=10000)")
	parser.add_option("--norm", action="store_true",dest="norm",default=False,
		help="Normalize particles")
	parser.add_option("-m",dest="onlymodel",type="int",metavar="#",
		help="only convert this model (optional, starts with 0)")
	parser.add_option("--shrink",type="int",metavar="#",
		help="meanshrink images")
	parser.add_option("--scale",type="float",metavar="FLOAT",
		help="multiply shifts by scaling factor")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
		help="debug")
	parser.add_option("--phaseflip", action="store_true", dest="phaseflip", default=False,
		help="Only phaseflip particles, creates spider file")

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
def checkConflicts(params):
	if not params['stack']:
		print "\nWarning: no stack specified\n"
	elif not os.path.exists(params['stack']):
		print "\nError: stack file '%s' does not exist\n" % params['stack']
		sys.exit()
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
	## find number of models included in reconstruction
	f=open(params['param'])
	mods = []
	for line in f:
		l = line.split()
		model=float(l[-1])
		if model == 888 :
			continue
		if model not in mods:
			mods.append(model)
	f.close()
	return len(mods)

#=========================
def Eman2Freali(az,alt,phi):
    t1 = Transform({"type":"eman","az":az,"alt":alt,"phi":phi,"mirror":False})
    #t_conv = Transform({"type":"eman","alt":31.717474411458415,"az":90,"phi":-90,"mirror":False})
    #t2 = t1*t_conv.inverse()
    d = t1.get_params("eman")
    psi = d["phi"]+90
    if psi >360:
	psi = psi-360
    theta= d["alt"]
    phi = d["az"]-90
    return psi,theta,phi

#=========================
def createFiles(params):
	parm=params['param']
	numMods = params['num']
	mag = params['mag']
	stack = params['stack']
	debug = params['debug']

	# open EMAN2 param file
	f=open(parm,'r')

	# for each model, create an output file
	mout=[]
	mtxt=[]
	count=[]
	for m in range(numMods):
		mout.append(open("%s_%02i_frealign"%(parm,m),'w'))
		mtxt.append(open("%s_%02i.txt"%(parm,m),'w'))
		count.append(1)

	print "Calculating euler angle conversion..."

	pcount=1
	for line in f:
		l = line.split()
	 	
		parmPSI = float(l[0])
		parmTHETA = float(l[1])
		parmPHI = float(l[2])
		sx =(float(l[3]))
		sy =(float(l[4]))
		if params['shrink'] > 1:
			sx/=params['shrink']
			sy/=params['shrink']
		if params['scale']:
			sx*=params['scale']
			sy*=params['scale']
		model = int(float(l[5]))

		psi,theta,phi = Eman2Freali(parmPSI,parmTHETA,parmPHI)	

		if model != 888:
			if model == 999: model = 0
			if debug is True:
				print 'Particle %s is included' %(pcount-1)

			mtxt[model].write("%s\n" %(pcount-1))
			ctf = linecache.getline(params['ctf'],pcount)
			if debug is True:
				print 'Reading line %s in ctf file' %(pcount)
				print ctf
			c = ctf.split()

			micro = float(c[7]) 
			df1 = float(c[8])
			df2 = float(c[9])
			astig = float(c[10])

			mout[model].write("%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.1f%6d%9.1f%9.1f%8.2f%7.2f%8.2f\n" %(count[model],psi,theta,phi,sx,sy,mag,micro,df1,df2,astig,0,0))
			count[model] += 1

		pcount+=1

	# close files
	f.close()
	for m in range(numMods):
		mout[m].close()
		mtxt[m].close()
		if params['onlymodel'] is not None:
			if m!=params['onlymodel']:
				os.remove("%s_%02i_frealign"%(parm,m))
				os.remove("%s_%02i.txt"%(parm,m))


	# exit if not converting stack
	if stack is None:
		return

	# get box size
	im=EMData.read_images(stack,[0])
	nx = im[0].get_xsize()
	if params['shrink'] > 1:
		nx/=params['shrink']
	del im

	for m in range(numMods):
		if params['onlymodel'] is not None:
			if m!=params['onlymodel']: continue

		text='%s_%02i.txt' %(parm,m)
		parts = open(text).readlines()
		nimg = len(parts)

		imstack = "%s_model%02i"%(os.path.splitext(stack)[0],m)

		if params['phaseflip'] is True:
			for i in xrange(nimg):
				p = int(float(parts[i]))
				d = EMData()
				d.read_image(stack,p)
				ctf_params = data[im].get_attr("ctf")
				d = filt_ctf(data[im], ctf_params, sign=-1, binary=1)
				d.write_image(imstack+".spi",i,EMUtil.ImageType.IMAGE_SPIDER)
				progress = int(float(i)/nimg*100)
				if progress%2==0:
					print "%3i%% complete\t\r"%progress,
		else:
			print "Generating %i particle stack for Model %i..."%(nimg,m)
			for i in xrange(nimg):
				p = int(float(parts[i]))
				d = EMData()
				d.read_image(stack, p)
				if params['shrink'] > 1:
					d.process_inplace("math.meanshrink",{"n":params['shrink']})
				if params['norm'] is True:
					d.process_inplace("normalize")
				# mirror image
				d.process_inplace("xform.mirror",{"axis":"y"})
				# write image
				d.write_image(imstack+".img",i,EMUtil.ImageType.IMAGE_IMAGIC)
				progress = int(float(i)/nimg*100)
				if progress%2==0:
					print "%3i%% complete\t\r"%progress,
		print "100% complete\t"
		os.remove(text)

		# convert imagic to 3D MRC stack
		print "Converting img file to mrc stack..."

		header = mrcHeader(nx,nimg)
		tmphf = "tmpmrcheaderfile.bin"
		h = open(tmphf,'wb')
		h.write(header)
		h.close()

		cmd = "cat %s %s.img > %s.mrc"%(tmphf,imstack,imstack)
		subprocess.Popen(cmd,shell=True).wait()

#		with open(imstack+".img",'rb') as old:
#			with open(imstack+".mrc",'wb') as new:
#				new.write(header)
#				for chunk in iter(lambda: old.read(1024),b""):
#					new.write(chunk)
#		old.close()
#		new.close()

		os.remove(tmphf)
		os.remove(imstack+".img")
		os.remove(imstack+".hed")

#=========================
def mrcHeader(box,nump):
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
#=========================
if __name__ == "__main__":
	params=setupParserOptions()

	getEMANPath()
	from EMAN2 import *
	from sparx  import *

	checkConflicts(params)
	params['num']=getNumModels(params)
	print "EMAN2 parameter file contains %s models"%params['num']

	createFiles(params)

