# goldhelix
Software for cryo-EM structural analysis of helical filaments

Dependencies: EMAN2, hsearch_lorentz of Egelman, Appion libraries (for conversion to mrc stack, this could be modified as current versions of EMAN2 have support for MRC stacks), EMAN1, Frealign v9.11


Contents:


EMAN2:

up_head.py: python script to convert image stacks into EMAN2 format


EMAN2/EMAN2_actin:

refine.py: main procedure script for EMAN2/SPARX projection matching and IHRSR: launched by run_single.sh, run_multi.sh, and RefineHalves_SingleModel.py

functions.py, hfunctions.py, reconstruction_rj.py: python modules containing functions used by refine.py


EMAN2/submission_scripts:

run_single.sh: SLURM submission script for single-model IHRSR using EMAN2/SPARX and hsearch_lorentz

run_multi.sh: SLURM submission script for multi-model IHRSR using EMAN2/SPARX and hsearch_lorentz


goldhelix:

HalfStacksFromMulti.py: python script to create two random half-data sets from one full data set including particle stacks, parameters, and ctf values

RefineHalves_SingleModel.py: python script for gold-standard FSC refinement of independent half datasets using EMAN2/SPARX and hsearch_lorentz


convert_to_frealign:

EMAN2toFREALIGN.py: python script to convert parameters and image stacks from EMAN2 format into FREALIGN format

combineStacks.py: python script for combining multiple stacks into a single: useful for combining output from, e.g. John Rubinstein's alignparts

CombineHalfParamFilesFre.py: python script for combining half-dataset parameter files

RegenHalfStacksfromCTF.py: python script for regenerating half stacks from a previous round of HalfStacksFromMulti.  Useful for maintaining half-dataset assignments between different versions of the data.


example_files:

example_ctf: Example of CTF parameters output by appion in FREALIGN format: input for HalfStacksFromMulti.py

example_paramout: Example of alignment parameters output by multi-model EMAN2/SPARX refinement: input for HalfStacksFromMulti.py

example_mparameters: Example parameter file for half-dataset refinement in FREALIGN.

example_frealign.par: Example frealign parameter file for first round of refinement


Basic protocol:

-Convert phase-flipped stack of particles to .hdf format using up_head.py.  A reference in HDF format is also required
-Perform initial alignment and reconstruction using the contents of EMAN2 directory.  Single and multi-model refinements are supported.  This stage will sort out multiple states and/or bad segments by cross-correlation cutoff.
-Using parameter file from round of refinement chosen based on reasonable FSC curve (falls to zero), particle stack, and per-particle ctf file, run HalfStacksFromMulti.py to generate half-stacks and corresponding CTF parameters.
-Refine independent half reconstructions using RefineHalves_SingleModel.py.  A strongly low-pass filtered version (>=35 angstroms) of the output from EMAN2 refinement can be used as the initial model at this stage.
-Use EMAN2toFREALIGN.py to generate frealign parameter files for each half-dataset utilizing the output alignment parameters from goldhelix, as well as outputs from HalfStacksFromMulti.py.  
-Combine half params using CombineParamFilesFre.py.  This parameter file, as well as the combined stack, can be used for FREALIGN v9.11 refinement, following instructions on the Grigorieff lab website.

Note that scripts will need to be modified to point to the correct paths for your workstation / cluster environment.

Please contact Greg Alushin (galushin [at] rockefeller [dot] edu) with questions, comments, or concerns.
