#!/usr/bin/env python

"""
Markov Chains Traffic Assignment (MCTA) Command Line Interface (CLI)

Syntax:
	ipython mcta_cli.py -- [options] {BaseName}

options:
	--verbose	|	-v		display verbose messages and detailed results
	--showfigs	|	-f		generate figures
	--figs2png	|	-s		save figures into PNG files
	--memdump	|	-m		save mem dump variables into JSON file
							({BaseName}-Results.json)

Input data files to load (must be in the script folder):
	mcta_vis_cfg.json				visualization settings JSON file
	{BaseName}-Map.png				Base map PNG image, used in road map rendering
	{BaseName}-Map.json				Road network JSON file
	{BaseName}-Map.geojson			Road network GeoJSON file (for road attributes)

Optional output files:
	{BaseName}-Results-Fig#.png		Figures png files
	{BaseName}-Results.json			mem dump variables into JSON file

Complete Code run-through test:
	ipython ./mcta.py -- --verbose --showfigs --figs2png --memdump {BaseName}
	ipython ./mcta.py -- -v -f -s -m {BaseName}

Code by Sinan Salman, 2017-2019
"""

HELP_MSG = "Options:\n\
\t--verbose  | -v  display verbose messages and detailed results\n\
\t--showfigs | -f  generate figures\n\
\t--figs2png | -s  save figures into PNG files\n\
\t--memdump  | -m  save mem dump variables into JSON file\n\
\t                   ({BaseName}-Results.json)\n"

__author__ = 	"Sinan Salman (sinan.salman@zu.ac.ae)"
__version__ = 	"Revision: 0.17"
__date__ = 		"Date: 2019/07/08"
__copyright__ = "Copyright (c)2017-2019 Sinan Salman"
__license__ = 	"GPLv3"

### Initialization #######################################################################

import sys
import mcta
import mcta_vis
import mcta_rw
import scipy as sp

# options
verbose = False
showfigs = False
figs2png = False
savememdump = False

# global variables
base = None

sp.set_printoptions(suppress=True,precision=3,linewidth=140)

### Data processing ######################################################################

def ProcessCLI():
	"""Process CLI parameters"""
	global base
	global verbose
	global showfigs
	global figs2png
	global savememdump

	if len(sys.argv) == 1:
		print("Missing argument\n\nSyntax:\n\tipython mcta_cli.py -- [options] {BaseName}")
		print(HELP_MSG)
		sys.exit(0)
	if '--verbose' in sys.argv or '-v' in sys.argv:
		print ("*** option: verbose mode")
		verbose = True
	if '--showfigs' in sys.argv or '-f' in sys.argv:
		print ("*** option: generate figures")
		showfigs = True
	if '--figs2png' in sys.argv or '-s' in sys.argv:
		print ("*** option: figures will be saved to PNG files")
		figs2png = True
	if '--memdump' in sys.argv or '-m' in sys.argv:
		print ("*** option: save MemDump to JSON file")
		savememdump = True
	base = sys.argv[len(sys.argv)-1]


### Main #################################################################################

if __name__ == "__main__":
	ProcessCLI()
	(Settings, JSON_Map, GeoJSON_Map) = mcta_rw.LoadDataFiles(base)
	lengths, lanes, FFS, P = mcta.SetupMCTA(JSON_Map, GeoJSON_Map, Verbose=verbose)
	Results = mcta.SolveMCTA(lengths=lengths, lanes=lanes, P=P, 
							VehiclesCount=JSON_Map['VehiclesCountEst'], 
							Objectives=['D','K','C','PI','E'], 
							FreeFlowSpeeds=FFS, SkipChecksForSpeed=False)
	if savememdump:
		mcta_rw.SaveResults(Results)
	if showfigs or figs2png:
		mcta_vis.Initialize(JSON_Map,Settings,base,ShowFigs=showfigs, SaveFigs2PNG=figs2png)
	for x in Results.keys():
		if x not in ['linkIDs','P_org','P_updated','eigenvalues']:
			print(f'\n{x}: {Results[x]}')
			if showfigs or figs2png:
				mcta_vis.Generate_Figure(Results,x)
	if showfigs:
		import matplotlib.pyplot as plt
		plt.show(block=True)