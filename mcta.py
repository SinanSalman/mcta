#!/usr/bin/env python

"""
Markov Chains Traffic Assignment (MCTA)

Analyze road network traffic using Descrete Time Markov Chain (DTMC)
based on mdoel by Crisostomi, et. al (2011). "A Google-like model of road network dynamics
and its application to regulation and control". International Journal of Control, 84(3),
633–651. http://doi.org/10.1080/00207179.2011.568005

Syntax:
	ipython mcta.py -- [options] {BaseName}

options:
	--verbose	|	-v		display verbose messages and detailed results
	--showfigs	|	-f		generate figures
	--figs2png	|	-s		save figures into PNG files
	--geojson	|	-g		save mcta results as GeoJSON file for analysis in qGIS
							({BaseName}-Results.geojson)
	--memdump	|	-m		save mem dump variables into JSON file
							({BaseName}-Results-MemDump.json)

Input data files to load (must be in the script folder):
	mcta.json						Rules and Assumptions JSON file
	{BaseName}-Map.png				Base map PNG image, used in road map rendering
	{BaseName}-Map.json				Road network JSON file
	{BaseName}-Map.geojson			Road network GeoJSON file (for road attributes)

Optional output files:
	{BaseName}-Results-Fig#.png		Figures png files
	{BaseName}-Results-MemDump.json	save mem dump variables into JSON file
	{BaseName}-Results.Geojson		Save results back into geoSJSON file for analysis in
									qGIS (same format as "{BaseName}-Map.geojson" but
									with added link attributes for results analysis)

Complete Code run-through test:
	ipython ./mcta.py --  --verbose --showfigs --figs2png --geojson --memdump {BaseName}
	ipython ./mcta.py --  -v -f -s -g -m {BaseName}

Code by Sinan Salman, 2017

ToDo:
# adjust speed acording to congestion = loop runs?
# add extra state (for optimization purposes)
Later:
# need more accurate distance estimation function
"""

HELP_MSG = "Options:\n\
\t--verbose  | -v  display verbose messages and detailed results\n\
\t--showfigs | -f  generate figures\n\
\t--figs2png | -s  save figures into PNG files\n\
\t--geojson  | -g  save mcta results as GeoJSON file for analysis in qGIS\n\
\t                   ({BaseName}-Results.geojson)\n\
\t--memdump  | -m  save mem dump variables into JSON file\n\
\t                   ({BaseName}-Results-MemDump.json)\n"

__author__ = 	"Sinan Salman (sinan.salman@zu.ac.ae)"
__version__ = 	"Revision: 0.16"
__date__ = 		"Date: 2017/02/12"
__copyright__ = "Copyright (c)2017 Sinan Salman"
__license__ = 	"GPLv3"

### Initialization #######################################################################

import os
import sys
import json
import math
import scipy as sp
import scipy.sparse as sprs
import matplotlib.pyplot as plt

# options
verbose = False
showfigs = False
figs2png = False
savegeojson = False
savememdump = False

# global variables
path = None
base = None
BASE = None
CMAP = None
FIGURE_SIZE = None
MemDump = {}
GeoJSON = {}
Settings = {}
VehiclesCount = 0

nodes = {}
links = {}
bounds = {}
mapbounds = {}
revlinks = []

P = None
results = {}

sp.set_printoptions(suppress=True,precision=3,linewidth=140)

### Data processing ######################################################################

def ProcessCLI():
	"""Process CLI parameters"""
	global base
	global verbose
	global showfigs
	global figs2png
	global savegeojson
	global savememdump

	if len(sys.argv) == 1:
		print("Missing argument\n\nSyntax:\n\tipython mcta.py -- [options] {BaseName}")
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
		if not showfigs:
			print ("      note: --showfigs automatically set to creat figures before saving")
			showfigs = True
	if '--geojson' in sys.argv or '-g' in sys.argv:
		print ("*** option: save roadnetwork and mcta results as GeoJSON file")
		savegeojson = True
	if '--memdump' in sys.argv or '-m' in sys.argv:
		print ("*** option: save MemDump to JSON file")
		savememdump = True
	base = sys.argv[len(sys.argv)-1]

def LoadDataFiles():
	"""load data and configuration files"""
	global path
	global BASE
	global GeoJSON
	global Settings
	global CMAP
	global FIGURE_SIZE
	global nodes
	global links
	global revlinks
	global bounds
	global mapbounds
	global VehiclesCount

	path = os.path.split(os.path.realpath(__file__))[0] + os.path.sep
	BASE = path + base
	if verbose:
		print ('data loading and processing:')
		print ("\tlooking in: " + path + "")
		print ('\treading configurations and rules from mcta.json file...')

	if os.path.isfile('mcta.json'):
		with open('mcta.json') as mctajsonfile:
			Settings = 		json.load(mctajsonfile)
			CMAP = 			plt.get_cmap(Settings['CMAP'])
			FIGURE_SIZE = 	Settings['FIGURE_SIZE']
	else:
		print ("\nerror - can't find 'mcta.json'.")
		sys.exit(0)

	if os.path.isfile(BASE + '-Map.json'):
		if verbose: print ('\treading ' + base + '-Map.json')
		with open(BASE + '-Map.json') as Mapjsonfile:
			mapjson = json.load(Mapjsonfile)
			VehiclesCount = mapjson['VehiclesCountEst']
			bounds = 		mapjson['bbox_roads']
			mapbounds = 	mapjson['bbox_map']
			nodes = 		mapjson['nodes']
			links = 		mapjson['links']
			revlinks =		mapjson['LinkReverseDirection']
			nodes = {int(k):v for k,v in nodes.items()} # convert dictionary keys to int since JSON require them to be string
			links = {int(k):v for k,v in links.items()} # convert dictionary keys to int since JSON require them to be string
	else:
		print ("\nerror - can't find JSON file for the specified BaseName: " + base)
		sys.exit(0)

	if os.path.isfile(BASE + '-Map.geojson'):
		if verbose: print ('\treading ' + base + '-Map.geojson')
		with open(BASE + '-Map.geojson') as Mapgeojsonfile:
			GeoJSON = json.load(Mapgeojsonfile)
	else:
		print ("\nerror - can't find JSON file for the specified BaseName: " + base)
		sys.exit(0)

### Module functions: for using MCTA to evaluate networks ################################

def test_AD():

	# original network
	# mg=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

	# solution of single object GA optimization run
	mg=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 0, 2, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

	# solution of single object GA optimization run - Crisis example
	# mg=[0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 2, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 2, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

	SetEnviroment(env_base='AD',
				  env_verbose=False,
				  env_showfigs=True,
				  env_figs2png=False,
				  env_CrisisRoads=[62,63])
	D, K, ResultMsg = ModifyNetworkAndSolveMC(mg)
	ExplainModification(mg)
	print ("D=",D,"\tK=",K,"\n",ResultMsg,'\n','{:} changes'.format(len([x for x in mg if x>0])))
	ShowModificationOnMap(mg)
	plt.show()

def test():
	SetEnviroment(env_base='Ex',
				  env_verbose=True,
				  env_showfigs=True,
				  env_figs2png=False,
				  env_CrisisRoads=[8,9])
	if False:
		D, K = ModifyNetworkAndSolveMC([0,0,0,0,0,0,0,0]) # reverse none
		print ("-\t",D,"\t",K)
		D, K = ModifyNetworkAndSolveMC([0,0,0,1,0,0,0,0]) # reverse 6
		print ("6\t",D,"\t",K)
		D, K = ModifyNetworkAndSolveMC([0,0,0,1,2,0,0,0]) # reverse 6 & 9
		print ("6,9\t",D,"\t",K)
	else:
		#     0  2  4  6  8 10 12 14
		# mg = [2, 0, 2, 1, 1, 2, 2, 1]
		mg = [0, 0, 0, 0, 0, 0, 0, 0]
		D, K, ResultMsg = ModifyNetworkAndSolveMC(mg)
		ExplainModification(mg)
		print ("D=",D,"\tK=",K,"\n",ResultMsg,'\n','{:} changes'.format(len([x for x in mg if x>0])))
		plt.show()

def ExplainModification(modification_gene):
	"""explain network modification using plain english"""
	msg="the following links direction are reversed: "
	for i in range(len(modification_gene)):
		if modification_gene[i] == 1:
			msg += str(revlinks[i][0]) + ", "
		elif modification_gene[i] == 2:
			msg += str(revlinks[i][1]) + ", "
	return msg

def ShowModificationOnMap(modification_gene):
	sol=sp.matrix([14]*len(modification_gene)*2).T
	for i in range(len(modification_gene)):
		if modification_gene[i] == 1:
			sol[revlinks[i][0]]=27
		elif modification_gene[i] == 2:
			sol[revlinks[i][1]]=27
	RenderMap(sol, 'Road Network Changes',LinkEndMark = '', LimitDensity_28=True)

def SetEnviroment(env_base='AD',
				  env_verbose=False,
				  env_showfigs=False,
				  env_figs2png=False,
				  env_savememdump=False,
				  env_savegeojson=False,
				  env_CrisisRoads=[]):
	"""set enviroment for mcta; used when imported as a module"""
	global base
	global verbose
	global showfigs
	global figs2png
	global savememdump
	global savegeojson

	base = env_base
	verbose = env_verbose
	showfigs = env_showfigs
	figs2png = env_figs2png
	savegeojson = env_savegeojson
	savememdump = env_savememdump

	LoadDataFiles()
	SetupTransitionMatrix()
	ModifyForLinkTravelTime()

	# print ('before crisis closure')
	# for j in feeding_CR:
	# 	for i in range(len(links)):
	# 		if P[j,i] > 0: print ('P[{:},{:}]={:}'.format(j,i,P[j,i]))
	if len(env_CrisisRoads):
		CloseCrisisRoads(env_CrisisRoads)
	# print ('after crisis closure')
	# for j in feeding_CR:
	# 	for i in range(len(links)):
	# 		if P[j,i] > 0: print ('P[{:},{:}]={:}'.format(j,i,P[j,i]))

def CloseCrisisRoads(CR):
	"""modify network by closing RD roads"""
	global P
	nl = len(links)
	ONES = sp.matrix(sp.ones((nl,1)))

	X = sp.matrix(sp.eye(nl))
	for i in CR: X[i,i] = 0	 # 0 = no flow AND 1 = full flow

	# close crisis Roads
	if verbose: print('closing crisis roads:' + str(CR) + '\n'+'*'*30)
	if verbose: print('P\n',P.todense())
	if verbose: print('X\n',X)
	P_d = sp.diag(sp.diag(P.todense()))
	if verbose: print('P_d\n',P_d)
	C1 = X*(P-P_d)*X
	if verbose: print('X*(P-P_d)*X\n',C1)
	P_d = X*P_d*X
	if verbose: print('new P_D = X*P_d*X\n',P_d)
	C2 = C1*ONES
	if verbose: print('Row sum of X*(P-P_d)*X\n',C2)
	for i in range(len(C2)):
		if C2[i] == 0:
			C2[i] = 1
	if verbose: print('Row sum of X*(P-P_d)*X, w/ zeros replaced with ones\n',C2)
	C3 = sp.matrix(sp.diag(C2.A1)).I
	if verbose: print('Factor\n',C3)
	C3 = C3 * (sp.eye(nl)-P_d)
	if verbose: print('Factor, w/ diagonal values taken into consideration\n',C3)
	C4 = C3*C1
	if verbose: print('multiply by factor matrix\n',C4)
	P = C4 + P_d
	if verbose: print('New P\n',P)
	if verbose: print('check P row sum\n',P*ONES)

def ModifyNetworkAndSolveMC (modification_gene):
	"""modify network using provided modification genes and solve the resulting MC"""

	nl = len(links)
	ONES = sp.matrix(sp.ones((nl,1)))

	closed = []
	r = [-1] * nl # link reverse list
	for i in range(len(modification_gene)):
		l1 = revlinks[i][0]
		l2 = revlinks[i][1]
		r[l1] = l2
		r[l2] = l1
		assert modification_gene[i] in [0,1,2] ,"modification gene is out of bound (0,1,2)"
		if   modification_gene[i] == 1:
			closed.append(l1)
		elif modification_gene[i] == 2:
			closed.append(l2)

	X = sp.matrix(sp.eye(nl))
	for i in closed: X[i,i] = 0	 # 0 = no flow AND 1 = full flow

	# Non-Ergodic Network Handling; recursively close all orphan links (no inbound or outbound probabilities) (after closure (setting a row and column i to zeros, more links may become orphan)
	C1 = X*P*X
	REMOVE = closed.copy()
	KEEP   = [i for i in range(nl) if i not in REMOVE]
	Zr = False
	Zc = False
	while True:
		# remove all zero columns
		tmp = ONES.T*C1
		rem = [i for i in KEEP if round(tmp[0,i],6) == 0]
		REMOVE += rem
		KEEP = [i for i in range(nl) if i not in REMOVE]
		if rem != []: Zc = True
		else: Zc = False
		for i in REMOVE: X[i,i] = 0
		C1 = X*P*X

		# remove all zero rows
		tmp = C1*ONES
		rem = [i for i in KEEP if round(tmp[i,0],6) == 0]
		REMOVE += rem
		KEEP = [i for i in range(nl) if i not in REMOVE]
		if rem != []: Zr = True
		else: Zr = False
		for i in REMOVE: X[i,i] = 0
		C1 = X*P*X

		if not Zc and not Zr: break

	Delete_States = True   # best results and faster computing time using 'True'

	# calculate vehicle/m density
	length = sp.array([links[i]['length'] for i in links.keys()])
	# lanes = sp.array([(2-X[r[id],r[id]])*links[id]['lanes'] for id in links.keys()])
	# replaced with new constraint to allow different number of lanes values for opposing roads
	lanes = sp.array([ X[i,i] * links[i]['lanes'] + ( 1 - X[r[i],r[i]] ) * links[r[i]]['lanes'] for i in links.keys() ])

	# Generate Q permutation matrix
	j=0
	Q=sp.matrix(sp.eye(nl))
	for i in range(nl):
		if X[i,i]==1:
			Q[[i,j]]=Q[[j,i]]
			j+=1

	if Delete_States:
		nl     = nl-len(REMOVE)
		ONES   = sp.delete(ONES,REMOVE,0)
		lanes  = sp.delete(lanes,REMOVE,0)
		length = sp.delete(length,REMOVE,0)
		C1     = sp.delete(C1,REMOVE,1)
		C1     = sp.delete(C1,REMOVE,0)
		C2     = (C1*ONES).A1
		C      = sp.matrix(sp.diag(C2)).I * C1
	else:
		C1 = X*P*X
		B  = Q*C1*Q.T
		C2 = [x if x>0 else 1 for x in (B*ONES).A1]
		C  = sp.matrix(sp.diag(C2)).I * B
		assert sp.allclose(Q*ONES,ONES), "not all rows in Q include a value of 1, not a proper permutation (orthogonal matrix), please check!\n{:}".format(Q)
		assert sp.allclose(Q.T*ONES,ONES), "not all columns in Q include a value of 1; not a proper permutation (orthogonal matrix), please check!\n{:}".format(Q)
		assert sp.allclose(Q*Q.T,sp.eye(nl)), "Q*Q' is not an identity matrix; not a proper permutation (orthogonal matrix), please check!\n{:}".format(Q)

		length = (sp.matrix(length)*Q.T).A1
		lanes  = (sp.matrix(lanes)*Q.T).A1

	#find the left-hand Perron eigenvector
	from scipy.linalg import eig
	eigenvalues, eigenvectors = eig(C.T)

	# Get Perron vectors and root
	idx = eigenvalues.argmax()
	assert round(eigenvalues[idx].real,6) == 1, 'Perron root Error\nNetwork has the following links\n\tGene: \t\t{:}\n\tClosed: \t{:}\n\tRemoved: \t{:}\n\tRemaining: \t{:}\n\tEigenvalues: {:}'.format(modification_gene,closed,[x for x in REMOVE if x not in closed],[x for x in range(len(links)) if x not in REMOVE],eigenvalues)
	PI = sp.matrix(eigenvectors[:,idx].real).A1
	SUM = PI.sum()
	if SUM != 1: # scaling an eigenvector by scalar value "a" does not change the eigenvector: P (a*v) = Lamda (a*v), where a*v is the eigen vector
		PI /= PI.sum()

	# Kemeny constant; where K = Ki, i = 1..n (using eigenvalues)
	if 1 in [round(x,6) for x in sp.delete(eigenvalues,idx)]:
		K = float("inf")
		ResultMsg = 'reducible network with multiple communicating classes (K=inf)'
	else:
		idx = eigenvalues[:].argsort()[::-1] #eigenvalue and in reversed order
		eigenvalues = eigenvalues[idx]
		K = (1/(1 - eigenvalues[1:])).sum()
		if not Delete_States and K != 0: K -=  len(REMOVE)
		ResultMsg = 'success'

	D = VehiclesCount * PI/(length*lanes)

	if verbose and ResultMsg != 'success':
		print (ResultMsg)
		print ('Eigenvalues: {:}'.format(eigenvalues))
		print ('Network has the following links\n\tGene: \t\t{:}\n\tClosed: \t{:}\n\tRemoved: \t{:}\n\tRemaining: \t{:}\n'.format(modification_gene,closed,[x for x in REMOVE if x not in closed],[x for x in range(len(links)) if x not in REMOVE]))

	if showfigs:
		print ('generating figures...')
		streets = sp.array(list(links.keys()))[sp.newaxis] + 1 # to make it into ndarray (scipy) and start from 1
		if Delete_States:
			ADD_BACK = sp.array([0]*len(REMOVE))
			PI = sp.append(PI,ADD_BACK)
			PI = (sp.matrix(PI)*Q).A1
			D = sp.append(D,ADD_BACK)
			D = (sp.matrix(D)*Q).A1
		else:
			streets = (sp.matrix(streets)*Q.T).A1
		# Plot2D(streets, PI, 'Roads', 'Probability density', 'Stationary distribution of vehicles')
		Plot2D(streets, D, 'Roads', 'Density (vehicle/km/lane)', 'Road traffic density',LimitDensity_28=True, Yaxis_lim=[-3,60])
		# RenderMap(PI, 'Stationary distribution of vehicles heatmap','')
		RenderMap(D, 'Road traffic density (vehicle/km/lane)','',LimitDensity_28=True)

	return round(max(D),6), round(K,6), ResultMsg

### Generate turning probabilities #######################################################

def list2string(lst): # to strings of the format (#:1,2,3,4,...)
	"""not used yet, but good for saving lists into qGIS string format"""
	if lst == []:	return ['(0:)']
	else: lst = "({:}:{:})".format(len(lst),str(lst).trim('[]').replace(' ',''))

def string2list(str, TYPE='int'): # for strings of the format (#:1,2,3,4,...)
	lst = str.split(':',1)[1].strip(')').split(',')
	if lst == ['']:
		return []
	elif TYPE=='int':
		return list(map(int,lst))
	elif TYPE=='float':
		return list(map(float,lst))
	else:
		raise ValueError('error - unknown type passed to function string2list: ' + TYPE)

def SetupTransitionMatrix():
	"""populate transition matrix from GeoJSON link attributes."""
	global P
	global links
	if verbose: print ('\n\tConstructing Transition Matrix...')
	linkcount = len(links)
	P = sprs.lil_matrix((linkcount,linkcount))
	turns = ['g.right','g.left','g.forward','g.u-turn']
	FoundErrors = ""
	corrections1 = 0
	corrections2 = 0
	for item in GeoJSON['features']:
		id = item['properties']['l.id']
		if verbose: print ('\tlink: {:}'.format(id))
		links[id]['coordinates'] = item['geometry']['coordinates']
		totalprob = 0
		for t in turns:
			if isinstance(item['properties'][t],list):
				conn = item['properties'][t]
				prob = item['properties'][t+'P']
			else:
				conn = string2list(item['properties'][t],'int')
				prob = string2list(item['properties'][t+'P'],'float')
			num = len(conn)
			if verbose: print ('\t{:>12}: {:20} Prob: {:}'.format(str(t),str(conn),str(prob)))
			for i in range(num):
				P[id,conn[i]] = prob[i]
				totalprob += prob[i]
		if totalprob != 1:
			if round(totalprob,5) == 1:
				if verbose: print ('\t\tcorrecting ({:}) rounding error at link: {:}'.format(1-totalprob,id))
				corrections1 += 1
				for i in range(linkcount):
					P[id,i] /= totalprob
				#double check!
				totalprob = 0
				for i in range(linkcount):
					totalprob += P[id,i]
				remainder = 1-totalprob
				if remainder == 0:
					if verbose: print ('\t\tcorrected.')
				elif -0.00001 < remainder and remainder < 0.00001:
					corrections2 += 1
					P[id,id] += remainder
					if verbose: print ('\t\tcorrecting ({:}) rounding error at link: {:} AGAIN! remainder added to P[i,i]'.format(remainder,id))
				else:
					assert False,"problem while correcting rounding errors at link {:}. TotalProb = {:}".format(id,totalprob)
			else:
				print ("\t*** probabilities don't add up to 1 for link {:}, totalprob = {:}".format(id,totalprob))
				FoundErrors += str(id) + ', '
	if verbose and corrections1 + corrections2 > 0:
		print ('\tcorrected probability rounding errors in:\n\t\t{:} links from the 1st try\n\t\t{:} links from the 2nd try\n'.format(corrections1, corrections2))
	assert FoundErrors == "",'error - please correct the above issues befor continuing for links: ' + FoundErrors[:-2]
	Psum = sp.sum(P.todense(),axis=1)
	assert sp.allclose(Psum,sp.ones(Psum.shape[0])),"error - The generated turning probabilities matrix 'P' is not a stochastic matrix"
	if savememdump:
		global MemDump
		MemDump['P_original'] = [ str([", ".join('{:.3f}'.format(i) for i in P.todense().tolist()[j])]).replace("'","").replace("[","").replace("]","") for j in range(P.shape[0]) ]

def ModifyForLinkTravelTime():
	"""normalize link travel times (ltt) and reflict them into the P probability matrix
	returns P and StepTime (time/step conversion factor)"""
	global P
	global results
	P_original = P
	if verbose: print ('normalize link travel times (ltt):')
	ltt = sp.array([v['length']/v['speed'] for k,v in sorted(links.items())]) # time = len / spd
	if verbose: print ('ltt=',ltt)
	ltt_min = min(ltt)
	assert ltt_min>0,"error - minimum link travel time = zero"
	if verbose: print ('ltt_min=',ltt_min)
	ltt /= ltt_min
	if verbose: print ('ltt/ltt_min=',ltt)
	lttv = (ltt-1)/ltt 		# probability of staying; Pii
	if verbose: print ('lttv=',lttv)
	lttd = sp.diag(lttv)	# lttv as a square matrix with a diagonal of its values
	if verbose: print ('lttd=',lttd)
	f = sp.transpose((1-lttv)[sp.newaxis]) # use sp.newaxis to create a 2D array that can be transposed
	if verbose: print ('f=',f)
	P = sprs.lil_matrix((P.A*f)+lttd)	# need to use array element-wise operations
	Psum = sp.sum(P.todense(),axis=1)
	if verbose: print ('Psum=',Psum)
	assert sp.allclose(Psum,sp.ones(Psum.shape[0])),"error - the modified 'P' matrix, for link travel time, is not a stochastic matrix"
	if not sp.allclose(P_original.todense(),P.todense()):
		if verbose:
			print('\tmodified P matrix for links travel time based on legth and speed.\n')
	else:
		print('warning - link travel time did not modify P! check link lengths, speeds and results.')
	results['StepTime'] = ltt_min
	if savememdump:
		global MemDump
		MemDump['P_updated'] = [ str([", ".join('{:.3f}'.format(i) for i in P.todense().tolist()[j])]).replace("'","").replace("[","").replace("]","") for j in range(P.shape[0]) ]

### Markov chains ########################################################################

def FindEigen():
	"""find the left-hand Perron eigenvectors (top 2)"""
	global results
	from scipy.linalg import eig
	eigenvalues, eigenvectors = eig(P.T.todense())

	# sort eigenvalues/vectors
	idx = eigenvalues.argsort()[::-1]
	eigenvalues = eigenvalues[idx]
	eigenvectors = eigenvectors[:,idx]

	if verbose: print ('\t',eigenvalues.shape[0], ' eigenvalues: ', eigenvalues)
	if verbose: print ('\teigenvector: \n', eigenvectors.T)

	#check that Perron Root is real
	assert eigenvalues[0].imag == 0j, "error - Perron Root value is complex ({:}+{:}j)".format(eigenvalues[0].real,eigenvalues[0].imag)
	assert round(eigenvalues[0].real,6) == 1, "error - cannot find Perron Root (eigenvalue {:} != 1)".format(eigenvalues[0].real)
	for j in range(eigenvectors.shape[1]):
		assert eigenvectors[j,0].imag == 0j, "error - Perron Root eigen vector include complex value @ {:}: {:}+{:}j".format(j,eigenvectors[0,j].real,eigenvectors[0,j].imag)

	PI = sp.matrix(eigenvectors[:, 0].real).A1
	SUM = PI.sum()
	if SUM != 1: # scaling an eigenvector by scalar value "a" does not change the eigenvector: P (a*v) = Lamda (a*v), where a*v is the eigen vector
		if verbose: print ('\tNormalizing eigen vector for Perron Root by a scaling factor of: {:}'.format(SUM))
		PI /= PI.sum()
	C  = sp.matrix(eigenvectors[:, 1].real).A1

	results['eigenvalues'] = eigenvalues
	results['eigenvectors'] = eigenvectors
	results['StationaryProb'] = PI
	results['Clusters'] = C
	return PI, C, eigenvalues

def GenInv():
	"""calculate general inverse (Darzin)"""
	from numpy.linalg import matrix_rank
	Q = sp.matrix(sp.identity(P.shape[0])-P)
	Qk = Q
	Qk1 = sp.dot(Qk,Qk)
	# Determine index of A - when rank(A^k) = rank(A^(k+1))
	# The Drazin inverse of a matrix of index 0 or 1 is called the group inverse
	assert matrix_rank(Qk) == matrix_rank(Qk1)
	X = sp.linalg.solve(Qk1,Qk) # A^(k+1) * X = A^k
	Qi = Qk * sp.dot(X,X)        # Ad =  A^(k) * X^(k+1)
	assert sp.allclose(sp.dot(Q,Qi),sp.dot(Qi,Q)) ,'general inverse 1st property invalid (Q)(Q#) = (Q#)(Q)'
	assert sp.allclose(sp.dot(sp.dot(Q,Qi),Q),Q)  ,'general inverse 2nd property invalid (Q)(Q#)(Q) = (Q)'
	assert sp.allclose(sp.dot(sp.dot(Qi,Q),Qi),Qi),'general inverse 3rd property invalid (Q#)(Q)(Q#) = (Q#)'
	return Qi

def MFPT(Qi,Pi):
	"""find Mean first passage times"""
	global results
	m = sp.matrix(sp.zeros(P.shape))
	for i in range(P.shape[0]):
	    for j in range(P.shape[1]):
	        if i==j:
	            m[i,j] = 0
	        else:
	            m[i,j] = (Qi[j,j] - Qi[i,j]) / Pi[j]
	results['MFPT'] = m
	return m

def SolveMC():
	"""solve markov chain traffic assignment problem"""
	global results
	import time

	print ('\nMCTA results:')
	streets = sp.array(list(links.keys()))[sp.newaxis] + 1 # to make it into ndarray (scipy) and start from 1
	tmr = time.time()

	# solve MC
	pi, c, eigenvalues = FindEigen()
	# calculate mean first passage time matrix "m"
	Qi = GenInv()
	mfpt = MFPT(Qi,pi)
	# Kemeny constant; where K = Ki, i = 1..n (using mean first passage times)
	K1=0
	for j in range(P.shape[0]):
		K1 += mfpt[0,j] * pi[j]

	# Kemeny constant; where K = Ki, i = 1..n (using eigenvalues)
	K2 = (1/(1 - eigenvalues[1:])).sum()
	assert round(K1-K2,6)==0,"Kemeny constant calculation methods do not agree! ( {:} != {:} )".format(K1, K2)

	# calculate vehicle/m density
	length = sp.array([links[id]['length'] for id in links.keys()])
	lanes = sp.array([links[id]['lanes'] for id in links.keys()])
	D = VehiclesCount * pi/(length*lanes)

	print ('\n\tKemeny constant (MFPT method): {:.4f} '.format(K1) + ' steps ({:.4f} '.format(K1*results['StepTime']*60) + ' min)')
	print ('\tKemeny constant (eigV method): {:.4f} '.format(K2) + ' steps ({:.4f} '.format(K2*results['StepTime']*60) + ' min)')
	print ('\texpecting time to mixing: {:.4f} '.format(K1+1) + ' steps ({:.4f} '.format((K1+1)*results['StepTime']*60) + ' min)')
	print ('\tcompleted MCTA in {:.3f} seconds.'.format(time.time()-tmr))
	if verbose: print ('\n')

	results['VehicleDensity'] = D
	results['KemenyConst'] = K1

	if showfigs:
		print ('generating figures...')
		Plot2D(streets, pi, 'Roads', 'Probability density', 'Stationary distribution of vehicles')
		Plot2D(streets, D, 'Roads', 'Density (vehicle/km/lane)', 'Road traffic density',LimitDensity_28=True)
		Plot2D(streets, c, 'Roads', 'Entries of the second eigenvector', 'Road clusters')
		RenderMap(pi, 'Stationary distribution of vehicles heatmap','')
		RenderMap(c, 'Clusters in road network','')
		RenderMap(D, 'Road traffic density (vehicle/km/lane)','',LimitDensity_28=True)
		Plot3D(mfpt, 'Mean first passage times (steps)')

### Plotting #############################################################################

def range_f(start, stop, step):
	"""like the internal range() but for float type"""
	i = start
	while i < stop:
		yield i
		i += step

def Plot2D(X, Y, Xlabel, Ylabel, Title, LimitDensity_28=False, Yaxis_lim=[]):
	"""plot 2D points"""
	import matplotlib.colors as colors
	if LimitDensity_28:
		DotColor = [CMAP(int(Y[i]/28*254)) if Y[i]<28 else CMAP(254) for i in range(Y.shape[0])]
	else:
		DotColor = [CMAP(int((Y[i]-Y.min())/(Y.max()-Y.min())*254)) for i in range(Y.shape[0])]
	fg, ax = plt.subplots(figsize=FIGURE_SIZE)
	ax.grid(linestyle='dotted')
	ax.set_xlabel(Xlabel)
	ax.set_ylabel(Ylabel)
	ax.set_title(Title)
	# ax.set_xlim([X.min(),X.max()])
	if Yaxis_lim != []:
		ax.set_ylim(Yaxis_lim)
	ax.scatter(X,Y,color=DotColor);
	SaveFigurePNG()
	return fg, ax

def Plot3D(Z, Title):
	"""plot 3D surface and countour plots"""
	from mpl_toolkits.mplot3d.axes3d import Axes3D
	import matplotlib.colors as colors
	fg1 = plt.figure(figsize=FIGURE_SIZE)
	ax1 = Axes3D(fg1)
	ax1.set_ylabel('from') # rows
	ax1.set_xlabel('to')   # columns
	ax1.set_zlabel('time(steps)')
	ax1.set_title(Title)
	X, Y = sp.meshgrid(range(Z.shape[0]),range(Z.shape[0]))
	ax1.view_init(30, 240)
	ax1.plot_surface(X,Y,Z,cmap=CMAP);
	SaveFigurePNG()

	fg2, ax2 = plt.subplots(figsize=FIGURE_SIZE)
	ax2.grid(True)
	ax2.set_title(Title)
	ax2.set_ylabel('from') # rows
	ax2.set_xlabel('to')   # columns
	CF = ax2.contourf(Z.A,cmap=CMAP)#,norm=colors.LogNorm())
	plt.colorbar(CF)
	SaveFigurePNG()
	return fg1, ax1, fg2, ax2

### Render road network ##################################################################

def onpickMap(event):
	"""roadmap onpick event handler"""
	thisline = event.artist
	thisline.set_visible(not thisline.get_visible())
	thisline.figure.canvas.draw()
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	points = tuple(zip(xdata[ind], ydata[ind]))
	print('picked object:', thisline.name, ' onpick points:', points)

def RenderMap(X, Title='RoadNetwork', LinkEndMark = 's', LimitDensity_28=False):
	"""render road map"""
	if LimitDensity_28:
		LinksColor = [CMAP(int(X[i]/28*254)) if X[i]<28 else CMAP(254) for i in range(X.shape[0])]
	else:
		LinksColor = [CMAP(int((X[i]-X.min())/(X.max()-X.min())*254)) for i in range(X.shape[0])]

	# add a 0.05 margin to the plot
	bbounds={}
	lonRange = abs(bounds['maxlon'] - bounds['minlon'])
	latRange = abs(bounds['maxlat'] - bounds['minlat'])
	bbounds['minlon'] = bounds['minlon'] - lonRange * 0.05
	bbounds['maxlon'] = bounds['maxlon'] + lonRange * 0.05
	bbounds['minlat'] = bounds['minlat'] - latRange * 0.05
	bbounds['maxlat'] = bounds['maxlat'] + latRange * 0.05

	fig = plt.figure(figsize=FIGURE_SIZE)

	# Using contourf to provide colorbar info, then clear the figure
	if LimitDensity_28:
		Z = [[0,0],[0,0]]
		levels = [x for x in range_f(0,28*1.1,28/10)]
	else:
		step = (X.max()-X.min())/10
		Z = [[0,0],[0,0]]
		levels = [x for x in range_f(X.min(),X.max()+step,step)]
	cb = plt.contourf(Z, levels, cmap=CMAP)
	plt.clf()

	# by setting limits before hand, plotting is about 3 times faster
	ax = fig.add_subplot(111,autoscale_on=False,
						 xlim=(bbounds['minlon'],bbounds['maxlon']),
						 ylim=(bbounds['minlat'],bbounds['maxlat']))
	ax.set_title(Title)

	if os.path.isfile(BASE + '-Map.png'):
		import PIL
		img = PIL.Image.open(BASE + "-Map.png")
		plt.imshow(img, zorder=-100, extent=[ mapbounds['minlon'],mapbounds['maxlon'],mapbounds['minlat'],mapbounds['maxlat'] ])
	else:
		print ("warning - '" + base + "-Map.png' file not found")

	renderingRules = Settings['RenderRules']
	for ID in links.keys():
		if links[ID]['type'] in renderingRules.keys():
			thisRendering = renderingRules[links[ID]['type']]
		else:
			thisRendering = renderingRules['default']
		ln, = plt.plot(
			[ links[ID]['coordinates'][0][0], links[ID]['coordinates'][1][0] ],
			[ links[ID]['coordinates'][0][1], links[ID]['coordinates'][1][1] ],
			marker			= LinkEndMark,
			linestyle		= thisRendering['linestyle'],
			linewidth		= thisRendering['linewidth'],
			color			= LinksColor[ID],
			solid_capstyle	= 'round',
			solid_joinstyle = 'round',
			zorder			= thisRendering['zorder'],
			picker = 5)
		ln.name = 'line:' + str(ID)
	plt.colorbar(cb) #, orientation='horizontal')
	plt.axis('off')
	fig.canvas.mpl_connect('pick_event', onpickMap)
	SaveFigurePNG()

### Save/Load results ####################################################################

def SaveFigurePNG():
	"""save figures into png files"""
	if figs2png:
		if verbose: print ("\tsaving figure " + str(SaveFigurePNG.FigNum) + " as png...")
		plt.savefig(BASE + "-Results-Fig" + str(SaveFigurePNG.FigNum) + ".png",dpi=Settings['DPI_SETTING'], bbox_inches='tight')
		SaveFigurePNG.FigNum += 1
SaveFigurePNG.FigNum = 1

def SaveResults():
	"""save results to json/geojson files"""

	if savememdump:
		if verbose: print ("saving memdump to " + base + '-Results-MemDump.json')
		global MemDump
		MemDump['eigenvalues'] = 	results['eigenvalues'].real.tolist()
		MemDump['eigenvectors'] = 	results['eigenvectors'].real.T.tolist()
		MemDump['StationaryProb'] = results['StationaryProb'].tolist()
		MemDump['Clusters'] = 		results['Clusters'].tolist()
		MemDump['VehicleDensity'] = results['VehicleDensity'].tolist()
		MemDump['K'] = 				results['KemenyConst']
		MemDump['MFPT'] = [ str([", ".join('{:05.0f}'.format(i) for i in results['MFPT'].tolist()[j])]).replace("'","").replace("[","").replace("]","") for j in range(results['MFPT'].shape[0]) ]
		with open(BASE + '-Results-MemDump.json', 'w') as outfile:
			json.dump(MemDump, outfile, indent=4, sort_keys=True)

	if savegeojson:
		if verbose: print ("saving results to " + base + '-Results.geojson')
		for k in links.keys():
			GeoJSON['features'][k]['properties']['r.Pi']= results['StationaryProb'][k]
			GeoJSON['features'][k]['properties']['r.C'] = results['Clusters'][k]
			GeoJSON['features'][k]['properties']['r.D'] = results['VehicleDensity'][k]
		with open(BASE + '-Results.geojson', 'w') as outfile:
			json.dump(GeoJSON, outfile, indent=4, sort_keys=True)

### Script Management ####################################################################

def RunScript():
	"""standard script execusion logic"""
	LoadDataFiles()
	SetupTransitionMatrix()
	ModifyForLinkTravelTime()
	SolveMC()
	SaveResults()
	if showfigs: plt.show()

### Main #################################################################################

if __name__ == "__main__":
	print ("runing as an python script...")
	ProcessCLI()
	RunScript()
else:
	# print ("mcta: runing as a python module...")
	pass
