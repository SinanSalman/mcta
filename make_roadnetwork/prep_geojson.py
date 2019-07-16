"""
Create map PNG, JSON and GeoJSON files from line network representation GeoGJSON file. The 
model processes a simple one-way line network geoJSON export from qGIS. The output GeoJSON
can be further fine-tuned in qGIS, and ultimately use in MCTA. The model will:
	- add two way roads (with displacement)
	- add link attributes (for editing in qGIS and use in MCTA)
	- add nodes info (for MCTA)
	- doanload a basemap for the roadnetwork
	- apply edits from the configuration file

Syntax:
	ipython prep_geojson.py -- [options] {BaseName}

options:
	--verbose	|	-v				display verbose messages and detailed results
	--basemap	|	-b				create basemap

Input data files to load (must be in the script folder): 
	{BaseName}-config.json			configuration (settings and rules) JSON file
	{BaseName}-Lines.geojson		qGIS lines GeoJSON export for the Road network 
									 (single link features)
Output file:
	{BaseName}-Map.geojson			GeoSJSON file with added link attributes and two-way
									 roads for analysis in qGIS and input into mcta.py
	{BaseName}-Map.json				Road network data in JSON format
	{BaseName}-Map.png				Basemap for the road network
	{BaseName}-Map.points.geojson	Basemap for the road network

Complete Code run-through test:
	ipython ./prep_geojson.py --  --verbose --basemap {BaseName}
	ipython ./prep_geojson.py --  -v -b {BaseName}	
	
Code by Sinan Salman, 2016
"""

HELP_MSG = "Options:\n\
\t--verbose  | -v  display verbose messages and detailed results\n\
\t--basemap  | -b  create basemap\n"

__author__ = 	"Sinan Salman (sinan.salman@zu.ac.ae)"
__version__ = 	"Revision: 0.1"
__date__ = 		"Date: 2016/08/29"
__copyright__ = "Copyright (c)2016 Sinan Salman"
__license__ = 	"GPLv3"

### Initialization #######################################################################

import os
import sys	
import json
import math 

# options
verbose = False
basemap = False

# global variables
path = ""
base = ""
BASE = ""
nodes = {}
links = {}
bounds = {'minlon': 0,'maxlon': 0,'minlat': 0,'maxlat': 0}
mapbounds = {'minlon': 0,'maxlon': 0,'minlat': 0,'maxlat': 0}
GeoJSON = {}
Settings = {}
Edits = {}
connectedlinks = {}
connectedlinksProb = {}

### Data processing ######################################################################

def ProcessCLI():
	"""Process CLI parameters"""
	global base
	global verbose
	global basemap

	if len(sys.argv) == 1:
		print("Missing argument\n\nSyntax:\n\tipython prep_geojson.py -- [options] {BaseName}")
		print(HELP_MSG)
		sys.exit(0)

	if '--verbose' in sys.argv or '-v' in sys.argv:
		print ("*** option: verbose mode")
		verbose = True
	if '--basemap' in sys.argv or '-b' in sys.argv:
		print ("*** option: create basemap")
		basemap = True
	base = sys.argv[len(sys.argv)-1]
 
def LoadDataFiles():
	"""load data and configuration files"""
	global path
	global BASE
	global GeoJSON
	global Settings
	global Edits

	print ('data loading and processing:')
	path = os.path.split(os.path.realpath(__file__))[0] + os.path.sep
	BASE = path + base
	if verbose: print ("\tlooking in: " + path + "")
	
	if verbose: print ('\treading configurations (settings/rules) file:' + base + '-Config.json')
	if os.path.isfile(BASE + '-Config.json'):
		with open(BASE + '-Config.json') as configjsonfile:
			data = 		json.load(configjsonfile)
			Settings =	data['Settings']
			Edits =		data['Edits']
	else:
		print ("\nerror - can't find '" + base + "-Config.json' file.")
		sys.exit(0)

	if os.path.isfile(BASE + '-Lines.geojson'):
		if verbose: print ('\treading GeoJSON input file: ' + base + '-Lines.geojson')
		with open(BASE + '-Lines.geojson') as Mapgeojsonfile:
			GeoJSON = json.load(Mapgeojsonfile)
			process_GeoJSON(GeoJSON)
	else:
		print ("\nerror - can't find JSON file for the specified BaseName: " + base)
		sys.exit(0)	

def process_GeoJSON(geojson):
	"""convert geojson file to json, assuming all lines are made of exactly two points (no multi-segment lines)"""
	maxlon = -1000 # large numbers to initialize
	maxlat = -1000
	for line in geojson['features']:
		assert len(line['geometry']['coordinates'])==2, 'multi-segment lines are not implemented. feature#' + str(int(lcounter/2))
		mid_lon = (line['geometry']['coordinates'][0][0]+line['geometry']['coordinates'][1][0])/2
		mid_lat = (line['geometry']['coordinates'][0][1]+line['geometry']['coordinates'][1][1])/2
		line['geometry']['midpoint'] = [mid_lon,mid_lat]
		if mid_lon > maxlon:# and mid_lat > maxlat: 
			maxlon = mid_lon
			maxlat = mid_lat
	for line in geojson['features']:
		line['geometry']['SortVal'] = distance_on_earth_est1(maxlat, maxlon, line['geometry']['midpoint'][1], line['geometry']['midpoint'][0])
	global links
	lcounter = 0
	for line in sorted(geojson['features'], key=lambda k: k['geometry']['SortVal']):
		ptA = addnode(line['geometry']['coordinates'][0])
		ptB = addnode(line['geometry']['coordinates'][1])
		RoadType = line['geometry']['type']
		speed = Settings['Defaults']['RoadSpeed']
		lanes = Settings['Defaults']['RoadLanes']
		length = distance_on_earth_est1(nodes[ptA]['lat'],nodes[ptA]['lon'],nodes[ptB]['lat'],nodes[ptB]['lon'])
		assert length>0, 'Line with zero length @lat:' + str(nodes[ptA]['lat']) + ',lon:' + str(nodes[ptA]['lon']) + ' feature#' + str(int(lcounter/2))
		assert speed>0, 'link with speed=zero. feature#' + str(int(lcounter/2))
		links[lcounter] = {'type':RoadType, 'org':ptA, 'dst':ptB, 'length':length, 'speed': speed, 'lanes': lanes, 'edits': ""}
		links[lcounter+1] = {'type':RoadType, 'org':ptB, 'dst':ptA, 'length':length, 'speed': speed, 'lanes': lanes, 'edits': ""}
		lcounter += 2 

def addnode(pt):
	"""find existing node from coordinates or create a new one and return its referance"""
	global bounds
	global nodes
	for k,v in nodes.items():
		if v['lat'] == pt[1] and v['lon'] == pt[0]:
			return k
	if addnode.ncounter == 1:
		bounds['minlon'] = pt[0]
		bounds['maxlon'] = pt[0]
		bounds['minlat'] = pt[1]
		bounds['maxlat'] = pt[1]
	else:
		if bounds['minlon'] > pt[0]: bounds['minlon'] = pt[0]
		if bounds['maxlon'] < pt[0]: bounds['maxlon'] = pt[0]		
		if bounds['minlat'] > pt[1]: bounds['minlat'] = pt[1]
		if bounds['maxlat'] < pt[1]: bounds['maxlat'] = pt[1]
	nodes[addnode.ncounter] = {'lat':pt[1], 'lon':pt[0]}
	addnode.ncounter += 1
	return addnode.ncounter - 1
addnode.ncounter = 0

def distance_on_earth_est1(lat1, long1, lat2, long2):
	"""Convert latitude and longitude to km using the great circule equation
	adapted from http://www.johndcook.com/blog/python_longitude_latitude"""
	degrees_to_radians = math.pi/180.0
 
	# phi = 90 - latitude
	phi1 = (90.0 - lat1)*degrees_to_radians
	phi2 = (90.0 - lat2)*degrees_to_radians
 
	# theta = longitude
	theta1 = long1*degrees_to_radians
	theta2 = long2*degrees_to_radians
 
	# Compute spherical distance from spherical coordinates.
	# For two locations in spherical coordinates
	# (1, theta, phi) and (1, theta', phi')
	# cosine( arc length ) =
	# sin phi sin phi' cos(theta-theta') + cos phi cos phi'
	# distance = rho * arc length 
	cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) +
	math.cos(phi1)*math.cos(phi2))
	arc = math.acos( cos )
 
	# The distance returned is relative to Earth’s radius. 
	# To get the distance in miles, multiply by 3959. 
	# To get the distance in kilometers, multiply by 6371.
	return arc * 6371

def distance_on_earth_est2(lat1, lon1, lat2, lon2):
	"""Convert latitude and longitude to km using rough estimates
	based on: http://gis.stackexchange.com/questions/5821/calculating-latitude-longitude-x-miles-from-point"""
	lat2km = 111.111
	dy = (lat2-lat1) * lat2km
	dx = (lon2-lon1) * math.cos(lat1) * lat2km 
	return math.sqrt(dx**2 + dy**2)

def test_distance_on_earth(): # for comparison and testing purposes only
	"""Test and compare the two estimation methods for earth distance using a ~20KM range"""
	est1 = distance_on_earth_est1(54.3181,24.3956,54.4888,24.5193)
	est2 = distance_on_earth_est2(54.3181,24.3956,54.4888,24.5193)
	print ('est1=',est1,'\nest2=',est2,'\ndiff=',(est1-est2)/est1*100,'%')
			
### Generate turning probabilities #######################################################

def BuildNetwork():
	"""Build road network by finding connections and calculating turning probabilities based on turn rules assumptions."""
	global connectedlinks
	global connectedlinksProb
	if verbose: print ('\nNetwork parsing detail:')
	for l1 in links.keys():
		connectedlinks[l1] = {'all':[], 'right':[],'forward':[],'left':[],'u-turn':[]}
		for l2 in links.keys():
			if links[l1]['dst'] == links[l2]['org']:
				connectedlinks[l1][assignTurnType(l1,l2)].extend([l2])
				connectedlinks[l1]['all'].extend([l2])
	ApplyEdits() # apply network edits, including: no turns, adding/removing, speed, and # of lanes
	if verbose: print ('\tCalculating turn probabilities:')
	turnProb = 	Settings['TurnProbability']		 
	for l1 in links.keys():
		if verbose: print ('\t\tlink:', l1, ' connected to:')
		connectedlinksProb[l1] = {'right':[],'forward':[],'left':[],'u-turn':[]}
		totalprob = sum([v for k,v in turnProb.items() if connectedlinks[l1][k] != []])
		assert totalprob > 0,'error - unconnected link:{:}'.format(l1)
		assert totalprob <= 1,'error - total probability more than 1.0 for link:{:}'.format(l1)
		for t in turnProb.keys():
			turnlinks = connectedlinks[l1][t]
			num = len(turnlinks)
			if verbose and num>0: print ('\t\t\t{:>8}:{:10}\t(Tot.prob:{:.5f})'.format(t,str(turnlinks),turnProb[t]/totalprob))
			for tl in turnlinks:
				connectedlinksProb[l1][t].extend([turnProb[t]/num/totalprob])
		assert round(sum([Pr for (k,v) in connectedlinksProb[l1].items() for Pr in v] ), 8) == 1, 'error - total probability sum not equal to 1 (link:{:})'.format(l1) 

def assignTurnType(l1,l2):
	"""calculate the angle between two lines and assign a turn type accordingly"""
	alpha1 = math.atan2(nodes[links[l1]['dst']]['lat'] - nodes[links[l1]['org']]['lat'],
						nodes[links[l1]['dst']]['lon'] - nodes[links[l1]['org']]['lon'])
	alpha2 = math.atan2(nodes[links[l2]['dst']]['lat'] - nodes[links[l2]['org']]['lat'],
						nodes[links[l2]['dst']]['lon'] - nodes[links[l2]['org']]['lon'])
	alpha = int((alpha2 - alpha1) / math.pi * 180) # get difference between line slope angles and convert it to degrees
	alpha -= 180 # make u-turn a zero (360) angle for simplicity 
	while alpha <=1: # (1 degree tolerance)
		alpha += 360 # make sure alpha is in the range (5,365]
	if 	alpha > Settings['TurnAssignment']['right'][0] and \
		alpha <= Settings['TurnAssignment']['right'][1]:
		return 'right'
	if 	alpha > Settings['TurnAssignment']['forward'][0] and \
		alpha <= Settings['TurnAssignment']['forward'][1]:
		return 'forward'
	if 	alpha > Settings['TurnAssignment']['left'][0] and \
		alpha <= Settings['TurnAssignment']['left'][1]:
		return 'left'
	if 	alpha > Settings['TurnAssignment']['u-turn'][0] and \
		alpha <= Settings['TurnAssignment']['u-turn'][1]:
		return 'u-turn'
	assert False,"error - could not determin turnType for links: " + l1 + " and " + l2 #if execution gets here there is a problem

def ApplyEdits():
	"""Apply edits into road network."""
	global connectedlinks
	print ('\tApplying network edits from configuration file:')
	turns = Settings['TurnProbability'].keys()
	#check edits types for errors
	for E in Edits.keys():
		if E in ['no_' + t + '_junc' for t in turns] + ['no_left_T_junc','add','remove']: pass
		elif E[:8] == "maxspeed": pass
		elif E[-5:] == "lanes": pass
		else: assert False, "error - unknown edit type in edits json file: " + E
	EditsDone = {}
	for id in links.keys():
		EditsDone[id] = ""

	for E in Edits.keys(): # remove junctions turns
		if E in ['no_' + t + '_junc' for t in turns]:
			for id in links.keys():
				turn = E[3:-5]
				if links[id]['dst'] in Edits[E]:
					if connectedlinks[id][turn] != []:
						if verbose: print ('\t\tremove {:} turn @node:{:}. links{:}->{:}'.format(turn,links[id]['dst'],id,str(connectedlinks[id][turn])))
						connectedlinks[id][turn]=[]
						EditsDone[id] += ",-{:<8}".format(turn)
					else:
						print ('warning -  no {:} link found, check network data. @node:{:} link:{:}'.format(turn,links[id]['dst'],id))
		if E == 'no_left_T_junc':
			for id in links.keys():
				if links[id]['dst'] in Edits[E]:
					if connectedlinks[id]['forward'] == [] and connectedlinks[id]['right'] != []:
						if verbose: print ('\t\tremove left turn in (T-Junc:{:}) {:}->{:}'.format(links[id]['dst'],id,str(connectedlinks[id]['left'])))
						connectedlinks[id]['left']=[]
						EditsDone[id] += ",-{:<8}".format('left_T')
					elif connectedlinks[id]['forward'] != [] :
						if verbose: print ('\t\tremove    U turn in (T-Junc:{:}) {:}->{:}'.format(links[id]['dst'],id,str(connectedlinks[id]['u-turn'])))
						connectedlinks[id]['u-turn']=[]
						if connectedlinks[id]['left'] != []:
							if verbose: print ('\t\tremove left turn in (T-Junc:{:}) {:}->{:}'.format(links[id]['dst'],id,str(connectedlinks[id]['left'])))
							connectedlinks[id]['left']=[]
						EditsDone[id] += ",-{:<8}".format('left_T')
					else:
						print ('warning -  problem with T-junc:{:} in link {:}, check network data'.format(links[id]['dst'],id))
		if E[:8] == "maxspeed":
			speed = int(E[-3:])
			for path in Edits[E]:
				l = [k for k,v in links.items() if v['org'] in path and v['dst'] in path]
				for id in l: 
					links[id]['speed'] = speed
					EditsDone[id] += ",=Speed{:<3}".format(speed)
				if verbose: print ('\t\tadjust max speed to:{:} for links:{:}'.format(speed,str(l)))
		if E[-5:] == "lanes":
			lanes = int(E[0])
			for path in Edits[E]:
				l = [k for k,v in links.items() if v['org'] in path and v['dst'] in path]
				for id in l: 
					links[id]['lanes'] = lanes
					EditsDone[id] += ",={:}Lanes  ".format(lanes)
				if verbose: print ('\t\tadjust # of lanes to:{:} for links:{:}'.format(lanes,str(l)))
	if "remove" in Edits.keys(): # remove turns
		for seq in Edits["remove"]:
			found = False
			for l1 in links.keys():
				if links[l1]['org'] == seq[0] and links[l1]['dst'] == seq[1]:
					conn = connectedlinks[l1]['all']
					for l2 in conn:
						if links[l2]['org'] == seq[1] and links[l2]['dst'] == seq[2]:
							for t in turns:
								if l2 in connectedlinks[l1][t]:
									connectedlinks[l1][t].remove(l2)
									connectedlinks[l1]['all'].remove(l2)
									found = True
									if verbose: print ('\t\tremoved turn {:}->{:} ({:})'.format(l1,l2,seq))
									EditsDone[l1] += ",-link{:<8}".format(str(l2))
									break
							break
					break
			if not found: print ("warning - could not find turn in remove edit: " + str(seq))
	if "add" in Edits.keys(): # add turns
		for seq in Edits["add"]:
			found = None
			for l1 in links.keys():
				if links[l1]['org'] == seq[0] and links[l1]['dst'] == seq[1]:
					conn = connectedlinks[l1]['all']
					for l2 in conn:
						if links[l2]['org'] == seq[1] and links[l2]['dst'] == seq[2]:
							found = l2
							break
					if found != None:
						ttype = assignTurnType(l1,l2)
						connectedlinks[l1][ttype].append(l2)
						connectedlinks[l1]['all'].append(l2)
						if verbose: print ('\t\tadded {:} turn {:}->{:} ({:})'.format(ttype,l1,l2,seq))
						EditsDone[l1] += ",+link{:<8}".format(l2)
					else:
						print ("warning - could not add turn as it already exists: " + str(seq))
					break
	for id in links.keys():
		links[id]['edits'] = EditsDone[id].strip(", ")
	if verbose: print ('\n\tedits summary:\n\t\t{:}'.format(str(EditsDone).replace("', ","\n\t\t ").replace("'","")))
	
def CheckNetworkReducability(link=None, network=None):
	"""check network reducability; will error if reducable, and stay silent if not reducable"""
	if link != None and network == None:
		if links[link]['dst'] in CheckNetworkReducability.RemNodes:
			CheckNetworkReducability.log += '\tgoing down link (' + str(link) + ') to node (' + str(links[link]['dst']) + ')\n'
			CheckNetworkReducability.RemNodes.remove(links[link]['dst'])
			for l in CheckNetworkReducability.Network[link]['all']:
				CheckNetworkReducability(link=l)
		else: 
			CheckNetworkReducability.log += '\tignoring link (' + str(link) + ') as node (' + str(links[link]['dst']) + ') was already visited\n'
	elif link == None and network != None:
		if verbose: print('\nChecking network reducability:')
		CheckNetworkReducability.Network = network
		CheckNetworkReducability.RemNodes = list(nodes.keys())
		CheckNetworkReducability.RemNodes.remove(0)
		CheckNetworkReducability.log += '\tStarting at node (0)\n'
		startinglinks = [k for k,v in links.items() if v['org'] == 0] # start with a list of links that originate from node 0
		for l in startinglinks:
			CheckNetworkReducability(link=l)
		assert len(CheckNetworkReducability.RemNodes) == 0, "Error - network is reducable, the following nodes are not connected to the main network: " + str(CheckNetworkReducability.RemNodes) + "\n\twalking log:\n" + CheckNetworkReducability.log
#		if verbose: print ("\nwalking log:\n" + CheckNetworkReducability.log)
		print('\tWalked the entire network with no nodes remaining.\n\t\t==> success: network is not reducable!')
	else:
		assert False, "Error - CheckNetworkReducability function used incorrectly."
CheckNetworkReducability.Network = {}
CheckNetworkReducability.RemNodes = []
CheckNetworkReducability.log = ""
	 
### Get basemap for the resulting bounding box ###########################################

def GetBaseMap():
	import smopy
	global mapbounds
	print ("downloading basemap: " + base + '-map.png')
	lonRange = abs(bounds['maxlon'] - bounds['minlon'])
	latRange = abs(bounds['maxlat'] - bounds['minlat'])
	mapbounds['minlon'] = bounds['minlon'] - lonRange * 0.1
	mapbounds['maxlon'] = bounds['maxlon'] + lonRange * 0.1
	mapbounds['minlat'] = bounds['minlat'] - latRange * 0.1
	mapbounds['maxlat'] = bounds['maxlat'] + latRange * 0.1
	smopy.TILE_SERVER = "http://tile.basemaps.cartocdn.com/light_all/{z}/{x}/{y}@2x.png"
	#smopy.TILE_SERVER = "http://tile.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}@2x.png"
	smopy.TILE_SIZE = 512
	map = smopy.Map((mapbounds['minlat'],mapbounds['minlon'],mapbounds['maxlat'],mapbounds['maxlon']), z=12, margin=0)
	x1, y1 = map.to_pixels(mapbounds['maxlat'],mapbounds['minlon'])
	x2, y2 = map.to_pixels(mapbounds['minlat'],mapbounds['maxlon'])
	img = map.img.crop((x1,y1,x2,y2))
	img.save(BASE + '-Map.png')

### Save/Load results ####################################################################

def SaveResults():
	"""save results to json/geojson files"""
	print ("saving to " + base + '-Map.geojson')
	GeoJSON['features'] = []
	lonRange = abs(bounds['maxlon'] - bounds['minlon'])
	latRange = abs(bounds['maxlat'] - bounds['minlat'])
	skew = latRange/lonRange	#could need change if the map is too thin in one dimension, works for now!
	displacement = min([lonRange, latRange]) * Settings['twowaydisplacement']
	for (k,v) in links.items():
		Ax = nodes[v['org']]['lon']
		Ay = nodes[v['org']]['lat']
		Bx = nodes[v['dst']]['lon']
		By = nodes[v['dst']]['lat']
		alpha = math.atan2(By-Ay,Bx-Ax) - math.pi/2
		Dx = displacement * math.cos(alpha) * skew
		Dy = displacement * math.sin(alpha)
		Ax += Dx
		Ay += Dy
		Bx += Dx
		By += Dy
		
		GeoJSON['features'].append({ 
			"type": "Feature", 
			"properties": { 
				"l.id": 		k,
				"edits":		v['edits'],
				"l.length": 	v['length'],
				"l.lanes":		v['lanes'],
				"l.speed": 		v['speed'],
				"l.type": 		v['type'],
				"p.org":		v['org'],
				"p.dst":		v['dst'],
				"g.forward":	connectedlinks[k]['forward'],
				"g.forwardP":	connectedlinksProb[k]['forward'],
				"g.left":		connectedlinks[k]['left'],
				"g.leftP":		connectedlinksProb[k]['left'],
				"g.right":		connectedlinks[k]['right'],
				"g.rightP":		connectedlinksProb[k]['right'],
				"g.u-turn":		connectedlinks[k]['u-turn'],
				"g.u-turnP":	connectedlinksProb[k]['u-turn']
				}, 
			"geometry": { 
				"type": v['type'], 
				"coordinates": [[Ax,Ay],[Bx,By]] 
				} 
			})
	with open(BASE + '-Map.geojson', 'w') as outfile:
		json.dump(GeoJSON, outfile, indent=4, sort_keys=True)		

	print ("saving to " + base + '-Map.points.geojson')
	GeoJSON['features'] = []
	for (k,v) in nodes.items():
		GeoJSON['features'].append({ 
			"type": "Feature", 
			"properties": {
				"id": k }, 
			"geometry": { 
				"type": "Point", 
				"coordinates": [ nodes[k]['lon'], nodes[k]['lat'] ] 
				} 
			})
	with open(BASE + '-Map.points.geojson', 'w') as outfile:
		json.dump(GeoJSON, outfile, indent=4, sort_keys=True)		

	data = {'VehiclesCountEst':	Settings['VehiclesCountEst'],
	        'bbox_roads':		bounds,
	        'bbox_map':			mapbounds,
	        'nodes':			nodes,
	        'links':			links}
	if verbose: print ("saving to " + base + '-Map.json')
	with open(BASE + '-Map.json', 'w') as outfile:
		json.dump(data, outfile, indent=4, sort_keys=True)
	if verbose: print ("\t" + str(len(nodes)) + " nodes")
	if verbose: print ("\t" + str(len(links)) + " links")

### Main #################################################################################

if __name__ == "__main__":
	print ("Starting...")
	ProcessCLI()
	LoadDataFiles()
	BuildNetwork()
	CheckNetworkReducability(network=connectedlinks)
	if basemap: GetBaseMap()	
	SaveResults()
else:
	print ("prep_geojson.py is not designed to run as a module.")
