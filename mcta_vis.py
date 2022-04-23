#!/usr/bin/env python

"""
Markov Chains Traffic Assignment (MCTA) visualization module

Usage:
	Initialize(GeoJSON, Settings, Base = 'MCTA', ShowFigs=True, SaveFigs2PNG=False)
	Generate_Figure(Variable, Limit=0)
	ShowModificationOnMap(EditedLinks)
Results:
	visualization figure on screen or saved to PNG files

_author = 	"Sinan Salman (sinan.salman@zu.ac.ae)"
_version = 	"Revision: 0.17"
_date = 	"Date: 2019/06/11"
_copyright= "Copyright (c)2017-2019 Sinan Salman"
_license =	"GPLv3"
"""

### Initialization #######################################################################

import mcta as _mcta
import scipy as _sp
import matplotlib.pyplot as _plt
import os as _os

# global variables
_BASE = None
_CMAP = None
_Settings = {}
_showfigs = False
_figs2png = False
_bounds = None
_mapbounds = None

def Initialize(GeoJSON, Settings, Base = 'MCTA', ShowFigs=True, SaveFigs2PNG=False):
	"""Initialize can configure visualization module
	Initialize(GeoJSON, Settings, Base = 'MCTA', ShowFigs=True, SaveFigs2PNG=False)
		GeoJSON			roadnetwork dictionary loaded from GeoJSON file (for map bounding box)
		Settings		figure settings (CMAP, Figure_Size, RenderingRules, DPI_Setting)
		Base			base name for figure files
		ShowFigs		show figures on screen
		SaveFigs2PNG	save figures to PNG files
	"""
	global _BASE
	global _showfigs
	global _figs2png
	global _bounds
	global _mapbounds
	global _Settings
	global _CMAP
	_BASE = Base
	_showfigs = ShowFigs
	_figs2png = SaveFigs2PNG
	_bounds = GeoJSON['mcta_json']['bbox_roads']
	_mapbounds = GeoJSON['mcta_json']['bbox_map']
	_Settings = Settings
	_CMAP = _plt.get_cmap(Settings['CMAP'])
	_Show_Save_Fig.FigNum = 1


### Plotting #############################################################################

def _Plot2D(X, Y, Xlabel, Ylabel, Title, Limit=0, Yaxis_lim=[]):
	"""plot 2D points"""
	import matplotlib.colors as colors
	if Limit > 0:
		DotColor = [_CMAP(int(Y[i]/Limit*254)) if Y[i]<Limit else _CMAP(254) for i in range(Y.shape[0])]
	else:
		DotColor = [_CMAP(int((Y[i]-Y.min())/(Y.max()-Y.min())*254)) for i in range(Y.shape[0])]
	fg, ax = _plt.subplots(figsize=_Settings['FIGURE_SIZE'])
	ax.grid(linestyle='dotted')
	ax.set_xlabel(Xlabel)
	ax.set_ylabel(Ylabel)
	ax.set_title(Title)
	# ax.set_xlim([X.min(),X.max()])
	if Yaxis_lim != []:
		ax.set_ylim(Yaxis_lim)
	ax.scatter(X,Y,color=DotColor)
	_Show_Save_Fig()
	return fg, ax


def _Plot3D(Z, Title):
	"""plot 3D surface and countour plots"""
	from mpl_toolkits.mplot3d.axes3d import Axes3D
	import matplotlib.colors as colors
	fg1 = _plt.figure(figsize=_Settings['FIGURE_SIZE'])
	ax1 = Axes3D(fg1)
	ax1.set_ylabel('from') # rows
	ax1.set_xlabel('to')   # columns
	ax1.set_zlabel('time(steps)')
	ax1.set_title(Title)
	X, Y = _sp.meshgrid(range(Z.shape[0]),range(Z.shape[0]))
	ax1.view_init(30, 240)
	ax1.plot_surface(X,Y,Z,cmap=_CMAP)
	_Show_Save_Fig()

	fg2, ax2 = _plt.subplots(figsize=_Settings['FIGURE_SIZE'])
	ax2.grid(True)
	ax2.set_title(Title)
	ax2.set_ylabel('Starting Road ID') # rows
	ax2.set_xlabel('Ending Road ID')   # columns
	CF = ax2.contourf(Z.A,cmap=_CMAP)  #,norm=colors.LogNorm())
	_plt.colorbar(CF)
	_Show_Save_Fig()
	return fg1, ax1, fg2, ax2


### Render road network ##################################################################

def _onpickMap(event):
	"""roadmap onpick event handler"""
	thisline = event.artist
	thisline.set_visible(not thisline.get_visible())
	thisline.figure.canvas.draw()
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	points = tuple(zip(xdata[ind], ydata[ind]))
	print('picked object:', thisline.name, ' onpick points:', points)


def _RenderMap(X, Title='RoadNetwork', LinkEndMark='s', Limit=0):
	"""render road map"""
	if Limit > 0:
		LinksColor = [_CMAP(int(X[i]/Limit*254)) if X[i]<Limit else _CMAP(254) for i in range(X.shape[0])]
	else:
		LinksColor = [_CMAP(int((X[i]-X.min())/(X.max()-X.min())*254)) for i in range(X.shape[0])]

	# add a 0.05 margin to the plot
	bbounds={}
	lonRange = abs(_bounds['maxlon'] - _bounds['minlon'])
	latRange = abs(_bounds['maxlat'] - _bounds['minlat'])
	bbounds['minlon'] = _bounds['minlon'] - lonRange * 0.05
	bbounds['maxlon'] = _bounds['maxlon'] + lonRange * 0.05
	bbounds['minlat'] = _bounds['minlat'] - latRange * 0.05
	bbounds['maxlat'] = _bounds['maxlat'] + latRange * 0.05

	fig = _plt.figure(figsize=_Settings['FIGURE_SIZE'])

	# Using contourf to provide colorbar info, then clear the figure
	if Limit > 0:
		Z = [[0,0],[0,0]]
		levels = [x for x in _sp.arange(0,Limit*1.1,Limit/10)]
	else:
		step = (X.max()-X.min())/10
		Z = [[0,0],[0,0]]
		levels = [x for x in _sp.arange(X.min(),X.max()+step,step)]
	cb = _plt.contourf(Z, levels, cmap=_CMAP)
	_plt.clf()

	# by setting limits before hand, plotting is about 3 times faster
	ax = fig.add_subplot(111,autoscale_on=False,
						 xlim=(bbounds['minlon'],bbounds['maxlon']),
						 ylim=(bbounds['minlat'],bbounds['maxlat']))
	ax.set_title(Title)

	if _os.path.isfile(_BASE + '-Map.png'):
		from PIL import Image
		img = Image.open(_BASE + "-Map.png")
		_plt.imshow(img, zorder=-100, extent=[ _mapbounds['minlon'],_mapbounds['maxlon'],_mapbounds['minlat'],_mapbounds['maxlat'] ])
	else:
		print ("warning - '" + _BASE + "-Map.png' file not found")

	renderingRules = _Settings['RenderRules']
	links = _mcta._links
	for ID in links.keys():
		if links[ID]['type'] in renderingRules.keys():
			thisRendering = renderingRules[links[ID]['type']]
		else:
			thisRendering = renderingRules['default']
		ln, = _plt.plot(
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
	_plt.colorbar(cb) #, orientation='horizontal')
	_plt.axis('off')
	fig.canvas.mpl_connect('pick_event', _onpickMap)
	_Show_Save_Fig()


### Save figures #########################################################################

def _Show_Save_Fig():
	"""save figures into png files"""
	if _figs2png:
		print ("\tsaving figure " + str(_Show_Save_Fig.FigNum) + " as png...")
		_plt.savefig(_BASE + "-Results-Fig" + str(_Show_Save_Fig.FigNum) + ".png",dpi=_Settings['DPI_SETTING'], bbox_inches='tight')
		_Show_Save_Fig.FigNum += 1
	if _showfigs:
		_plt.show()
_Show_Save_Fig.FigNum = 1


### API ####################################################################

def Generate_Figure(Results, Variable, Limit=0):
	info = {'StationaryProb':{'x':'Roads IDs', 'y':'Probability density', 'title':'Stationary probability distribution of vehicles'},
			'Density':{'x':'Roads IDs', 'y':'Density', 'title':'Road traffic density $\\left(\\frac{vehicles}{km \\cdot lane}\\right)$'},
			'Clusters':{'x':'Roads IDs', 'y':'Entries of the second eigenvector', 'title':'Clusters in road network'},
			'Emission':{'x':'Roads IDs', 'y':"Road's average CO emissions (g/km)", 'title':'Average CO emissions (g/km)'},
         	'EmissionCost($/hr)': {'x': 'Roads IDs', 'y': "Road Emissions External Cost ($/hr)", 'title': "Road Emissions External Costs($/hr)"},
         	'EmissionCost($/km)': {'x': 'Roads IDs', 'y': "Road Emissions External Cost ($/km)", 'title': "Road Emissions External Costs($/km)"}}

	if Variable in ['linkIDs','Message']:
		return  # ignore
	if Variable not in Results.keys():
		print (f'Variable not in MCTA results: {Variable}')
		return
	if Variable == 'MFPT':
		_Plot3D(Results['MFPT'], 'Mean first passage times (steps)')
	elif Variable in info.keys():
		_Plot2D(Results['linkIDs'], Results[Variable], info[Variable]['x'], info[Variable]['y'], info[Variable]['title'],Limit)
		_RenderMap(Results[Variable], info[Variable]['title'],'',Limit)
	elif hasattr(Results[Variable],'shape'):
		if Results[Variable].shape == Results['linkIDs'].shape:
			print(f'attempting to plot {Variable}...')
			_Plot2D(Results['linkIDs'], Results[Variable], 'Roads IDs', Variable, Variable, Limit)
			_RenderMap(Results[Variable], Variable,'',Limit)
	else:
		print (f'Visualization is not supported yet for: {Variable}')
		return


def ShowModificationOnMap(EditedLinks): # 14 and 27 are arbitrary density values for color coding only
	revlinks = _mcta._ReverseLinks
	sol=_sp.matrix([14]*len(EditedLinks)*2).T
	for i in range(len(EditedLinks)):
		if EditedLinks[i] == 1:
			sol[revlinks[i][0]]=27
		elif EditedLinks[i] == 2:
			sol[revlinks[i][1]]=27
	_RenderMap(sol, 'Road Network Changes',LinkEndMark = '', Limit=28)


### Main #################################################################################

if __name__ == "__main__":
	print(f'{_os.path.basename(__file__)} is designed for use as module, please see the doc string for usage.')
