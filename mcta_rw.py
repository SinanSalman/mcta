#!/usr/bin/env python

"""
Markov Chains Traffic Assignment (MCTA) read/write module

Usage:
	(Settings, GeoJSON) = LoadDataFiles(base)
	SaveResults(results)

_author = 	"Sinan Salman (sinan.salman@zu.ac.ae)"
_version = 	"Revision: 0.18"
_date = 	"Date: 2019/07/17"
_copyright= "Copyright (c)2017-2019 Sinan Salman"
_license =	"GPLv3"
"""

### Initialization #######################################################################

import os as _os
import sys as _sys
import json as _json
import numpy as _np

# global variables
_config_file = 'mcta_vis.json'
_BASE = ''
_base = ''

### JSON Serialization of complex numbers ################################################

class ComplexEncoder(_json.JSONEncoder):
   def default(self, obj):
       if isinstance(obj, complex):
           return [obj.real, obj.imag, 'j']
       return _json.JSONEncoder.default(self, obj) # Let the base class default method raise the TypeError


### Data processing ######################################################################

def LoadDataFiles(base):
	"""load data and configuration files
	LoadDataFiles(base)
		base	input/output files prefix (.json .geojson .png)
	Results:
		Settings	visualization settings (CMAP, Figure_Size, RenderingRules, DPI_Setting)
		GeoJSON		road network dictionary loaded from GeoJSON file (for links, turning probabilities, lenght, lanes, etc.)
		"""
	global _BASE
	global _base

	path = _os.path.split(_os.path.realpath(__file__))[0] + _os.path.sep
	_base = base
	_BASE = path + base
	print ("\tlooking in: " + path + "")
	print (f'\treading configurations and rules from {_config_file} file...')

	if _os.path.isfile(_config_file):
		with open(_config_file) as mctajsonfile:
			Settings = _json.load(mctajsonfile)
	else:
		print (f"\nerror - can't find {_config_file}.")
		_sys.exit(0)

	if _os.path.isfile(_BASE + '-Map.geojson'):
		print ('\treading ' + base + '-Map.geojson')
		with open(_BASE + '-Map.geojson') as geojsonfile:
			GeoJSON = _json.load(geojsonfile)
	else:
		print ("\nerror - can't find JSON file for the specified BaseName: " + base)
		_sys.exit(0)
	
	return Settings, GeoJSON


### Save/Load results ####################################################################

def SaveResults(Results):
	"""save results to json file"""
	print ("\tsaving memdump to " + _base + '-Results.json')
	MemDump = {}
	for k,v in Results.items():
		if k in ['P_org','P_updated','StepTime','KemenyConst','Message','TotalNetworkEmissionCost']:
			MemDump[k] = v
		elif k in ['eigenvalues','StationaryProb','Clusters','Density','linkIDs','EmissionCost','Speeds']:
			MemDump[k] = v.real.tolist()
		elif k in ['eigenvectors']:
			MemDump[k] = v.real.T.tolist()
		elif k in ['MFPT']:
			MemDump[k] = [ str([", ".join('{:05.0f}'.format(i) for i in v.tolist()[j])]).replace("'","").replace("[","").replace("]","") for j in range(v.shape[0]) ]
		else:
			print(f'Unrecognized result {k}. it will not be included in the saved results.')	
	with open(_BASE + '-Results.json', 'w') as outfile:
		_json.dump(MemDump, outfile, indent=4, sort_keys=True, cls=ComplexEncoder)


### Main #################################################################################

if __name__ == "__main__":
	print(f'{_os.path.basename(__file__)} is designed for use as module, please see the doc string for usage.')