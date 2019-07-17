#!/usr/bin/env python

"""
Markov Chains Traffic Assignment (MCTA) module

Analyze road network traffic using Descrete Time Markov Chain (DTMC)

Usage:
	SetupMCTA(JSON_Map, GeoJSON_Map, Verbose = False)
	SolveMCTA(VehiclesCount = 1)
Results are stored in:
	Results		a dictionary object with mcta results 

_author = 	"Sinan Salman (sinan.salman@zu.ac.ae)"
_version = 	"Revision: 0.19"
_date = 	"Date: 2019/07/17"
_copyright= "Copyright (c)2017-2019 Sinan Salman"
_license =	"GPLv3"
"""

### Initialization #######################################################################

import os as _os
import scipy as _sp
import scipy.sparse as _sprs
import emissions as _em

# global variables
_links = {}
_ReverseLinks = []
_verbose = False
_results = {}

### Helper function ######################################################################

def _string2list(str, TYPE='int'): # for strings of the format (#:1,2,3,4,...)
	lst = str.split(':',1)[1].strip(')').split(',')
	if lst == ['']:
		return []
	elif TYPE=='int':
		return list(map(int,lst))
	elif TYPE=='float':
		return list(map(float,lst))
	else:
		raise ValueError('error - unknown type passed to function _string2list: ' + TYPE)


### Setup MCTA ###########################################################################

def _SetupTransitionMatrix_P(GeoJSON):
	"""populate transition matrix from GeoJSON link attributes."""
	if _verbose:
		print ('\n\tConstructing Transition Matrix...')
	linkcount = len(_links)
	P = _sprs.lil_matrix((linkcount,linkcount))
	turns = ['g.right','g.left','g.forward','g.u-turn']
	FoundErrors = ""
	corrections1 = 0
	corrections2 = 0
	for item in GeoJSON['features']:
		id = item['properties']['l.id']
		if _verbose:
			print ('\tlink: {:}'.format(id))
		_links[id]['coordinates'] = item['geometry']['coordinates']
		totalprob = 0
		for t in turns:
			if isinstance(item['properties'][t],list):
				conn = item['properties'][t]
				prob = item['properties'][t+'P']
			else:
				conn = _string2list(item['properties'][t],'int')
				prob = _string2list(item['properties'][t+'P'],'float')
			num = len(conn)
			if _verbose:
				print ('\t{:>12}: {:20} Prob: {:}'.format(str(t),str(conn),str(prob)))
			for i in range(num):
				P[id,conn[i]] = prob[i]
				totalprob += prob[i]
		if totalprob != 1:
			if round(totalprob,5) == 1:
				if _verbose:
					print ('\t\tcorrecting ({:}) rounding error at link: {:}'.format(1-totalprob,id))
				corrections1 += 1
				for i in range(linkcount):
					P[id,i] /= totalprob
				#double check!
				totalprob = 0
				for i in range(linkcount):
					totalprob += P[id,i]
				remainder = 1-totalprob
				if remainder == 0:
					if _verbose:
						print ('\t\tcorrected.')
				elif -0.00001 < remainder and remainder < 0.00001:
					corrections2 += 1
					P[id,id] += remainder
					if _verbose:
						print ('\t\tcorrecting ({:}) rounding error at link: {:} AGAIN! remainder added to P[i,i]'.format(remainder,id))
				else:
					assert False,"problem while correcting rounding errors at link {:}. TotalProb = {:}".format(id,totalprob)
			else:
				print ("\t*** probabilities don't add up to 1 for link {:}, totalprob = {:}".format(id,totalprob))
				FoundErrors += str(id) + ', '
	if _verbose and corrections1 + corrections2 > 0:
		print ('\tcorrected probability rounding errors in:\n\t\t{:} links from the 1st try\n\t\t{:} links from the 2nd try\n'.format(corrections1, corrections2))
	assert FoundErrors == "",'error - please correct the above issues befor continuing for links: ' + FoundErrors[:-2]
	Psum = _sp.sum(P.todense(),axis=1)
	assert _sp.allclose(Psum,_sp.ones(Psum.shape[0])),"error - The generated turning probabilities matrix 'P' is not a stochastic matrix"
	if _verbose:
		global _results
		_results['P_org'] = [ str([", ".join('{:.3f}'.format(i) for i in P.todense().tolist()[j])]).replace("'","").replace("[","").replace("]","") for j in range(P.shape[0]) ]
	return 	P


def _ModifyForLinkTravelTime(P):
	"""normalize link travel times (ltt) and reflict them into the P probability matrix
	returns P and StepTime (time/step conversion factor)"""
	P_org = P
	if _verbose:
		print ('\nnormalize link travel times (ltt):')
	ltt = _sp.array([v['length']/v['speed'] for k,v in sorted(_links.items())]) # time = len / spd
	if _verbose:
		print ('ltt=',ltt)
	ltt_min = min(ltt)
	assert ltt_min>0,"error - minimum link travel time = zero"
	if _verbose:
		print ('ltt_min=',ltt_min)
	ltt /= ltt_min
	if _verbose: print ('ltt/ltt_min=',ltt)
	lttv = (ltt-1)/ltt 		# probability of staying; Pii
	if _verbose:
		print ('lttv=',lttv)
	lttd = _sp.diag(lttv)	# lttv as a square matrix with a diagonal of its values
	if _verbose:
		print ('lttd=',lttd)
	f = _sp.transpose((1-lttv)[_sp.newaxis]) # use _sp.newaxis to create a 2D array that can be transposed
	if _verbose:
		print ('f=',f)
	P = _sprs.lil_matrix((P.A*f)+lttd)	# need to use array element-wise operations
	Psum = _sp.sum(P.todense(),axis=1)
	if _verbose:
		print ('Psum=',Psum)
	assert _sp.allclose(Psum,_sp.ones(Psum.shape[0])),"error - the modified 'P' matrix, for link travel time, is not a stochastic matrix"
	if not _sp.allclose(P_org.todense(),P.todense()):
		if _verbose:
			print('\tmodified P matrix for links travel time based on legth and speed.\n')
	else:
		print('warning - link travel time did not modify P! check link lengths, speeds and results.')
	if _verbose:
		global _results
		_results['P_updated'] = [ str([", ".join('{:.3f}'.format(i) for i in P.todense().tolist()[j])]).replace("'","").replace("[","").replace("]","") for j in range(P.shape[0]) ]
		_results['StepTime'] = ltt_min
	return 	P

def SetupMCTA(GeoJSON, Verbose = False):	
	"""Setup MCTA network
	SetupMCTA(GeoJSON, Verbose = False)
		GeoJSON		road network dictionary loaded from GeoJSON file (for links, turning probabilities, lenght, lanes, etc.)
		Verbose		show detailed outputs"""
	global _verbose
	global _links
	global _ReverseLinks
	global _results
	_verbose = Verbose
	_links = {int(k):v for k,v in GeoJSON['mcta_json']['links'].items()} # convert dictionary keys to int since JSON require them to be string
	_ReverseLinks = GeoJSON['mcta_json']['LinkReverseDirection']
	_results['linkIDs'] = _sp.array(list(_links.keys())) + 1 # streets; start from 1
	lengths = _sp.array([_links[id]['length'] for id in _links.keys()])
	lanes = _sp.array([_links[id]['lanes'] for id in _links.keys()])
	FreeFlowSpeeds = _sp.array([_links[id]['speed'] for id in _links.keys()])
	P = _ModifyForLinkTravelTime(_SetupTransitionMatrix_P(GeoJSON))
	return lengths, lanes, FreeFlowSpeeds,P


### Solve MCTA ########################################################################

def _FindEigen(P, Calc_C = True, SkipChecksForSpeed = False):
	"""find the left-hand Perron eigenvectors (top 2)"""
	from scipy.linalg import eig
	if SkipChecksForSpeed:
		eigenvalues, eigenvectors = eig(P.T.todense(),check_finite=False)
	else:
		eigenvalues, eigenvectors = eig(P.T.todense())

	# sort eigenvalues/vectors
	idx = eigenvalues.argsort()[::-1] # Perron Root is @ index 0
	eigenvalues = eigenvalues[idx]
	eigenvectors = eigenvectors[:,idx]
	if _verbose and not SkipChecksForSpeed:
		print ('\t',eigenvalues.shape[0], ' eigenvalues: ', eigenvalues)
		print ('\teigenvector: \n', eigenvectors.T)

	#check that Perron Root is real
	if not SkipChecksForSpeed:
		assert eigenvalues[0].imag == 0j, "error - Perron Root value is complex ({:}+{:}j)".format(eigenvalues[0].real,eigenvalues[0].imag)
		assert round(eigenvalues[0].real,6) == 1, "error - cannot find Perron Root (eigenvalue {:} != 1)".format(eigenvalues[0].real)
		for j in range(eigenvectors.shape[1]):
			assert eigenvectors[j,0].imag == 0j, "error - Perron Root eigen vector include complex value @ {:}: {:}+{:}j".format(j,eigenvectors[0,j].real,eigenvectors[0,j].imag)

	PI = _sp.matrix(eigenvectors[:, 0].real).A1
	SUM = PI.sum()
	if SUM != 1: # scaling an eigenvector by scalar value "a" does not change the eigenvector: P (a*v) = Lamda (a*v), where a*v is the eigen vector
		if _verbose and not SkipChecksForSpeed:
			print ('\tNormalizing eigen vector for Perron Root by a scaling factor of: {:}'.format(SUM))
		PI /= PI.sum()

	if _verbose and not SkipChecksForSpeed:
		global _results 
		_results['eigenvalues'] = eigenvalues
		_results['eigenvectors'] = eigenvectors

	if Calc_C:
		C  = _sp.matrix(eigenvectors[:, 1].real).A1
		return PI, C, eigenvalues
	else:
		return PI, None, eigenvalues


def _GeneralInverse(P):
	"""calculate general inverse (Darzin)"""
	from numpy.linalg import matrix_rank
	Q = _sp.matrix(_sp.identity(P.shape[0])-P)
	Qk = Q
	Qk1 = _sp.dot(Qk,Qk)
	# Determine index of A - when rank(A^k) = rank(A^(k+1))
	# The Drazin inverse of a matrix of index 0 or 1 is called the group inverse
	assert matrix_rank(Qk) == matrix_rank(Qk1)
	X = _sp.linalg.solve(Qk1,Qk)  # A^(k+1) * X = A^k
	Qi = Qk * _sp.dot(X,X)  # Ad =  A^(k) * X^(k+1)
	assert _sp.allclose(_sp.dot(Q,Qi),_sp.dot(Qi,Q)) ,'general inverse 1st property invalid (Q)(Q#) = (Q#)(Q)'
	assert _sp.allclose(_sp.dot(_sp.dot(Q,Qi),Q),Q)  ,'general inverse 2nd property invalid (Q)(Q#)(Q) = (Q)'
	assert _sp.allclose(_sp.dot(_sp.dot(Qi,Q),Qi),Qi),'general inverse 3rd property invalid (Q#)(Q)(Q#) = (Q#)'
	return Qi


def _MFPT(P, Pi):
	"""find Mean first passage times"""
	Qi = _GeneralInverse(P)
	m = _sp.matrix(_sp.zeros(P.shape))
	for i in range(P.shape[0]):
		for j in range(P.shape[1]):
			if i == j:
				m[i,j] = 0
			elif Pi[j] == 0: # added to account for when closed states are not deleted from the matrix (when mcta_edit.Delete_States = False)
				m[i,j] = 0
			else:
				m[i,j] = (Qi[j,j] - Qi[i,j]) / Pi[j]
	if _verbose:
		global _results
		print(f'MFPT = {m}')
		_results['MFPT'] = m
	return m
	
def SolveMCTA(lengths, lanes, P, VehiclesCount=1, Objectives=['D','K','C','PI','E'], FreeFlowSpeeds=[], SkipChecksForSpeed=False):
	"""Solve markov chain traffic assignment problem
	SolveMCTA(lengths, lanes, P, VehiclesCount=1, Objectives=['D','K','C','PI','E'], SkipChecksForSpeed=False)
		lengths			links' length in km
		lanes			links' number of lanes
		P				State transition propability matrix (turning probabilities)
		VehiclesCount	number of vehicles on the network at any time
		Objectives		list objectives to produce; 'D' = Density, 
													'K' = Kemeny constant, 
													'C' = clusters, 
													'PI' = probability distribution
													'E' = CO Emissions
		SkipChecksForSpeed	true to speed up the code by skipping verification and validation steps; it will also turn verbosity off"""
	global _results
	if _verbose and not SkipChecksForSpeed:
		print ('\nMCTA results:')
		import time
		tmr = time.time()

	pi, c, eigenvalues = _FindEigen(P, 'C' in Objectives, SkipChecksForSpeed)  # solve MC
	if 'C' in Objectives:
		_results['Clusters'] = c
	if 'PI' in Objectives:
		_results['StationaryProb'] = pi

	if 'K' in Objectives:	
		# Kemeny constant; where K = Ki, i = 1..n (using eigenvalues)
		# if 1 in [round(x,6) for x in eigenvalues[1:]]:
		if any([x>0.999999 for x in eigenvalues[1:]]):  # removed round call for speed
			K = float("inf")
			ResultMsg = 'reducible network with multiple communicating classes (K=inf)'
		else:
			K = (1/(1 - eigenvalues[1:])).sum()
			ClosedStatesCount = P.shape[0]-_sp.count_nonzero(P.sum(axis=0))
			if ClosedStatesCount > 0: # are there any closed but not deleted states? (when mcta_edit.Delete_States = False)
				K -= ClosedStatesCount
				if _verbose and not SkipChecksForSpeed:
					assert P.sum(axis=0) == P.sum(axis=1), "Error, number of closed states in P columns and rows do not match!"
			ResultMsg = 'success'
		if not SkipChecksForSpeed:
			mfpt = _MFPT(P,pi)  # calculate mean first passage time matrix
			K1 = 0
			# Kemeny constant; where K = Ki, i = 1..n (using mean first passage times); with ClosedStatesCount added to Kemeny to adjust for closed states
			for j in range(P.shape[0]):
				K1 += mfpt[0,j] * pi[j]
			assert round(K-K1,6)==0,"Kemeny constant calculation methods do not agree! ( {:} != {:} )".format(K, K1)
		_results['KemenyConst'] = K

	if 'D' in Objectives:
		D = VehiclesCount * pi/(lengths*lanes)  # calculate vehicle density - vehicles/(km.lane)
		# if D.min()<0:
		# 	print(f'found non-ergodic network with negative Densities\nD_min = {D.min()}\K={K}')
		D = _sp.clip(D, 0, None)	# clip Density values to the interval (0,inf), as non-ergodic networks will produce 
									# negative eigenvalues, and so densities. The negative values (PI and D) will also 
									# be acombinied by extreamly large values in other states, and these solutions will 
									# be eliminated by the optimizer. in these cases K = inf
		_results['Density'] = D

	if 'E' in Objectives:
		assert len(FreeFlowSpeeds) == len(lanes), f'Missing FreeFlowSpeeds for some or all lanes ({FreeFlowSpeeds})'
		Speeds, _ = _em.Speed_via_Density(D,'veh/(km.lane)',FreeFlowSpeeds,'km/h') # estimate vehicles speed on each link
		_results['Speeds']=Speeds
		EC, _ = _em.ExternalCost(Speeds, 'km/h', Method='TRANSYT7f') # Calculate average vehicle emissions cost $/(km.veh)
		_results['EmissionCost'] = EC * VehiclesCount * pi  # Calculate average link emissions $/km (in less operations)
		_results['TotalNetworkEmissionCost'] = sum(EC * lengths * D * Speeds)  # Calculate total network emissions $/hr (Emissions * flow, while flow = Density * Speed)

	if _verbose and not SkipChecksForSpeed:
		print ('\tMCTA Message: ' + ResultMsg)
		if 'K' in Objectives:	
			print ('\tKemeny constant: {:.4f} '.format(K) + ' steps ({:.4f} '.format(K*_results['StepTime']*60) + ' min)')
			print ('\texpecting time to mixing: {:.4f} '.format(K+1) + ' steps ({:.4f} '.format((K+1)*_results['StepTime']*60) + ' min)')
		if 'D' in Objectives:
			print ('\tMax density: {:.4f} '.format(max(D)) + ' vehicles/(km.lane)')
		if 'E' in Objectives:
			print(f'\tMin/Avg/Max network speeds (based on BPR): {min(Speeds):.2f}/{Speeds.mean():.2f}/{max(Speeds):.2f} km/h')
			print(f'\tMin/Avg/Max vehicle emission cost on a road (based on its speed): {min(EC):.5f}/{EC.mean():.5f}/{max(EC):.5f} g/(km.veh)')
			print(f"\tMin/Avg/Max link emissions cost: {min(_results['EmissionCost']):.2f}/{_results['EmissionCost'].mean():.2f}/{max(_results['EmissionCost']):.2f} $/km")
			print(f"\tTotal network emission cost: {_results['TotalNetworkEmissionCost']:.2f} $/hr")
		print ('\tcompleted MCTA in {:.3f} seconds.'.format(time.time()-tmr))

	_results['Message'] = ResultMsg
	return _results

### Main ########################################################################

if __name__ == "__main__":
	print(f'{_os.path.basename(__file__)} is designed for use as module, please see the doc string for usage.')