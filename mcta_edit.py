#!/usr/bin/env python

"""
Markov Chains Traffic Assignment (MCTA) network modification module

Usage:
	msg = ExplainModification(LinksEdits)
	P = CloseCrisisRoads(P, CR)
	Results = ModifyNetworkAndSolveMC(P, LinksEdits, VehiclesCount=1, SkipChecksForSpeed=False)


_author = 	"Sinan Salman (sinan.salman@zu.ac.ae)"
_version = 	"Revision: 0.18"
_date = 	"Date: 2019/07/17"
_copyright= "Copyright (c)2017-2019 Sinan Salman"
_license =	"GPLv3"
"""

import os as _os
import scipy as _sp
import scipy.sparse as _sprs
import mcta as _mcta


def ExplainModification(LinksEdits):
	"""explain network modification using plain english"""
	revlinks = _mcta._ReverseLinks
	msg="the following links direction are reversed: "
	for i in range(len(LinksEdits)):
		if LinksEdits[i] == 1:
			msg += str(revlinks[i][0]) + ", "
		elif LinksEdits[i] == 2:
			msg += str(revlinks[i][1]) + ", "
	return msg


# Does not yet work with Delete_States=False.
def CloseCrisisRoads(P, CR): 
	"""modify network by closing road IDs in CR"""
	verbose = _mcta._verbose
	links = _mcta._links
	nl = len(links)
	ONES = _sp.matrix(_sp.ones((nl,1)))

	X = _sp.matrix(_sp.eye(nl))
	for i in CR: X[i,i] = 0	 # 0 = no flow AND 1 = full flow

	# close crisis Roads
	if verbose: 
		print('closing crisis roads:' + str(CR) + '\n'+'*'*30)
	if verbose: 
		print('P\n',P.todense())
	if verbose: 
		print('X\n',X)
	P_d = _sp.diag(_sp.diag(P.todense()))
	if verbose: 
		print('P_d\n',P_d)
	C1 = X*(P-P_d)*X
	if verbose: 
		print('X*(P-P_d)*X\n',C1)
	P_d = X*P_d*X
	if verbose: 
		print('new P_D = X*P_d*X\n',P_d)
	C2 = C1*ONES
	if verbose: 
		print('Row sum of X*(P-P_d)*X\n',C2)
	for i in range(len(C2)):
		if C2[i] == 0:
			C2[i] = 1
	if verbose: 
		print('Row sum of X*(P-P_d)*X, w/ zeros replaced with ones\n',C2)
	C3 = _sp.matrix(_sp.diag(C2.A1)).I
	if verbose: 
		print('Factor\n',C3)
	C3 = C3 * (_sp.eye(nl)-P_d)
	if verbose: 
		print('Factor, w/ diagonal values taken into consideration\n',C3)
	C4 = C3*C1
	if verbose: 
		print('multiply by factor matrix\n',C4)
	P = _sprs.lil_matrix(C4 + P_d)
	if verbose: 
		print('New P\n',P)
	if verbose: 
		print('check P row sum\n',P*ONES)
	return P

def ModifyNetworkAndSolveMC(P, LinksEdits, VehiclesCount=1, Objectives=['D','K','E'], SkipChecksForSpeed=False):
	"""modify network using LinksEdits and solve the resulting MC
	for each link 0=no change, 1=close the link, 2=close the reverse side"""
	verbose = _mcta._verbose
	links = _mcta._links
	# revlinks but be formated as pairs of opposing direction links; 
	# total number of pairs will be half the notoal number of links 
	# (every link must have a reverse side)
	revlinks = _mcta._ReverseLinks 

	nl = len(links)
	ONES = _sp.matrix(_sp.ones((nl,1)))

	if not SkipChecksForSpeed:
		assert len(LinksEdits) == nl/2, "Size of LinksEdits must be half of the number of links in the network"
		assert all(x>=0 and x<=2 for x in LinksEdits), "LinksEdits include value(s) out of bound (0,1,2)"

	closed = []
	r = [-1] * nl # link reverse list
	for i in range(len(LinksEdits)):
		l1 = revlinks[i][0]
		l2 = revlinks[i][1]
		r[l1] = l2
		r[l2] = l1
		if LinksEdits[i] == 1:
			closed.append(l1)
		elif LinksEdits[i] == 2:
			closed.append(l2)

	X = _sp.matrix(_sp.eye(nl))
	for i in closed: 
		X[i,i] = 0	 # 0 = no flow AND 1 = full flow

	# Non-Ergodic Network Handling; recursively close all 
	# orphan links (no inbound or outbound probabilities) 
	# (after closure (setting a row and column i to zeros) more links may become orphan
	C1 = X*P*X
	REMOVE = closed.copy()
	KEEP = [i for i in range(nl) if i not in REMOVE]
	Zr = False
	Zc = False
	while True:
		# remove all zero columns
		tmp = ONES.T*C1
		# rem = [i for i in KEEP if round(tmp[0,i],6) == 0]
		rem = [i for i in KEEP if tmp[0,i] < 0.000001]  # removed round call for speed
		REMOVE += rem
		KEEP = [i for i in range(nl) if i not in REMOVE]
		if rem != []: Zc = True
		else: Zc = False
		for i in REMOVE: X[i,i] = 0
		C1 = X*P*X

		# remove all zero rows
		tmp = C1*ONES
		# rem = [i for i in KEEP if round(tmp[i,0],6) == 0]
		rem = [i for i in KEEP if tmp[i,0] < 0.000001]  # removed round call for speed
		REMOVE += rem
		KEEP = [i for i in range(nl) if i not in REMOVE]
		if rem != []: Zr = True
		else: Zr = False
		for i in REMOVE: X[i,i] = 0
		C1 = X*P*X

		if not Zc and not Zr: break

	lengths = _sp.array([links[i]['length'] for i in links.keys()])
	FFS = _sp.array([links[i]['speed'] for i in links.keys()])
	# constraint to allow different number of lanes values for opposing roads
	lanes = _sp.array([ X[i,i] * links[i]['lanes'] + ( 1 - X[r[i],r[i]] ) * links[r[i]]['lanes'] for i in links.keys() ]) 

	# Generate Q permutation matrix
	j=0
	Q=_sp.matrix(_sp.eye(nl))
	for i in range(nl):
		if X[i,i]==1:
			Q[[i,j]]=Q[[j,i]]
			j+=1

	Delete_States = True   # best results and faster computing time using 'True'
	if Delete_States:
		nl     = nl-len(REMOVE)
		ONES   = _sp.delete(ONES,REMOVE,0)
		lanes  = _sp.delete(lanes,REMOVE,0)
		FFS    = _sp.delete(FFS,REMOVE,0)
		lengths= _sp.delete(lengths,REMOVE,0)
		C1     = _sp.delete(C1,REMOVE,1)
		C1     = _sp.delete(C1,REMOVE,0)
		C2     = (C1*ONES).A1
		C      = _sprs.lil_matrix(_sp.matrix(_sp.diag(C2)).I * C1)
		if verbose:
			print('\t** state closure method: Remove States from P matrix')
	else:
		C1 = X*P*X
		B  = Q*C1*Q.T
		C2 = [x if x>0 else 1 for x in (B*ONES).A1]
		C  = _sprs.lil_matrix(_sp.matrix(_sp.diag(C2)).I * B)
		if not SkipChecksForSpeed:
			assert _sp.allclose(Q*ONES,ONES), "not all rows in Q include a value of 1, not a proper permutation (orthogonal matrix), please check!\n{:}".format(Q)
			assert _sp.allclose(Q.T*ONES,ONES), "not all columns in Q include a value of 1; not a proper permutation (orthogonal matrix), please check!\n{:}".format(Q)
			assert _sp.allclose(Q*Q.T,_sp.eye(nl)), "Q*Q' is not an identity matrix; not a proper permutation (orthogonal matrix), please check!\n{:}".format(Q)
		lengths= (_sp.matrix(lengths)*Q.T).A1
		FFS	   = (_sp.matrix(FFS)*Q.T).A1
		lanes  = (_sp.matrix(lanes)*Q.T).A1
		if verbose:
			print('\t** state closure method: Set closed states probabilities to zero in P matrix, and reorder P col/rows')

	Results = _mcta.SolveMCTA(lengths=lengths, lanes=lanes, P=C, VehiclesCount=VehiclesCount, 
                         	  Objectives=Objectives, FreeFlowSpeeds=FFS, 
							  SkipChecksForSpeed=SkipChecksForSpeed)
	
	# Post-solve edit
	NoPostSolveEditsNeeded = ['linkIDs','P_org','P_updated','StepTime','KemenyConst','Message']  # skip items not changed in pre-solve edits
	Removed_Links_Count = len(REMOVE)
	if Removed_Links_Count > 0:
		if Delete_States:
			# AB: Add_Back
			AB1 = _sp.array([0]*Removed_Links_Count)
			AB2 = _sp.zeros((nl,Removed_Links_Count))
			AB3 = _sp.zeros((Removed_Links_Count,nl))
			AB4 = _sp.zeros((Removed_Links_Count,Removed_Links_Count))
		for k in Results.keys():
			if k in NoPostSolveEditsNeeded:
				continue
			if hasattr(Results[k],'shape'):
				if len(Results[k].shape) == 1:
					if Delete_States:
						Results[k] = _sp.append(Results[k],AB1)
					Results[k] = (Results[k]*Q).A1
				elif len(Results[k].shape) == 2:
					if Delete_States:
						Results[k] = _sp.block([ [Results[k],AB2], [AB3,AB4] ])
					Results[k] = Q.T*Results[k]*Q

	if verbose and Results['Message'] != 'success':
		print (Results['Message'])
		print ('Network has the following links\n\tGene: \t\t{:}\n\tClosed: \t{:}\n\tRemoved: \t{:}\n\tRemaining: \t{:}\n'.format(LinksEdits,closed,[x for x in REMOVE if x not in closed],[x for x in range(len(links)) if x not in REMOVE]))

	return Results


### Main ########################################################################

if __name__ == "__main__":
	print(f'{_os.path.basename(__file__)} is designed for use as module, please see the doc string for usage.')