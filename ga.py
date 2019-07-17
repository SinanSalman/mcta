#!/usr/bin/env python

import time, numpy, random, multiprocessing, pickle
from deap import creator, base, tools

import mcta
import mcta_rw
import mcta_edit
# import mcta_vis

AllGenes = 			[0, 1, 2]
chromosome_size = 	184 # 180 for AD, 184 for AD2
population_size = 	192 #192
KeepBetweenGens = 	56 #56
num_generations = 	5 #1000
CXPB, MUTPB, IndpProb = 0.50, 0.18, 0.01  # crossover prob, mutation prob, individual gene mutations prob
nObjectives = 		2		 			  # multi-objective optimization?
MCTA_BaseFile =		'AD2'

nFitnessValues = None
ga_data = {}
ga_data['GA_Parameters'] = {'MCTA_BaseFile':MCTA_BaseFile,'AllGenes':AllGenes,'chromosome_size':chromosome_size,'population_size':population_size,'KeepBetweenGens':KeepBetweenGens,'num_generations':num_generations,'CXPB':CXPB,'MUTPB':MUTPB,'nObjectives':nObjectives}
LOG = ""

# MCTA initialize and load roadnetwork data
(Settings, GeoJSON) = mcta_rw.LoadDataFiles(MCTA_BaseFile)
_, _, _, P = mcta.SetupMCTA(GeoJSON, Verbose=False)
# mcta_vis.Initialize(GeoJSON,Settings,MCTA_BaseFile,ShowFigs=True, SaveFigs2PNG=False)


def print2(txt):
	global LOG
	print(txt)
	LOG += '\n' + txt


def evalObj(individual):
	R = mcta_edit.ModifyNetworkAndSolveMC(P, individual, GeoJSON['mcta_json']['VehiclesCountEst'], SkipChecksForSpeed = True)
	# with numpy.errstate(all='ignore'):
	return (R['Density'].max(), 
			# R['KemenyConst'].real, 
			R['TotalNetworkEmissionCost'])		# return the fitness values


def mutFlipNum(individual, indpb):
	"""Flip the value of the gene to one of the other two possible values in [0,1,2] \n:param individual: Individual to be mutated.\n:param indpb: Independent probability for each attribute to be flipped.\n:returns: A tuple of one individual."""
	for i in range(chromosome_size):
		if random.random() < indpb:
			genes = AllGenes.copy()
			genes.remove(individual[i])
			individual[i] = genes[random.randint(0, 1)]
	return individual,


def ObjStats(CurBest, OldBest):
	delta = []
	for i in range(nFitnessValues):
		if OldBest[i] != 0:	
			delta.append(f'{(OldBest[i]-CurBest[i])/OldBest[i] * 100: 7.2f}%')
		else:				
			delta.append('---.--%')
	fit = [f'{x: 8.3f}' for x in CurBest]
	return "Best: ({:}), Delta: ({:})".format(", ".join(fit),", ".join(delta))


# if nObjectives > 1:
creator.create("FitnessMin", base.Fitness, weights=(-1,) * 4) 
# else:
# 	creator.create("FitnessMin", base.Fitness, weights=(-1,-1e-10))
creator.create("Individual", list, fitness=creator.FitnessMin)

# Structure initializers
toolbox = base.Toolbox()									# initialize toolbox
toolbox.register("gene_generator", random.randint, 0, 2)    # RandomGene
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.gene_generator, chromosome_size)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("evaluate", evalObj)						# register the goal / fitness function
toolbox.register("mate", tools.cxTwoPoint)					# register the crossover operator
toolbox.register("mutate", mutFlipNum, indpb=IndpProb)		# register a mutation operator with a probability to flip each attribute/gene

logbook = tools.Logbook()		# create logbook for storing statistics
if nObjectives > 1:
	toolbox.register("select", tools.selNSGA2)  # operator for selecting individuals for breeding the next generation: non dominant sorting using NSGA2 algorithm
	HoF = tools.ParetoFront()
else:
	toolbox.register("select", tools.selBest)		# operator for selecting best K individuals for breeding the next generation
	HoF = tools.HallOfFame(population_size)


def main():
	global logbook

	pool = multiprocessing.Pool()		# Process Pool of workers for parallel processing
	ga_data['CPUcores'] = pool._processes
	print2("using {:} parallel processes".format(pool._processes))
	toolbox.register("map", pool.map)

	tic = time.time()

	pop = toolbox.population(n=population_size)		# create an initial population
	for j in range(round(population_size*0.05)):
		for i in range(chromosome_size):			# set 5% of individuals to current road network state (all genes 0)
			pop[j][i] = 0

	fitnesses = list(toolbox.map(toolbox.evaluate, pop))		# Evaluate the entire population
	for ind, fit in zip(pop, fitnesses):
		ind.fitness.values = fit

	CurrentSolution = pop[0].fitness.values
	global nFitnessValues
	nFitnessValues = len(CurrentSolution)
	print2(f'Optimizing using {nObjectives} objective(s)\nTracking {nFitnessValues} fitness value(s)\n')

	print2(f"\nCurrent Solution: {str(pop[0])}\nwith fitness: {str(CurrentSolution)}")
	print2("\nGA Parameters:\n\tchromosome_size: %s\n\tpopulation_size: %s\n\tKeepBetweenGens: %s\n\tnum_generations: %s\n\t           CXPB: %s\n\t          MUTPB: %s\n\t       IndpProb: %s" % (chromosome_size,population_size,KeepBetweenGens,num_generations,CXPB, MUTPB, IndpProb))
	ga_data['CurrentSolution'] = CurrentSolution

	# stats setup
	
	stats = tools.Statistics(key=lambda ind: ind.fitness.values)
	stats.register("min", numpy.min, axis=0)
	stats.register("max", numpy.max, axis=0)
	stats.register("std", numpy.std, axis=0)
	stats.register("mean", numpy.mean, axis=0)

	print2("\nStart of evolution\n" + "*"*18)

	best_ind = tools.selBest(pop, 1)[0]
	Cur_best = best_ind.fitness.values
	Old_best = best_ind.fitness.values
	print2("  Generation {: 4d} - Evaluated {: 4d} individuals. {:}. time:{: 8.3f} min.".format(0,len(pop),ObjStats(Cur_best,Old_best),(time.time()-tic)/60))
	# with numpy.errstate(all='ignore'):
	RECORD = stats.compile(pop)
	logbook.record(gen=0, evals=len(pop), best_ind=Cur_best, **RECORD)

	for g in range(num_generations):			# Begin the evolution

		offspring = list(toolbox.map(toolbox.clone, pop))		# Clone the selected individuals

		for child1, child2 in zip(offspring[::2], offspring[1::2]):		# Apply crossover and mutation on the offspring
			if random.random() < CXPB:			# cross two individuals with probability CXPB
				toolbox.mate(child1, child2)
				del child1.fitness.values		# fitness values of the children must be recalculated later
				del child2.fitness.values

		for mutant in offspring:				# mutate an individual with probability MUTPB
			if random.random() < MUTPB:
				toolbox.mutate(mutant)
				del mutant.fitness.values

		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit

		tmp = toolbox.select(offspring, KeepBetweenGens)					# Select the next generation individuals
		pop = list(toolbox.map(toolbox.clone, tmp)) * (len(pop)//KeepBetweenGens)		# Clone the selected individuals
		rem = len(pop) % KeepBetweenGens
		if rem > 0:
			tmp = toolbox.select(offspring, rem)
			pop += list(toolbox.map(toolbox.clone, tmp))

		best_ind = tools.selBest(pop, 1)[0]
		Old_best = Cur_best
		Cur_best = best_ind.fitness.values
		print2("  Generation {: 4d} - Evaluated {: 4d} individuals. {:}. time:{: 8.3f} min.".format(g+1,len(invalid_ind),ObjStats(Cur_best,Old_best),(time.time()-tic)/60))

		RECORD = stats.compile(pop)
		logbook.record(gen=g+1, evals=len(invalid_ind), best_ind=Cur_best, **RECORD)

		HoF.update(pop)

	pool.close()

	toc = time.time()

	ga_data['logbook'] = logbook
	print2("-- End of (successful) evolution --")
	print2("\nrun statistics:\n" + "*"*18)
	print2(logbook.__str__())
	
	best_ind = tools.selBest(pop, 1)[0]
	ga_data['best_ind'] = {'best_ind_solution':best_ind,'fitness':best_ind.fitness.values}
	print2("Best individual is %s\n%s\n" % (best_ind, ObjStats(Cur_best,CurrentSolution)))
	print2(mcta_edit.ExplainModification(best_ind))
	print2("\nRun time = {: 8.3f} min".format((toc-tic)/60))

	if nObjectives > 1:
		ga_data['ParetoFront'] = list(zip(*[x.fitness.values for x in HoF]))
	
	print2("\nSaving data to 'ga_data' pickle...")
	pickle.dump(ga_data, open("ga_data.pkl", "wb"))

	# write log to text file
	print2("Saving screen log to 'ga.log'...")
	with open('ga.log','w') as outputfile:
		outputfile.write(LOG)


def ga_batchrun(MaxTime=600, nGn=999999, Psz=200, Ktp=0.1, Cpr=0.5, Mpr=0.30, gMpr=0.01):

	# # debug
	# import inspect
	# frame = inspect.currentframe()
	# args, _, _, values = inspect.getargvalues(frame)
	# print('function name "%s"' % inspect.getframeinfo(frame)[2])
	# for i in args:
	# 	print("    %s = %s" % (i, values[i]))

	if nObjectives > 1:
		raise ValueError('must use single objective for batch runs. check nObjectives.')

	# batch runs are always single objective runs. for use in finetuning
	global num_generations
	global population_size
	global KeepBetweenGens
	global CXPB
	global MUTPB
	global IndpProb

	num_generations = nGn
	population_size = Psz
	KeepBetweenGens = round(Ktp*Psz)
	CXPB = Cpr
	MUTPB = Mpr
	IndpProb = gMpr

	pool = multiprocessing.Pool()		# Process Pool of workers for parallel processing
	# print("using {:} parallel processes".format(pool._processes))
	toolbox.register("map", pool.map)

	tic = time.time()

	pop = toolbox.population(n=population_size)		# create an initial population
	for j in range(round(population_size*0.05)):
		for i in range(chromosome_size):			# set 5% of individuals to current road network state (all genes 0)
			pop[j][i] = 0

	obj_evals = len(pop)
	fitnesses = list(toolbox.map(toolbox.evaluate, pop))		# Evaluate the entire population
	for ind, fit in zip(pop, fitnesses):
		ind.fitness.values = fit

	# CurrentSolution = pop[0].fitness.values

	# print("\nStart of evolution\n" + "*"*18)

	# best_ind = tools.selBest(pop, 1)[0]
	# Cur_best = best_ind.fitness.values
	# Old_best = best_ind.fitness.values
	# print("  Generation {:4d} - Evaluated {:4d} individuals. {:}. time:{:4.4f} min.".format(0,len(pop),ObjStats(Cur_best,Old_best),(time.time()-tic)/60))

	for g in range(num_generations):			# Begin the evolution

		offspring = list(toolbox.map(toolbox.clone, pop))		# Clone the selected individuals

		for child1, child2 in zip(offspring[::2], offspring[1::2]):		# Apply crossover and mutation on the offspring
			if random.random() < CXPB:			# cross two individuals with probability CXPB
				toolbox.mate(child1, child2)
				del child1.fitness.values		# fitness values of the children must be recalculated later
				del child2.fitness.values

		for mutant in offspring:				# mutate an individual with probability MUTPB
			if random.random() < MUTPB:
				toolbox.mutate(mutant)
				del mutant.fitness.values

		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		obj_evals += len(invalid_ind)
		fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit

		tmp = toolbox.select(offspring, KeepBetweenGens)					# Select the next generation individuals
		pop = list(toolbox.map(toolbox.clone, tmp)) * (len(pop)//KeepBetweenGens)		# Clone the selected individuals
		rem = len(pop) % KeepBetweenGens
		if rem > 0:
			tmp = toolbox.select(offspring, rem)
			pop += list(toolbox.map(toolbox.clone, tmp))

		# best_ind = tools.selBest(pop, 1)[0]
		# Old_best = Cur_best
		# Cur_best = best_ind.fitness.values
		# print("  Generation {:4d} - Evaluated {:4d} individuals. {:}. time:{:4.4f} min.".format(g+1,len(invalid_ind),ObjStats(Cur_best,Old_best),(time.time()-tic)/60))

		if time.time()-tic > MaxTime:
			break

	pool.close()

	toc = time.time()

	# print("-- End of (successful) evolution --")
	# print("\nrun statistics:\n" + "*"*18)

	best_ind = tools.selBest(pop, 1)[0]
	# print("Best individual is %s, %s\n" % (best_ind, ObjStats(Cur_best,CurrentSolution)))
	# print("\nRun time = {:4.4f} min".format((toc-tic)/60))

	return {'fitness':best_ind.fitness.values, 'runtime': toc-tic, 'obj_evals': obj_evals}


# Main ########################################################################
if __name__ == "__main__":
	print2("\nGA runing as an python script...")
	main()

# Keep for future use #########################################################

# probability distribution for genegeneration
# RandGeneDist = { 0:0.98, 1:0.01, 2:0.01 } # slow conversion to best (near optimal) solution.
# RandGeneDist = { 0:0.90, 1:0.05, 2:0.05 }
# AllGenes = list(RandGeneDist.keys())
# GenesNum = len(AllGenes)
# RandGeneCumDist = {}
# CumDist=0
# for i in range(GenesNum):
#     CumDist += RandGeneDist[i]
#     RandGeneCumDist[i] = CumDist
# del CumDist
# assert RandGeneCumDist[AllGenes[-1]]==1,"Random Genes Distribution (RandGeneDist) does not add to one."
#
# def RandomGene():
# 	p = random.random()
# 	for i in range(GenesNum):
# 		if p < RandGeneCumDist[i]:
# 			return i
