import time, numpy, random, multiprocessing
from deap import creator, base, tools

AllGenes = 			[0, 1, 2]
chromosome_size = 	180
population_size = 	200
KeepBetweenGens = 	20
num_generations = 	1000
CXPB, MUTPB, IndpProb = 0.50, 0.30, 0.01  # crossover prob, mutation prob, individual gene mutations prob
PlotResults = 		False	 			  # plot results?
nObjectives = 		2		 			  # multi-objective optimization?
MCTA_BaseFile =		'AD'

ga_data = {}
ga_data['GA_Parameters'] = {'MCTA_BaseFile':MCTA_BaseFile,'AllGenes':AllGenes,'chromosome_size':chromosome_size,'population_size':population_size,'KeepBetweenGens':KeepBetweenGens,'num_generations':num_generations,'CXPB':CXPB,'MUTPB':MUTPB,'nObjectives':nObjectives}

# MCTA initialize and load roadnetwork data
import mcta
mcta.SetEnviroment(env_base=MCTA_BaseFile,env_verbose=False)

LOG = ""


def print2(txt):
	global LOG
	print(txt)
	LOG += '\n' + txt


def evalObj(individual):
	tmp = mcta.ModifyNetworkAndSolveMC(individual)
	# with numpy.errstate(all='ignore'):
	return (tmp[0],tmp[1].real)		# return the first two fitness values


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
	for i in range(nObjectives):
		if OldBest[i] != 0:	delta.append("{:03.3f}%".format((OldBest[i]-CurBest[i])/OldBest[i] * 100))
		else:				delta.append("---.---%")
	fit = ["{:03.3f}".format(CurBest[i]) for i in range(nObjectives)]
	return "Best: ({:}), Delta: ({:})".format(", ".join(fit),", ".join(delta))


if nObjectives == 2:	creator.create("FitnessMin", base.Fitness, weights=(-1,-1))  # (-65,-1)		(-1e-10,-1)
else:					creator.create("FitnessMin", base.Fitness, weights=(-1,-1e-10))
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
if nObjectives == 2:
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
	for i in range(chromosome_size):				# set individual 0 to current road network state (all genes 0)
		pop[0][i] = 0

	fitnesses = list(toolbox.map(toolbox.evaluate, pop))		# Evaluate the entire population
	for ind, fit in zip(pop, fitnesses):
		ind.fitness.values = fit

	CurrentSolution = pop[0].fitness.values
	print2("\nCurrent Solution: %s with fitness: %s" % (str(pop[0]),str(CurrentSolution)))
	print2("\nGA Parameters:\n\tchromosome_size: %s\n\tpopulation_size: %s\n\tKeepBetweenGens: %s\n\tnum_generations: %s\n\t           CXPB: %s\n\t          MUTPB: %s\n\t       IndpProb: %s" % (chromosome_size,population_size,KeepBetweenGens,num_generations,CXPB, MUTPB, IndpProb))
	ga_data['CurrentSolution'] = CurrentSolution

	# stats setup
	stats_d = tools.Statistics(key=lambda ind: ind.fitness.values[0])
	stats_d.register("min", numpy.min)
	stats_d.register("max", numpy.max)
	stats_d.register("std", numpy.std)
	stats_d.register("avg", numpy.mean)
	stats_kd = tools.Statistics()
	if nObjectives == 2:
		stats_k = tools.Statistics(key=lambda ind: ind.fitness.values[1])
		stats_k.register("min", numpy.min)
		stats_k.register("max", numpy.max)
		stats_k.register("std", numpy.std)
		stats_k.register("avg", numpy.mean)
		mstats = tools.MultiStatistics(D=stats_d, K_Dmin=stats_kd, K=stats_k)
	else:
		mstats = tools.MultiStatistics(D=stats_d, K_Dmin=stats_kd)

	print2("\nStart of evolution\n" + "*"*18)

	best_ind = tools.selBest(pop, 1)[0]
	Cur_best = best_ind.fitness.values
	Old_best = best_ind.fitness.values
	print2("  Generation {:4d} - Evaluated {:4d} individuals. {:}. time:{:4.4f} min.".format(0,len(pop),ObjStats(Cur_best,Old_best),(time.time()-tic)/60))
	with numpy.errstate(all='ignore'):
		RECORD = mstats.compile(pop)
	RECORD.update({'K_Dmin':{'K_Dmin':best_ind.fitness.values[1].real}})
	logbook.record(gen=0, evals=len(pop), record=RECORD)

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
		print2("  Generation {:4d} - Evaluated {:4d} individuals. {:}. time:{:4.4f} min.".format(g+1,len(invalid_ind),ObjStats(Cur_best,Old_best),(time.time()-tic)/60))

		with numpy.errstate(all='ignore'):
			RECORD = mstats.compile(pop)
		RECORD.update({'K_Dmin':{'K_Dmin':best_ind.fitness.values[1].real}})
		logbook.record(gen=g+1, evals=len(invalid_ind), record=RECORD)

		HoF.update(pop)

	pool.close()

	toc = time.time()

	print2("-- End of (successful) evolution --")
	print2("\nrun statistics:\n" + "*"*18)
	import io
	from contextlib import redirect_stdout
	output = ""
	with io.StringIO() as buf, redirect_stdout(buf):
		print(logbook)
		output = buf.getvalue()
	print2(output)
	ga_data['logbook'] = logbook

	best_ind = tools.selBest(pop, 1)[0]
	ga_data['best_ind'] = {'best_ind_solution':best_ind,'fitness':best_ind.fitness.values}
	print2("Best individual is %s, %s\n" % (best_ind, ObjStats(Cur_best,CurrentSolution)))
	print2(mcta.ExplainModification(best_ind))
	print2("\nRun time = {:4.4f} min".format((toc-tic)/60))

	# data for fine-tuning
	print2("\nSaving HoF_GeneStat to 'ga_HoF_GeneStat' as pickle...")
	import pickle
	HoF_GeneStat = {}
	for p in range(chromosome_size):
		for g in AllGenes:
			HoF_GeneStat[(p,g)] = [0,0]
	for ind in HoF:
		for p in range(chromosome_size):
			HoF_GeneStat[(p,ind[p])][0] += 1
			HoF_GeneStat[(p,ind[p])][1] += CurrentSolution[0]/ind.fitness.values[0]
	# print(str(HoF_GeneStat).replace(', ((','\n(('))
	pickle.dump(HoF_GeneStat, open("ga_HoF_GeneStat.pkl", "wb"))

	if nObjectives == 2:
		HoF_D, HoF_K = zip(*[x.fitness.values for x in HoF])
		ga_data['ParetoFront'] = {'Dmax':HoF_D, 'K':HoF_K}
	print2("\nSaving data to 'ga_data' as json and pickle...")
	pickle.dump(ga_data, open("ga_data.pkl", "wb"))
	import json
	with open("ga_data.json", 'w') as outfile:
		json.dump(ga_data, outfile, indent=4, sort_keys=True)

	# write log to text file
	print2("Saving screen log to 'ga.log'...")
	with open('ga.log','w') as outputfile:
		outputfile.write(LOG)
		outputfile.close()

	# plotting
	if PlotResults:
		print("Generating plots...")
		import matplotlib.pyplot as plt
		gen = logbook.select("gen")

		D_min = logbook.chapters['record'].chapters["D"].select("min")
		fig, ax1 = plt.subplots()
		line1 = ax1.plot(gen[:], D_min[:], "b-", label="min Dmax")
		ax1.set_xlabel("Generation")
		ax1.set_ylabel("Vehicle Density")
		lns = line1

		K_Dmin = logbook.chapters['record'].chapters["K_Dmin"].select("K_Dmin")
		ax2 = ax1.twinx()
		line2 = ax2.plot(gen[:], K_Dmin[:], "g-", label="Kemeny(min Dmax)")
		ax2.set_ylabel("Kemeny")
		lns += line2

		if nObjectives == 2:
			K_min = logbook.chapters['record'].chapters["K"].select("min")
			line3 = ax2.plot(gen[:], K_min[:], "r-", label="min Kemeny")
			lns += line3

		labs = [l.get_label() for l in lns]
		ax1.legend(lns, labs, loc="upper right")
		plt.title("GA Run Results")

		if nObjectives == 2:
			fig, ax = plt.subplots()
			ax.plot(HoF_D, HoF_K, "o", label="Pareto Front")
			ax.set_xlabel("Dmax")
			ax.set_ylabel("Kemeny")
			plt.title("Pareto Front")

		plt.show()


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
