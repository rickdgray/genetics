import random
from tabulate import tabulate

#variables
num_chromos = 20
pco = [0.7] #prob of crossover
target_genes = 682 #decimal of 1010101010
mask = 992 #decimal of 1111100000
inverse_mask = ~mask #decimal of 0000011111

random.seed()

#class contains genes and the fitness value
class chromosome:
	genes = 0
	fit = 0
	fit_cumul = 0

	def __init__(self, genes):
		self.genes = genes

	#less than function for making comparisons
	def __lt__(self, other):
		return self.fit < other.fit

	#bitwise operations to calculate fitness
	def calc_fit(self):
		self.fit = 1 + (10 - (bin(target_genes ^ self.genes).count("1")))

#function loops through all chromosomes and updates fitness. returns boolean for if target gene is found
def set_fit(chromosomes):
	not_found = True
	for i in chromosomes:
		i.calc_fit()
		if (i.fit == 11):
			not_found = False
	return not_found

#function loops through all chromosomes and updates cumulative fitness. this is used for probability calculation
def set_cumul_fit(chromosomes):
	current_fit_cumul = 0
	for i in chromosomes:
		current_fit_cumul += i.fit
		i.fit_cumul = current_fit_cumul
	return current_fit_cumul

#bitwise operation to crossbreed genes
def cross_breed(genes1, genes2):
	temp_genes = genes1
	genes1 = ((genes1 & mask) | (genes2 & inverse_mask))
	genes2 = ((genes2 & mask) | (temp_genes & inverse_mask))
	return (genes1, genes2)

all_avg_gens = []
init_pop_chromosomes = []

#generate random population
for i in range(num_chromos):
	genes = random.randint(0, 1023) #generate random bits from 0 to 1111111111
	init_pop_chromosomes.append(chromosome(genes))

#run entire algorithm for each pco value
for h in pco:

	avg_gens = 0
	all_num_gens = []

	#run evolution algorithm 20 times to get average per pco value
	for j in range(20):

		chromosomes = init_pop_chromosomes
		not_found = set_fit(chromosomes)
		curr_fit_cumul = set_cumul_fit(chromosomes)

		gens = 1

		#evolve chromosomes to find target gene
		while not_found:
			evolved = []

			#replication
			for i in range(int(((1 - h) * num_chromos) + 0.5)):
				rand = random.randint(1, curr_fit_cumul)
				chosen = 0
				for k in chromosomes:
					if (rand > k.fit_cumul):
						chosen += 1
				evolved.append(chromosomes[chosen])

			#crossover
			for i in range(int(((h * num_chromos) / 2) + 0.5)):
				same_chromosome = True
				while same_chromosome:
					rand1 = random.randint(1, curr_fit_cumul)
					rand2 = random.randint(1, curr_fit_cumul)
					chosen1 = 0
					chosen2 = 0
					for k in chromosomes:
						if (rand1 > k.fit_cumul):
							chosen1 += 1
						if (rand2 > k.fit_cumul):
							chosen2 += 1
					if (chosen1 != chosen2):
						same_chromosome = False
				evolved_genes1, evolved_genes2 = cross_breed(chromosomes[chosen1].genes, chromosomes[chosen2].genes)
				evolved.append(chromosome(evolved_genes1))
				evolved.append(chromosome(evolved_genes2))

			#bit mutation
			rand1 = random.randint(0, 19)
			rand2 = (1 << random.randint(0, 9))
			if ((evolved[rand1].genes & rand2) == rand2):
				evolved[rand1].genes -= rand2
			else:
				evolved[rand1].genes += rand2

			chromosomes = evolved
			not_found = set_fit(chromosomes)
			curr_fit_cumul = set_cumul_fit(chromosomes)
			gens += 1

			#calculate average fitness
			avg_fit = 0
			for i in chromosomes:
				avg_fit += i.fit
			avg_fit /= 20

			#print genes for each generation in a table
			print("\n\nGeneration ", gens, "\n")
			table = []
			for current in chromosomes:
				table.append([current.genes, (current.fit - 1)])
			print(tabulate(table, headers=['Genes', 'Fitness']) + '\n')
			print(avg_fit)

		all_num_gens.append(gens)
		avg_gens += gens

	#print table with number of generations for each trial and average generations across all trials
	print("Number of Generations for ", h, ":")
	table = []
	numslol = 1
	for i in all_num_gens:
		table.append([numslol, i])
		numslol += 1
	print(tabulate(table, headers=['Trial', 'Number of Generations']) + '\n')
	all_avg_gens.append(avg_gens / 20)
	print("\nAverage generations for pco ", h, ": ", (avg_gens / 20))

	input("Press Enter")

print("\n\n")

#finally, print initial population
print("Initial Population")
table = []
for current in init_pop_chromosomes:
	table.append([current.genes, (current.fit - 1)])
print(tabulate(table, headers=['Genes', 'Fitness']) + '\n')