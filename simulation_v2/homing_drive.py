#homing drive
'''
'''
import random
import numpy as np 
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
from matplotlib.pyplot import MultipleLocator

#The fitness cost of drive element and loss-of-function allele of target gene can be set by global variable.
#The gender of every individuals are initially set to hermaphroditic, i.e. bisexual, for most of plant species have
# both male and female reproductive organs. The gender will be changed to male or female before they mate.
DRIVE_FITNESS=1 #0.95
R2_FITNESS=1 #0.99
CAPACITY=10000

class Individual_homing:
    def __init__(self,sex,cutrate,chr1,fitCoef):
        self.sex=sex
        self.cutrate=cutrate
        self.chr1=chr1
        self.fitCoef=fitCoef
    
    def homing(self):
        if self.chr1 in [("dr","wt"),("wt","dr")]:
            if random.random() < self.cutrate:
                self.chr1=("dr","dr")

    def getFitness(self):
        c1_value=self.fitCoef if "dr"==self.chr1[0] else 1
        c2_value=self.fitCoef if "dr"==self.chr1[1] else 1

        return c1_value+c2_value
    
class Population_homing:
    def __init__(self,capacity,size,cutrate,chr1,fitCoef):
        self.size = size #population size, an integer
        self.individuals = [] #A list storing the individual objects.
        self.cutrate=cutrate
        self.capacity=capacity
        self.fitCoef=fitCoef

        for i in range(size):
            sex="hermaphroditic"
            self.individuals.append(Individual_homing(sex, cutrate, chr1, fitCoef)) #Initialize the list with size of individual objects.


    def add_pop(self,subpop): #Copy the individual objects from other Population object.
        for ind in subpop.individuals:
            self.individuals.append(ind)
        self.size=len(self.individuals)
    
    def get_size(self):
        self.size=len(self.individuals)
        return self.size
     
    def next_generation(self):
        offspring = []
        weights=[]
        for i in self.individuals:
            weights.append(i.getFitness())
        self.size=self.get_size()

        for i in range(int(self.size/2)):
            parent1 = random.choice(self.individuals) #If the drive and KO has fitness cost, use random.choices(self.individuals,weights=weights)[0] to set different probability to each individuals.
            parent1.sex="male"     #Set the individual's gender when crossing.
            parent2 = random.choice(self.individuals) #random.choices(self.individuals,weights=weights)[0]
            parent2.sex="female"   #Set the individual's gender when crossing.
            offspring.extend(self.produce_offspring(parent1, parent2))

            parent1.sex="hermaphroditic" #Change the gender to the bisexual state.
            parent2.sex="hermaphroditic" #Change the gender to the bisexual state.

        self.individuals = offspring   #Cover the parent generation using offspring generation
        self.size=len(self.individuals)
        random.shuffle(self.individuals) #Shuffle the individual list

    def produce_offspring(self, male, female):
        offspring = []
        capacity_fitness_scalling=10/((9*self.size/self.capacity)+1)
        p=capacity_fitness_scalling*0.04
        num_offspring=np.random.binomial(50, p, 1)[0]
        for i in range(num_offspring):
            male.homing()
            genome1=random.choice(male.chr1)
            female.homing()
            genome2=random.choice(female.chr1)

            chr1=(genome1,genome2)
            sex="hermaphroditic"
            offspring.append(Individual_homing(sex, self.cutrate,chr1, self.fitCoef))
        return offspring


    def count_gene_drive_alleles(self,allele):
        count = 0
        if allele=="True": #If the allele is True, count the drive allele number, else, count the drive carriers number.
            for individual in self.individuals:
                for i in range(2):
                    if 'dr' == individual.chr1[i]:
                        count += 1
        else:
            for individual in self.individuals:
                if 'dr' in individual.chr1:
                    count += 1
        return count
       