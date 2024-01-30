'''
Author: Bingke Jiao, supervised by Professor Wenfeng Qian.
Ph.D. Yang Liu and Weiguang Wang gave valuable advice in the design of this script.

This script implements the dynamic simulation of CRISPR-mediated toxin-antidote gene drive.
The main features of this scipt include:
1. It is a forward genetic simulation. Set parameters to the initial population and the population will propagate. You can monitor the frequency of drive carriers with the increasing of generations.
2. It is an individual based, stocastic model, which using Wright-Fisher model.
3  It uses two classes: the individual class and population class.

The design of this script referenced to the SLiM software written by Benjamin C. Haller and Philipp W. Messer. and the Toxin-antidote scripts written by Samuel E. Champer and Jackson Champer
'''

import random
import copy
import numpy as np 
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
from matplotlib.pyplot import MultipleLocator

#The fitness cost of drive element and loss-of-function allele of target gene can be set by global variable.
#The gender of every individuals are initially set to hermaphroditic, i.e. bisexual, for most of plant species 
#have both male and female reproductive organs. The gender will be changed to male or female before they mate.
# DRIVE_FITNESS=1 #0.95
# KO_FITNESS=1 #0.99
# CAPACITY=10000

class Individual:
    def __init__(self, sex, embryo_cutrate,germline_cutrate,chr1,chr2,incomPene,varyFit):
        self.sex = sex
        self.chr1=chr1
        self.chr2=chr2
        self.embryo_cutrate=embryo_cutrate
        self.germline_cutrate=germline_cutrate

        self.germ_chr1=self.chr1   #copy the zygotic genome to the germline genome. The parent will transmit the germ_chr1 and germ_chr2 to the children.
        self.germ_chr2=self.chr2
        self.incomPene=incomPene
        self.varyFit=varyFit 
        self.fitness=self.getFitness()
    
    def embryo_cut(self):  #This method cut the "wt" in chr2 when the "dr" is present. Store the cutted genotype in germ_chr1 and germ_chr2.
        l=[0,0]
        if "dr" in self.chr1:
            if "wt" in self.chr2:
                for i in range(2):
                    if self.chr2[i]=="wt":
                        if random.random() < self.embryo_cutrate:
                            l[i]="ko"
                        else:
                            l[i]=self.chr2[i]
                    else:
                        l[i]=self.chr2[i]
        if(l!=[0,0]):
            self.germ_chr2=tuple(l) #If the list l has been changed, copy the cutted zygotic genome to the germline genome

    def germline_cut(self): #This method do the same thing as the embryo_cut() but function in male germline. Store the cutted genotype in germ_chr1 and germ_chr2.
        l=[0,0]
        if "dr" in self.chr1:
            if "wt" in self.chr2:
                for i in range(2):
                    if self.chr2[i]=="wt":
                        if random.random() < self.germline_cutrate:
                            l[i]="ko"
                        else:
                            l[i]=self.chr2[i]
                    else:
                        l[i]=self.chr2[i]
        if(l!=[0,0]):
            self.germ_chr2=tuple(l) #If the list l has been changed, copy the cutted zygotic genome to the germline genome
        
    def recom(self):    #This recom() method recombine the chr1 and chr2 and transmit haplotype genome to the offsprings.
        rand=random.random()
        if self.sex=="male":
            self.germline_cut() #Cut target sites in male germline
            if self.germ_chr1 in [("dr","wt"),("wt","dr")]: #If the "dr" is heterozygous in chr1.
                if self.germ_chr2 in [("wt","ko"),("ko","wt")]:
                    if rand<100/(300+self.incomPene):
                        return ("dr","ko")
                    elif rand >= 100/(300+self.incomPene) and rand < 200/(300+self.incomPene):
                        return ("dr","wt")
                    elif rand >= 200/(300+self.incomPene) and rand < 300/(300+self.incomPene): #If there are incomplete penetrance, the "wt,wt" haplotype will be viable.
                        return ("wt","wt")
                    else:
                        return ("wt","ko")
                elif self.germ_chr2==("ko","ko"):
                    if rand<100/(100+self.incomPene):
                        return ("dr","ko")
                    else:
                        return ("wt","ko")
                elif self.germ_chr2==("wt","wt"):
                    return (random.choice(self.germ_chr1),random.choice(self.germ_chr2))
            elif self.germ_chr1 == ("wt","wt"):   #If the "dr" is not present
                if self.germ_chr2 in [("wt","ko"),("ko","wt")]:
                    if rand < 100/(100+self.incomPene):
                        return ("wt","wt")
                    else:
                        return ("wt","ko")
                else:
                    return (random.choice(self.germ_chr1),random.choice(self.germ_chr2))
            else:  #If the "dr" is homozygous in chr1.
                s1=random.choice(self.germ_chr1)
                s2=random.choice(self.germ_chr2)
                return (s1,s2)
        else:
            self.embryo_cut() #Cut target site in embryo.
            s1=random.choice(self.germ_chr1)
            s2=random.choice(self.germ_chr2)
            return (s1,s2)
        
    def getFitness(self): #Modify the fitness of individual according to the dosage of "dr" and "ko".
        fitness=1
        if self.sex=="male":
            if self.germ_chr1 in [("dr","wt"),("wt","dr")]: #If the individual is drive heterozygous and both NPG1 are knock out, the pollen grains can be half of that of wild type.
                if self.germ_chr2 in [("ko","ko")]:
                    fitness=0.5*self.varyFit
                elif self.germ_chr2 in [("ko","wt"),("wt","ko")]:
                    fitness=0.75*self.varyFit
        return fitness


class Population:
    def __init__(self,size,capacity,embryo_cutrate,germline_cutrate,chr1,chr2,incomPene,varyFit):
        self.size = size #population size, an integer
        self.capacity=capacity #environment capacity
        self.individuals = [] #A list storing the individual objects.
        self.embryo_cutrate=embryo_cutrate
        self.germline_cutrate=germline_cutrate
        self.incomPene=incomPene #incomplete penetrance
        self.varyFit=varyFit

        for i in range(size):
            sex="hermaphroditic"
            self.individuals.append(Individual(sex, self.embryo_cutrate,self.germline_cutrate, chr1, chr2, self.incomPene, self.varyFit)) #Initialize the list with size of individual objects.

    def get_size(self):
        self.size=len(self.individuals)
        return self.size
    
    def add_pop(self,subpop): #Copy the individual objects from other Population object.
        for ind in subpop.individuals:
            self.individuals.append(ind)
        
    def next_generation(self):
        offspring = []
        weights=[]
        for i in self.individuals:
            weights.append(i.fitness)

        parents_list=random.choices(self.individuals,k=self.size)
        male_list=parents_list[0:int(len(parents_list)/2)]
        female_list=parents_list[int(len(parents_list)/2):]
        self.size=len(self.individuals)
        for i in range(int(self.size/2)):
            [parent1,parent2] =[copy.deepcopy(male_list[i]),copy.deepcopy(female_list[i])]  #If the drive and KO has fitness cost, use random.choices(self.individuals,weights=weights)[0] to set different probability to each individuals.
            parent1.sex="male"     #Set the individual's gender when crossing.
            parent2.sex="female"   #Set the individual's gender when crossing.    
            offspring.extend(self.produce_offspring(parent1, parent2))

        self.individuals = offspring #Cover the parents generation with offsprings
        random.shuffle(self.individuals) #Shuffle the individual list

    def produce_offspring(self, male, female):
        offspring = []
        capacity_fitness_scalling=10/((9*self.get_size()/self.capacity)+1)
        p=capacity_fitness_scalling*0.04
        num_offspring=np.random.binomial(50, p, 1)[0]

        genome1=male.recom()
        genome2=female.recom()
        male_fitness=male.getFitness()
        if male_fitness > 1:
            male_fitness=1
        fit_dict={0.5:(1,0),0.6:(0.8,0.2),0.7:(0.6,0.4),0.8:(0.4,0.6),0.9:(0.2,0.8),1:(0,1)} #Use a dictionary to map the fertility of male drive heterozygote to the probability of producing either one or two offspring.
        if male_fitness in fit_dict:
            if random.random() < fit_dict[male_fitness][0]:
                num_offspring=num_offspring/2
            else:
                num_offspring=num_offspring
        
        for i in range(int(num_offspring)):
            chr1=(genome1[0],genome2[0])
            chr2=(genome1[1],genome2[1])
            sex="hermaphroditic"
            offspring.append(Individual(sex, self.embryo_cutrate, self.germline_cutrate,chr1, chr2, self.incomPene, self.varyFit))
        return offspring


    def count_gene_drive_alleles(self):
        count = 0
        for individual in self.individuals:
            if 'dr' in individual.chr1:
                count += 1
        return count
    
