'''
Author: Bingke Jiao, supervised by Professor Wenfeng Qian.
Ph.D. Yang Liu and Weiguang Wang gave valuable advice in the design of this script.

This script implements the dynamics simulation of CRISPR-mediated toxin-antidote gene drive.
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
'''
The fitness cost of drive element and loss-of-function allele of target gene can be set by global variable.
The gender of every individuals are initially set to hermaphroditic, i.e. bisexual, for most of plant species have both male 
and female reproductive organs.The gender will be changed to male or female before they mate.
'''
class Individual:
    def __init__(self, sex, embryo_cutrate,germline_cutrate,chr1,chr2,incomPene,birthsCoef,fitCoef,drive="Yes",linkage="No"):
        self.sex = sex
        self.chr1=chr1
        self.chr2=chr2
        self.embryo_cutrate=embryo_cutrate
        self.germline_cutrate=germline_cutrate

        self.germ_chr1=self.chr1  #copy the zygotic genome to the germline genome. The parent will transmit the germ_chr1 and germ_chr2 to the children.
        self.germ_chr2=self.chr2
        self.incomPene=incomPene
        self.birthsCoef=birthsCoef
        self.fitCoef=fitCoef
        self.drive=drive
        self.linkage=linkage 
    
    def embryo_cut(self):  #This method used cut the "wt" in chr2 when the "dr" is present. Store the cutted genotype in germ_chr1 and germ_chr2.
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
            self.germ_chr2=tuple(l) #If the list l has been changed, store the cutted zygotic genome to the germline genome
            

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
            self.germ_chr2=tuple(l) #If the list l has been changed, store the cutted zygotic genome to the germline genome
    
    def getGamete(self): #This method is applied when the driver follows Mendelian inheritance laws.
        return (random.choice(self.chr1),random.choice(self.chr2))

    def recom(self):    #This recom() method recombine the chr1 and chr2 and transmit haplotype genome to the offspring.
        rand=random.random()
        
        if self.linkage=="No": #If the drive is not in linkage with the target gene
            if self.sex=="male":
                #self.germline_cut() #Cut target sites in male germline
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
                #self.embryo_cut() #Cut target site in embryo.
                s1=random.choice(self.germ_chr1)
                s2=random.choice(self.germ_chr2)
                return (s1,s2)
            
        elif self.linkage=="Yes": #If the drive is in linkage with the target gene.
            if self.sex=="male":
                if self.germ_chr1==("dr","dr"):
                    return ("dr",random.choice(self.germ_chr2))
                if self.germ_chr1 in [("dr","wt"),("wt","dr")]:
                    if self.germ_chr2==("wt","wt"):
                        return random.choice([("dr","wt"),("wt","wt")])
                if self.germ_chr1==("dr","wt"):
                    if self.germ_chr2==("ko","wt"):
                        return random.choice([("dr","ko"),("wt","wt")])
                    elif self.germ_chr2==("wt","ko"):
                        return ("dr","wt")
                    elif self.germ_chr2==("ko","ko"):
                        return ("dr","ko")
                elif self.germ_chr1==("wt","dr"):
                    if self.germ_chr2==("ko","wt"):
                        return ("dr","wt")
                    elif self.germ_chr2==("wt","ko"):
                        return random.choice([("wt","wt"),("dr","ko")])
                    elif self.germ_chr2==("ko","ko"):
                        return ("dr","ko")
                elif self.germ_chr1==("wt","wt"):
                    if self.germ_chr2 in [("ko","wt"),("wt","ko")]:
                        return ("wt","wt")
                    elif self.germ_chr2==("wt","wt"):
                        return ("wt","wt")
                    elif self.germ_chr2==("ko","ko"):
                        print("Male is sterile!")
                        return None
            if self.sex=="female":
                return random.choice([(self.germ_chr1[0],self.germ_chr2[0]),(self.germ_chr1[1],self.germ_chr2[1])])
        
    def getBirthsCoef(self): #Modify the birth coefficient of individual according to its' genotype.
        coef=1
        if self.sex=="male":
            if self.germ_chr1 in [("dr","wt"),("wt","dr")]:
                if self.germ_chr2 in [("ko","ko")]:
                    coef=0.5*self.birthsCoef
                elif self.germ_chr2 in [("ko","wt"),("wt","ko")]:
                    coef=0.75*self.birthsCoef
        return coef
    
    def getViableCoef(self): #Calculate the fitness value of individual according to the dosage of dr
        c1_value=self.fitCoef if "dr"==self.chr1[0] else 1
        c2_value=self.fitCoef if "dr"==self.chr1[1] else 1
        self.viableCoef=c1_value+c2_value #additive effect
        return self.viableCoef



class Population:
    def __init__(self,size,capacity,embryo_cutrate,germline_cutrate,chr1,chr2,incomPene,birthsCoef,fitCoef,inbCoef,drive="Yes",linkage="No"):
        self.size = size #population size, an integer
        self.capacity=capacity #environment capacity
        self.individuals = [] #A list storing the individual objects.
        self.embryo_cutrate=embryo_cutrate #embryo cut rate of Cas9, used in female individual
        self.germline_cutrate=germline_cutrate #germline cut rate of Cas9, used in male individual
        self.incomPene=incomPene #incomplete penetrance
        self.birthsCoef=birthsCoef #births coefficient for drive carrier
        self.fitCoef=fitCoef #fitness value of one drive allele
        self.inbCoef=inbCoef #inbreeding coefficient
        self.drive=drive #The boolean value indicates whether CAIN functions as a driver or is simply inherited according to Mendelian laws. 
        self.linkage=linkage #The boolean value indicates whether the CAIN is in linkage with the target gene or not.

        for i in range(size):
            sex="hermaphroditic"
            self.individuals.append(Individual(sex, self.embryo_cutrate,self.germline_cutrate, chr1, chr2, self.incomPene, self.birthsCoef, self.fitCoef, self.drive, self.linkage)) #Initialize the list with size of individual objects.

    def get_size(self):
        self.size=len(self.individuals)
        return self.size
    
    def add_pop(self,subpop): #Copy the individual objects from other Population object.
        for ind in subpop.individuals:
            self.individuals.append(ind)
        self.size=self.get_size()
        
    def next_generation(self):
        offspring = []
        self.size=len(self.individuals)

        weights=[]
        for i in self.individuals:
            weights.append(i.getViableCoef())

        parents_list=random.choices(self.individuals,k=self.size,weights=weights)
        male_list=parents_list[0:int(len(parents_list)/2)]
        female_list=parents_list[int(len(parents_list)/2):]

        for i in range(int(self.size/2)):
            [parent1,parent2] =[copy.deepcopy(male_list[i]),copy.deepcopy(female_list[i])]  #If the drive and KO has fitness cost, use random.choices(self.individuals,weights=weights)[0] to set different probability to each individuals.
            parent1.sex="male"     #Set the individual's gender when crossing.

            if random.random()<=self.inbCoef:
                parent2=copy.deepcopy(parent1)
                
            parent2.sex="female"   #Set the individual's gender when crossing.    
            offspring.extend(self.produce_offspring(parent1, parent2))
        
        self.individuals =offspring  #Cover the parent generation using offspring generation
        
        random.shuffle(self.individuals) #Shuffle the individual list

    def produce_offspring(self, male, female):
        offspring = []
        capacity_fitness_scalling=10/((9*self.get_size()/self.capacity)+1) #Use the density regulation production strategy
        p=capacity_fitness_scalling*0.04  #0.4
        num_offspring=np.random.binomial(50, p, 1)[0] #the mean of offspring number is 50*p=2

        if self.drive=="Yes": #If the driver can induce the transmission ratio distortion
            male.germline_cut() #Cas9 cuts the NPG1 in the male individual if possible
            female.embryo_cut() #Cas9 cuts the NPG1 in the female individual if possible

            male_birthsCoef=male.getBirthsCoef()
            if male_birthsCoef > 1:
                male_birthsCoef=1
            num_offspring=num_offspring*male_birthsCoef  #The male_birthsCoef affects its number of offsprings.

            for i in range(int(num_offspring)):
                genome1=male.recom()
                genome2=female.recom()

                chr1=(genome1[0],genome2[0])
                chr2=(genome1[1],genome2[1])
                sex="hermaphroditic"
                offspring.append(Individual(sex, self.embryo_cutrate, self.germline_cutrate,chr1, chr2, self.incomPene, self.birthsCoef,self.fitCoef,self.drive,self.linkage))
        elif self.drive=="No": #If the driver behaves according to the Mendelian law.
            for i in range(int(num_offspring)):
                genome1=male.getGamete()
                genome2=female.getGamete()
                chr1=(genome1[0],genome2[0])
                chr2=(genome1[1],genome2[1])                
                sex="hermaphroditic"
                offspring.append(Individual(sex, self.embryo_cutrate, self.germline_cutrate,chr1, chr2, self.incomPene, self.birthsCoef,self.fitCoef,self.drive,self.linkage))
        
        return offspring


    def count_gene_drive_alleles(self):
        count = 0

        for individual in self.individuals:
            if 'dr' in individual.chr1:
                count += 1
        return count
    