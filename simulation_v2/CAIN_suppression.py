#TADS suppression use chr only, remove the embro_chr.
'''
Author: Bingke Jiao, supervised by Professor Wenfeng Qian.
Ph.D. Yang Liu and Weiguang Wang gave valuable advice in the design of this script.

This script implements the dynamic simulation of CRISPR-mediated toxin-antidote gene drive.
The main features of this scipt include:
1. It is a forward genetic simulation. Set parameters to the initial population and the population will propagate. You can monitor the frequency of drive carriers with the increasing of generations.
2. It is an individual based, stocastic model, which using Wright-Fisher model.
3  It uses two classes: the individual class and population class.

The design of this script draws inspiration from the SLiM software developed by Benjamin C. Haller and Philipp W. Messer. and the Toxin-antidote scripts written by Samuel E. Champer and Jackson Champer
'''
import random
import copy
import numpy as np 

#The fitness cost of drive element and loss-of-function allele of target gene can be set by global variable.
#The gender of every individuals are initially set to hermaphroditic, i.e. bisexual, for most of plant species have 
#both male and female reproductive organs. The gender will be changed to male or female before they mate.


class Individual:
    def __init__(self, sex, embryo_cutrate,germline_cutrate,chr1,chr2,incomPene,sup,fitCoef):
        self.sex = sex
        self.chr1=chr1
        self.chr2=chr2
        self.embryo_cutrate=embryo_cutrate
        self.germline_cutrate=germline_cutrate

        self.germ_chr1=self.chr1   #copy the zygotic genome to the germline genome. The parents will transmit the germ_chr1 and germ_chr2 to the children.
        self.germ_chr2=self.chr2   
        
        self.incomPene=incomPene
        self.sup=sup           #If sup is female, the female of drive homozygote is sterile. If the sup is male, the male of drive homozygote is sterile. If the sup is both, the homozygotes will die.
        self.fitCoef=fitCoef
        self.fitness=self.getFitness()
    
    def embryo_cut(self):  #This method cut the "wt" in chr2 when the "dr" is present, which functions in the female. Store the cutted genotype in germ_chr1 and germ_chr2.
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
                            print("The self.chr2 is "+str(l[i]))
                    else:
                        l[i]=self.chr2[i]
        if(l!=[0,0]):
            self.germ_chr2=tuple(l) #If the list l has been changed, copy the cutted zygotic genome to the germline genome


    def recom(self):    #This recom() method recombine the chr1 and chr2 and transmit haplotype genome to the offsprings.
        rand=random.random()
        if self.sex=="male":
            self.germline_cut() #Cut target sites in germline
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
                if self.sup == "male":
                    #Male homozygous is sterile
                    return None
                else:
                    #Male homozygous is fertile
                    s1=random.choice(self.germ_chr1)
                    s2=random.choice(self.germ_chr2)
                    return (s1,s2)
    
        elif self.sex=="female":
            self.embryo_cut() #Cut target site in embryo.
            if self.sup == "female":
                #Female homozygous is sterile
                if self.germ_chr1 == ("dr","dr"):
                    return None
                else:
                    s1=random.choice(self.germ_chr1)
                    s2=random.choice(self.germ_chr2)
                    return (s1,s2)    
            else:
                #Female homozygous is fertile
                s1=random.choice(self.germ_chr1)
                s2=random.choice(self.germ_chr2)
                return (s1,s2)
        else:
            print("This individual is herma")

        
    def getFitness(self): #Modify the fitness of individual according to the dosage of "dr" and "ko".
        c1_value=self.fitCoef if "dr"==self.chr1[0] else 1
        c2_value=self.fitCoef if "dr"==self.chr1[1] else 1

        # c1_value=c1_value*KO_FITNESS if "ko"==self.chr2[0] else c1_value*1
        # c2_value=c2_value*KO_FITNESS if "ko"==self.chr2[1] else c2_value*1

        return c1_value+c2_value



class Population:
    def __init__(self,size,capacity,embryo_cutrate,germline_cutrate,chr1,chr2,incomPene,sup,fitCoef):
        self.size = size #population size, an integer
        self.capacity=capacity #environment capacity
        self.individuals = [] #A list storing the individual objects.
        self.embryo_cutrate=embryo_cutrate 
        self.germline_cutrate=germline_cutrate
        self.incomPene=incomPene #incomplete penetrance
        self.sup=sup #homozygotes of this are sterile. Set "both" to make the homozygotes non-viable.
        self.fitCoef=fitCoef

        for i in range(size):
            sex="hermaphroditic"
            self.individuals.append(Individual(sex, self.embryo_cutrate,self.germline_cutrate, chr1, chr2, self.incomPene,self.sup,self.fitCoef)) #Initialize the list with size of individual objects.


    def add_pop(self,subpop): #Copy the individual objects from other Population object.
        for ind in subpop.individuals:
            self.individuals.append(ind)
        self.size=len(self.individuals)

    def next_generation(self):
        offspring = []
        weights=[]
        for i in self.individuals:
            weights.append(i.fitness)
        self.size=len(self.individuals)
        if self.size==0:
            return None
        else:
            parents_list=random.choices(self.individuals,k=self.size,weights=weights) 
            male_list=parents_list[0:int(len(parents_list)/2)]
            female_list=parents_list[int(len(parents_list)/2):]

            for i in range(int(self.size/2)):
                [parent1,parent2] =[copy.deepcopy(male_list[i]),copy.deepcopy(female_list[i])]  #use copy.deepcopy to get a real copy of an individual.
                parent1.sex="male"     #Set the individual's gender when crossing.    
                parent2.sex="female"   #Set the individual's gender when crossing.    
                offspring.extend(self.produce_offspring(parent1, parent2))

            self.individuals = offspring   #Cover the parent generation using offspring generation
            self.size=len(self.individuals)
            random.shuffle(self.individuals) #Shuffle the individual list

    def produce_offspring(self, male, female):
        offspring = []
        capacity_fitness_scalling=10/((9*self.size/self.capacity)+1)
        p=capacity_fitness_scalling*0.04
        num_offspring=np.random.binomial(50, p, 1)[0]
        for i in range(num_offspring):
            genome1=male.recom()
            genome2=female.recom()
            if genome1 != None and genome2 != None:
                chr1=(genome1[0],genome2[0])
                chr2=(genome1[1],genome2[1])
                if self.sup=="both":
                    if chr1==("dr","dr"):
                        pass
                    else:
                        sex="hermaphroditic"
                        ind=Individual(sex, self.embryo_cutrate, self.germline_cutrate,chr1, chr2, self.incomPene,self.sup,self.fitCoef)
                        offspring.append(ind)
                else:
                        sex="hermaphroditic"
                        ind=Individual(sex, self.embryo_cutrate, self.germline_cutrate,chr1, chr2, self.incomPene,self.sup,self.fitCoef)
                        offspring.append(ind)
        return offspring


    def count_gene_drive_alleles(self,allele):
        count = 0
        heter=0
        homo=0
        if allele=="True": #If the allele is True, count the drive allele number, else, count the drive carriers number.
            for individual in self.individuals:
                for i in range(2):
                    if 'dr' == individual.chr1[i]:
                        count += 1
                if 'dr' in individual.chr1:
                    if 'dr' == individual.chr1[0] and 'dr' == individual.chr1[1]:
                        homo+=1
                    else:
                        heter+=1
        else:
            for individual in self.individuals:
                if 'dr' in individual.chr1:
                    count += 1
                    if 'dr' == individual.chr1[0] and 'dr' == individual.chr1[1]:
                        homo+=1
                    else:
                        heter+=1

        return (count,homo,heter)
    