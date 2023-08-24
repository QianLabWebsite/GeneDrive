#CAIN suppression
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

import argparse
import sys
import pandas as pd
import random
import copy
import numpy as np 
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.pyplot import MultipleLocator

#The fitness cost of drive element and loss-of-function allele of target gene can be set by global variable.
#The gender of every individuals are initially set to hermaphroditic, i.e. bisexual, for most of plant species have both male and female reproductive organs. The gender will be changed to male or female before they mate.
global DRIVE_FITNESS
DRIVE_FITNESS=1 #0.95
global KO_FITNESS
KO_FITNESS=1 #0.99
global CAPACITY
CAPACITY=100000

class Individual:
    def __init__(self, sex, embryo_cutrate,germline_cutrate,chr1,chr2,incomPene,sup):
        self.sex = sex
        self.chr1=chr1
        self.chr2=chr2
        self.embryo_cutrate=embryo_cutrate
        self.germline_cutrate=germline_cutrate

        self.germ_chr1=copy.copy(self.chr1)   #copy the zygotic genome to the germline genome. The male parent will transmit the germ_chr1 and germ_chr2 to the children.
        self.germ_chr2=copy.copy(self.chr2)
        self.embryo_chr1=copy.copy(self.chr1) #copy the zygotic genome to the embryo genome. The female parent will transmit the embryo_chr1 and embryo_chr2 to the children.
        #Attention, the germ_chr and embryo_chr are independent to each other.
        self.embryo_chr2=copy.copy(self.chr2)
        self.incomPene=incomPene
        self.sup=sup #If sup is female, the female of drive homozygote is sterile. If the sup is male, the male of drive homozygote is sterile. 
        self.fitness=self.getFitness()
    
    def embryo_cut(self):  #This method cut the "wt" in chr2 when the "dr" is present. Store the cutted genotype in embryo_chr1 and embryo_chr2.
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
            self.embryo_chr2=tuple(l) #If the list l has been changed, copy the cutted zygotic genome to the embryo genome

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
                if self.embryo_chr1 == ("dr","dr"):
                    #print("female return none")
                    return None
                else:
                    s1=random.choice(self.embryo_chr1)
                    s2=random.choice(self.embryo_chr2)
                    return (s1,s2)    
            else:
                #Female homozygous is fertile
                s1=random.choice(self.embryo_chr1)
                s2=random.choice(self.embryo_chr2)
                return (s1,s2)
        
    def getFitness(self): #Modify the fitness of individual according to the dosage of "dr" and "ko".
        c1_value=DRIVE_FITNESS if "dr"==self.chr1[0] else 1
        c2_value=DRIVE_FITNESS if "dr"==self.chr1[1] else 1

        c1_value=c1_value*KO_FITNESS if "ko"==self.chr2[0] else c1_value*1
        c2_value=c2_value*KO_FITNESS if "ko"==self.chr2[1] else c2_value*1

        return c1_value*c2_value



class Population:
    def __init__(self,size,embryo_cutrate,germline_cutrate,chr1,chr2,incomPene,sup):
        self.size = size #population size, an integer
        self.individuals = [] #A list storing the individual objects.
        self.embryo_cutrate=embryo_cutrate
        self.germline_cutrate=germline_cutrate
        self.incomPene=incomPene #incomplete penetrance
        self.sup=sup

        for i in range(size):
            sex="hermaphroditic"
            # sex=""
            # rand=random.random()
            # if rand<0.5:
            #     sex="male"
            # else:
            #     sex="female"
            self.individuals.append(Individual(sex, self.embryo_cutrate,self.germline_cutrate, chr1, chr2, self.incomPene,self.sup)) #Initialize the list with size of individual objects.


    def add_pop(self,subpop): #Copy the individual objects from other Population object.
        for ind in subpop.individuals:
            self.individuals.append(ind)
        self.size=len(self.individuals)

    def next_generation(self):
        offspring = []
        weights=[]
        for i in self.individuals:
            weights.append(i.fitness)
        
        for i in range(self.size):
            parent1 = copy.deepcopy(random.choice(self.individuals)) #If the drive and KO has fitness cost, use random.choices(self.individuals,weights=weights)[0] to set different probability to each individuals.
            parent1.sex="male"     #Set the individual's gender when crossing.
            
            parent2 = copy.deepcopy(random.choice(self.individuals)) #random.choices(self.individuals,weights=weights)[0]
            parent2.sex="female"   #Set the individual's gender when crossing.
            #print([parent2.embryo_chr1,parent2.embryo_chr2])

            offspring.extend(self.produce_offspring(parent1, parent2))

        self.individuals = offspring   #Cover the parent generation using offspring generation
        self.size=len(self.individuals)
        random.shuffle(self.individuals) #Shuffle the individual list

    def produce_offspring(self, male, female):
        offspring = []
        capacity_fitness_scalling=10/((9*self.size/CAPACITY)+1)
        p=capacity_fitness_scalling*0.02
        num_offspring=np.random.binomial(50, p, 1)[0]
        for i in range(num_offspring):
        #for i in range(1):
            genome1=male.recom()
            genome2=female.recom()
            #if genome2 == None:
            #    print("No")
            if genome1 != None and genome2 != None:
                #print("Yes")
                chr1=(genome1[0],genome2[0])
                chr2=(genome1[1],genome2[1])
                sex="hermaphroditic"
                #sex=random.choice(["male","female"])
                offspring.append(Individual(sex, self.embryo_cutrate, self.germline_cutrate,chr1, chr2, self.incomPene,self.sup))
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
        else:
            for individual in self.individuals:
                if 'dr' in individual.chr1:
                    count += 1
                    if 'dr' == individual.chr1[0] and 'dr' == individual.chr1[1]:
                        homo+=1
                    else:
                        heter+=1

        return (count,homo,heter)


__doc__='''
Dynamics simulation of the CAIN(TADS) suppression drive.

Set the parameters to the initial populations and the script will generate frequency and individual number of drive carriers along with the increasing of generations
'''
__author__="Bingke Jiao"
__mail__="bkjiao@genetics.ac.cn"
__date__="2023.8.23"
__version__="1.0"

def main():
    parser=argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawTextHelpFormatter,
                                   epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
    parser.add_argument('-w',help='Size of wild population.[default 9900]',dest='wild',type=int,default='9900',required=False)
    parser.add_argument('-d',help='Size of heterozygous individuals carring drive.[default 100]',dest='drive',type=int,default='100',required=False)
    parser.add_argument('-t',help='Generation of propagation.[default 50]',dest='generation',type=int,default='50',required=False)
    parser.add_argument('-e',help='Embryo DNA cleavage efficiency, float number, between 0-1',dest='embryo_cutrate',type=float,required=True)
    parser.add_argument('-g',help='Germline DNA cleavage efficiency,float number, between 0-1',dest='germline_cutrate',type=float,required=True)
    parser.add_argument('-s',help='Set the fertility gene where drive located in, character, male or female',dest='sex',type=str,default="male",required=False)
    parser.add_argument('-o',help='Output file prefix, character',dest='output',type=str,default="CAIN_suppression",required=False)
    args=parser.parse_args()
    #pop.size=args.wild+args.drive
    wild_num=args.wild
    drive_num=args.drive
    CAPACITY=args.wild+args.drive

    #for i in [(args.embryo_cutrate,args.germline_cutrate,args.sex)]: #[(0,0.5,0),(0.5,0.5,0),(0,0.984,4.0),(0.941,0.984,4.0),(1,1,0)]

    a,b,c=args.embryo_cutrate,args.germline_cutrate,args.sex
    pop=Population(size=wild_num,embryo_cutrate=a,germline_cutrate=b,incomPene=0,chr1=("wt","wt"),chr2=("wt","wt"),sup=c) #
    pop.add_pop(Population(size=drive_num,embryo_cutrate=a,germline_cutrate=b,incomPene=0,chr1=("dr","wt"),chr2=("wt","wt"),sup=c)) #

    freq=[drive_num/pop.size]
    freq_homo=[0]
    freq_heter=[drive_num/pop.size]
    
    popsize=[pop.size]
    drive_carrier=[drive_num]
    homo_carrier=[0]
    heter_carrier=[drive_num]

    for n in range(args.generation):
        pop.next_generation()
        (drive_ind,homo_ind,heter_ind)=pop.count_gene_drive_alleles(allele="False")

        if pop.size == 0:
            freq.append(0)
            freq_homo.append(0)
            freq_heter.append(0)
            popsize.append(0)
            drive_carrier.append(0)
            homo_carrier.append(0)
            heter_carrier.append(0)
        else:
            freq.append(drive_ind/pop.size)
            freq_homo.append(homo_ind/pop.size)
            freq_heter.append(heter_ind/pop.size)
            popsize.append(pop.size)
            drive_carrier.append(drive_ind)
            homo_carrier.append(homo_ind)
            heter_carrier.append(heter_ind)

    # print(freq)
    # print(freq_homo)
    # print(freq_heter)
    # print(drive_carrier)
    # print(homo_carrier)
    # print(heter_carrier)
    # print(popsize)

    dataframe = pd.DataFrame({'freq': freq,'freq_homo':freq_homo,'freq_heter':freq_heter,'drive_carrier': drive_carrier,'homo_carrier':homo_carrier,'heter_carrier':heter_carrier,'popsize': popsize})
    #将DataFrame存储为csv,index表示是否显示行名，default=True
    dataframe.to_csv("CAIN_suppresion.density_regulate_production.hermaphroditic.csv",index=True,sep=',') #,mode="a"

    end=0 
    i=0
    while i < len(dataframe):
        if dataframe['freq'][i]==1.0:
            end=i
        i+=1
    if end==0:
        print("The population did not collapse within 50 generations. Please increase the number of generations.")
        sys.exit(0)
    mpl.rcParams['pdf.fonttype']=42
    mpl.rcParams['ps.fonttype']=42
  
    fig, axs = plt.subplots(2, 1, figsize=(10,5), sharex=True)
    freq=dataframe['freq'][:end+1] #27
    freq_homo=dataframe['freq_homo'][:end+1]
    freq_heter=dataframe['freq_heter'][:end+1]
    popsize=dataframe['popsize'][:end+2] #28
    drive_carrier=dataframe['drive_carrier'][:end+2] #28
    homo_carrier=dataframe['homo_carrier'][:end+2]
    heter_carrier=dataframe['heter_carrier'][:end+2]
    x=list(range(end+1))
    #axs[0].hlines(1,0,25,color="dodgerblue")
    axs[0].plot(freq,linestyle="-",label="Frequency of CAIN carriers",color="orange")
    #axs[0].plot(freq_homo,linestyle="-",label="Frequency of CAIN homozygotes",color="darkgoldenrod")
    axs[0].plot(freq_heter,linestyle="-",label="Frequency of CAIN heterozygotes",color="yellow")
    axs[0].fill_between(x,freq,1,alpha=1,facecolor="dodgerblue")
    axs[0].fill_between(x,freq_heter,freq,alpha=1,facecolor="orange")
    axs[0].fill_between(x,freq_heter,alpha=1,facecolor="yellow")

    axs[1].plot(popsize,linestyle="-",label="Number of wild-type individuals",color="dodgerblue")
    axs[1].plot(drive_carrier,linestyle="-",label="Number of CAIN carriers",color="orange")
    #axs[1].plot(homo_carrier,linestyle="-",label="Number of CAIN homozygotes",color="darkgoldenrod")
    axs[1].plot(heter_carrier,linestyle="-",label="Number of CAIN heterozygotes",color="yellow")
    x=list(range(end+2))
    axs[1].fill_between(x,drive_carrier,popsize,alpha=1,facecolor="dodgerblue")
    # axs[1].plot(popsize2,linestyle="-",label="Female of drive homozygote is sterile")
    axs[1].fill_between(x,heter_carrier,drive_carrier,alpha=1,facecolor="orange")
    axs[1].fill_between(x,heter_carrier,alpha=1,facecolor="yellow")
    axs[0].set_ylim(0,1.1)
    axs[1].set_ylim(0,CAPACITY*1.05)

    plt.margins(x=0)
    plt.xticks(range(0,end+4,2))
    plt.xlabel("Generation",fontsize=12)
    axs[0].set_ylabel("Frequency",fontsize=10)
    axs[1].set_ylabel("Number of individuals in population",fontsize=10)


    axs[0].axvline(end+1,linestyle="dotted")
    axs[1].axvline(end+1,linestyle="dotted")
    yellow_patch = mpatches.Patch(color='yellow', label='CAIN heterozygotes')
    orange_patch = mpatches.Patch(color='orange', label='CAIN homozygotes')
    blue_patch = mpatches.Patch(color='dodgerblue', label='wild-type individuals')
    plt.legend(loc="upper left",bbox_to_anchor=(1.05, 1.0),handles=[yellow_patch,orange_patch,blue_patch])
    axs[0].set_title("CAIN suppresion(with density regulation production)")
    plt.tight_layout()
    plt.show()

    plt.savefig(args.output+".density_regulation_production.hermaphroditic.pdf",format="pdf")

if __name__=="__main__":
	main()
