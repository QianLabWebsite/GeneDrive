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
import argparse
import copy
import numpy as np 
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
from matplotlib.pyplot import MultipleLocator

#The fitness cost of drive element and loss-of-function allele of target gene can be set by global variable.
#The gender of every individuals are initially set to hermaphroditic, i.e. bisexual, for most of plant species have both male and female reproductive organs. The gender will be changed to male or female before they mate.
global DRIVE_FITNESS
DRIVE_FITNESS=1 #0.95
global KO_FITNESS
KO_FITNESS=1 #0.99

class Individual:
	def __init__(self, sex, embryo_cutrate,germline_cutrate,chr1,chr2,incomPene):
		self.sex = sex
		self.chr1=chr1
		self.chr2=chr2
		self.embryo_cutrate=embryo_cutrate
		self.germline_cutrate=germline_cutrate

		self.germ_chr1=copy.copy(self.chr1)   #copy the zygotic genome to the germline genome. The male parent will transmit the germ_chr1 and germ_chr2 to the children.
		self.germ_chr2=copy.copy(self.chr2)
		self.embryo_chr1=copy.copy(self.chr1) #copy the zygotic genome to the embryo genome. The female parent will transmit the embryo_chr1 and embryo_chr2 to the children. Attention, the germ_chr and embryo_chr are independent to each other.
		self.embryo_chr2=copy.copy(self.chr2)
		self.incomPene=incomPene*100
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

	def recom(self):	#This recom() method recombine the chr1 and chr2 and transmit haplotype genome to the offsprings.
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
	def __init__(self,size,embryo_cutrate,germline_cutrate,chr1,chr2,incomPene):
		self.size = size #population size, an integer
		self.individuals = [] #A list storing the individual objects.
		self.embryo_cutrate=embryo_cutrate
		self.germline_cutrate=germline_cutrate
		self.incomPene=incomPene #incomplete penetrance

		for i in range(size):
			sex="hermaphroditic"
			self.individuals.append(Individual(sex, self.embryo_cutrate,self.germline_cutrate, chr1, chr2, self.incomPene)) #Initialize the list with size of individual objects.


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
			parent1 = random.choice(self.individuals) #If the drive and KO has fitness cost, use random.choices(self.individuals,weights=weights)[0] to set different probability to each individuals.
			parent1.sex="male"	 #Set the individual's gender when crossing.
			parent2 = random.choice(self.individuals) #random.choices(self.individuals,weights=weights)[0]
			parent2.sex="female"   #Set the individual's gender when crossing.
			offspring.extend(self.produce_offspring(parent1, parent2))

			parent1.sex="hermaphroditic" #Change the gender to the bisexual state.
			parent2.sex="hermaphroditic" #Change the gender to the bisexual state.

		self.individuals = offspring   #Cover the parent generation using offspring generation
		random.shuffle(self.individuals) #Shuffle the individual list

	def produce_offspring(self, male, female):
		offspring = []
		for i in range(1):
			genome1=male.recom()
			genome2=female.recom()

			chr1=(genome1[0],genome2[0])
			chr2=(genome1[1],genome2[1])
			sex="hermaphroditic"
			offspring.append(Individual(sex, self.embryo_cutrate, self.germline_cutrate,chr1, chr2, self.incomPene))
		return offspring


	def count_gene_drive_alleles(self):
		count = 0

		for individual in self.individuals:
			if 'dr' in individual.chr1:
				count += 1
		return count

__doc__='''
	Dynamics simulation of the CRISPR-mediated toxin-antidote gene drive.

	Set the parameters to the initial populations and the script will generate frequency of drive carriers along with the increasing of generations
'''
__author__="Bingke Jiao"
__mail__="bkjiao@genetics.ac.cn"
__date__="2023.7.11"
__version__="1.0"

def main():
	parser=argparse.ArgumentParser(description=__doc__,
								   formatter_class=argparse.RawTextHelpFormatter,
								   epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-w',help='Size of wild population.[default 9900]',dest='wild',type=int,default='9900',required=False)
	parser.add_argument('-d',help='Size of heterozygous individuals carring drive.[default 100]',dest='drive',type=int,default='100',required=False)
	parser.add_argument('-t',help='Generation of propagation.[default 50]',dest='generation',type=int,default='50',required=False)
	parser.add_argument('-e',help='Embryo DNA cleavage efficiency, float number, between 0-1',dest='embryo_cutrate',type=float,required=True)
	parser.add_argument('-g',help='Male germline cleavage efficiency,float number, between 0-1',dest='germline_cutrate',type=float,required=True)
	parser.add_argument('-i',help='Incomplete penetrance rate, float number, between 0-1',dest='incomp_pene',type=float,required=True)
	args=parser.parse_args()
	total_num=args.wild+args.drive

	mpl.rcParams['pdf.fonttype']=42
	mpl.rcParams['ps.fonttype']=42
	plt.figure(figsize=(10,5))
	for i in [(args.embryo_cutrate,args.germline_cutrate,args.incomp_pene)]: #[(0,0.5,0),(0.5,0.5,0),(0,0.984,4.0),(0.941,0.984,4.0),(1,1,0)]
		a,b,c=i[0],i[1],i[2]
		pop=Population(size=args.wild,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),incomPene=c) #
		pop.add_pop(Population(size=args.drive,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),incomPene=c)) #
		num=[args.drive]
		for n in range(args.generation):
			pop.next_generation()
			count=pop.count_gene_drive_alleles()
			num.append(count)
		
		ratio=[x/total_num for x in num]
		print(ratio)
		f=open("drive_carriers_freq.embryoRate{0}_germRate{1}_incompene{2}.txt".format(args.embryo_cutrate,args.germline_cutrate,args.incomp_pene),"w")
		f.writelines("Generation\tDriveCarriersFreq\n")
		for r in range(len(ratio)):
			f.writelines(str(r)+'\t'+str(ratio[r])+'\n')
		f.close()
		#linestyle="--"
		plt.plot(ratio,linestyle="-",label="Embryonic DNA cleavage efficiency: "+str(format(a*100,'.1f'))+"%;\n"+"Gross male germline cleavage efficiency: "+str(format(b*100,'.1f'))+"%;\n"+"Incomplete penetrance: "+str(format(c*100,'.1f'))+"%") 
		#plt.plot(ratio,linestyle="-",label="Gross male germline cleavage efficiency: "+str(format(b*100,'.1f'))+"%;\n"+"Incomplete penetrance rate: "+str(c)+"%") 
		plt.axis([0, 50, 0, 1.1]);#define the plotting range
		plt.xlabel("Generation",fontsize=12)
		plt.ylabel("Frequency of drive carriers",fontsize=12)
		plt.legend(loc="lower right")

	plt.savefig("drive_carriers_freq.embryoRate{0}_germRate{1}_incompene{2}.pdf".format(args.embryo_cutrate,args.germline_cutrate,args.incomp_pene),format="pdf")
	# array = np.ones([1,100], dtype = int)
	# row=0
	# for i in [(args.embryo_cutrate,args.germline_cutrate,args.incomp_pene)]: #range(10)
	#	 a,b,c=i[0],i[1],i[2]
	#	 num=[]
	#	 for j in range(100):
	#		 pop=Population(size=args.wild,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),incomPene=c) 
	#		 pop.add_pop(Population(size=args.drive,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),incomPene=c)) 
	#		 n=0
	#		 while 1:
	#			 n+=1
	#			 pop.next_generation()
	#			 count=pop.count_gene_drive_alleles()
	#			 ratio=count/total_num
	#			 if ratio>0.99:
	#				 num.append(n)
	#				 break
	#			 if n > 100:
	#				 num.append(n)
	#				 break

	#	 array[row]=num
	#	 row+=1
	# np.savetxt("fixed_generation.txt",array,delimiter=',', fmt='%.0f')

if __name__=="__main__":
	main()
