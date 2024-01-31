'''
Author: Bingke Jiao, supervised by Professor Wenfeng Qian.
Ph.D. Yang Liu and Weiguang Wang gave valuable advice in the design of this script.

This script implements the dynamic simulation of CRISPR-mediated toxin-antidote gene drive.
The main features of this scipt include:
1. It is a forward genetic simulation. Set parameters to the initial population and the population will propagate. You can monitor the frequency of drive carriers with the increasing of generations.
2. It is an individual based, stocastic model, which using Wright-Fisher model.
3  It uses two classes: the individual class and population class.

The design of this script draws inspiration from the SLiM software written by Benjamin C. Haller and Philipp W. Messer. and the Toxin-antidote scripts written by Samuel E. Champer and Jackson Champer
'''

import random
import argparse
import copy
import numpy as np 
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
from matplotlib.pyplot import MultipleLocator
from CAIN_viability_and_inbreeding import Population

__doc__='''
	Dynamics simulation of the CRISPR-mediated toxin-antidote gene drive.

	Set the parameters to the initial populations and the script will generate frequency of drive carriers along with the increasing of generations
'''
__author__="Bingke Jiao"
__mail__="bkjiao@genetics.ac.cn"
__date__="2024.1.31"
__version__="1.0"

def main():
	parser=argparse.ArgumentParser(description=__doc__,
								   formatter_class=argparse.RawTextHelpFormatter,
								   epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-w',help='Size of wild population.[default 9900]',dest='wild',type=int,default='9900',required=False)
	parser.add_argument('-d',help='Size of heterozygous individuals carring drive.[default 100]',dest='drive',type=int,default='100',required=False)
	parser.add_argument('-c',help='Environment capacity.[default 10000]',dest='capa',type=int,default='10000',required=False)
	parser.add_argument('-t',help='Generation of propagation.[default 50]',dest='generation',type=int,default='50',required=False)
	parser.add_argument('-f',help='Female cleavage efficiency, float number, between 0-1',dest='embryo_cutrate',type=float,required=True)
	parser.add_argument('-g',help='Male germline cleavage efficiency,float number, between 0-1',dest='germline_cutrate',type=float,required=True)
	parser.add_argument('-i',help='Incomplete penetrance rate, float number, between 0-1',dest='incomp_pene',type=float,required=True)
	parser.add_argument('-b',help='Births coefficient affecting the fertility of CAIN heterozygotes, float number, between 1-2',dest='birth_fit',type=float,default='2',required=False)
	parser.add_argument('-r',help='Relative fitness value comparing of CAIN comparing with the wild-type, float number, between 0-1',dest='drive_fitness',type=float,default='1',required=False)
	parser.add_argument('-in',help='Inbreeding level, float number, between 0-1',dest='inbreeding',type=float,default='0',required=False)
	parser.add_argument('-dr',help='Value indicating whether CAIN is inherited as a drive or obey the rule of Mendelian law, string, Yes or No', dest='supMend',type=str,default='Yes',required=False)
	parser.add_argument('-li',help='Value indicating whether CAIN is in linkage with the target gene or not, string, Yes or No', dest='linkage',type=str,default='No',required=False)
	args=parser.parse_args()
	#total_num=args.wild+args.drive

	mpl.rcParams['pdf.fonttype']=42
	mpl.rcParams['ps.fonttype']=42
	plt.figure(figsize=(10,5))
	x=np.empty(shape=[0,51],dtype=float)
	for i in [(args.embryo_cutrate,args.germline_cutrate,args.incomp_pene)]: #[(0,0.5,0),(0.5,0.5,0),(0,0.984,4.0),(0.941,0.984,4.0),(1,1,0)]
		a,b,c=i[0],i[1],i[2]
		pop=Population(size=args.wild,capacity=args.capa,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),
			 incomPene=c,birthsCoef=args.birth_fit,fitCoef=args.drive_fitness,inbCoef=args.inbreeding,drive=args.supMend,linkage=args.linkage)
		pop.add_pop(Population(size=args.drive,capacity=args.capa,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),
			 incomPene=c,birthsCoef=args.birth_fit,fitCoef=args.drive_fitness,inbCoef=args.inbreeding,drive=args.supMend,linkage=args.linkage)) #
		#num=[args.drive]
		ratio=[args.drive/pop.get_size()]
		for n in range(args.generation):
			pop.next_generation()
			count=pop.count_gene_drive_alleles()
			ratio.append(count/pop.get_size())
		
		print(ratio)
		x=np.append(x,[ratio],axis=0)
		
		f=open("drive_carriers_freq.femaleRate{0}_maleRate{1}_incompene{2}.txt".format(args.embryo_cutrate,args.germline_cutrate,args.incomp_pene),"w")
		f.writelines("Generation\tDriveCarriersFreq\n")
		for r in range(len(ratio)):
			f.writelines(str(r)+'\t'+str(ratio[r])+'\n')
		f.close()
		#linestyle="--"
		plt.plot(ratio,linestyle="-",label="Female DNA cleavage efficiency: "+str(format(a*100,'.1f'))+"%;\n"+"Male germline cleavage efficiency: "+str(format(b*100,'.1f'))+"%;\n"+"Incomplete penetrance: "+str(format(c*100,'.1f'))+"%") 
		#plt.plot(ratio,linestyle="-",label="Gross male germline cleavage efficiency: "+str(format(b*100,'.1f'))+"%;\n"+"Incomplete penetrance rate: "+str(c)+"%") 
		plt.axis([0, 50, 0, 1.1]);#define the plotting range
		plt.xlabel("Generation",fontsize=12)
		plt.ylabel("Frequency of drive carriers",fontsize=12)
		plt.legend(loc="lower right")

	plt.savefig("drive_carriers_freq.femaleRate{0}_maleRate{1}_incompene{2}.pdf".format(args.embryo_cutrate,args.germline_cutrate,args.incomp_pene),format="pdf")
	plt.savefig("drive_carriers_freq.femaleRate{0}_maleRate{1}_incompene{2}.png".format(args.embryo_cutrate,args.germline_cutrate,args.incomp_pene),format="png")

if __name__=="__main__":
	main()
