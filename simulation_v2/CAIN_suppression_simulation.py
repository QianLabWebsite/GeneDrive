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
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.pyplot import MultipleLocator
from CAIN_suppression import Population

__doc__='''
Dynamics simulation of the CAIN(TADS) suppression drive.

Set the parameters to the initial populations and the script will generate frequency and individual number of drive carriers along with the increasing of generations
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
    parser.add_argument('-f',help='female germline DNA cleavage efficiency, float number, between 0-1',dest='embryo_cutrate',type=float,required=True)
    parser.add_argument('-g',help='male germline DNA cleavage efficiency, float number, between 0-1',dest='germline_cutrate',type=float,required=True)
    parser.add_argument('-i',help='Incomplete penetrance rate, float number, between 0-1',dest='incomp_pene',type=float,default=0,required=False)
    parser.add_argument('-s',help='Set the fertility gene where drive located in, character, male, female or both. "both" means homozygotes are not viable.',dest='sex',type=str,default="male",required=False)
    parser.add_argument('-r',help='Relative fitness value comparing of CAIN comparing with the wild-type, float number, between 0-1',dest='drive_fitness',type=float,default='1',required=False)
    parser.add_argument('-o',help='Output file prefix, character',dest='output',type=str,default="CAIN_suppression",required=False)
    args=parser.parse_args()
    #pop.size=args.wild+args.drive
    wild_num=args.wild
    drive_num=args.drive
    capacity=args.capa

    #for i in [(args.embryo_cutrate,args.germline_cutrate,args.sex)]: #[(0,0.5,0),(0.5,0.5,0),(0,0.984,4.0),(0.941,0.984,4.0),(1,1,0)]

    a,b,c,d,e=args.embryo_cutrate,args.germline_cutrate,args.incomp_pene,args.sex,args.drive_fitness
    pop=Population(size=wild_num,capacity=capacity,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),incomPene=c,sup=d,fitCoef=e) #
    pop.add_pop(Population(size=drive_num,capacity=capacity,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),incomPene=c,sup=d,fitCoef=e)) #

    freq=[drive_num/(pop.size*2)]
    freq_homo=[0]
    freq_heter=[drive_num/pop.size]
    
    popsize=[pop.size]
    drive_carrier=[drive_num]
    homo_carrier=[0]
    heter_carrier=[drive_num]

    for n in range(args.generation):
        pop.next_generation()
        (drive_ind,homo_ind,heter_ind)=pop.count_gene_drive_alleles(allele="True")

        if pop.size == 0:
            freq.append(0)
            freq_homo.append(0)
            freq_heter.append(0)
            popsize.append(0)
            drive_carrier.append(0)
            homo_carrier.append(0)
            heter_carrier.append(0)
        else:
            freq.append(drive_ind/(pop.size*2))
            freq_homo.append(homo_ind/pop.size)
            freq_heter.append(heter_ind/pop.size)
            popsize.append(pop.size)
            drive_carrier.append(homo_ind+heter_ind)
            homo_carrier.append(homo_ind)
            heter_carrier.append(heter_ind)


    #Save the data
    dataframe = pd.DataFrame({'freq': freq,'freq_homo':freq_homo,'freq_heter':freq_heter,'drive_carrier': drive_carrier,'homo_carrier':homo_carrier,'heter_carrier':heter_carrier,'popsize': popsize})
    dataframe.to_csv(args.output+"."+args.sex+".csv",index=True,sep=',') #,mode="a"

    #Draw the picture
    gender=args.sex
    end=0 
    if gender=="male":
        i=0
        while i < len(dataframe):
            if dataframe['freq'][i]>0.99:
                end=i
            i+=1
        if end==0:
            print("The population did not collapse within 50 generations. Please increase the number of generations.")
            sys.exit(0)
    elif gender=="female" or gender=="both":
        end=48

    mpl.rcParams['pdf.fonttype']=42
    mpl.rcParams['ps.fonttype']=42
    fig, axs = plt.subplots(2, 1, figsize=(10,5), sharex=True)
    #freq=dataframe['freq'][:end+1] #27 [:end+1] for male fertility [:end+2] for female fertility
    freq_homo=dataframe['freq_homo'][:end+1]
    freq_heter=dataframe['freq_heter'][:end+1]
    freq=0
    if gender=="male":
        freq=freq_homo+freq_heter #frequency of CAIN carriers
    elif gender=="female" or gender =="both": 
        freq=dataframe['freq'][:end+1] #frequency of CAIN allele

    popsize=dataframe['popsize'][:end+2] #28
    drive_carrier=dataframe['drive_carrier'][:end+2] #28
    homo_carrier=dataframe['homo_carrier'][:end+2]
    heter_carrier=dataframe['heter_carrier'][:end+2]
    x=list(range(end+1)) #end+1 for male fertility

    if gender=="male":
        axs[0].plot(freq,linestyle="-",label="Frequency of CAIN carriers",color="orange")
        axs[0].plot(freq_heter,linestyle="-",label="Frequency of CAIN heterozygotes",color="yellow")
        axs[0].fill_between(x,freq,1,alpha=1,facecolor="dodgerblue")
        axs[0].fill_between(x,freq_heter,freq,alpha=1,facecolor="orange")
        axs[0].fill_between(x,freq_heter,alpha=1,facecolor="yellow")
    elif gender in ["female","both"]:
        axs[0].plot(freq,linestyle="-",label="Frequency of CAIN alleles",color="orange") #"Frequency of CAIN carriers

    axs[1].plot(popsize,linestyle="-",label="Number of wild-type individuals",color="dodgerblue")
    axs[1].plot(drive_carrier,linestyle="-",label="Number of CAIN carriers",color="orange")
    axs[1].plot(heter_carrier,linestyle="-",label="Number of CAIN heterozygotes",color="yellow")
    x=list(range(end+2))
    axs[1].fill_between(x,drive_carrier,popsize,alpha=1,facecolor="dodgerblue")
    axs[1].fill_between(x,heter_carrier,drive_carrier,alpha=1,facecolor="orange")
    axs[1].fill_between(x,heter_carrier,alpha=1,facecolor="yellow")
    axs[0].set_ylim(0,1.1)
    axs[1].set_ylim(0,args.capa*1.05)

    if gender=="female":
        extraticks=[0.75]
        axs[0].set_yticks(list(axs[0].get_yticks()) + extraticks)
        axs[0].axhline(0.75,linestyle="dotted")
    elif gender=="both":
        extraticks=[0.5]
        axs[0].set_yticks(list(axs[0].get_yticks()) + extraticks)
        axs[0].axhline(0.5,linestyle="dotted")
    plt.margins(x=0)
    plt.xticks(range(0,end+2,2)) #end+4 for male fertility
    plt.xlabel("Generation",fontsize=12)

    if gender=="male":
        axs[0].set_ylabel("Frequency of CAIN carriers",fontsize=10)
    elif gender in ["female","both"]:
        axs[0].set_ylabel("Frequency of CAIN allele",fontsize=10)
    axs[1].set_ylabel("Number of individuals in population",fontsize=10)


    yellow_patch = mpatches.Patch(color='yellow', label='CAIN heterozygotes')
    orange_patch = mpatches.Patch(color='orange', label='CAIN homozygotes')
    blue_patch = mpatches.Patch(color='dodgerblue', label='wild-type individuals')
    plt.legend(loc="upper left",bbox_to_anchor=(1.05, 1.0),handles=[yellow_patch,orange_patch,blue_patch])
    if gender=="male":
        axs[0].set_title("CAIN suppresion(with density regulation production, male driver homozygote sterile)")
    elif gender=="female":
        axs[0].set_title("CAIN suppresion(with density regulation production, female driver homozygote sterile)")
    elif gender=="both":
        axs[0].set_title("CAIN suppresion(with density regulation production, driver homozygote not viable)")

    plt.tight_layout()

    plt.savefig(args.output+"."+args.sex+".pdf",format="pdf")


if __name__=="__main__":
	main()
