#TADS suppression
import pandas as pd
import numpy as np
import sys
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
from matplotlib.pyplot import MultipleLocator
from CAIN_suppression import Population

mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
CAPACITY=100000 #10000
drive_num=1000 #1000
wild_num=CAPACITY-drive_num

output=sys.argv[1] #CAIN_suppression_male_homo_sterile
i=(0,1,0)
a,b,c=i[0],i[1],i[2]
pop=Population(size=wild_num,capacity=CAPACITY,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),incomPene=c,sup="male",fitCoef=1) #
pop.add_pop(Population(size=drive_num,capacity=CAPACITY,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),incomPene=c,sup="male",fitCoef=1)) #

freq=[drive_num/(pop.size*2)]
freq_homo=[0]
freq_heter=[drive_num/pop.size]
popsize=[pop.size]
drive_carrier=[drive_num]
homo_carrier=[0]
heter_carrier=[drive_num]
for n in range(50):
    pop.next_generation()
    #print("This is "+str(n)+" generation")
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

print(freq)
print(freq_homo)
print(freq_heter)
print(drive_carrier)
print(homo_carrier)
print(heter_carrier)
print(popsize)

dataframe = pd.DataFrame({'freq': freq,'freq_homo':freq_homo,'freq_heter':freq_heter,'drive_carrier': drive_carrier,'homo_carrier':homo_carrier,'heter_carrier':heter_carrier,'popsize': popsize})
dataframe.to_csv(output+".csv",index=True,sep=',')


