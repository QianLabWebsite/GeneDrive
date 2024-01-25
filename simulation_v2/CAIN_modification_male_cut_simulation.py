############## TADS varying male cleavage efficiency and penetrance, Fig6a, 2024.1.25
############## Inbreeding coefficient=0, Fitness of driver=1, Fertility of driver=1.
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
from matplotlib.pyplot import MultipleLocator
from CAIN_viability_and_inbreeding import Population

output=sys.argv[1] #"CAIN.Fig6a"

CAPACITY=10000 #The environment capacity of population
drive_num=100 #The number of CAIN carriers initially release to the wild population
wild_num=CAPACITY-drive_num #The number of individuals in wild population

end=50 #The reproduction generation
x=np.empty(shape=[0,end+1],dtype=float)

DRIVE_FITNESS=1 #The fithess of driver, which affects the weights of chosing individuals to be as parent.
INBREEDING=0    #The inbreeding coefcient
BIRTH_FIT=2     #The births coefcient which affects the number of offsping when the CAIN heterouzygous individual is the male parent.
for i in [(0,0.5,50),(0,0.984,4.0),(0,0.984,0),(0,1,4.0),(0,1,0)]: #Test the different combination of embryo cut rate, germline cut rate and penetrance
    print(i)
    a,b,c=i[0],i[1],i[2]
    pop=Population(size=wild_num,capacity=CAPACITY,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),incomPene=c,birthsCoef=BIRTH_FIT,fitCoef=DRIVE_FITNESS,inbCoef=INBREEDING,drive="Yes") #
    '''The parameters in the Population class include:
    1. the population size, 
    2. the environment capacity, 
    3. the embryo cut rate of Cas9, which was used as the cut rate in female, 
    4. the germline cut rate of Cas9, which was used as the cut rate in male, 
    5. the chr1 genotype, which bears the CAIN, 
    6. the chr2 genotype, which bears the NPG1, 
    7. the incomplete penetrance,
    8. the births coefficient,
    9. the fitness of CAIN allele,
    10.the inbreeding coefficent,
    11.the boolean value indicates whether CAIN functions as a driver or is simply inherited according to Mendelian laws.
    '''
    pop.add_pop(Population(size=drive_num,capacity=CAPACITY,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),incomPene=c,birthsCoef=BIRTH_FIT,fitCoef=DRIVE_FITNESS,inbCoef=INBREEDING,drive="Yes")) #
    #Add the drive carriers to the wild population.
    ratio=[drive_num/pop.get_size()]
    popSize=[]
    popSize.append(pop.get_size())
    for n in range(50): #reproduce 50 generations
        pop.next_generation()
        count=pop.count_gene_drive_alleles() 
        ratio.append(count/pop.get_size())
        popSize.append(pop.get_size())
    print(ratio)
    print(popSize)
    x=np.append(x,[ratio],axis=0)

dataframe = pd.DataFrame(x)
dataframe.to_csv(output+".csv",index=True,sep=',')


##################
#Draw the picture
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
plt.figure(figsize=(15,5))
x=pd.read_csv(output+".csv",index_col=0)
cList=["royalblue","darkorange","seagreen","crimson","darkorchid","sienna","hotpink","gray","#d62728","#9467bd","firebrick"]
n=0
for i in [(0,0.5,50),(0,0.984,4.0),(0,0.984,0),(0,1,4.0),(0,1,0)]:
    a,b,c=i[0],i[1],i[2]
    #print(i)
    #print(x.iloc[i])
    #plt.plot(x.iloc[n],linestyle="-",label="Male germline cleavage efficiency:{:.1f}%;Penetrance rate:{:.1f}%".format(b*100,100-c))
    plt.plot(x.iloc[n],linestyle="-",label="$\it{M_e}$:"+"{:.1f}%;".format(b*100)+"$\it{P_e}$:"+"{:.1f}%;".format(100-c),color=cList[n])
    plt.scatter(range(0,end+1),x.iloc[n],color=cList[n],s=10)
    n+=1
plt.axis([0, 50, 0, 1.1]);#define the plotting range
plt.xlabel("Generation",fontsize=12)
plt.ylabel("Frequency of CAIN carriers",fontsize=12)

plt.legend(loc="center right")
plt.tight_layout()
plt.savefig(output+".pdf",format="pdf")

