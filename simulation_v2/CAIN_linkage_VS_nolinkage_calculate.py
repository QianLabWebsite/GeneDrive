############## Compare the spread of CAIN when CAIN is in linkage with NPG1 gene and when they are not.
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from CAIN_viability_and_inbreeding import Population

output=sys.argv[1]

mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
plt.figure(figsize=(15,5))

CAPACITY=10000
drive_num=100
wild_num=CAPACITY-drive_num

x=np.empty(shape=[0,51],dtype=float) #51
x2=np.empty(shape=[0,51],dtype=float)

#TADS

BIRTH_FIT=2
INBREEDING=0
DRIVE_FITNESS=1
cList="rgbcmyk" #color list
cList=["royalblue","darkorange","seagreen","crimson","darkorchid","sienna","hotpink","gray","#d62728","#9467bd","firebrick"]
#inbreList= [i/10 for i in range(0, 11)]
inbreList=[0]
patch_list=[]
MaleCleave=[i/10 for i in range(0,10,2)]
cIndex=0
for i in MaleCleave: #range(10) [(1,0)]   [(0,0.5,50),(0,0.5,0),(0,0.984,4.0),(0,0.984,0),(0,1,0)] (0,0.5,4.0),(0.5,0.5,4.0),(0,0.984,4.0)    |(0,1,0)
    print(i)
    a,c=0,0
    b=i
    pop=Population(size=wild_num,capacity=CAPACITY,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),incomPene=c,
                   birthsCoef=BIRTH_FIT,fitCoef=DRIVE_FITNESS,inbCoef=INBREEDING,drive="Yes",linkage="No") 
    pop.add_pop(Population(size=drive_num,capacity=CAPACITY,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),incomPene=c,
                           birthsCoef=BIRTH_FIT,fitCoef=DRIVE_FITNESS,inbCoef=INBREEDING,drive="Yes",linkage="No")) #
    pop2=Population(size=wild_num,capacity=CAPACITY,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),incomPene=c,
                    birthsCoef=BIRTH_FIT,fitCoef=DRIVE_FITNESS,inbCoef=INBREEDING,drive="Yes",linkage="Yes") #
    pop2.add_pop(Population(size=drive_num,capacity=CAPACITY,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),incomPene=c,
                            birthsCoef=BIRTH_FIT,fitCoef=DRIVE_FITNESS,inbCoef=INBREEDING,drive="Yes",linkage="Yes")) #

    ratio=[drive_num/pop.get_size()]
    ratio2=[drive_num/pop2.get_size()]
    
    popSize=[]
    popSize2=[]
    popSize.append(pop.get_size())
    popSize2.append(pop2.get_size())
    for n in range(50): #50
        pop.next_generation()
        count=pop.count_gene_drive_alleles() #allele="True"
        ratio.append(count/pop.get_size())
        popSize.append(pop.get_size())
        
        pop2.next_generation()
        count=pop2.count_gene_drive_alleles() #allele="True"
        ratio2.append(count/pop2.get_size())
        popSize2.append(pop2.get_size())

    print(ratio)
    print(ratio2)
    print(popSize)
    print(popSize2)
    x=np.append(x,[ratio],axis=0)
    x2=np.append(x2,[ratio2],axis=0)

    i=i*10
    plt.plot(ratio,linestyle="-.",label="$\it{M_e}$:"+"{:.1f}%;No linkage".format(i/10*100),color=cList[int(i)])
    #plt.plot(x2.iloc[i],linestyle="-",label="Male germline cleavage efficiency:{:.1f}%;Linkage".format(i/10*100),color=cList[int(i/2)])
    plt.plot(ratio2,linestyle="-",label="$\it{M_e}$:"+"{:.1f}%;Linkage".format(i/10*100),color=cList[int(i)])
    patch = mpatches.Patch(color=cList[int(i)], label="$\it{M_e}$:"+"{:.1f}%".format(i/10*100))
    patch_list.append(patch)

first_legend=plt.legend(loc="upper left",bbox_to_anchor=(1.05, 0.8),handles=patch_list)
ax=plt.gca().add_artist(first_legend)

custom_lines = [Line2D([0], [0], color="black", lw=2, linestyle="-"),
                Line2D([0], [0], color="black", lw=2, linestyle="-.")]
legend2=plt.legend(custom_lines,["CAIN and target gene in linkage","CAIN and target gene not in linkage"],loc="lower left",bbox_to_anchor=(1.05, 0.8))

plt.axis([0, 50, 0, 1.1]);#define the plotting range
plt.xlabel("Generation",fontsize=12)
plt.ylabel("Frequency of CAIN carriers",fontsize=12)

plt.title("The dynamics compare of drive and target gene linked or not linked(stochastic model)\n"+
          "Penetrance rate:{}%;Fecundity of drive carriers:{}%; Viability of drive carriers:{}%\n".format(100,100,100)+
          r"$\it{M_e}$ denotes the $\bf{m}$ale germline cleavage $\bf{e}$fficiency")

plt.tight_layout()
#plt.show()
plt.savefig(output+".pdf",format="pdf") #CAIN.10-100_release_ratio.stochastic_model.not_linked_VS_linked.0108

dataframe = pd.DataFrame(x)
dataframe.to_csv(output+".nolinkage.csv",index=True,sep=',') #CAIN.1-100_release_ratio.stochastic_model.CAIN_NPG1_not_linkage.0108

dataframe = pd.DataFrame(x2)
dataframe.to_csv(output+".linkage.csv",index=True,sep=',') #CAIN.1-100_release_ratio.stochastic_model.CAIN_NPG1_linkage.0108