############## TADS varing male germline cleavage efficiency and fixed penetrance rate, to see the dynamics of drive spread
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
from matplotlib.pyplot import MultipleLocator
from CAIN_male_pollen_half import Population


output=sys.argv[1]

CAPACITY=10000
drive_num=100
wild_num=CAPACITY-drive_num

for j in range(1):
    x=np.empty(shape=[0,51],dtype=float) #51
    fitList=[1,1.2,1.4,1.6,1.8,2]
    for varyFit in fitList:
        print(varyFit)
        male_cut=[1]
        for i in male_cut: #range(10) [(1,0)]   [(0,0.5,50),(0,0.984,4.0),(0,0.984,0),(0,1,4.0),(0,1,0)]    |(0,1,0)
            #print(i)
            b=i
            a,c=0,0
            pop=Population(size=wild_num,capacity=CAPACITY,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),incomPene=c,varyFit=varyFit) #
            pop.add_pop(Population(size=drive_num,capacity=drive_num,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),incomPene=c,varyFit=varyFit)) #
            ratio=[drive_num/pop.get_size()]
            #ratio=[drive_num/(pop.size*2)]
            popSize=[]
            popSize.append(pop.get_size())
            for n in range(50): #50
                pop.next_generation()
                count=pop.count_gene_drive_alleles() #allele="True"
                ratio.append(count/pop.get_size())
                popSize.append(pop.get_size())
            print(ratio)
            print(popSize)
            x=np.append(x,[ratio],axis=0)

    dataframe = pd.DataFrame(x)
    #将DataFrame存储为csv,index表示是否显示行名，default=True
    dataframe.to_csv(output+".csv",index=True,sep=',')

mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
plt.figure(figsize=(15,5))

#x=pd.read_csv("CAIN.1-100_release_ratio.male_half_pollen.0121.0.csv",index_col=0)
x=pd.read_csv(output+".csv",index_col=0)
cList=["royalblue","darkorange","seagreen","crimson","darkorchid","sienna","hotpink","gray","#d62728","#9467bd","firebrick"]
n=0
for varyFit in [1,1.2,1.4,1.6,1.8,2]:
    print(varyFit)
    #print(x.iloc[i])
    plt.plot(x.iloc[n],linestyle="-",label="Male driver fertility:{:.1f}".format(varyFit*0.5),color=cList[n])
    n+=1
    
plt.axis([0, 50, 0, 1.1]);#define the plotting range
plt.xlabel("Generation",fontsize=12)
plt.ylabel("Frequency of CAIN carriers",fontsize=12)

b=1
a,c=0,0
plt.title("Male germline cleavage efficiency:{:.1f}%;Penetrance rate:{:.1f}%".format(b*100,100-c)) 
plt.legend(loc="upper left",bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()
#plt.show()
plt.savefig(output+".pdf",format="pdf")