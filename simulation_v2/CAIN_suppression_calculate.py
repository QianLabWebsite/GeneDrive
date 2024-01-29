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
CAPACITY=10000 #10000
drive_num=100 #1000
wild_num=CAPACITY-drive_num

output=sys.argv[1] #CAIN_suppression_male_homo_sterile
i=(0,1,0)
a,b,c=i[0],i[1],i[2]
pop=Population(size=wild_num,capacity=10000,embryo_cutrate=a,germline_cutrate=b,chr1=("wt","wt"),chr2=("wt","wt"),incomPene=c,sup="male",fitCoef=1) #
pop.add_pop(Population(size=drive_num,capacity=10000,embryo_cutrate=a,germline_cutrate=b,chr1=("dr","wt"),chr2=("wt","wt"),incomPene=c,sup="male",fitCoef=1)) #

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


dataframe=pd.read_csv(output+'.csv',header=[0],index_col=[0])

end=0 
i=0
while i < len(dataframe):
    #if dataframe['freq'][i]==1.0:
    if dataframe['freq'][i]>0.99:
        end=i
    i+=1
# if end==0:
#     print("The population did not collapse within 50 generations. Please increase the number of generations.")
#     sys.exit(0)

#end=48 #for female sterile
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
#plt.figure(figsize=(10,5))
fig, axs = plt.subplots(2, 1, figsize=(10,5), sharex=True)
freq=dataframe['freq'][:end+1] #27 [:end+1] for male fertility [:end+2] for female fertility
freq_homo=dataframe['freq_homo'][:end+1]
freq_heter=dataframe['freq_heter'][:end+1]
freq=freq_homo+freq_heter
popsize=dataframe['popsize'][:end+2] #28
drive_carrier=dataframe['drive_carrier'][:end+2] #28
homo_carrier=dataframe['homo_carrier'][:end+2]
heter_carrier=dataframe['heter_carrier'][:end+2]
x=list(range(end+1)) #end+1 for male fertility
axs[0].hlines(1,0,25,color="dodgerblue")
axs[0].plot(freq,linestyle="-",label="Frequency of CAIN alleles",color="orange") #"Frequency of CAIN carriers
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
#axs[0].tick_params(axis='y', labelsize=6)
axs[1].set_ylim(0,CAPACITY*1.05)

plt.margins(x=0)
plt.xticks(range(0,end+2,2)) #end+4 for male fertility
plt.xlabel("Generation",fontsize=12)
axs[0].set_ylabel("Frequency of CAIN alleles",fontsize=10)
axs[1].set_ylabel("Number of individuals in population",fontsize=10)


yellow_patch = mpatches.Patch(color='yellow', label='CAIN heterozygotes')
orange_patch = mpatches.Patch(color='orange', label='CAIN homozygotes')
blue_patch = mpatches.Patch(color='dodgerblue', label='wild-type individuals')
plt.legend(loc="upper left",bbox_to_anchor=(1.05, 1.0),handles=[yellow_patch,orange_patch,blue_patch])
axs[0].set_title("CAIN suppresion(with density regulation production, male driver homozygote sterile)")

plt.tight_layout()

plt.savefig(output+".pdf",format="pdf")