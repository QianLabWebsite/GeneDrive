import pandas as pd
import numpy as np
import sys
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt #matplotlib is plotting library
import matplotlib as mpl
from matplotlib.pyplot import MultipleLocator
output=sys.argv[1]
gender=sys.argv[2]
dataframe=pd.read_csv(output+'.csv',header=[0],index_col=[0])

CAPACITY=100000
end=0 
if gender=="male":
    i=0
    while i < len(dataframe):
        if dataframe['freq'][i]>0.99:
            end=i
        i+=1
elif gender=="female" or gender=="both":
    end=48

mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
#plt.figure(figsize=(10,5))
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
#axs[0].hlines(1,0,25,color="dodgerblue")

#axs[0].plot(freq_homo,linestyle="-",label="Frequency of CAIN homozygotes",color="darkgoldenrod")
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

plt.savefig(output+".pdf",format="pdf")
