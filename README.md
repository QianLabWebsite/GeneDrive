# GeneDrive
Selfish genetic elements are commonly found in nature and can be transmitted to progeny at super-Mendelian (>50%) frequencies, although it could cause potertial fitness or fecundity costs to the host. Inspired by the selfish genetic elements, we constructed CAIN (CRISPR-Assisted Inheritance utilizing NPG1), a synthetic toxin-antidote gene drive specifically developed for plants, making a red fluorescent marker transmitting at frequency above 88%. 

This repository contains scripts related to this CRISPR-mediated toxin-antidote gene drive.
## 1. Calling variant from  mpileup file and output vcf file
To determine the types and efficiency of gRNA editing , we first mapped Illumina reads from the PCR product of the target region to the reference using bwa. We then used samtools to call variants from the resulting bam file, which generated mpileup files. To facilitate easy analysis of variant information, we developed a script that converts the mpileup file to a vcf-like format, for the vcf format is easy for get a detailed information of variants for eyes.

The scipt and example files are in mpileup_to_vcf directory. Use `perl mpileup_to_vcf.pl` and type enter to see the usage.
```
Description
        Calling variant from a single mpileup file and output vcf.

Usage
        perl mpileup_to_vcf.pl -i file.mpileup -o file.vcf -s sample_name -t threshuold_for_snp

Parameters
        -i <s> input file in mpileup format
        -o <s> output file in vcf-like format
        -s <s> sample name used in vcf header line
        -t <s> threshold for snp, about e-3 for illumina reads.[default 0.005]
```
## 2. Dynamics simulation of gene drive
### 1. CAIN modification drive    
The script `toxin_antidote_dynamics_simulation.py` in gene_drive_simulation directory simulate the propagation of the CRISPR-mediated toxin-antidote gene drive under different parameters such as embryo DNA cleavage efficiency, male germline cleavage efficiency, incomplete penetrance rate, and initial release ratio of drive heterozygotes. The model is stochastic, and therefore, multiple simulations are necessary to produce reliable results(uncommentting the code in the end of script will run the same parameters 100 times). The script outputs the frequency of drive carriers at each generation and can also be used to estimate the number of generations required for the gene drive to reach fixation. 
The main features of this scipt include:  
1. It is a forward genetic simulation. Set parameters to the initial population and the population will propagate. You can monitor the frequency of drive carriers along with the increasing of generations.  
2. The model is individual-based and stochastic, with population properties based on the Wright-Fisher model. This model is characterized by finite individuals, random mating, and non-overlapping generations.    
3. It uses two classes: the individual class and population class.  
```
usage: toxin_antidote_dynamics_simulation.py [-h] [-w WILD] [-d DRIVE] [-t GENERATION] -e EMBRYO_CUTRATE -g GERMLINE_CUTRATE -i INCOMP_PENE

Dynamics simulation of the CRISPR-mediated toxin-antidote gene drive.

Set parameters to the initial populations and the script will generate frequency of drive carriers along with the increasing of generations

options:
  -h, --help           show this help message and exit
  -w WILD              Size of wild population.[default 9900]
  -d DRIVE             Size of heterozygous individuals carring drive.[default 100]
  -t GENERATION        Generation of propagation.[default 50]
  -e EMBRYO_CUTRATE    Embryo DNA cleavage efficiency, float number, between 0-1
  -g GERMLINE_CUTRATE  Male germline cleavage efficiency,float number, between 0-1
  -i INCOMP_PENE       Incomplete penetrance rate, float number, between 0-1

author: Bingke Jiao
mail:   bkjiao@genetics.ac.cn
date:   2023.7.11
version:        1.0

```
The output of simulation are the dynamics curve picture of drive carriers and the frequency of drive carriers in txt file.
![gene drive simulation](https://github.com/QianLabWebsite/GeneDrive/blob/main/gene_drive_simulation/drive_carriers_freq.embryoRate0.941_germRate0.984_incompene0.04.png)


|Generation      |DriveCarriersFreq|
|----------------|-----------------|
|0       |0.01|
|1       |0.0151|
|2       |0.0193|
|3       |0.0285|
|4       |0.0422|
|5       |0.0624|
|6       |0.0947|
|7       |0.1377|
|8       |0.1954|
|9       |0.2654|
|10      |0.3625|
|11      |0.4796|
|12      |0.6082|
|13      |0.7397|
|14      |0.8505|
|15      |0.9295|
|16      |0.9749|
|17      |0.993|
|18      |0.9981|
|19      |1.0|
### 2. CAIN suppression drive
The script CAIN_suppression.py is designed to simulate the dynamics of a suppression version of the CAIN drive. The CAIN drive can be located in and therefore disrupt a fertility gene in either males or females. Individuals who are homozygous for the CAIN drive will become sterile. As the number of individuals with homozygous CAIN drive increases, the population will eventually collapse.
```
usage: CAIN_suppression.py [-h] [-w WILD] [-d DRIVE] [-t GENERATION] -e EMBRYO_CUTRATE -g GERMLINE_CUTRATE [-s SEX] [-o OUTPUT]

Dynamics simulation of the CAIN(TADS) suppression drive.

Set the parameters to the initial populations and the script will generate frequency and individual number of drive carriers along with the increasing of generation              s

options:
  -h, --help           show this help message and exit
  -w WILD              Size of wild population.[default 9900]
  -d DRIVE             Size of heterozygous individuals carring drive.[default 100]
  -t GENERATION        Generation of propagation.[default 50]
  -e EMBRYO_CUTRATE    Embryo DNA cleavage efficiency, float number, between 0-1
  -g GERMLINE_CUTRATE  Germline DNA cleavage efficiency,float number, between 0-1
  -s SEX               Set the fertility gene where drive located in, character, male or female
  -o OUTPUT            Output file prefix, character

author: Bingke Jiao
mail:   bkjiao@genetics.ac.cn
date:   2023.8.23
version:        1.0
```
![CAIN_suppression simulation](https://github.com/QianLabWebsite/GeneDrive/blob/main/gene_drive_simulation/CAIN_suppression.density_regulation_production.hermaphroditic.png)
