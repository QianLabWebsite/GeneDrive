# GeneDrive
Selfish genetic elements are commonly found in nature and can be transmitted to progeny at super-Mendelian (>50%) frequencies, although it could cause potertial fitness or fecundity costs to the host. Inspired by the selfish genetic elements, we constructed CAIN (CRISPR-Assisted Inheritance utilizing NPG1), a synthetic toxin-antidote gene drive specifically developed for plants, making a red fluorescent marker transmitting at frequency above 88%. 

This repository contains scripts related to this CRISPR-mediated toxin-antidote gene drive.
## 1. Calling variant from a single mpileup file and output vcf
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
This script in gene_drive_simulation implements the dynamic simulation of CRISPR-mediated toxin-antidote gene drive.  
The main features of this scipt include:  
1. It is a forward genetic simulation. Set parameters to the initial population and the population will propagate. You can monitor the frequency of drive carriers with the increasing of generations.  
2. It is an individual-based, stocastic model.The population is in Wright-Fisher model.    
3. It uses two classes: the individual class and population class.  
```
usage: simulate.py [-h] [-w WILD] [-d DRIVE] [-t GENERATION] -e EMBRYO_CUTRATE -g GERMLINE_CUTRATE -i INCOMP_PENE

Dynamics simulation of the CRISPR-mediated toxin-antidote gene drive.

Set the parameters to the initial populations and the script will generate frequency of drive carriers along with the increasing of generations

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
gene_drive_simulation
![gene drive simulation](https://github.com/QianLabWebsite/GeneDrive/blob/main/gene_drive_simulation/drive_carriers_freq.embryoRate0.941_germRate0.984_incompene0.04.pdf)https://github.com/QianLabWebsite/GeneDrive/blob/main/gene_drive_simulation/drive_carriers_freq.embryoRate0.941_germRate0.984_incompene0.04.png)
