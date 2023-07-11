# GeneDrive
Selfish genetic elements are commonly found in nature and can be transmitted to progeny at super-Mendelian (>50%) frequencies, although it could cause potertial fitness or fecundity costs to the host. Inspired by the selfish genetic elements, we constructed CAIN (CRISPR-Assisted Inheritance utilizing NPG1), a synthetic toxin-antidote gene drive specifically developed for plants, making a red fluorescent marker transmitting at frequency above 88%. 

This repository contains scripts related to this CRISPR-mediated toxin-antidote gene drive.
## 1. Converting the mpileup to vcf-like file
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

