
The classes are stored in `CAIN_viability_and_inbreeding.py`, `homing_drive.py`, `CAIN_suppression.py`, `CAIN_male_pollen_half.py`. 

Other files whose names include "calculate", "simulation" or "draw" are the implementation of the classes files.  

The files whose suffixes are ".csv" or ".pdf" are the output files.  
##### 1. To simulate the impact of male germline cleavage efficiency and incomplete penetrance on the CAIN spread(Fig 6a), put the `CAIN_modification_male_cut_simulation.py` and the `CAIN_viability_and_inbreeding.py` in one directory, and then run   
```
python CAIN_modification_male_cut_simulation.py CAIN.Fig6a #To change the parameters, such as environment capacity, inbreeding coefficient, you can edit script before running
```
The output files are `CAIN.Fig6a.csv` and `CAIN.Fig6a.pdf`
##### 2. To simulate the impact of addition of female germline cleavage efficiency on the CAIN spread(Fig 6b), put the `CAIN_viability_and_inbreeding.py`, `homing_drive.py` and the `CAIN_homing_calculate.py` in one directory, and then run
```
python CAIN_homing_calculate.py CAIN_Homing.0124.10k
```
The output files are `CAIN_Homing.0124.100k.csv` and `CAIN_Homing.0124.100k.pdf`
##### 3. To simulate the CAIN suppression, put the `CAIN_suppression.py` and the `CAIN_suppression_calculate.py` in one directory, and then run
If the CAIN is inserted to a male fertility gene(Fig 6c), run
```
python CAIN_suppression_calculate.py CAIN_suppression_male_homo_sterile #The output file is CAIN_suppression_male_homo_sterile.csv
python CAIN_suppression_draw.py CAIN_suppression_male_homo_sterile male #The output file is CAIN_suppression_male_homo_sterile.pdf
```
If the CAIN is inserted to a female fertility gene(Extended Data Figure 8), run
```
python CAIN_suppression_calculate.py CAIN_suppression_female_homo_sterile female #The output file is CAIN_suppression_female_homo_sterile.csv
python CAIN_suppression_draw.py CAIN_suppression_female_homo_sterile female #The output file is CAIN_suppression_female_homo_sterile.pdf
```
If the CAIN is inserted to a viable gene, where the CAIN homozygotes are not viable(Extended Data Figure 8), run
```
python CAIN_suppression_calculate.py CAIN_suppression_homo_nonviable both #The output file is CAIN_suppression_homo_nonviable.csv
python CAIN_suppression_draw.py CAIN_suppression_homo_nonviable both  #The output file is CAIN_suppression_homo_nonviable.pdf
```

##### 4. To simulate the effect of half pollen of CAIN heterozygous male parent on the fertility(Extended Data Figure 8), put the `CAIN_male_pollen_half.py` and `CAIN_male_pollen_half_calculate.py` in one directory and run the command below.
```
python CAIN_male_pollen_half_calculate.py CAIN.male_half_pollen.summary
```
The output files are `CAIN.male_half_pollen.summary.csv` and `CAIN.male_half_pollen.summary.pdf`
##### 5. To simulate the scenario where the CAIN is in linkage with the target gene NPG1 and compare it with the unlinkage one(Extended Data Figure 8), put the `CAIN_viability_and_inbreeding.py` and `CAIN_linkage_VS_nolinkage_calculate.py` and then run
```
python CAIN_linkage_VS_nolinkage_calculate.py CAIN_spread_linkage_VS_nolinkage
```
The output files are `CAIN_spread_linkage_VS_nolinkage.linkage.csv`, `CAIN_spread_linkage_VS_nolinkage.nolinkage.csv`,`CAIN_spread_linkage_VS_nolinkage.pdf`.
