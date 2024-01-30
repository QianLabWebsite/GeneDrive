##### 1. To produce the Fig 6a, put the `CAIN_modification_male_cut_simulation.py` and the `CAIN_viability_and_inbreeding.py` in one directory, and then run   
```
python CAIN_modification_male_cut_simulation.py CAIN.Fig6a #To change the parameters, such as environment capacity, inbreeding coefficient, you can edit script before running
```
The output files are `CAIN.Fig6a.csv` and `CAIN.Fig6a.pdf`
##### 2. To produce the Fig 6b, put the `CAIN_modification_male_cut_simulation.py`, `homing_drive.py` and the `CAIN_homing_calculate.py` in one directory, and then run
```
python CAIN_homing_calculate.py CAIN_Homing.0124.10k
```
The output files are `CAIN_Homing.0124.100k.csv` and `CAIN_Homing.0124.100k.pdf`
##### 3. To produce the Fig 6c, put the `CAIN_suppression.py` and the `CAIN_suppression_calculate.py` in one directory, and then run
```
python CAIN_suppression_calculate.py CAIN_suppression_male_homo_sterile #The output file is CAIN_suppression_male_homo_sterile.csv
python CAIN_suppression_draw.py CAIN_suppression_male_homo_sterile male #The output file is CAIN_suppression_male_homo_sterile.pdf
```
If the CAIN is inserted to a female fertility gene(Extended Data Figure 8), run
```
python CAIN_suppression_calculate.py CAIN_suppression_female_homo_sterile female #The output file is CAIN_suppression_female_homo_sterile.csv
python CAIN_suppression_draw.py CAIN_suppression_female_homo_sterile female #The output file is CAIN_suppression_female_homo_sterile.pdf
```
##### 4. To simulate the effect of half pollen of CAIN heterozygous male parent on the fertilty, run the command below.
```
python CAIN_male_pollen_half_calculate.py CAIN.male_half_pollen.summary
```
The output files are `CAIN.male_half_pollen.summary.csv` and `CAIN.male_half_pollen.summary.pdf`
