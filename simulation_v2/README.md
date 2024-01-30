##### 1. To produce the Fig 6a, put the `CAIN_modification_male_cut_simulation.py` and the `CAIN_viability_and_inbreeding.py` in one directory, and then run   
```
python CAIN_modification_male_cut_simulation.py CAIN.Fig6a #To change the parameters, such as environment capacity, inbreeding coefficient, you can edit script before running
```
The output files are `CAIN.Fig7a.0124.csv` and `CAIN.Fig7a.0124.pdf`
##### 2. To produce the Fig 6b, put the `CAIN_modification_male_cut_simulation.py`, `homing_drive.py` and the `CAIN_homing_calculate.py` in one directory, and then run
```
python CAIN_homing_calculate.py CAIN_Homing.0124.10k
```
The output files are `CAIN_Homing.0124.100k.csv` and `CAIN_Homing.0124.100k.pdf`
##### 3. To produce the Fig 6c, put the `CAIN_suppression.py` and the `CAIN_suppression_calculate.py` in one directory, and then run
```
python CAIN_suppression_calculate.py CAIN_suppression_male_homo_sterile
```
The output files are `CAIN_suppression_male_homo_sterile.csv` and `CAIN_suppression_male_homo_sterile.pdf`
