#!/usr/bin/env python3

import subprocess
import csv
import numpy as np

ks = [0.135, 1.35, 13.5]

for k in ks: 
  print("running pSpatiocyte irreversible reaction model with k:%f" %k)
  dirname = 'k_'+str(k)
  result = subprocess.run(['mpirun', '-np', '8', 'irreversible', dirname,
    str(k)], stdout=subprocess.PIPE)
  time = float(result.stdout.decode('utf-8').split('\n')[-2])
  subprocess.run(['python3', '../../../scripts/gather_timecourse.py', dirname])
  print('\telapsed time:',time)




