#!/usr/bin/env python3

import subprocess
import csv
import numpy as np

Ds = [1e-12, 5e-12, 10e-12]
iterations = 10000

for D in Ds: 
  for i in range(iterations): 
    print("running pSpatiocyte dilute diffusion model with D:%e"
        %D, "iteration:%d/%d" %(i+1, iterations))
    dirname = 'D_'+str(D)+'__'+str(i)
    result = subprocess.run(['mpirun', '-np', '8', 'dilute', dirname,
      str(D), str(i)], stdout=subprocess.PIPE)
    time = float(result.stdout.decode('utf-8').split('\n')[-2])
    subprocess.run(['python3', '../../../scripts/gather_coordinates.py',
      dirname])
    print('\telapsed time:',time)



