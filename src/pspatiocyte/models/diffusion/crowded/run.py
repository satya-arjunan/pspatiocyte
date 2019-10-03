#!/usr/bin/env python3

import subprocess
import csv
import numpy as np

fractions = [0.0, 0.3, 0.5]
iterations = 1000

for fraction in fractions: 
  for i in range(iterations): 
    print("running pSpatiocyte crowded diffusion model with volume fraction:%e"
        %fraction, "iteration:%d/%d" %(i+1, iterations))
    dirname = 'fraction_'+str(fraction)+'__'+str(i)
    result = subprocess.run(['mpirun', '-np', '8', 'crowded', dirname,
      str(fraction), str(i)], stdout=subprocess.PIPE)
    time = float(result.stdout.decode('utf-8').split('\n')[-2])
    subprocess.run(['python3', '../../../scripts/gather_coordinates.py',
      dirname])
    print('\telapsed time:',time)



