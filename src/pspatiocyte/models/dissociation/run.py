#!/usr/bin/env python3

import subprocess
import csv
import numpy as np

searches = [0, 1]
fractions = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
iterations = 100

for search in searches:
  for fraction in fractions: 
    for i in range(iterations): 
      print("running pSpatiocyte dissociation model with forced search:%d"
          %search, "crowding fraction:%.1f" %fraction,
          "iteration:%d/%d" %(i+1, iterations))
      fraction_str = '{:.2f}'.format(fraction)
      dirname = 'forced_'+str(search)+'__fraction_'+fraction_str+'_'+str(i)
      result = subprocess.run(['mpirun', '-np', '8', 'dissociation', dirname,
        str(search), str(fraction)], stdout=subprocess.PIPE)
      time = float(result.stdout.decode('utf-8').split('\n')[-2])
      print('\telapsed time:',time)




