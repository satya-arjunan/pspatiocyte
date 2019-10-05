#!/usr/bin/env python3

import subprocess
import csv
import numpy as np

print("running pSpatiocyte reversible reaction model")
result = subprocess.run(['mpirun', '-np', '8', 'reversible'
  ], stdout=subprocess.PIPE)
time = float(result.stdout.decode('utf-8').split('\n')[-2])
subprocess.run(['python3', '../../../scripts/gather_timecourse.py'])
print('\telapsed time:',time)



