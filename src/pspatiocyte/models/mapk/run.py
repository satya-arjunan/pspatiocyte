#!/usr/bin/env python3

import subprocess
import csv
import numpy as np

ratios = np.logspace(-1.5,1.5,12)
Ds = [0.06, 4.0]
data = []

for D in Ds:
  for ratio in ratios: 
    print("running pSpatiocyte MAPK model with diffusion coefficient, D:", D,
        "ratio:", ratio)
    ratio_str = '{:.4f}'.format(ratio)
    dirname = 'D_'+str(D)+'__ratio_'+ratio_str
    result = subprocess.run(['mpirun', '-np', '8', 'mapk', dirname, str(D),
      str(ratio)], stdout=subprocess.PIPE)
    time = float(result.stdout.decode('utf-8').split('\n')[-2])
    data.append([D, ratio, time])
    print('\telapsed time:',time)

with open("elapsed_time.txt", "w+") as f:
  writer = csv.writer(f)
  writer.writerow(['D', 'ratio', 'elapsed_time'])
  for row in data:
    writer.writerow(row)

